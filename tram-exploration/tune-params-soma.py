import tconsumption
from tramparse import tram_csv_parser
import ruptures as rpt
from fastdtw import fastdtw
import random
from copy import deepcopy
import numpy as np

class ConsumptionFunction:
    def __init__(
            self, bounds_and_steps,
            csv_path, geojson_path,
            mass_locomotive=34500, mass_wagon=0, power_limit=600*1000,
            elevation_smoothing=100, curve_smoothing=10
        ):
        self.bounds_and_steps = bounds_and_steps
        self.csv_path = csv_path
        self.geojson_path = geojson_path
        self.mass_locomotive = mass_locomotive
        self.mass_wagon = mass_wagon
        self.power_limit = power_limit
        self.elevation_smoothing = elevation_smoothing
        self.curve_smoothing = curve_smoothing

        self.velocity_segment_cache = None
        self.station_cache = None
        
    def get_vals(self, params):
        # Measured data
        df = tram_csv_parser(self.csv_path)

        # Calculated data
        c = tconsumption.Consumption()

        # Set variables
        c.params["mass_locomotive"] = self.mass_locomotive
        c.params["mass_wagon"] = self.mass_wagon
        c.params["power_limit"] = self.power_limit

        # Variables to tune
        c.variable_params = params

        # Fixed variables
        c.variable_params["Elevation smoothing"] = self.elevation_smoothing
        c.variable_params["Curve smoothing"] = self.curve_smoothing

        c.load_from_file(self.geojson_path)

        # Inferring stops
        if self.station_cache is None:
            stations_stopped = []
            last_row_zero = 0
            for i,row in enumerate(df["tm_rychlost_3_napravy"]):
                if row < 1:
                    if last_row_zero < i-1:
                        stations_stopped.append(i)
                    last_row_zero = i
            c.stations = [c.stations[0]] + stations_stopped + [c.stations[-1]]
            self.station_cache = c.stations.copy()
        else:
            c.stations = self.station_cache.copy()

        # Getting velocity segments
        if self.velocity_segment_cache is None:
            algo = rpt.Pelt(model="rbf").fit(df["tm_rychlost_3_napravy"].values)
            result = algo.predict(pen=100)
            self.velocity_segment_cache = result.copy()
        else:
            result = self.velocity_segment_cache.copy()

        to_remove = []
        for i,r in enumerate(result):
            if i == 0:
                continue
            slice = df["tm_rychlost_3_napravy"][result[i-1]:r]
            if len(slice)/4 < len([x for x in slice if x <= 0]):
                to_remove.append(i-1)
        to_remove.reverse()
        for i in to_remove:
            del result[i]

        max_velocities = []
        for i,r in enumerate(result):
            if i == 0:
                segment = df["tm_rychlost_3_napravy"][:r]
            else:
                segment = df["tm_rychlost_3_napravy"][result[i-1]:r]
            max_velocities += [max(segment)] * len(segment)

        c.max_velocities_in_mps = max_velocities

        # Running the simulation
        c.run()

        # Compare two sets
        energy_calculated = [x/3600000 for x in c.series["energy_from_exerted_force"]]
        energy_real = df["cela_energie_vyrovnana_bez_spotrebicu"]
        return fastdtw(energy_calculated, energy_real)[0]
    
    def get_random_params(self):
        random_params = {}
        for b in self.bounds_and_steps:
            random_params[b["name"]] = random.uniform(b["min"], b["max"])
        return random_params

class Individual:
    def __init__(self, function: ConsumptionFunction, params: dict):
        self.function = function
        self.params = params
        self.find_eval()
        
    def find_eval(self):
        self.eval = self.function.get_vals(self.params)
        
    def snap_coords_back(self):
        for p in self.function.bounds_and_steps:
            if self.params[p["name"]] < p["min"]:
                self.params[p["name"]] = p["min"]
            if self.params[p["name"]] > p["max"]:
                self.params[p["name"]] = p["max"]

        self.find_eval()

class RandomIndividual(Individual):
    def __init__(self, function: ConsumptionFunction):
        Individual.__init__(self, function, function.get_random_params())

class SOMA:
    def __init__(self, migrations, pop_size, prt, step, path_length):
        self.migrations = migrations
        self.pop_size = pop_size
        self.prt = prt
        self.step = step
        self.path_length = path_length

        # Introduction
        print("=============================================")
        print("Picked settings:")
        print("=============================================")
        print("Migrations:", migrations)
        print("Population size:", pop_size)
        print("PRT:", prt)
        print("Path length:", path_length)
        print("Step:", step)
        print("=============================================\n")

    def run(self, function: ConsumptionFunction):
        population = [RandomIndividual(function) for _ in range(self.pop_size)]
        best_ind = min(population, key=lambda ind: ind.eval)

        for i in range(self.migrations):
            print(f"⏲️️  Migration {i+1}")
            new_population = []
            for ind_idx, ind in enumerate(population, start=1):
                print(f"  🏃 Individual {ind_idx}")
                t = 0
                best_pos = ind
                while t <= self.path_length:
                    prt_vector = []
                    for j in range(len(ind.params)):
                        if np.random.uniform() < self.prt:
                            prt_vector.append(1)
                        else:
                            prt_vector.append(0)
                    next_pos_params = {}
                    for j,p in enumerate(function.bounds_and_steps):
                        param_name = p["name"]
                        next_pos_params[param_name] = ind.params[param_name]+(best_ind.params[param_name]-ind.params[param_name])*t*prt_vector[j]
                    for j,p in enumerate(ind.params.keys()):
                        next_pos_params[p] = ind.params[p]+(best_ind.params[p]-ind.params[p])*t*prt_vector[j]
                    next_pos = Individual(function, next_pos_params)
                    next_pos.snap_coords_back()
                    if next_pos.eval < best_pos.eval:
                        best_pos = next_pos
                    
                    t += self.step
                new_population.append(best_pos)
            population = deepcopy(new_population) 
            best_ind = min(population, key=lambda ind: ind.eval)
            print("✊  Currently best individual:", best_ind.eval, "\n")
        print("🫅  Best params:")
        print(best_ind.params)
        return best_ind

if __name__ == "__main__":
    bounds_and_steps = [
        {
            "name": "Recuperation coefficient",
            "min": 0,
            "max": 1,
            "step": 0.01,
        },
        {
            "name": "Curve A",
            "min": 1,
            "max": 1000,
            "step": 1,
        },
        {
            "name": "Curve B",
            "min": 0,
            "max": 100,
            "step": 1,
        },
        {
            "name": "Running a",
            "min": 0,
            "max": 2,
            "step": 0.01,
        },
        {
            "name": "Running b",
            "min": 0,
            "max": 0.1,
            "step": 0.0001,
        },
        {
            "name": "Running c",
            "min": 0,
            "max": 0.01,
            "step": 0.00001,
        },
        {
            "name": "Comfortable acceleration",
            "min": 0,
            "max": 1,
            "step": 0.01,
        }
    ]
    soma = SOMA(migrations=20, pop_size=20, prt=1, step=1, path_length=10)
    csf = ConsumptionFunction(
        bounds_and_steps=bounds_and_steps,
        csv_path="../testing-data/DPO/Jízda_Poruba_Zátiší_20.04.2022/01-Vozovna-Vresinska_2022-04-20.csv",
        geojson_path="../testing-data/DPO/Jízda_Poruba_Zátiší_20.04.2022/01-Vozovna-Vresinska_2022-04-20.geojson"
    )
    soma.run(csf)