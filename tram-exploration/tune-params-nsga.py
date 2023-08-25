import tconsumption
from tramparse import tram_csv_parser
import ruptures as rpt
from fastdtw import fastdtw
import random
from tabulate import tabulate

# defines
BETTER = 0
WORSE = 1
NOTHING = 2

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
    def __init__(self, functions: list, params: dict):
        self.functions = functions
        self.params = params
        self.S_dominated = []
        self.n_dominating = 0
        self.new_zero = False
        self.snap_coords_back()
        
    def find_evals(self):
        self.evals = []
        for f in self.functions:
            self.evals.append(f.get_vals(self.params))
        
    def snap_coords_back(self):
        for p in self.functions[0].bounds_and_steps:
            if self.params[p["name"]] < p["min"]:
                self.params[p["name"]] = p["min"]
            if self.params[p["name"]] > p["max"]:
                self.params[p["name"]] = p["max"]

        self.find_evals()

    def is_new_zero(self):
        to_return = self.new_zero
        self.new_zero = False
        return to_return

class RandomIndividual(Individual):
    def __init__(self, functions: list):
        Individual.__init__(self, functions, functions[0].get_random_params())

class Optimizer:
    def __init__(self, functions: list, individual_count=10, rounds = 100):
        # Introduction
        print("⚙️  Picked settings:")
        table = table = [["Individual count", individual_count], ["Rounds", rounds], ["Number of functions", len(functions)]]
        print(tabulate(table, tablefmt="rounded_grid"))        

        print("⏳ Initializing...", end="", flush=True)
        self.individual_count = individual_count
        self.rounds = rounds
        self.functions = functions
        self.individuals = self.generate_first_population()
        self.pareto_ranks = []
        print("Done")

    def generate_first_population(self):
        return [RandomIndividual(self.functions) for _ in range(self.individual_count)]

    def compare_evals(self, ind_a, ind_b):
        better_count = 0
        worse_count = 0
        for i,eval_a in enumerate(ind_a.evals):
            eval_b = ind_b.evals[i]
            if eval_a < eval_b:
                better_count += 1
            elif eval_a > eval_b:
                worse_count += 1

        if better_count > worse_count:
            return BETTER
        elif better_count < worse_count:
            return WORSE
        else:
            return NOTHING

    # def compare_evals(self, ev_a, ev_b):
    #     if ev_a[0] < ev_b[0] and ev_a[1] > ev_b[1] or ev_a[0] > ev_b[0] and ev_a[1] < ev_b[1]:
    #         return NOTHING
    #     if ev_a[0] < ev_b[0] or ev_a[1] < ev_b[1]:
    #         return BETTER
    #     elif ev_a[0] > ev_b[0] or ev_a[1] > ev_b[1]:
    #         return WORSE
    #     else:
    #         return NOTHING

    def get_ranks(self):
        q1_rank = []
        for ind in self.individuals:
            if ind.n_dominating == 0:
                q1_rank.append(ind)
        self.pareto_ranks.append(q1_rank)

        for ind in self.pareto_ranks[-1]:
            for dominated in ind.S_dominated:
                dominated.n_dominating -= 1
                if dominated.n_dominating == 0:
                    dominated.new_zero = True
        self.pareto_ranks.append([x for x in self.individuals if x.is_new_zero()])
        self.pareto_ranks.append([x for x in self.individuals if x not in self.pareto_ranks[0] and x not in self.pareto_ranks[1]])

        # Clear metadata
        for ind in self.individuals:
            ind.n_dominating = 0
            ind.S_dominated = []

        # print("pareto thick:", [len(x) for x in self.pareto_ranks])

    def crossover_individuals(self, ind_a: Individual, ind_b: Individual):
        new_params = {}
        if random.uniform(0,1) < .5:
            for p in ind_a.params:
                new_params[p] = (ind_a.params[p] + ind_b.params[p])/2
        else:
            for p in ind_a.params:
                new_params[p] = (ind_a.params[p] - ind_b.params[p])/2
        return Individual(self.functions, new_params)

    def mutate_individual(self, ind: Individual):
        new_params = {}
        for p in ind.params:
            new_params[p] = ind.params[p] + random.uniform(0,1)
        return Individual(self.functions, new_params)

    def run(self):
        for r in range(self.rounds):
            # Create new individuals
            new_individuals = []
            for _ in range(0, self.individual_count):
                random_choice = random.randint(0,2)
                # Create random new
                if random_choice == 0:
                    new_individuals.append(RandomIndividual(self.functions))
                # Crossover two existing
                elif random_choice == 1:
                    random_parent_a = self.individuals[random.randint(0, self.individual_count-1)]
                    while True:
                        random_parent_b = self.individuals[random.randint(0, self.individual_count-1)]
                        if random_parent_a != random_parent_b:
                            break
                    new_individuals.append(self.crossover_individuals(random_parent_a, random_parent_b))
                # Mutate existing
                else:
                    random_parent = self.individuals[random.randint(0, self.individual_count-1)]
                    new_individuals.append(self.mutate_individual(random_parent))
            self.individuals += new_individuals
            
            # # Evaluate parameters in individuals
            # for ind in self.individuals:
            #     if ind.evals is None:
            #         ind.evals = (self.f1.evalutation(ind.r, ind.h), self.f2.evalutation(ind.r, ind.h))

            # Compare evaluations
            for ind_a in self.individuals:
                for ind_b in self.individuals:
                    if ind_a == ind_b:
                        continue
                    comparison = self.compare_evals(ind_a, ind_b)
                    if comparison == BETTER:
                        ind_a.S_dominated.append(ind_b)
                    elif comparison == WORSE:
                        ind_a.n_dominating += 1
            # Create pareto ranks
            self.get_ranks()
            best_individuals = []
            counter = 0
            for rank in self.pareto_ranks:
                for ind in rank:
                    best_individuals.append(ind)
                    counter += 1
                    if counter == self.individual_count:
                        break
            self.individuals = best_individuals
            self.pareto_ranks.clear()
            
            # Best individual
            best_ind = min(self.individuals, key=lambda x: sum(x.evals))

            # Console info
            print(f"⏲️️  Round {r+1}")
            print("✊ Currently best individual:")
            print("  Params:", best_ind.params)
            print("  Evals:")
            for i,e in enumerate(best_ind.evals):
                print(f"    {i}: {e}")
            print("")
        print("🫅 Best evals:")
        for i,e in enumerate(best_ind.evals):
            print(f"  {i}: {e}")

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
    functions = [
        ConsumptionFunction(
            bounds_and_steps,
            "../testing-data/DPO/Jízda_Poruba_Zátiší_20.04.2022/01-Vozovna-Vresinska_2022-04-20.csv",
            "../testing-data/DPO/Jízda_Poruba_Zátiší_20.04.2022/01-Vozovna-Vresinska_2022-04-20.geojson"
        ),
        ConsumptionFunction(
            bounds_and_steps,
            "../testing-data/DPO/Jízda_Poruba_Zátiší_20.04.2022/02-Vresinska-Zatisi_2022-04-20.csv",
            "../testing-data/DPO/Jízda_Poruba_Zátiší_20.04.2022/02-Vresinska-Zatisi_2022-04-20.geojson"
        ),
        ConsumptionFunction(
            bounds_and_steps,
            "../testing-data/DPO/Jízda_Poruba_Zátiší_20.04.2022/03-Zatisi-Vresinska_2022-04-20.csv",
            "../testing-data/DPO/Jízda_Poruba_Zátiší_20.04.2022/03-Zatisi-Vresinska_2022-04-20.geojson"
        ),
        ConsumptionFunction(
            bounds_and_steps,
            "../testing-data/DPO/Jizdy_Centrum_07-08-12_07.2022/1710_02.csv",
            "../testing-data/DPO/Jizdy_Centrum_07-08-12_07.2022/1710_02.geojson"
        )
    ]
    opti = Optimizer(functions, 10, 1000)
    opti.run()