import tconsumption
from tramparse import tram_csv_parser
import ruptures as rpt
import os
from fastdtw import fastdtw
import multiprocessing

DATA_PATH="../../enet-sz-data/real_rides/"

REF_PAIRS = [
    (
        DATA_PATH + "DPO/Jízda_Poruba_Zátiší_20.04.2022/01-Vozovna-Vresinska_2022-04-20.csv",
        DATA_PATH + "DPO/Jízda_Poruba_Zátiší_20.04.2022/01-Vozovna-Vresinska_2022-04-20.geojson"
    ),
    (
        DATA_PATH + "DPO/Jízda_Poruba_Zátiší_20.04.2022/02-Vresinska-Zatisi_2022-04-20.csv",
        DATA_PATH + "DPO/Jízda_Poruba_Zátiší_20.04.2022/02-Vresinska-Zatisi_2022-04-20.geojson"
    ),
    (
        DATA_PATH + "DPO/Jízda_Poruba_Zátiší_20.04.2022/03-Zatisi-Vresinska_2022-04-20.csv",
        DATA_PATH + "DPO/Jízda_Poruba_Zátiší_20.04.2022/03-Zatisi-Vresinska_2022-04-20.geojson"
    ),
    (
        DATA_PATH + "DPO/Jizdy_Centrum_07-08-12_07.2022/1710_02.csv",
        DATA_PATH + "DPO/Jizdy_Centrum_07-08-12_07.2022/1710_02.geojson"
    )
]

class MPCaching(multiprocessing.Process):
    def __init__(self, csv_path, geojson_path):
        super(MPCaching, self).__init__()
        self.csv_path = csv_path
        self.geojson_path = geojson_path
        self.velocity_segment_cache = multiprocessing.Queue()
        self.station_cache = multiprocessing.Queue()
                 
    def run(self):
        # Measured data
        df = tram_csv_parser(self.csv_path)

        # Calculated data
        c = tconsumption.Consumption()
        c.load_from_file(self.geojson_path)

        # Inferring stops
        stations_stopped = []
        last_row_zero = 0
        for i,row in enumerate(df["tm_rychlost_3_napravy"]):
            if row < 1:
                if last_row_zero < i-1:
                    stations_stopped.append(i)
                last_row_zero = i
        c.stations = [c.stations[0]] + stations_stopped + [c.stations[-1]]
        self.station_cache.put(c.stations.copy())

        # Getting velocity segments
        algo = rpt.Pelt(model="rbf").fit(df["tm_rychlost_3_napravy"].values)
        result = algo.predict(pen=100)
        self.velocity_segment_cache.put(result.copy())

class MPTuning(multiprocessing.Process):
    def __init__(self, csv_path, geojson_path):
        super(MPTuning, self).__init__()
        self.csv_path = csv_path
        self.geojson_path = geojson_path
        self.working_params = None
        self.velocity_segment_cache = None
        self.station_cache = None
        self.distance_queue = multiprocessing.Queue()
                 
    def run(self):
        # Measured data
        df = tram_csv_parser(self.csv_path)

        # Calculated data
        c = tconsumption.Consumption()

        # Set variables
        c.params["mass_locomotive"] = 34500
        c.params["power_limit"] = 600*1000

        # Variables to tune
        c.variable_params = self.working_params

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
        distance = fastdtw(energy_calculated, energy_real)[0]
        self.distance_queue.put(distance)

cur_tp_idx = 0
new_tp = True
def tweak_params(to_tweak, working_params, best_params):
    global cur_tp_idx, new_tp

    t = to_tweak[cur_tp_idx]
    if new_tp:
        new_tp = False
        working_params[t["name"]] = 0
    
    if working_params[t["name"]] >= t["max"]:
        working_params[t["name"]] = best_params[t["name"]]
        cur_tp_idx += 1
        new_tp = True
        if cur_tp_idx >= len(to_tweak):
            cur_tp_idx = 0
            return None
        return tweak_params(to_tweak, working_params, best_params)

    working_params[t["name"]] += t["step"]

    return t["name"]
    # for k in to_tweak:
    #     if not k["done"]:
    #         if working_params[k["name"]] < k["max"]:  
    #             working_params[k["name"]] += k["step"]
    #             return k
    #         working_params[k["name"]] = best_params[k["name"]]
    #         k["done"] = True
    # return None

if __name__ == "__main__":
    station_cache = {}
    velocity_segment_cache = {}
    working_params = {
        "Elevation smoothing": 100,
        "Curve smoothing": 10,
        "Curve A": 650,
        "Curve B": 55,
        "Running a": 1.35,
        "Running b": 0.0008,
        "Running c": 0.00033,
        "Recuperation coefficient": 0.5,
        "Comfortable acceleration": 0.89,
        "Compensation polynomial": None
    }
    to_tweak = [
        {
            "name": "Recuperation coefficient",
            "max": 1,
            "step": 0.01,
            "done": False
        },
        {
            "name": "Curve A",
            "max": 1000,
            "step": 1,
            "done": False
        },
        {
            "name": "Curve B",
            "max": 100,
            "step": 1,
            "done": False
        },
        {
            "name": "Running a",
            "max": 2,
            "step": 0.01,
            "done": False
        },
        {
            "name": "Running b",
            "max": 0.1,
            "step": 0.0001,
            "done": False
        },
        {
            "name": "Running c",
            "max": 0.01,
            "step": 0.00001,
            "done": False
        },
        {
            "name": "Comfortable acceleration",
            "max": 1,
            "step": 0.01,
            "done": False
        }
    ]
    best_distance = None
    best_params = None
    prev_best_distance = None

    # Generate caches
    mpcs = [MPCaching(x[0], x[1]) for x in REF_PAIRS]
    # mpcs = [MPCaching(REF_PAIRS[0][0], REF_PAIRS[0][1])]
    for m in mpcs:
        m.start()
    for m in mpcs:
        m.join()
        m.velocity_segment_cache = m.velocity_segment_cache.get()
        m.station_cache = m.station_cache.get()

    # RUN STARTS HERE
    opt_round = 1
    while True:
        distances = []
        tweak_name = tweak_params(to_tweak, working_params, best_params)
        if tweak_name is None:
            if prev_best_distance is not None and best_distance >= prev_best_distance:
                break
            prev_best_distance = best_distance
            continue
        mps = []
        for mc in mpcs:
            m = MPTuning(mc.csv_path, mc.geojson_path)
            m.velocity_segment_cache = mc.velocity_segment_cache
            m.station_cache = mc.station_cache
            m.working_params = working_params
            m.start()
            mps.append(m)
        for m in mps:
            m.join()
            distances.append(m.distance_queue.get())

        print(f"Round {opt_round} (tweaking '{tweak_name}')")
        opt_round += 1
        print("Params:", working_params)
        print(f"😎 Best e. distance (previously): {best_distance}")
        avg_distance = sum(distances)/len(distances)
        print(f"⚙️ Current e. distance: {avg_distance}\n")
        if best_distance is None or avg_distance < best_distance:
            best_distance = avg_distance
            best_params = working_params.copy()
    
    print(f"Best params:")
    for key in best_params.keys():
        print(f"\t{key}: {best_params[key]}")
