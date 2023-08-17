import pandas as pd
import numpy as np
import json

def remove_leader(path: str):
    return path.split(".")[-1]

def convert_sig_unsig(number):
    if number < 32767:
        return number
    else:
        return number-65535

def tram_csv_parser(csv_path):
    # Load from csv
    df = pd.read_csv(csv_path, delimiter=";", skiprows=73, decimal=",")
    default_cols = list(df.columns)

    # Rename cols
    rename_dict = {}
    for col in df.columns:
        new_name = remove_leader(col)
        rename_dict[col] = new_name
    df.rename(columns=rename_dict, inplace=True)

    # New cols (calculated)
    # Time
    df.insert(2, "DateTime", pd.to_datetime(df["Date"] + ' ' + df["Time"], format="%d.%m.%Y %H:%M:%S:%f"))
    del df["Date"]
    del df["Time"]
    df["delta_T"] = df["DateTime"].diff(periods=-1).abs()

    # Trakcni menic 1
    df["tm_stejnoproud_trakce"] = df["AM_A_C1_14_IwTCU_ILF1"].apply(convert_sig_unsig)/10
    df["tm_stejnoproud_spotrebice"] = df["AM_A_C1_14_IwAUX_I_DClink"]*0.1
    df["tm_stridaproud_vystupu"] = df["AM_A_C1_14_IwAUX_I_rms"]*0.01
    df["tm_vstupni_napeti"] = df["AM_A_C1_14_IwTCU_UD1Voltage"]/10
    df["tm_tah_vyvijeny_motory"] = df["AM_B_C1_15_IwTCU_E_MC1"].apply(convert_sig_unsig)/100
    df["tm_rychlost_3_napravy"] = df["VrTCU_A_A1_Velocity"]/3.6
    df["tm_prikon_trakce"] = df["tm_stejnoproud_trakce"]*df["tm_vstupni_napeti"]/1000
    df["tm_vykon_trakce"] = df["tm_tah_vyvijeny_motory"]*df["tm_rychlost_3_napravy"]
    df["tm_prikon_spotrebice"] = df["AM_A_C1_14_IwTCU_P_AUX"]/100
    df["tm_prikon_topeni_klima"] = 1.73*400*df["tm_stridaproud_vystupu"]/1000
    df["tm_vykon_zmareny"] = np.where((df["tm_vykon_trakce"] > 0) & (df["tm_vykon_trakce"] < df["tm_prikon_trakce"]), df["tm_prikon_trakce"]-df["tm_vykon_trakce"], 0)

    # Trakcni menic 2
    df["tm2_stejnoproud_trakce"] = df["AM_B_C1_15_IwTCU_ILF1"].apply(convert_sig_unsig)/10
    df["tm2_vstupni_napeti"] = df["AM_B_C1_15_IwTCU_UD1Voltage"]/10
    df["tm2_tah_vyvijeny_motory"] = df["AM_B_C1_15_IwTCU_E_MC1"].apply(convert_sig_unsig)/100
    df["tm2_rychlost_5_napravy"] = df["VrTCU_B_A1_Velocity"]/3.6
    df["tm2_prikon_trakce"] = df["tm2_stejnoproud_trakce"]*df["tm2_vstupni_napeti"]/1000
    df["tm2_vykon_trakce"] = df["tm2_tah_vyvijeny_motory"]*df["tm2_rychlost_5_napravy"]
    df["tm2_vykon_zmareny"] = np.where((df["tm2_vykon_trakce"] > 0) & (df["tm2_vykon_trakce"] < df["tm2_prikon_trakce"]), df["tm2_prikon_trakce"]-df["tm2_vykon_trakce"], 0)

    # Cela tramvaj
    df["cela_vykon_spotrebovany"] = np.where(df["tm_prikon_trakce"]+df["tm_prikon_spotrebice"]+df["tm2_prikon_trakce"] > 0, df["tm_prikon_trakce"]+df["tm_prikon_spotrebice"]+df["tm2_prikon_trakce"], 0)
    df["cela_vykon_spotrebovany_bez_spotrebicu"] = np.where(df["tm_prikon_trakce"]+df["tm2_prikon_trakce"] > 0, df["tm_prikon_trakce"]+df["tm2_prikon_trakce"], 0)
    df["cela_vykon_rekuperovany"] = np.where(df["tm_prikon_trakce"]+df["tm_prikon_spotrebice"]+df["tm2_prikon_trakce"] < 0, -(df["tm_prikon_trakce"]+df["tm_prikon_spotrebice"]+df["tm2_prikon_trakce"]), 0)
    df["cela_vykon_rekuperovany_bez_spotrebicu"] = np.where(df["tm_prikon_trakce"]+df["tm2_prikon_trakce"] < 0, -(df["tm_prikon_trakce"]+df["tm2_prikon_trakce"]), 0)
    df["cela_vykon_zmareny"] = df["tm_vykon_zmareny"]+df["tm2_vykon_zmareny"]
    df["cela_vykon_trakce"] = df["tm_vykon_trakce"]+df["tm2_vykon_trakce"]
    df["cela_prikon_topeni_klima"] = df["tm_prikon_topeni_klima"]

    df["cela_energie_spotrebovana"] = (df["cela_vykon_spotrebovany"]*df["delta_T"].dt.total_seconds()/3600).cumsum().shift(1)
    df.loc[0, "cela_energie_spotrebovana"] = 0.0

    df["cela_energie_spotrebovana_bez_spotrebicu"] = (df["cela_vykon_spotrebovany_bez_spotrebicu"]*df["delta_T"].dt.total_seconds()/3600).cumsum().shift(1)
    df.loc[0, "cela_energie_spotrebovana_bez_spotrebicu"] = 0.0

    df["cela_energie_rekuperovana"] = (df["cela_vykon_rekuperovany"]*df["delta_T"].dt.total_seconds()/3600).cumsum().shift(1)
    df.loc[0, "cela_energie_rekuperovana"] = 0.0

    df["cela_energie_rekuperovana_bez_spotrebicu"] = (df["cela_vykon_rekuperovany_bez_spotrebicu"]*df["delta_T"].dt.total_seconds()/3600).cumsum().shift(1)
    df.loc[0, "cela_energie_rekuperovana_bez_spotrebicu"] = 0.0

    df["cela_energie_zmarena"] = (df["cela_vykon_zmareny"]*df["delta_T"].dt.total_seconds()/3600).cumsum().shift(1)
    df.loc[0, "cela_energie_zmarena"] = 0.0

    df["cela_energie_topeni_klima"] = (df["cela_prikon_topeni_klima"]*df["delta_T"].dt.total_seconds()/3600).cumsum().shift(1)
    df.loc[0, "cela_energie_topeni_klima"] = 0.0

    df["cela_vykon_vyrovnany_bez_spotrebicu"] = df["cela_vykon_spotrebovany_bez_spotrebicu"]-df["cela_vykon_rekuperovany_bez_spotrebicu"]
    df["cela_energie_vyrovnana_bez_spotrebicu"] = df["cela_energie_spotrebovana_bez_spotrebicu"]-df["cela_energie_rekuperovana_bez_spotrebicu"]

    df["km_total"] = df["QdwMesitTotalDistance"]*0.1

    return df
    return default_cols

def make_geojson(lat, long):
    geojson = {
        "type": "LineString",
    }
    geojson["coordinates"] = list(zip(long, lat))
    return json.dumps(geojson)