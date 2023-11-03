import tconsumption, tconsumption_old
from time import sleep

if __name__ == "__main__":
    cn = tconsumption.Consumption()
    cn.load_from_file("../testing-data/olo-opava.geojson")
    cn.run()

    co = tconsumption_old.Consumption()
    co.load_from_file("../testing-data/olo-opava.geojson")
    co.run()

    # exit(0)
    all_good = True

    all_keys = list(cn.series.keys())
    for tested_key in all_keys:
        print("old:", type(co.series[tested_key]), "new:", type(cn.series[tested_key]))
        for i in range(len(cn.series[tested_key])):
            not_good = co.series[tested_key][i] != cn.series[tested_key][i]
            all_good &= not not_good
            if not_good:
                print("Error: Not matching")
                print(f"{i}:", co.series[tested_key][i], "->", cn.series[tested_key][i])

    if all_good:
        print("Everything is good")
    else:
        print("Encountered some differences")