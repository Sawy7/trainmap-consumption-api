from tconsumption import Consumption

c = Consumption()
c.load_from_file("../testing-data/dpo-vresina-zatisi.geojson")
c.run()
