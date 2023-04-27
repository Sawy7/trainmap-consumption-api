from tconsumption import Consumption
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

class PlotData:
    def __init__(self, name, pretty_name):
        self.name = name
        self.pretty_name = pretty_name
        self.axes_obj = None
        self.lines = []

class ConsumptionGUI(Consumption):
    def update_plot_data(self, plot_data: PlotData):
        if plot_data.name not in self.series:
            return
        plot_data.lines[0].set_ydata(self.series[plot_data.name])
        if plot_data.name in self.comparison_series:
            if len(plot_data.lines) > 1:
                plot_data.lines[1].set_ydata(self.comparison_series[plot_data.name])
            else:
                line, = plot_data.axes_obj.plot(self.comparison_series["dist_values"], self.comparison_series[plot_data.name], label="(cmp)")
                plot_data.lines.append(line)

            self.get_dtw(plot_data.name)

    def render_plot_window(self):
        self.plots = [
            PlotData("acceleration_values", "Zrychlení (m/s2)"),
            PlotData("elevation_values", "Výška (m)"), # This has different length (len(points) != len(point2point))
            PlotData("velocity_values", "Rychlost (m/s)"),
            PlotData("force_values", "Síla (N)"),
            PlotData("energy_from_exerted_force", "Energie z vydané síly (J)")
        ]

        self.fig, axs = plt.subplots(len(self.plots), figsize=(20,4), dpi=100, sharex=True)
        for i,p in enumerate(self.plots):
            self.plots[i].axes_obj = axs[i]
            line, = axs[i].plot(c.series["dist_values"], c.series[p.name], label=p.pretty_name)
            self.plots[i].lines.append(line)
            axs[i].legend(loc="center left", bbox_to_anchor=(1, 0.5))
            if i == len(self.plots)-1:
                axs[i].set_xlabel("Vzdálenost (m)")
    
        self.params_to_sliders() # ENABLE SLIDERS HERE
        plt.show()

    def params_to_sliders(self):
        self.fig.subplots_adjust(left=0.40)
        pos_step = .05
        pos_left = .03
        for vp in self.variable_params.keys():
            if self.variable_params[vp] is not None:
                c.create_slider(vp, 0, 200, pos_left),
                pos_left += pos_step

    def create_slider(self, key_label, min, max, pos_left):
        axslider_window = self.fig.add_axes([pos_left, 0.18, 0.0225, 0.63])
        self.sliders.append(
            Slider(
                ax=axslider_window,
                label=key_label.replace(" ", "\n"),
                valmin=min,
                valmax=max,
                valinit=self.variable_params[key_label],
                orientation="vertical"
        ))
        self.update_slider(self.sliders[-1], key_label)

    def update_slider(self, slider: Slider, key_label: str):
        def update(val):
            self.variable_params[key_label] = val
            self.clean()
            self.run()
            for p in self.plots:
                self.update_plot_data(p)
            self.fig.canvas.draw_idle()

        slider.on_changed(update)

if __name__ == "__main__":
    c = ConsumptionGUI()
    # c.load_from_file("../testing-data/olo-opava.geojson")
    c.load_from_file("/tmp/49002.json")
    c.run()

    # # Testing comparison
    # acceleration_test_cmp = [x-0.1 for x in c.series["acceleration_values"]]
    # c.insert_comparsion("dist_values", c.series["dist_values"]) # This is very neccesary
    # c.insert_comparsion("acceleration_values", acceleration_test_cmp)

    c.render_plot_window()
