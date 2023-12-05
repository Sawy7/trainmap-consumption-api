# Source: https://github.com/yuma-m/matplotlib-draggable-plot/blob/master/draggable_plot.py
# Source: https://stackoverflow.com/a/19829987

import math
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.backend_bases import MouseEvent


class DraggablePlotExample(object):
    u""" An example of plot with draggable markers """

    def __init__(self, bg_data = []):
        self._figure, self._axes, self._line = None, None, None
        self._dragging_point = None
        self._points = {}
        self._bg_data = bg_data
        self._ctrl_down = 0
        self._mouse_press = None
        self._cur_xlim = None
        self._cur_ylim = None
        self.x0 = None
        self.y0 = None
        self.x1 = None
        self.y1 = None
        self.xpress = None
        self.ypress = None
        self.redraw_delayer_def = 10
        self.redraw_delayer = self.redraw_delayer_def

        self._init_plot()

    def _init_plot(self):
        self._figure = plt.figure("Example plot")
        self._figure.set_figwidth(20)
        axes = plt.subplot(1, 1, 1)
        axes.set_xlim(0, 100)
        axes.set_ylim(0, 100)
        axes.grid(which="both")
        self._axes = axes
        if len(self._bg_data) > 0:
            self._axes.plot(range(len(self._bg_data)), self._bg_data, "r")
            axes.set_xlim(0, len(self._bg_data))
        self.background = self._axes.figure.canvas.copy_from_bbox(self._axes.bbox)

        self._figure.canvas.mpl_connect('button_press_event', self._on_click)
        self._figure.canvas.mpl_connect('scroll_event', self._on_scroll)
        self._figure.canvas.mpl_connect('key_press_event', self._on_key_press)
        self._figure.canvas.mpl_connect('key_release_event', self._on_key_release)
        self._figure.canvas.mpl_connect('button_release_event', self._on_release)
        self._figure.canvas.mpl_connect('motion_notify_event', self._on_motion)
        plt.show()

    def _update_plot(self):
        if not self._points:
            self._line.set_data([], [])
        else:
            x, y = zip(*sorted(self._points.items()))
            # Add new plot
            if not self._line:
                self._line, = self._axes.plot(x, y, "b", marker="o", markersize=5)
            # Update current plot
            else:
                self._line.set_data(x, y)
        self._figure.canvas.draw()

    def _add_point(self, x, y=None):
        if isinstance(x, MouseEvent):
            x, y = x.xdata, x.ydata
        self._points[x] = y
        # print(self._points)
        return x, y

    def _remove_point(self, x, _):
        if x in self._points:
            self._points.pop(x)

    def _find_neighbor_point(self, event):
        u""" Find point around mouse position

        :rtype: ((int, int)|None)
        :return: (x, y) if there are any point around mouse else None
        """
        distance_threshold = 3.0
        nearest_point = None
        min_distance = math.sqrt(2 * (100 ** 2))
        for x, y in self._points.items():
            distance = math.hypot(event.xdata - x, event.ydata - y)
            if distance < min_distance:
                min_distance = distance
                nearest_point = (x, y)
        if min_distance < distance_threshold:
            return nearest_point
        return None

    def _on_click(self, event):
        u""" callback method for mouse click event

        :type event: MouseEvent
        """

        # Don't do modify points when ctrl not down
        if self._ctrl_down == 0:
            if event.button == 1 and event.inaxes in [self._axes]:
                self._cur_xlim = self._axes.get_xlim()
                self._cur_ylim = self._axes.get_ylim()
                self._mouse_press = self.x0, self.y0, event.xdata, event.ydata
                self.x0, self.y0, self.xpress, self.ypress = self._mouse_press
            return

        # left click
        if event.button == 1 and event.inaxes in [self._axes]:
            point = self._find_neighbor_point(event)
            if point:
                self._dragging_point = point
            else:
                self._add_point(event)
            self._update_plot()
        # right click
        elif event.button == 3 and event.inaxes in [self._axes]:
            point = self._find_neighbor_point(event)
            if point:
                self._remove_point(*point)
                self._update_plot()

    def _on_scroll(self, event):
        u""" callback method for mouse scroll event

        :type event: MouseEvent
        """
        base_scale = 2.
        cur_xlim = self._axes.get_xlim()
        cur_ylim = self._axes.get_ylim()

        xdata = event.xdata # get event x location
        ydata = event.ydata # get event y location

        if event.button == "up":
            # deal with zoom in
            scale_factor = 1 / base_scale
        elif event.button == "down":
            # deal with zoom out
            scale_factor = base_scale
        else:
            # deal with something that should never happen
            scale_factor = 1
            print(event.button)

        new_width = (cur_xlim[1] - cur_xlim[0]) * scale_factor
        new_height = (cur_ylim[1] - cur_ylim[0]) * scale_factor

        relx = (cur_xlim[1] - xdata)/(cur_xlim[1] - cur_xlim[0])
        rely = (cur_ylim[1] - ydata)/(cur_ylim[1] - cur_ylim[0])

        self._axes.set_xlim([xdata - new_width * (1-relx), xdata + new_width * (relx)])
        self._axes.set_ylim([ydata - new_height * (1-rely), ydata + new_height * (rely)])
        self._axes.figure.canvas.draw()

    def _on_key_press(self, event):
        u""" callback method for key press event

        :type event: KeyEvent
        """
        if event.key == "control":
            self._ctrl_down += 1

    def _on_key_release(self, event):
        u""" callback method for key press event

        :type event: KeyEvent
        """
        if event.key == "control":
            self._ctrl_down -= 1

    def _on_release(self, event):
        u""" callback method for mouse release event

        :type event: MouseEvent
        """
        if event.button == 1 and event.inaxes in [self._axes] and self._dragging_point:
            self._dragging_point = None
            self._update_plot()
        else:
            self._mouse_press = None
            self._axes.figure.canvas.draw()

    def _on_motion(self, event):
        u""" callback method for mouse motion event

        :type event: MouseEvent
        """
        # Don't do modify points when ctrl not down
        if self._mouse_press is not None:
            self._cur_xlim -= event.xdata - self.xpress
            self._cur_ylim -= event.ydata - self.ypress
            self._axes.set_xlim(self._cur_xlim)
            self._axes.set_ylim(self._cur_ylim)

            if self.redraw_delayer == 0:
                # Instead of redrawing the entire canvas, just blit the updated region
                self._axes.figure.canvas.restore_region(self.background)
                self._axes.draw_artist(self._axes)

                # Blit the updated region to the screen
                self._axes.figure.canvas.blit(self._axes.bbox)

                # Store the background for the next event
                self.background = self._axes.figure.canvas.copy_from_bbox(self._axes.bbox)
                self.redraw_delayer = self.redraw_delayer_def
            else:
                self.redraw_delayer -= 1

            return

        if not self._dragging_point:
            return
        if event.xdata is None or event.ydata is None:
            return
        self._remove_point(*self._dragging_point)
        self._dragging_point = self._add_point(event)
        self._update_plot()


if __name__ == "__main__":
    DATA_PATH="../../enet-sz-data/real_rides/"
    df = pd.read_csv(DATA_PATH + "Suchdol-Svinov/RegioPanter/Location.csv", delimiter=",")
    plot = DraggablePlotExample(df["speed"])