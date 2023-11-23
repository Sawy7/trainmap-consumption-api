import math
import numpy as np
from scipy.signal import savgol_filter
from fastdtw import fastdtw
import json
from pyproj import Proj


# Constants: ##################################################################
G_TO_MS2 = 9.80665
TRAIN_ACC_G = 0.1
TRAIN_ACC_MS2 = TRAIN_ACC_G * G_TO_MS2
TRAIN_DEC_G = 0.12
TRAIN_DEC_MS2 = TRAIN_DEC_G * G_TO_MS2
###############################################################################


# Helper functions: ###########################################################
# Source: https://www.movable-type.co.uk/scripts/latlong.html
def calc_distance_two_points(point_a, point_b):
    R = 6371e3  # meters
    φ1 = point_a[1] * math.pi/180  # φ, λ in radians
    φ2 = point_b[1] * math.pi/180
    Δφ = (point_b[1]-point_a[1]) * math.pi/180
    Δλ = (point_b[0]-point_a[0]) * math.pi/180

    a = math.sin(Δφ/2) * math.sin(Δφ/2)  \
        + math.cos(φ1) * math.cos(φ2) * math.sin(Δλ/2) * math.sin(Δλ/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))

    return R * c  # meters


def calc_radius_curve(point_a, point_b, point_c):
    x1 = float(point_a[0])
    y1 = float(point_a[1])
    x2 = float(point_b[0])
    y2 = float(point_b[1])
    x3 = float(point_c[0])
    y3 = float(point_c[1])

    xCoefficientArray = [x2 - x1, y2 - y1]
    yCoefficientArray = [x3 - x1, y3 - y1]

    coefficientArray = np.array([xCoefficientArray, yCoefficientArray])
    constantArray = np.array([
        (pow(x2, 2) + pow(y2, 2) - pow(x1, 2) - pow(y1, 2))/2,
        (pow(x3, 2) + pow(y3, 2) - pow(x1, 2) - pow(y1, 2))/2])

    try:
        center = np.linalg.solve(coefficientArray, constantArray)
        return calc_distance_two_points(point_a, center)
    except:
        return None


def calc_curve_resistance(radius, numerator, denominator):
    return numerator/(radius-denominator)  # [N/kN]


def calc_curve_resistance_force(point_a, point_b, point_c,
                                mass, numerator=650, denominator=55):
    radius = calc_radius_curve(point_a, point_b, point_c)
    if radius is None:
        return 0
    resistance = calc_curve_resistance(radius, numerator, denominator)
    mass_in_tons = mass/1000
    return resistance * mass_in_tons * G_TO_MS2


def get_elevation_slope_cos(point_a, point_b, dist):
    elevation_delta = point_b[2] - point_a[2]
    slope_distance = math.sqrt(elevation_delta**2 + dist**2)
    if slope_distance == 0:
        return 1, slope_distance
    else:
        return dist/slope_distance, slope_distance


def get_elevation_slope_sin(point_a, point_b, dist):
    elevation_delta = abs(point_b[2] - point_a[2])
    slope_distance = math.sqrt(elevation_delta**2 + dist**2)
    if slope_distance == 0:
        return 0, slope_distance
    else:
        return elevation_delta/slope_distance, slope_distance


def calc_normal_force(mass, angle_cos):
    return mass*G_TO_MS2*angle_cos


def calc_adhesion_ck(velocity):
    return 7500/(velocity+44)+161  # μ_max


def calc_tangential_force(mass, angle_cos, velocity):
    velocity_in_kph = velocity*3.6
    normal_force_in_kn = calc_normal_force(mass, angle_cos)/1000
    return normal_force_in_kn*calc_adhesion_ck(velocity_in_kph)


def calc_parallel_g_force(mass, angle_sin):
    return mass*G_TO_MS2*angle_sin


def calc_running_resistance(velocity, a, b, c):
    if velocity == 0:
        a = 0
    velocity_in_kph = velocity*3.6
    return a+b*velocity_in_kph+c*velocity_in_kph**2


def calc_running_resistance_force(velocity, mass, a=1.35, b=0.0008, c=0.00033):
    resistance = calc_running_resistance(velocity, a, b, c)
    mass_in_tons = mass/1000
    return resistance * mass_in_tons * G_TO_MS2


def calc_acceleration(force, mass):
    return force/mass


def calc_velocity(acceleration, distance, init_velocity=0):
    try:
        # Changing parameters can break this
        return math.sqrt(init_velocity**2 + 2*acceleration*distance)
    except:
        return 0


def calc_reverse_acceleration(velocity, distance, init_velocity):
    if distance == 0:
        return 0
    return (velocity**2-init_velocity**2)/(2*distance)


def calc_force(mass, acceleration):
    return mass*acceleration


def parse_points_from_geojson(geojson_raw):
    geojson = json.loads(geojson_raw)
    return geojson["coordinates"]


def parse_stations_from_geojson(geojson_raw):
    geojson = json.loads(geojson_raw)
    return geojson["station_orders"]


def velocity_ways_to_max_velocity(velocity_ways):
    max_velocities = []
    for v in velocity_ways:
        max_velocities += [v["velocity"]]*(v["end"]-v["start"]+1)
    return max_velocities


def parse_velocity_ways_from_geojson(geojson_raw):
    geojson = json.loads(geojson_raw)
    velocity_ways = geojson["velocity_ways"]
    return velocity_ways


def get_power_raw(force_values, velocity_values):
    return [force_values[i]*velocity_values[i] for i in range(len(force_values))]


def get_energy_from_force(force_values, dist_values):
    to_return = []
    running_sum = 0
    for i in range(len(force_values)):
        if i > 0:
            cur_dist = dist_values[i] - dist_values[i-1]
            # print(dist_values[i], dist_values[i-1], cur_dist)
        else:
            cur_dist = dist_values[i]
            # print(dist_values[i], cur_dist)
        running_sum += force_values[i]*cur_dist
        to_return.append(running_sum)
    return to_return


def calc_power_from_acceleration(mass, acceleration, velocity):
    # P = mav
    return mass*acceleration*velocity


def calc_acceleration_from_power(power, mass, velocity):
    # a = P/(m*v)
    return power/(mass*velocity)
###############################################################################


###############################################################################
class ConsumptionPart:
    def __init__(
            self, series, mass_locomotive, mass_wagon, points,
            max_velocities, filter_window_elev, filter_window_curve,
            curve_res_p: tuple, running_res_p: tuple,
            power_limit, recuperation_coefficient,
            comfortable_acceleration, comp_poly):
        # Input parameters
        self.mass_locomotive = mass_locomotive
        self.mass_wagon = mass_wagon
        self.points = points
        self.max_velocities = max_velocities
        self.filter_window_elev = filter_window_elev
        self.filter_window_curve = filter_window_curve
        self.curve_res_p = curve_res_p
        self.running_res_p = running_res_p
        self.power_limit = power_limit
        self.recuperation_coefficient = recuperation_coefficient
        self.comfortable_acceleration = comfortable_acceleration

        # Compensation polynomial (velocity)
        self.comp_poly = comp_poly

        # Working lists
        self.series = series
        self.curve_res_force_all_l = [0]
        self.curve_res_force_all_w = [0]

    def get_curve_resistance_force(self):
        # # Curve forces
        # for i in range(len(self.points)):
        #     if i % 3 == 0:
        #         if i+2 < len(self.points):
        #             curve_res_force_l = calc_curve_resistance_force(
        #                     self.points[i], self.points[i+1], self.points[i+2],
        #                     self.mass_locomotive, *self.curve_res_p)/3
        #             curve_res_force_w = calc_curve_resistance_force(
        #                     self.points[i], self.points[i+1], self.points[i+2],
        #                     self.mass_wagon, *self.curve_res_p)/3
        #         else:
        #             curve_res_force_l = 0
        #             curve_res_force_w = 0
        #     self.curve_res_force_all_l.append(curve_res_force_l)
        #     self.curve_res_force_all_w.append(curve_res_force_w)

        # # Don't filter if windows == 0
        # if self.filter_window_curve <= 0:
        #     return

        # self.curve_res_force_all_l = savgol_filter(
        #         self.curve_res_force_all_l, self.filter_window_curve,
        #         0, mode="nearest")
        # self.curve_res_force_all_w = savgol_filter(
        #         self.curve_res_force_all_w, self.filter_window_curve,
        #         0, mode="nearest")

        # return

        x = [p[0] for p in self.points]
        y = [p[1] for p in self.points]

        czech_proj = Proj("epsg:5514")  # Hardcoded for CZ/SK
        czech_x, czech_y = czech_proj(x, y)

        dx = np.diff(czech_x)
        dy = np.diff(czech_y)

        distances = np.sqrt(dx**2 + dy**2)
        distances[distances == 0] = np.inf

        normalized_dx = dx / distances
        normalized_dy = dy / distances

        d2x = np.gradient(normalized_dx)
        d2y = np.gradient(normalized_dy)

        denominator = (normalized_dx**2 + normalized_dy**2)**1.5
        denominator[denominator == 0] = np.inf
        radius_of_curvature = np.abs((normalized_dx * d2y - normalized_dy * d2x) / denominator)

        for radius in radius_of_curvature:
            res = calc_curve_resistance(radius, *self.curve_res_p)
            # NOTE: resistance is negative, so we flip the sign (and subtract it later in the calcs)
            self.curve_res_force_all_l.append(-res * (self.mass_locomotive/1000) * G_TO_MS2)
            self.curve_res_force_all_w.append(-res * (self.mass_wagon/1000) * G_TO_MS2)

        # Don't filter if windows == 0
        if self.filter_window_curve <= 0:
            return

        self.curve_res_force_all_l = savgol_filter(
                self.curve_res_force_all_l, self.filter_window_curve,
                0, mode="nearest")
        self.curve_res_force_all_w = savgol_filter(
                self.curve_res_force_all_w, self.filter_window_curve,
                0, mode="nearest")


    def cap_acceleration(self, mass, acceleration, velocity):
        if self.power_limit is None:
            return acceleration
        else:
            uncapped_power = calc_power_from_acceleration(
                    mass, acceleration, velocity)
            if uncapped_power > self.power_limit:
                return calc_acceleration_from_power(
                        self.power_limit, mass, velocity)
            else:
                return acceleration

    def get_compensation(self, dist):
        if self.comp_poly is None:
            return 0
        return max(0, self.comp_poly(dist))

    def slow_down_to_max_limit_six(self, max_velocity, slow_point_index):
        end_force = []
        end_exerted_force = []
        end_velocity = [max_velocity]
        deceleration_values = []
        curve_res_force_l = 0
        curve_res_force_w = 0

        for i in range(slow_point_index, 0, -1):
            immediate_distance = calc_distance_two_points(self.points[i], self.points[i-1])
            angle_cos, slope_distance = get_elevation_slope_cos(self.points[i-1], self.points[i], immediate_distance)
            angle_sin = get_elevation_slope_sin(self.points[i-1], self.points[i], immediate_distance)[0]

            # Forces on locomotive
            tangential_force_l = calc_tangential_force(self.mass_locomotive, angle_cos, end_velocity[-1])
            parallel_g_force_l = calc_parallel_g_force(self.mass_locomotive, angle_sin)
            running_res_force_l = calc_running_resistance_force(end_velocity[-1], self.mass_locomotive, *self.running_res_p)

            # Forces on wagon
            parallel_g_force_w = calc_parallel_g_force(self.mass_wagon, angle_sin)
            running_res_force_w = calc_running_resistance_force(end_velocity[-1], self.mass_wagon, *self.running_res_p)

            # Curve forces
            curve_res_force_l = self.curve_res_force_all_l[i]
            curve_res_force_w = self.curve_res_force_all_w[i]

            # Is it incline/decline?
            if self.points[i-1][2] - self.points[i][2] > 0:  # Incline
                final_force = tangential_force_l - parallel_g_force_l - parallel_g_force_w - running_res_force_l - running_res_force_w
            else:  # Decline
                final_force = tangential_force_l + parallel_g_force_l + parallel_g_force_w - running_res_force_l - running_res_force_w
            final_force += - curve_res_force_l - curve_res_force_w
            # Exerted force doesn't 100 % translate to gains - multiply by coefficient
            exerted_force = tangential_force_l * self.recuperation_coefficient
            acceleration = calc_acceleration(final_force, self.mass_locomotive+self.mass_wagon)
            # Make acceleration when slowing down less
            # Acceleration capping must be implemented differently (slowing down is not motor-capped)
            prev_acc = acceleration
            acceleration = max(min(acceleration, self.comfortable_acceleration), -self.comfortable_acceleration)
            if prev_acc > acceleration:
                external_force = final_force - tangential_force_l
                final_force = calc_force(self.mass_locomotive+self.mass_wagon, acceleration)
                exerted_force = max(final_force - external_force, 0)
                if exerted_force > 0:
                    exerted_force *= self.recuperation_coefficient
            new_velocity = calc_velocity(acceleration, slope_distance, end_velocity[-1])
            if len(end_force) > 0:
                end_velocity.append(new_velocity)
            end_force.append(-final_force) # NOTE: Braking force has opposite direction
            end_exerted_force.append(-exerted_force)
            deceleration_values.append(-acceleration) # NOTE: Deceleration has opposite direction

            if new_velocity >= self.series["velocity_values"][i]:
                if len(end_velocity) == 1:
                    new_velocity = end_velocity[0]
                else:
                    new_velocity = self.series["velocity_values"][i]
                acceleration = calc_reverse_acceleration(end_velocity[-1], slope_distance, new_velocity) # Velocities are in 'reverse' order here (compared to next function)
                reverse_force = calc_force(self.mass_locomotive+self.mass_wagon, -acceleration)
                # Coeficient again
                exerted_force -= (final_force-reverse_force) * self.recuperation_coefficient
                final_force = reverse_force
                end_velocity[-1] = new_velocity
                end_force[-1] = final_force  # NOTE: Force here already has opposite direction (from reverse_force)
                end_exerted_force[-1] = -exerted_force
                deceleration_values[-1] = -acceleration
                break
        return end_force, end_exerted_force, end_velocity, deceleration_values

    def get_ramp_up_six(self):
        self.series["force_values"][0] = 0
        self.series["exerted_force_values"][0] = 0
        self.series["velocity_values"][0] = 0
        self.series["dist_values"][0] = 0
        self.series["acceleration_values"][0] = 0
        velocity_reached = False
        start_compensated = False

        compensated_part = []
        compensated_dist_offset = 0

        for i in range(len(self.points)-1):
            series_i = i+1
            immediate_distance = calc_distance_two_points(self.points[i], self.points[i+1])
            angle_cos, slope_distance = get_elevation_slope_cos(self.points[i], self.points[i+1], immediate_distance)
            angle_sin = get_elevation_slope_sin(self.points[i], self.points[i+1], immediate_distance)[0]
            last_velocity = self.series["velocity_values"][series_i-1]
            last_dist = self.series["dist_values"][series_i-1]

            # Forces on locomotive
            tangential_force_l = calc_tangential_force(self.mass_locomotive, angle_cos, last_velocity)
            parallel_g_force_l = calc_parallel_g_force(self.mass_locomotive, angle_sin)
            running_res_force_l = calc_running_resistance_force(last_velocity, self.mass_locomotive, *self.running_res_p)

            # Forces on wagon
            parallel_g_force_w = calc_parallel_g_force(self.mass_wagon, angle_sin)
            running_res_force_w = calc_running_resistance_force(last_velocity, self.mass_wagon, *self.running_res_p)

            # Curve forces
            curve_res_force_l = self.curve_res_force_all_l[i]
            curve_res_force_w = self.curve_res_force_all_w[i]

            if velocity_reached and not start_compensated:
                compensated_part = [max(0, x-compensated_part[0]) for x in compensated_part]
                compensated_part.reverse()
                for j, cp in enumerate(compensated_part):
                    self.series["velocity_values"][series_i-j-1] = self.series["velocity_values"][series_i-j-1] - cp
                start_compensated = True
                compensated_part = []

            if velocity_reached and (self.max_velocities[i] > self.max_velocities[i-1] or last_velocity < self.max_velocities[i]):
                velocity_reached = False

            if not velocity_reached:
                ##########################################################################################
                if len(compensated_part) == 0:
                    end_force_slow = self.slow_down_to_max_limit_six(0, i)[0]
                    slowdown_point_count = len(end_force_slow)
                    if last_dist == 0:
                        compensated_dist_offset = 0
                    else:
                        compensated_dist_offset = last_dist - self.series["dist_values"][series_i-slowdown_point_count]
                        compensated_dist_offset = last_dist - compensated_dist_offset
                ##########################################################################################
                # Is it incline/decline?
                if self.points[i+1][2] - self.points[i][2] > 0:  # Incline
                    final_force = tangential_force_l - parallel_g_force_l - parallel_g_force_w - running_res_force_l - running_res_force_w
                else:  # Decline
                    final_force = tangential_force_l + parallel_g_force_l + parallel_g_force_w - running_res_force_l - running_res_force_w
                final_force += - curve_res_force_l - curve_res_force_w
                acceleration = calc_acceleration(final_force, self.mass_locomotive+self.mass_wagon)
                prev_acc = acceleration
                acceleration = self.cap_acceleration(self.mass_locomotive+self.mass_wagon, acceleration, last_velocity)
                # Try capping it to sane values - this should combat imprecise data
                acceleration = max(min(acceleration, self.comfortable_acceleration), -self.comfortable_acceleration)

                # NOTE: No acceleration capping - REMOVED (power capping)
                new_velocity = calc_velocity(acceleration, slope_distance, last_velocity)
                ###############################################################################
                # Decrease velocity by polynom (found by comparing /w real data)
                compensated_dist = last_dist-compensated_dist_offset+slope_distance
                # print(last_dist, compensated_dist_offset, slope_distance, compensated_dist)
                compensated_velocity = self.get_compensation(compensated_dist)
                compensated_part.append(compensated_velocity)
                start_compensated = False
                ###############################################################################
                exerted_force = tangential_force_l
                if prev_acc > acceleration:
                    external_force = final_force - tangential_force_l
                    final_force = calc_force(self.mass_locomotive+self.mass_wagon, acceleration)
                    exerted_force = final_force - external_force
                    if exerted_force < 0:
                        exerted_force *= self.recuperation_coefficient

                # Clamp it down, but only if the limit is the same (otherwise we could be clamping by A LOT)
                if new_velocity > self.max_velocities[i] and self.max_velocities[i] == self.max_velocities[i-1]:
                    new_velocity = self.max_velocities[i]
                    # TODO: compensate velocity here also
                    acceleration = calc_reverse_acceleration(new_velocity, slope_distance, last_velocity)
                    reverse_force = calc_force(self.mass_locomotive+self.mass_wagon, acceleration)
                    exerted_force -= final_force-reverse_force
                    final_force = reverse_force
                    velocity_reached = True

                if velocity_reached and new_velocity-compensated_velocity < self.max_velocities[i]:
                    velocity_reached = False

            else:
                new_velocity = last_velocity
                # NOTE: when surface is not flat, we need to exert force to keep our speed
                # Is it incline/decline?
                if self.points[i+1][2] - self.points[i][2] > 0:  # Incline
                    final_force = parallel_g_force_l + parallel_g_force_w + running_res_force_l + running_res_force_w + curve_res_force_l + curve_res_force_w
                else:  # Decline
                    final_force = -parallel_g_force_l - parallel_g_force_w + running_res_force_l + running_res_force_w + curve_res_force_l + curve_res_force_w
                acceleration = 0  # NOTE: If no change in speed, acceleration is ZERO
                exerted_force = final_force
                if exerted_force < 0:
                    exerted_force *= self.recuperation_coefficient

            self.series["force_values"][series_i] = final_force
            self.series["exerted_force_values"][series_i] = exerted_force
            self.series["velocity_values"][series_i] = new_velocity
            self.series["acceleration_values"][series_i] = acceleration
            self.series["dist_values"][series_i] = last_dist+slope_distance

            if i > 0 and new_velocity > self.max_velocities[i]:
                end_force_slow, end_exerted_force_slow, end_velocity_slow, deceleration_values_slow = self.slow_down_to_max_limit_six(self.max_velocities[i], i)
                for j in range(len(end_velocity_slow)):
                    self.series["force_values"][series_i-j] = end_force_slow[j]
                    self.series["exerted_force_values"][series_i-j] = end_exerted_force_slow[j]
                    self.series["velocity_values"][series_i-j] = end_velocity_slow[j]
                    self.series["acceleration_values"][series_i-j] = deceleration_values_slow[j]

    def run(self):
        # Manually filter elevation values
        if self.filter_window_elev > 0:
            elevation = [x[2] for x in self.points]
            filtered = savgol_filter(elevation, self.filter_window_elev, 0, mode="nearest")
            self.points = [[x[0], x[1], filtered[i]] for i,x in enumerate(self.points)]

        # Calculate all curve resistance force ahead of time
        self.get_curve_resistance_force()

        # The master calculation
        self.get_ramp_up_six()

        # Braking (going backwards)
        end_force, end_exerted_force, end_velocity, deceleration_values = self.slow_down_to_max_limit_six(0, len(self.points)-1)
        for i in range(len(end_velocity)):
            self.series["force_values"][-i-1] = end_force[i]
            self.series["exerted_force_values"][-i-1] = end_exerted_force[i]
            self.series["velocity_values"][-i-1] = end_velocity[i]
            self.series["acceleration_values"][-i-1] = deceleration_values[i]


class Consumption:
    def __init__(self):
        # Params
        self.params = {
            "mass_locomotive": 30000,   # kg
            "mass_wagon": 0,            # kg
            "power_limit": 4500*1000    # 4500 kW
        }
        self.variable_params = {
            "Elevation smoothing": 100,
            "Curve smoothing": 10,
            "Curve A": 650,
            "Curve B": 55,
            "Running a": 1.35,
            "Running b": 0.0008,
            "Running c": 0.00033,
            "Recuperation coefficient": 1,
            "Comfortable acceleration": 0.89,
            # "Compensation polynomial": np.poly1d(
            #     [-1.9494156333280808e-14, 4.9173340030234155e-11,
            #      -4.531261473352924e-08, 1.5457227860221452e-05,
            #      0.00332072288476179, 0.7453216365390641])
            "Compensation polynomial": np.poly1d(
                [-1.3169756479815293e-14, 4.271539912516026e-11,
                 -3.874542069136512e-08, 2.677139343735735e-08,
                 0.00962245960144532, 0.619862570600036])
        }

        # Internal data
        self.points = None
        self.stations = None
        self.max_velocities_in_mps = None
        self.sliders = []
        self.comparison_series = {}

    def init_series(self, first_station, last_station):
        length = last_station + 1 - first_station
        self.series = {
            "force_values": np.empty(length),
            "exerted_force_values": np.empty(length),
            "dist_values": np.empty(length),
            "velocity_values": np.empty(length),
            "acceleration_values": np.empty(length),
            "energy_from_exerted_force": [],
        }

    def insert_comparsion(self, name, data):
        self.comparison_series[name] = data

    def load_from_file(self, file_path):
        with open(file_path) as f:
            lines = f.readlines()
            geojson_raw = "".join(lines)
            self.load(
                parse_points_from_geojson(geojson_raw),
                parse_stations_from_geojson(geojson_raw),
                parse_velocity_ways_from_geojson(geojson_raw)
            )

    def load(self, points, stations, velocity_ways):
        # First initialize working series
        self.init_series(stations[0], stations[-1])
        # Now load data
        self.points = points
        self.stations = [x-stations[0] for x in stations]
        self.stations.sort()
        # Max velocities get loaded and converted to m/s
        max_velocities = velocity_ways_to_max_velocity(velocity_ways)
        self.max_velocities_in_mps = [x/3.6 for x in max_velocities]
        # Elevation is cropped from first to last station (to make everything same length) # TODO: make into numpy array
        self.series["elevation_values"] = [e[2] for e in self.points][self.stations[0]:self.stations[-1]+1]

    def run(self):
        station_offset = 1
        for i in range(len(self.stations)-1):
            if i == len(self.stations)-2:
                station_offset = 0

            # Get offsets for ConsumptionPart
            start_offset = self.stations[i]
            end_offset = self.stations[i+1]-station_offset+1

            split_points = self.points[start_offset:end_offset]
            split_max_velocities_in_mps = self.max_velocities_in_mps[start_offset:end_offset]
            split_series = {}
            for s in self.series.keys():
                split_series[s] = self.series[s][start_offset:end_offset]

            consumption_part = ConsumptionPart(
                split_series,
                self.params["mass_locomotive"],
                self.params["mass_wagon"],
                split_points,
                split_max_velocities_in_mps,
                self.variable_params["Elevation smoothing"],
                self.variable_params["Curve smoothing"],
                (self.variable_params["Curve A"], self.variable_params["Curve B"]),
                (self.variable_params["Running a"], self.variable_params["Running b"], self.variable_params["Running c"]),
                self.params["power_limit"],
                self.variable_params["Recuperation coefficient"],
                self.variable_params["Comfortable acceleration"],
                self.variable_params["Compensation polynomial"]
            )
            consumption_part.run()

            if start_offset > 0:
                prev_dist = self.series["dist_values"][start_offset-1]
                self.series["dist_values"][start_offset:end_offset] += prev_dist

        self.series["energy_from_exerted_force"] = get_energy_from_force(self.series["exerted_force_values"], self.series["dist_values"])

    def get_dtw(self, name):
        distance, path = fastdtw(self.series[name], self.comparison_series[name])
        print(f"{name} distance: {distance}")

# Testing
if __name__ == "__main__":
    c = Consumption()
    c.load_from_file("/home/sawy/Dev/enet-sz-consumption-api/testing-data/olo-opava.geojson")
    c.run()