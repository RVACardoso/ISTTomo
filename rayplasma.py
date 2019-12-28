from raytomo import RayTomo
from plasmaemit import SinglePlasmaEmit, ModePlasmaEmit, NpySinglePlasmaEmit
from cmap import get_cmap

import numpy as np
import matplotlib.pyplot as plt
import sys
from math import cos, sin, radians, sqrt

from raysect.primitive import Cylinder, Sphere
from raysect.optical import translate, rotate
from raysect.optical.material import AbsorbingSurface

from daniel.simulateSignal import simulate_signal


def rebin(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)


class RayPlasma(RayTomo):
    csv_file_name = None

    def __init__(self, coord, radius_div=20, angle_div=40, pixel_nr=30, count=1):
        super().__init__()
        if coord != 'polar' and coord != 'cartesian':
            print("ERROR: coordinates must be polar or cartesians.")
            sys.exit()

        self.coord = coord
        self.radius_div = radius_div
        self.angle_div = angle_div
        self.pixel_nr = pixel_nr

        self.square_values = []

        self.radius_list = []
        self.angle_list = []
        self.xx_list = []
        self.yy_list = []

        self.plasma_count = count

    def place_source(self):

        # self.source = Cylinder(radius=self.vessel_in_rad-0.002, height=self.vessel_width,
        #                material=self.PlasmaEmit(shape=self.plasma_shape, count=self.plasma_count, center=self.plasma_center,
        #                                         width=self.plasma_width, power_gain=self.plasma_power_gain, flatlen=self.plasma_flat),
        #                parent=self.world, transform=translate(0, 0, 0) * rotate(0, 0, 0), name='PlasmaSource')

        self.source = Cylinder(radius=self.vessel_in_rad - 0.0015, height=self.vessel_width,
                               parent=self.world, transform=translate(0, 0, 0) * rotate(0, 0, 0), name='PlasmaSource')

        self.lid_back = Cylinder(radius=self.vessel_in_rad+0.002, height=0.001,
                       material=AbsorbingSurface(),
                       parent=self.world, transform=translate(0, 0, self.vessel_width+0.0005) * rotate(0, 0, 0))

        self.lid_front = Cylinder(radius=self.vessel_in_rad + 0.002, height=0.001,
                                 material=AbsorbingSurface(),
                                 parent=self.world, transform=translate(0, 0, -0.0015) * rotate(0, 0, 0))

        # print_scenegraph(self.source)

        # maybe add random noise??

    def get_tomogram_data(self, sinth_plasma, pixel_div=5):

        if self.coord == 'polar':
            total_radius = self.vessel_in_rad - 0.005
            radius_step = total_radius/self.radius_div
            angle_step = 360.0/self.angle_div
            for base_radius in np.arange(start=0.0, stop=total_radius, step=radius_step):
                for base_angle in np.arange(start=0.0, stop=360, step=angle_step):

                    square_points = 0
                    square_sum = 0
                    for delta_radius in np.arange(start=0, stop=radius_step, step=radius_step/pixel_div):
                        for delta_angle in np.arange(start=0, stop=angle_step, step=angle_step/pixel_div):

                            point_x = -(base_radius+delta_radius)*cos(radians(base_angle+delta_angle))
                            point_y = (base_radius+delta_radius)*sin(radians(base_angle+delta_angle))

                            for i in range(self.plasma_count):
                                # radius = sqrt((point_x - self.plasma_center[i][0]) ** 2 + (point_y - self.plasma_center[i][1]) ** 2)
                                #square_sum += self.plasma_power_gain * self.plasma_shape(radius, 0, self.plasma_width[i])  # cos((shift + 5) * radius)**4
                                square_sum += sinth_plasma.power_gain * sinth_plasma.compute_shape(x=point_x, y=point_y)
                            square_points += 1

                    self.square_values.append(square_sum/square_points)
                    self.radius_list.append(base_radius)
                    self.angle_list.append(base_angle)

        elif self.coord == 'cartesian':
            self.total_side = 0.10  # divertor at 8.5cm
            xx = -np.linspace(start=-self.total_side, stop=self.total_side, num=self.pixel_nr*pixel_div)
            yy = -np.linspace(start=-self.total_side, stop=self.total_side, num=self.pixel_nr*pixel_div)
            xv, yv = np.meshgrid(xx, yy)
            computed = np.zeros(xv.shape)
            for i in range(self.plasma_count):
                computed += (sinth_plasma.compute_shape(np.ravel(xv), np.ravel(yv)).reshape((self.pixel_nr*pixel_div, self.pixel_nr*pixel_div)))
            self.square_values = np.ravel(rebin(computed, shape=(self.pixel_nr, self.pixel_nr)))

        else:
            print("ERROR: coordinates must be either polar or cartesian!")

    def plot_tomogram(self):
        font_size=25
        plt.rcParams.update({'font.size': font_size})

        if self.coord == 'polar':
            delta_r = 1 / self.radius_div
            delta_ang = 6.28 / self.angle_div

            self.radius_list = self.radius_list / max(self.radius_list)
            self.angle_list = (self.angle_list / max(self.angle_list)) * (6.28 - delta_ang)
            # print(self.angle_list)
            self.square_values = np.array(self.square_values) / max(self.square_values)

            plt.figure()
            ax = plt.subplot(projection='polar')
            for i in range(len(self.radius_list)):
                plt.axvspan(ymin=self.radius_list[i], ymax=self.radius_list[i] + delta_r,
                            xmin=self.angle_list[i], xmax=self.angle_list[i] + delta_ang,
                            # facecolor=(self.square_values[i], self.square_values[i], self.square_values[i]), alpha=0.5)
                            facecolor='r', alpha=self.square_values[i])
            ax.grid(False)
            ax.set_yticklabels([round(0.2 * self.vessel_in_rad, 2), round(0.4 * self.vessel_in_rad, 2),
                                round(0.6 * self.vessel_in_rad, 2), round(0.8 * self.vessel_in_rad, 2),
                                round(1.0 * self.vessel_in_rad, 2)])

        elif self.coord == 'cartesian':
            import matplotlib
            # self.square_values = self.square_values / np.max(self.square_values)
            self.square_values = self.square_values.reshape((self.pixel_nr, self.pixel_nr))
            plt.figure()
            plt.imshow(self.square_values, cmap=matplotlib.cm.get_cmap("plasma"))
            plt.xlabel("x (m)", fontsize=font_size)
            plt.ylabel("z (m)", fontsize=font_size)
            # plt.xticks(np.linspace(0, 49, num=11), np.round(np.arange(start=-self.total_side, stop=self.total_side+0.02, step=0.02), 2))
            plt.xticks(np.linspace(0, 59, num=5),
                       np.round(np.arange(start=-self.total_side, stop=self.total_side + 0.02, step=0.05), 2))
            plt.yticks(np.linspace(0, 59, num=5), np.flipud(np.round(np.arange(start=-self.total_side, stop=self.total_side+0.02, step=0.05), 2)))
            cbar = plt.colorbar() #format='%.0e')
            cbar.ax.set_ylabel('Emissivity ($kW/m^{-3}$)', rotation=90, fontsize=font_size)

    def create_save(self):

        with open(RayPlasma.csv_file_name, "a+") as file:
            file.write("top0, top1, top2, top3, top4, top5, top6, top7, top8, top9, top10, top11, top12,"
                "top13, top14, top15, out0, out1, out2, out3, out4, out5, out6, out7, out8, out9, out10, out11, out12,"
                "out13, out14, out15")

            if self.coord == 'polar':
                pixel_count = self.radius_div*self.angle_div
            else:
                pixel_count = self.pixel_nr*self.pixel_nr

            for tomo_pixel in range(pixel_count):
                file.write(", px_" + str(tomo_pixel))
            file.write('\n')

    def load_data_to_save(self):
        with open(RayPlasma.csv_file_name, "a") as file:
            for pixel in range(16):
                file.write(str(self.top_data[pixel]) + ", ")
            for pixel in range(15):
                file.write(str(self.out_data[pixel]) + ", ")
            file.write(str(self.out_data[15]))

            for square in self.square_values:
                file.write(", " + str(square))
            file.write('\n')

    def plot_compare_los(self):

        # emissivity = self.square_values.reshape((self.pixel_nr, self.pixel_nr))
        # signals = simulate_signal(emissivity, "projections")
        los_top_data = np.array([2.27633276e-10, 1.80745613e-09, 8.90145992e-05, 2.19800454e-03,
                         6.23114228e-03, 1.42721551e-02, 2.19316503e-02, 5.22197529e-02,
                         1.26506393e-01, 2.05487226e-01, 2.13189247e-01, 1.70844403e-01,
                         1.02238920e-01, 5.10939790e-02, 2.35229398e-02, 1.01751719e-02])  # np.zeros(16)  # signals[:16]
        los_out_data = np.array([6.00351266e-03, 8.98786383e-03, 2.04519392e-02, 6.78003553e-02,
                         1.64823107e-01, 2.53159468e-01, 2.54643611e-01, 1.41516615e-01,
                         4.91891950e-02, 2.15590879e-02, 8.78117741e-03, 3.00394924e-03,
                         8.01107897e-05, 5.77734040e-09, 1.43561229e-09, 5.72060751e-13])  #np.zeros(16)  # signals[16:]

        los_top_data = los_top_data / sum(los_top_data)
        los_out_data = los_out_data / sum(los_out_data)
        self.top_data = np.array(self.top_data) / sum(self.top_data)
        self.out_data = np.array(self.out_data) / sum(self.out_data)

        cam_nr = np.arange(16)

        width = 0.4
        leg_size = 25
        # font_size = 30

        # plt.rc('xtick', labelsize=font_size)
        # plt.rc('ytick', labelsize=font_size)

        plt.figure()
        plt.subplot(211)
        plt.bar(cam_nr - width / 2, los_top_data, width, color='r', bottom=0, label='Without reflections')
        plt.bar(cam_nr + width / 2, self.top_data, width, color='b', bottom=0, label='With reflections')
        plt.xticks(cam_nr, cam_nr)
        # plt.text(3.0, .145, 'TOP camera', horizontalalignment='center')
        plt.title("Top camera")
        plt.xlabel("Detector Number")
        plt.ylabel("Fraction of total\ncollected power")
        plt.legend(prop={'size': leg_size-8}, bbox_to_anchor=(0.5, 0.95))

        plt.subplot(212)
        plt.bar(cam_nr - width / 2, los_out_data, width, color='r', bottom=0, label='Without reflections')
        plt.bar(cam_nr + width / 2, self.out_data, width, color='b', bottom=0, label='With reflections')
        # plt.text(12.5, .17, 'OUTER camera', horizontalalignment='center')
        plt.title("Outer camera")
        plt.xticks(cam_nr, cam_nr)
        plt.xlabel("Detector Number")
        plt.ylabel("Fraction of total\ncollected power")
        # plt.legend(prop={'size': leg_size-8}, bbox_to_anchor=(0.3, 0.95))

        plt.subplots_adjust(left=0.5, hspace=0.4)