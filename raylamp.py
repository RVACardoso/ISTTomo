from raytomo import RayTomo
from reallamp import RealLamp

import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin, radians

from raysect.optical import World, translate, rotate
from raysect.optical.material import UniformVolumeEmitter, AbsorbingSurface
from raysect.optical.spectralfunction import ConstantSF
from raysect.primitive import Cylinder


class RayLamp(RayTomo):
    csv_file_name = None

    def __init__(self):
        super().__init__()
        self.lamp_radius = None
        self.lamp_angle = None

    def place_source(self, radius, angle, legs):
        self.lamp_radius = radius
        self.lamp_angle = angle
        pos_x = -self.lamp_radius * cos(radians(self.lamp_angle))
        pos_y = self.lamp_radius * sin(radians(self.lamp_angle))
        self.source = Cylinder(radius=0.00205, height=0.15, parent=self.world,
                                    transform=translate(pos_x, pos_y, -0.0454), # scale=(1/0.007)*(699/(26*5))*1e8* (1/(1.175*1e6)))
                                    material=UniformVolumeEmitter(ConstantSF(1.0), scale=1.0), #65.372925*1e3),
                                    name="light cylinder")

        if legs:
            leg_radius = 0.07705
            first_leg_angle = angle + 22.5
            if radius != 0.0:
                for i in range(8):
                    leg_angle = first_leg_angle + i * 45
                    leg_n = Cylinder(radius=0.00455, height=0.15, parent=self.world, material=AbsorbingSurface(),
                                     transform=translate(-leg_radius * cos(radians(leg_angle)),
                                                         leg_radius * sin(radians(leg_angle)),
                                                         -0.0454))

    def calculate_top_loss(self):
        real = RealLamp()
        real.get_data(lamp_radius=self.lamp_radius, lamp_angle=self.lamp_angle)
        # real.convert_to_watt()
        real_top, _ = real.provide_data()
        real_top_frac = real_top/np.sum(real_top)
        top_data_frac = self.top_data/np.sum(self.top_data)
        # top_loss = np.sum(np.abs(top_data_frac - real_top_frac)/np.sum(real_top_frac))*100
        return real_top_frac, top_data_frac  # top_loss

    def calculate_out_loss(self):
        real = RealLamp()
        real.get_data(lamp_radius=self.lamp_radius, lamp_angle=self.lamp_angle)
        # real.convert_to_watt()
        _, real_out = real.provide_data()
        real_out_frac = real_out / np.sum(real_out)
        out_data_frac = self.out_data / np.sum(self.out_data)
        # out_loss = np.mean(np.abs(out_data_frac - real_out_frac)/np.sum(real_out_frac))*100
        return real_out_frac, out_data_frac  # out_loss

    def plot_compare(self):
        real = RealLamp()
        real.get_data(lamp_radius=self.lamp_radius, lamp_angle=self.lamp_angle)
        real.convert_to_watt()
        real_top, real_out = real.provide_data()

        # plt.figure()
        # plt.subplot(337)
        # angles = np.linspace(0, 360, 360)
        # vessel_x = self.vessel_in_rad * np.cos(angles)
        # vessel_y = self.vessel_in_rad * np.sin(angles)
        # plt.scatter(x=vessel_x, y=vessel_y)
        # plt.scatter(x=self.lamp_radius * np.cos(radians(self.lamp_angle)),
        #             y=self.lamp_radius * np.sin(radians(self.lamp_angle)), c="r",
        #             s=100)
        # plt.scatter(x=0.100, y=0, c='g', s=50)
        # plt.scatter(x=0.005, y=0.100, c='g', s=50)
        # plt.gca().set_aspect('equal', adjustable='box')
        #
        # cam_nr = np.arange(16)
        #
        # plt.subplot(334)
        # plt.bar(cam_nr, self.top_data)
        # plt.xticks(cam_nr, cam_nr)
        # plt.xlabel("SIMULATED TOP total power (P = " + str(self.top_data_sum) + ")")
        # plt.ylabel("Pixel number")
        #
        # plt.subplot(338)
        # plt.barh(cam_nr, np.flip(self.out_data))
        # plt.yticks(cam_nr, np.flip(cam_nr))
        # plt.xlabel("SIMULATED OUTER total power (P = " + str(self.out_data_sum) + ")")
        # plt.ylabel("Pixel number")
        #
        # plt.subplot(331)
        # cam_nr = np.arange(16)
        # sum_top_pixels = sum(real_top)
        # plt.bar(cam_nr, real_top)
        # plt.xticks(cam_nr, cam_nr)
        # plt.xlabel("REAL TOP total power (P ~~ " + str(sum_top_pixels) + ")")
        # plt.ylabel("Pixel number")
        #
        # plt.subplot(339)
        # sum_out_pixels = sum(real_out)
        # real_out = real_out
        # plt.barh(cam_nr, np.flip(real_out))
        # plt.yticks(cam_nr, np.flip(cam_nr))
        # plt.xlabel("REAL OUTER (P ~~ " + str(sum_out_pixels) + ")")
        # plt.ylabel("Pixel number")
        #
        # plt.suptitle("Radius: {}, angle: {}".format(self.lamp_radius, self.lamp_angle))
        # plt.savefig("./comparison_plots/r" + str(int(self.lamp_radius * 100)) + "_ang" + str(int(self.lamp_angle)) + ".png")

        print("THIS!!!!")

        real_top = real_top/sum(real_top)
        real_out = real_out/sum(real_out)
        self.top_data = np.array(self.top_data)/sum(self.top_data)
        self.out_data = np.array(self.out_data)/sum(self.out_data)

        cam_nr = np.arange(16)

        width=0.4
        leg_size = 15
        font_size=15

        plt.rc('xtick', labelsize=font_size)
        plt.rc('ytick', labelsize=font_size)

        plt.figure()
        plt.bar(cam_nr - width/2, real_top, width, color='r', bottom=0, label='Real')
        plt.bar(cam_nr + width/2, self.top_data, width, color='b', bottom=0, label='Simulated')
        plt.xticks(cam_nr, cam_nr)
        plt.text(3.0, .8, 'TOP camera', horizontalalignment='center', fontsize=20)
        plt.xlabel("Detector Number", fontsize=font_size)
        plt.ylabel("Fraction of total power received", fontsize=font_size)
        plt.legend(prop={'size': leg_size}, bbox_to_anchor=(0.45,0.85))

        plt.figure()
        plt.bar(cam_nr - width/2, real_out, width, color='r', bottom=0, label='Real')
        plt.bar(cam_nr + width/2, self.out_data, width, color='b', bottom=0, label='Simulated')
        plt.text(3.5, .66, 'OUTER camera', horizontalalignment='center', fontsize=20)
        plt.xticks(cam_nr, cam_nr)
        plt.xlabel("Detector Number", fontsize=font_size)
        plt.ylabel("Fraction of total power received", fontsize=font_size)
        plt.legend(prop={'size': leg_size}, bbox_to_anchor=(0.42,0.84))


    @classmethod
    def create_save(cls):
        with open(cls.csv_file_name, "a+") as file:
            file.write("radius, angle, top0, top1, top2, top3, top4, top5, top6, top7, top8, top9, top10, top11, top12,"
                "top13, top14, top15, out0, out1, out2, out3, out4, out5, out6, out7, out8, out9, out10, out11, out12,"
                "out13, out14, out15\n")

    def load_data_to_save(self):
        with open(RayLamp.csv_file_name, "a") as file:
            file.write(str(self.lamp_radius) + "," + str(self.lamp_angle) + ",")
            for pixel in range(16):
                file.write(str(self.top_data[pixel]) + ",")
            for pixel in range(15):
                file.write(str(self.out_data[pixel]) + ",")
            file.write(str(self.out_data[15]) + ",")
            file.write('\n')
