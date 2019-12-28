import numpy as np
import matplotlib.pyplot as plt
from math import radians


class RealTomo:
    def __init__(self):
        self.top_pixels, self.out_pixels = [], []
        self.shot_nr = None
        self.lamp_radius = None
        self.lamp_angle = None

    def get_data(self):
        raise NotImplementedError

    def convert_to_watt(self):
        self.top_pixels = self.top_pixels / (1.175 * 1e6)
        self.out_pixels = self.out_pixels / (1.175 * 1e6)

    def plot_data(self):
        plt.figure()
        plt.subplot(223)
        angles = np.linspace(0, 360, 360)
        vessel_x = 0.102 * np.cos(angles)
        vessel_y = 0.102 * np.sin(angles)
        plt.scatter(x=vessel_x, y=vessel_y)
        if self.lamp_radius != None:
            plt.scatter(x=self.lamp_radius * np.cos(radians(self.lamp_angle)),
                        y=self.lamp_radius * np.sin(radians(self.lamp_angle)),
                        c="r", s=100)
        else:
            pass
        plt.scatter(x=0.102, y=0, c='g', s=50)
        plt.scatter(x=0.005, y=0.1018, c='g', s=50)
        plt.gca().set_aspect('equal', adjustable='box')

        cam_nr = np.arange(16)

        plt.subplot(221)
        sum_top_pixels = sum(self.top_pixels)
        plt.bar(cam_nr, self.top_pixels)
        plt.xticks(cam_nr, cam_nr)
        plt.ylabel("Fraction of TOP total power (P ~~ " + str(sum_top_pixels) + ")")
        plt.xlabel("Pixel number")

        plt.subplot(224)
        sum_out_pixels = sum(self.out_pixels)
        plt.barh(cam_nr, np.flip(self.out_pixels))
        plt.yticks(cam_nr, np.flip(cam_nr))
        plt.xlabel("Fraction of OUTER total power (P ~~ " + str(sum_out_pixels) + ")")
        plt.ylabel("Pixel number")

        print("Out/top power ratio")
        print(sum_out_pixels/sum_top_pixels)

        if self.lamp_radius != None:
            plt.suptitle("REAL DATA" + " - r=" + str(int(self.lamp_radius * 100)) + ", ang=" + str(self.lamp_angle))
        else:
            plt.suptitle("REAL DATA - " + str(self.shot_nr))

        # plt.savefig("./withlegs_real_data_plots/" + "real_r"+str(int(lamp_radius*100))+"_ang"+str(lamp_angle) + ".png")

    def new_plot_data(self):

        import matplotlib
        font = {'family': 'normal',
                'size': 20}
        matplotlib.rc('font', **font)

        plt.figure()
        cam_nr = np.arange(16)

        plt.subplot(121)
        sum_top_pixels = sum(self.top_pixels)
        plt.bar(cam_nr, self.top_pixels)
        plt.xticks(cam_nr, cam_nr)
        # plt.ylabel("Fraction of TOP total power (P ~~ " + str(sum_top_pixels) + ")")
        plt.ylabel("Detector measurements (V)")
        plt.xlabel("#detector")
        plt.title("Top camera")

        plt.subplot(122)
        sum_out_pixels = sum(self.out_pixels)
        plt.bar(cam_nr, self.out_pixels)
        plt.xticks(cam_nr, cam_nr)
        plt.ylabel("Detector measurements (V)")
        plt.xlabel("#detector")
        plt.title("Outer camera")

        # print("Out/top power ratio")
        # print(sum_out_pixels / sum_top_pixels)
        plt.suptitle("Average detector measurements for shot 47258 between t=300ms and t=305ms")

        # if self.lamp_radius != None:
        #     plt.suptitle("REAL DATA" + " - r=" + str(int(self.lamp_radius * 100)) + ", ang=" + str(self.lamp_angle))
        # else:
        #     plt.suptitle("REAL DATA - " + str(self.shot_nr))

        # plt.savefig("./withlegs_real_data_plots/" + "real_r"+str(int(lamp_radius*100))+"_ang"+str(lamp_angle) + ".png")

    def provide_data(self):
        return self.top_pixels, self.out_pixels
