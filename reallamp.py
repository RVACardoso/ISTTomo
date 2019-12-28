from realtomo import RealTomo
import numpy as np


class RealLamp(RealTomo):
    def __init__(self):
        super().__init__()

        self.real_data = []
        with open("resources/withlegs_july_data.csv") as f:
        # with open("legless_july_data.csv") as f:
            content = f.readlines()
        for single_line in content[1:]:
            self.real_data.append(single_line.split(', ')[:-1])

    def get_data(self, lamp_radius, lamp_angle):
        self.lamp_radius = lamp_radius
        self.lamp_angle = lamp_angle

        self.top_pixels = np.zeros(16)
        self.out_pixels = np.zeros(16)
        count = 0
        for line in self.real_data:
            if int(line[1]) == int(self.lamp_radius * 100) and float(line[2]) == self.lamp_angle:
                self.top_pixels += np.array([float(number) for number in line[3:19]])
                self.out_pixels += np.array([float(number) for number in line[19:35]])
                count += 1
            elif count != 0:
                break
        if count == 0:
            print("Requested radius and angle values not found.")
        else:
            integ_time = 0.007
            self.top_pixels = self.top_pixels / (count*integ_time)
            self.out_pixels = self.out_pixels / (count*integ_time)
