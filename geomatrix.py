from rayplasma import RayPlasma
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from cmap import get_cmap


class GeoMatrix(RayPlasma):

    def __init__(self, matrix_file=None, matrix_data=None):
        self.pre_pixel_nr = 60
        super().__init__(coord='cartesian', pixel_nr=self.pre_pixel_nr, count=1)
        if matrix_file != None:
            self.geo_matrix = np.load(matrix_file)
        else:
            self.geo_matrix = matrix_data

    # get tomogram data

    def get_matrix_sensors(self):
        sensors = np.dot(self.geo_matrix, self.square_values)
        self.top_data = sensors[:16]
        self.top_data_sum = np.sum(self.top_data)
        self.out_data = sensors[16:]
        self.out_data_sum = np.sum(self.out_data)


    def plot_matrix(self, detector_nr):
        plt.rcParams.update({'font.size': 22})
        detector_view = self.geo_matrix[detector_nr].reshape((self.pre_pixel_nr, self.pre_pixel_nr))
        print(detector_view[0, 0].shape)
        print(detector_view.flatten()[25])
        # plt.imshow(detector_view, cmap=get_cmap())
        plt.imshow(detector_view, norm=matplotlib.colors.SymLogNorm(1e-14), cmap=get_cmap())
        cbar = plt.colorbar()  # format='%.0e')
        cbar.ax.set_ylabel('Influence ($m^3$)', rotation=90, fontsize=22)
        plt.xlabel("#pixel x")
        plt.ylabel("#pixel z")



    # plot tomogram

    # plot original simulations with raysect