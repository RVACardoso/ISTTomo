from rayplasma import RayPlasma
from raysect.optical.material import InhomogeneousVolumeEmitter
from math import floor
import numpy as np
import time
import matplotlib.pyplot as plt

from raysect.optical.material import UniformSurfaceEmitter
from raysect.optical import ConstantSF

from plasmaemit import GenericSinglePlasma
from raysect.optical.material.emitter.inhomogeneous import NumericalIntegrator


class PixelMaterial(InhomogeneousVolumeEmitter):

    def __init__(self, pixel_side, pixel_number):

        self.INTEG_STEP = 0.0005
        super().__init__(integrator=NumericalIntegrator(step=self.INTEG_STEP, min_samples=1))

        self.pixel_number = pixel_number
        self.pixel_side = pixel_side

        self.pixel_size = 0.2/self.pixel_side
        # print(floor(self.pixel_number/self.pixel_side))
        # print()
        self.pos_x = -(self.pixel_number % self.pixel_side)*self.pixel_size + 0.1
        self.pos_y = -floor(self.pixel_number/self.pixel_side)*self.pixel_size + 0.1

        # print(self.pixel_side)
        # print(self.pixel_number)
        # print(self.pos_x)
        # print(self.pos_y)

    def compute_shape(self, x_array, y_array):
        out = []
        for x, y in zip(x_array, y_array):
            if x < self.pos_x and x > self.pos_x-self.pixel_size and y < self.pos_y and y > self.pos_y-self.pixel_size:
                out.append(1.0)
            else:
                out.append(0.0)
        return np.array(out)

    def emission_function(self, point, direction, spectrum, world, ray, primitive, to_local, to_world):

        if point.x < self.pos_x and point.x > self.pos_x-self.pixel_size and point.y < self.pos_y and point.y > self.pos_y-self.pixel_size:
            spectrum.samples += 1.0
        else:
            spectrum.samples += 0.0
        return spectrum


## Generate geometric matrix

# RayPlasma.csv_file_name = 'csv_test1.csv'
pixel_side = 60

ray_plasma = RayPlasma(coord='cartesian', pixel_nr=pixel_side, count=1)
ray_plasma.add_isttok(pixel_samples=5000)
ray_plasma.place_source()


start = time.time()

geo_matrix = np.empty((32, 0))

for pixel in range(pixel_side*pixel_side):
    print(pixel)
    one_pixel = PixelMaterial(pixel_side=pixel_side, pixel_number=pixel)
    ray_plasma.source.material = one_pixel
    # ray_plasma.check_scene()

    ray_plasma.simulate_top_rays()
    ray_plasma.top_data = np.array(ray_plasma.top_data)
    ray_plasma.simulate_out_rays()
    ray_plasma.out_data = np.array(ray_plasma.out_data)

    # print(ray_plasma.top_data.shape)
    # print(ray_plasma.out_data.shape)
    # print()
    # print(ray_plasma.top_data)
    # print(ray_plasma.out_data)
    column = np.hstack((ray_plasma.top_data, ray_plasma.out_data)).reshape((32, 1))
    # print(column.shape)
    # print(column)
    geo_matrix = np.hstack((geo_matrix, column))
    # ray_plasma.get_tomogram_data(pixel_div=10, sinth_plasma=one_pixel)
    # ray_plasma.plot_tomogram()
    # ray_plasma.plot_data(fig_name='ahh')
    # print(geo_matrix.shape)
    # print(geo_matrix)
    # print()
    # plt.show()

print(time.time() - start)

print(geo_matrix.shape)
np.save("geo_mat60-newpin-ref-px5000", geo_matrix)


## Plot geometric matrix

from matplotlib.colors import SymLogNorm
import matplotlib
# from cmap import get_cmap

loaded_matrix = np.load("geo_mat60-newpin-ref-px5000.npy")
# for i in range(32):

# foto capa
plt.imshow(loaded_matrix[28].reshape(60, 60), cmap=matplotlib.cm.get_cmap("plasma"), norm=SymLogNorm(75e-15))
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    left=False,         # ticks along the top edge are off
    labelleft=False) # labels along the bottom edge are off

plt.show()



import matplotlib
font = {'family' : 'normal',
        'size'   : 25}
matplotlib.rc('font', **font)

matplotlib.rcParams['axes.titlepad'] = 20

plt.subplot(121)
plt.imshow(loaded_matrix[5].reshape(60, 60), cmap=matplotlib.cm.get_cmap("plasma"), norm=SymLogNorm(100e-15))
plt.title("Viewing geometry for\nTop detector #5", fontsize=24)
cbar = plt.colorbar(fraction=0.046, pad=0.04)
cbar.ax.set_ylabel('Influence ($m^3$)', rotation=90)
plt.xticks(np.linspace(0, 59, num=3),
           np.round(np.arange(start=-0.1, stop=0.1 + 0.1, step=0.1), 2))
plt.yticks(np.linspace(0, 59, num=5),
           np.flipud(np.round(np.arange(start=-0.1, stop=0.1 + 0.02, step=0.05), 2)))
plt.xlabel("x (m)")
plt.ylabel("z (m)")


plt.subplot(122)
plt.imshow(loaded_matrix[31].reshape(60, 60), cmap=matplotlib.cm.get_cmap("plasma"), norm=SymLogNorm(100e-15))
plt.title("Viewing geometry for\nOuter detector #15", fontsize=24)
cbar = plt.colorbar(fraction=0.046, pad=0.04)
cbar.ax.set_ylabel('Influence ($m^3$)', rotation=90)
plt.xticks(np.linspace(0, 59, num=3),
           np.round(np.arange(start=-0.1, stop=0.1 + 0.1, step=0.1), 2))
plt.yticks(np.linspace(0, 59, num=5),
           np.flipud(np.round(np.arange(start=-0.1, stop=0.1 + 0.02, step=0.05), 2)))
plt.xlabel("x (m)")
plt.ylabel("z (m)")

plt.subplots_adjust(left=None, bottom=0.2, right=None, top=0.8, wspace=0.65, hspace=0.5)

plt.show()

# plt.subplot(121)
# plt.imshow(loaded_matrix[0].reshape(60,60))
# plt.subplot(122)
# plt.imshow(loaded_matrix[22].reshape(60,60))
# plt.show()

# matrix = np.zeros((60,60))
# for i in range(32):
#     matrix += loaded_matrix[i].reshape(60,60)
# # matrix += loaded_matrix[17].reshape(60, 60)
#
# plt.imshow(matrix)
# plt.show()

# print(loaded_matrix.shape)
# print(sum(loaded_matrix.flatten()))
# print(loaded_matrix[8][1250])
# print(loaded_matrix)

# ray_plasma.get_tomogram_data(pixel_div=10, sinth_plasma=one_pixel)
# ray_plasma.load_data_to_save()
# print(time.time() - start)

# ray_plasma.plot_data(fig_name='fig name')
# ray_plasma.plot_tomogram()
# plt.show()