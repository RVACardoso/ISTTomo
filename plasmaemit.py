import numpy as np
from math import pi, sqrt, log, exp, pow
import sys
import warnings
from scipy import interpolate
from scipy.spatial import cKDTree
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

from raysect.optical.material import InhomogeneousVolumeEmitter
from raysect.optical.material.emitter.inhomogeneous import NumericalIntegrator

from plasmainteg import MyTestNumIntegrator


class ModePlasmaEmit(InhomogeneousVolumeEmitter):

    def __init__(self, shape, center, width, count, power_gain):
        super().__init__(integrator=NumericalIntegrator(step=0.001, min_samples=5))
        self.count = count
        self.center = center
        self.width = width
        self.power_gain = power_gain

        if shape == "gaussian":
            self.plaama_shape = self.gaussian
        elif shape == "lorentzian":
            self.plasma_shape = self.lorentzian
        elif shape == "quadratic":
            self.plasma_shape = self.quadratic
        elif shape == "flat_top":
            self.plasma_shape = self.flat_top
        elif shape == "beta":
            self.plasma_shape = self.beta_dist
        else:
            print("Plasma shape must be gaussian, lorentzian, quadratic, flat_top or beta.")
            sys.exit()

    @staticmethod
    def gaussian(x, y, mu_x, mu_y, fwhm_x, fwhm_y):
        gauss_amp = 1.0
        sig_x = fwhm_x / 2.3548
        sig_y = fwhm_y / 2.3548
        return gauss_amp * np.exp(
            -np.power(x - mu_x, 2.) / (2 * np.power(sig_x, 2.)) - np.power(y - mu_y, 2.) / (2 * np.power(sig_y, 2.)))

    @staticmethod
    def lorentzian(x, y, mu_x, mu_y, fwhm_x, fwhm_y):
        lorentz_amp = 1200.0  # 1.5708660015462401
        wid_x = 2 / fwhm_x  # /2
        wid_y = 2 / fwhm_y  # /2
        return lorentz_amp * (1 / (1 + ((y - mu_y) * wid_y) ** 2)) * (1 / (1 + ((x - mu_x) * wid_x) ** 2))

    @staticmethod
    def quadratic(x, y, mu_x, mu_y, fwhm_x, fwhm_y):
        quad_amp = 1.0  # 1.5708660015462401
        wid_x = 2 / (fwhm_x * fwhm_x)
        wid_y = 2 / (fwhm_y * fwhm_y)
        wid2 = 1.0
        return (quad_amp * (-wid_x * (x - mu_x) ** 2 - wid_y * (y - mu_y) ** 2 + wid2 * (x - mu_x) * (y - mu_y) + 1.0)).clip(min=0)

    @staticmethod
    def flat_top(x, y, mu_x, mu_y, fwhm_x, fwhm_y):
        flat_amp = 1.0
        wid_x = fwhm_x / 1.8248
        wid_y = fwhm_y / 1.8248
        n = 2
        return flat_amp * np.exp(- ((x - mu_x) / wid_x) ** (2 * n) - ((y - mu_y) / wid_y) ** (2 * n))

    @staticmethod
    def beta_dist(x, y, mu_x, mu_y, fwhm_x, fwhm_y):
        x = 10.0 * x #.copy()
        y = 10.0 * y #.copy()
        fwhm_x *= 10.0
        fwhm_y *= 10.0

        beta_amp = 1.0
        wid_x = -1.0 / (log(1.0 - fwhm_x * fwhm_x, 2))  # 75.0 # fwhm_x +1.2#- 1.0
        wid_y = -1.0 / (log(1.0 - fwhm_y * fwhm_y, 2))  # 75.0 # fwhm_y +1.2 #- 1.0
        x = x + 0.5 - mu_x * 10.0
        y = y + 0.5 - mu_y * 10.0
        # print(4 ** (wid_x) * 1e-300 * 4 ** (wid_x))
        # print(sum(beta_amp * (x*(1.0-x))**wid_x * (y*(1.0-y))**wid_y))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            prenorm_result = beta_amp * (x * (1.0 - x)) ** wid_x * (y * (1.0 - y)) ** wid_y * 4 ** wid_x * 4 ** wid_y
        prenorm_result = np.nan_to_num(prenorm_result)
        # norm = np.max(prenorm_result)
        norm = beta_amp * 0.25 ** wid_x * 0.25 ** wid_y * 4 ** wid_x * 4 ** wid_y
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return np.array(beta_amp * prenorm_result / norm, np.float64)

    def emission_function(self, point, direction, spectrum, world, ray, primitive, to_local, to_world):
        for i in range(self.count):
            spectrum.samples += self.power_gain*self.plasma_shape(point.x, point.y, mu_x=self.center[i][0], mu_y=self.center[i][1],
                                                           fwhm_x=self.width[i][0], fwhm_y=self.width[i][1])

        return spectrum


class SinglePlasmaEmit(InhomogeneousVolumeEmitter):

    def __init__(self, shape, center, width, power_gain=1.0, flatlen=None):
        super().__init__(integrator=NumericalIntegrator(step=0.005, min_samples=1))
        # super().__init__(integrator=MyTestNumIntegrator(step=0.001, min_samples=5))
        self.center = np.array(center[0])
        self.width = np.array(width[0])*2.0
        self.power_gain = power_gain
        self.wid_x = None
        self.wid_y = None

        if shape == "gaussian":
            self.plasma_shape = self.gaussian
            self.wid_x_up = self.width[0][0] / 2.3548
            self.wid_x_dn = self.width[0][1] / 2.3548
            self.wid_y_up = self.width[1][1] / 2.3548
            self.wid_y_dn = self.width[1][0] / 2.3548
            self.mu_x = self.center[0]
            self.mu_y = self.center[1]
        elif shape == "lorentzian":
            self.plasma_shape = self.lorentzian
            self.wid_x_up = 2.0/self.width[0][0]
            self.wid_x_dn = 2.0/self.width[0][1]
            self.wid_y_up = 2.0/self.width[1][1]
            self.wid_y_dn = 2.0/self.width[1][0]
            self.mu_x = self.center[0]
            self.mu_y = self.center[1]
        elif shape == "quadratic":
            self.plasma_shape = self.quadratic
            self.wid_x_up = 2.0 / (self.width[0][0] * self.width[0][0])
            self.wid_x_dn = 2.0 / (self.width[0][1] * self.width[0][1])
            self.wid_y_up = 2.0 / (self.width[1][1] * self.width[1][1])
            self.wid_y_dn = 2.0 / (self.width[1][0] * self.width[1][0])
            self.mu_x = self.center[0]
            self.mu_y = self.center[1]
        elif shape == "flat_top":
            self.plasma_shape = self.flat_top
            self.wid_x_up = self.width[0][0] / 1.8248
            self.wid_x_dn = self.width[0][1] / 1.8248
            self.wid_y_up = self.width[1][1] / 1.8248
            self.wid_y_dn = self.width[1][0] / 1.8248
            self.mu_x = self.center[0]
            self.mu_y = self.center[1]
        # elif shape == "beta":
        #     self.plasma_shape = self.beta
        #     self.scaling = 5.0
        #     self.wid_x_up = -1.0 / (log(1.0 - 100.0*self.width[0][0] * self.width[0][0], 2))
        #     self.wid_x_dn = -1.0 / (log(1.0 - 100.0*self.width[0][1] * self.width[0][1], 2))
        #     self.wid_y_up = -1.0 / (log(1.0 - 100.0*self.width[1][0] * self.width[1][0], 2))
        #     self.wid_y_dn = -1.0 / (log(1.0 - 100.0*self.width[1][1] * self.width[1][1], 2))
        #     self.mu_x = 0.5 - self.center[0] * self.scaling
        #     self.mu_y = 0.5 - self.center[1] * self.scaling
        elif shape == "bell":
            self.plasma_shape = self.bell
            self.wid_x_up = self.width[0][0]/2.0
            self.wid_x_dn = self.width[0][1]/2.0
            self.wid_y_up = self.width[1][1]/2.0
            self.wid_y_dn = self.width[1][0]/2.0
            self.mu_x = self.center[0]
            self.mu_y = self.center[1]
            self.flat = flatlen
        else:
            print("Plasma shape must be gaussian, lorentzian, quadratic, flat_top or beta (nope).")
            sys.exit()

    def set_skewness(self, x, y):
        if x >= self.mu_x and y >= self.mu_y:
            self.wid_x = self.wid_x_up
            self.wid_y = self.wid_y_up
            # self.norm = (0.25 ** self.wid_x_up) * (0.25 ** self.wid_y_up)
        elif x > self.mu_x and y < self.mu_y:
            self.wid_x = self.wid_x_up
            self.wid_y = self.wid_y_dn
            # self.norm = (0.25 ** self.wid_x_up) * (0.25 ** self.wid_y_dn)
        elif x < self.mu_x and y > self.mu_y:
            self.wid_x = self.wid_x_dn
            self.wid_y = self.wid_y_up
            # self.norm = (0.25 ** self.wid_x_dn) * (0.25 ** self.wid_y_up)
        else:
            self.wid_x = self.wid_x_dn
            self.wid_y = self.wid_y_dn
            # self.norm = (0.25 ** self.wid_x_dn) * (0.25 ** self.wid_y_dn)

    def gaussian(self, x, y):
        self.set_skewness(x=x, y=y)
        gauss_amp = 1.0
        return gauss_amp * exp(
                    -((x - self.mu_x) ** 2 / (2 * self.wid_x * self.wid_x))
                    -((y - self.mu_y) ** 2 / (2 * self.wid_y * self.wid_y)))

    def lorentzian(self, x, y):
        self.set_skewness(x=x, y=y)
        lorentz_amp = 1.0
        return lorentz_amp * (1 / (1 + ((y - self.mu_y) * self.wid_y) * ((y - self.mu_y) * self.wid_y))) * \
               (1 / (1 + ((x - self.mu_x) * self.wid_x) * ((x - self.mu_x) * self.wid_x) ))

    def quadratic(self, x, y):
        self.set_skewness(x=x, y=y)
        quad_amp = 1.0
        result = quad_amp * (-self.wid_x * (x - self.mu_x)**2 - self.wid_y * (y - self.mu_y)**2 + 1.0)
        if result >=0:
            return result
        else:
            return 0

    def flat_top(self, x, y):
        self.set_skewness(x=x, y=y)
        flat_amp = 1.0
        return flat_amp * exp(- ((x - self.mu_x) / self.wid_x) ** 4 - ((y - self.mu_y) / self.wid_y) ** 4)

    # def beta(self, x, y):
    #     x = self.scaling * x + self.mu_x #.copy()
    #     y = self.scaling * y + self.mu_y #.copy()
    #     self.set_skewness(x=x, y=y)
    #     beta_amp = 1.0
    #     try:
    #         return beta_amp * (pow((x * (1.0 - x)), self.wid_x)) * (pow((y * (1.0 - y)), self.wid_y)) / self.norm
    #     except ValueError:
    #         # print('excep')
    #         return 0

    def bell(self, x, y):
        self.set_skewness(x=x, y=y)
        bell_amp = 1.0
        return bell_amp * (1.0/(1.0+(abs((x-self.mu_x)/self.wid_x)) ** self.flat)) * (1.0/(1.0+abs((y-self.mu_y)/self.wid_y) ** self.flat))


    def compute_shape(self, x, y):
        return self.power_gain * self.plasma_shape(x=x, y=y)


    def emission_function(self, point, direction, spectrum, world, ray, primitive, to_local, to_world):
        spectrum.samples += self.power_gain*self.plasma_shape(x=point.x, y=point.y)
        return spectrum


class NpySinglePlasmaEmit(InhomogeneousVolumeEmitter):

    def __init__(self, shape, center, width, power_gain=1.0, flatlen=None):
        # super().__init__(integrator=NumericalIntegrator(step=0.001, min_samples=5))
        super().__init__(integrator=MyTestNumIntegrator(step=0.01, min_samples=5))
        self.center = np.array(center[0])
        self.width = np.array(width[0])*2.0
        self.power_gain = power_gain
        self.wid_x = None
        self.wid_y = None

        if shape == "gaussian":
            self.plasma_shape = self.gaussian
            self.wid_x_up = self.width[0][0] / 2.3548
            self.wid_x_dn = self.width[0][1] / 2.3548
            self.wid_y_up = self.width[1][1] / 2.3548
            self.wid_y_dn = self.width[1][0] / 2.3548
            self.mu_x = self.center[0]
            self.mu_y = self.center[1]
        elif shape == "lorentzian":
            self.plasma_shape = self.lorentzian
            self.wid_x_up = 2.0/self.width[0][0]
            self.wid_x_dn = 2.0/self.width[0][1]
            self.wid_y_up = 2.0/self.width[1][1]
            self.wid_y_dn = 2.0/self.width[1][0]
            self.mu_x = self.center[0]
            self.mu_y = self.center[1]
        elif shape == "quadratic":
            self.plasma_shape = self.quadratic
            self.wid_x_up = 2.0 / (self.width[0][0] * self.width[0][0])
            self.wid_x_dn = 2.0 / (self.width[0][1] * self.width[0][1])
            self.wid_y_up = 2.0 / (self.width[1][1] * self.width[1][1])
            self.wid_y_dn = 2.0 / (self.width[1][0] * self.width[1][0])
            self.mu_x = self.center[0]
            self.mu_y = self.center[1]
        elif shape == "flat_top":
            self.plasma_shape = self.flat_top
            self.wid_x_up = self.width[0][0] / 1.8248
            self.wid_x_dn = self.width[0][1] / 1.8248
            self.wid_y_up = self.width[1][1] / 1.8248
            self.wid_y_dn = self.width[1][0] / 1.8248
            self.mu_x = self.center[0]
            self.mu_y = self.center[1]
        # elif shape == "beta":
        #     self.plasma_shape = self.beta
        #     self.scaling = 5.0
        #     self.wid_x_up = -1.0 / (log(1.0 - 100.0*self.width[0][0] * self.width[0][0], 2))
        #     self.wid_x_dn = -1.0 / (log(1.0 - 100.0*self.width[0][1] * self.width[0][1], 2))
        #     self.wid_y_up = -1.0 / (log(1.0 - 100.0*self.width[1][0] * self.width[1][0], 2))
        #     self.wid_y_dn = -1.0 / (log(1.0 - 100.0*self.width[1][1] * self.width[1][1], 2))
        #     self.mu_x = 0.5 - self.center[0] * self.scaling
        #     self.mu_y = 0.5 - self.center[1] * self.scaling
        elif shape == "bell":
            self.plasma_shape = self.bell
            self.wid_x_up = self.width[0][0]/2.0
            self.wid_x_dn = self.width[0][1]/2.0
            self.wid_y_up = self.width[1][1]/2.0
            self.wid_y_dn = self.width[1][0]/2.0
            self.mu_x = self.center[0]
            self.mu_y = self.center[1]
            self.flat = flatlen
        else:
            print("Plasma shape must be gaussian, lorentzian, quadratic, flat_top or beta (nope).")
            sys.exit()

    def compute_branches(self, x, y):
        upup_idx = np.intersect1d(np.where(x >= self.mu_x)[0], np.where(y > self.mu_y)[0])
        updn_idx = np.intersect1d(np.where(x > self.mu_x)[0], np.where(y <= self.mu_y)[0])
        dnup_idx = np.intersect1d(np.where(x < self.mu_x)[0], np.where(y >= self.mu_y)[0])
        dndn_idx = np.intersect1d(np.where(x <= self.mu_x)[0], np.where(y < self.mu_y)[0])

        out = np.zeros(len(x))
        out[upup_idx] = self.plasma_shape(x=np.take(x, upup_idx), y=np.take(y, upup_idx), wid_x=self.wid_x_up, wid_y=self.wid_y_up)
        out[updn_idx] = self.plasma_shape(x=np.take(x, updn_idx), y=np.take(y, updn_idx), wid_x=self.wid_x_up, wid_y=self.wid_y_dn)
        out[dnup_idx] = self.plasma_shape(x=np.take(x, dnup_idx), y=np.take(y, dnup_idx), wid_x=self.wid_x_dn, wid_y=self.wid_y_up)
        out[dndn_idx] = self.plasma_shape(x=np.take(x, dndn_idx), y=np.take(y, dndn_idx), wid_x=self.wid_x_dn, wid_y=self.wid_y_dn)

        return out

    def gaussian(self, x, y, wid_x, wid_y):
        gauss_amp = 1.0
        return gauss_amp * np.exp(
                    -((x - self.mu_x) ** 2 / (2 * wid_x * wid_x))
                    -((y - self.mu_y) ** 2 / (2 * wid_y * wid_y)))

    def lorentzian(self, x, y, wid_x, wid_y):
        lorentz_amp = 1.0
        return lorentz_amp * (1 / (1 + ((y - self.mu_y) * wid_y) * ((y - self.mu_y) * wid_y))) * \
               (1 / (1 + ((x - self.mu_x) * wid_x) * ((x - self.mu_x) * wid_x) ))

    def quadratic(self, x, y, wid_x, wid_y):
        quad_amp = 1.0
        result = quad_amp * (-wid_x * (x - self.mu_x)**2 - wid_y * (y - self.mu_y)**2 + 1.0)
        return result.clip(min=0)

    def flat_top(self, x, y, wid_x, wid_y):
        flat_amp = 1.0
        return flat_amp * np.exp(- ((x - self.mu_x) / wid_x) ** 4 - ((y - self.mu_y) / wid_y) ** 4)

    # def beta(self, x, y):
    #     x = self.scaling * x + self.mu_x #.copy()
    #     y = self.scaling * y + self.mu_y #.copy()
    #     self.set_skewness(x=x, y=y)
    #     beta_amp = 1.0
    #     try:
    #         return beta_amp * (pow((x * (1.0 - x)), self.wid_x)) * (pow((y * (1.0 - y)), self.wid_y)) / self.norm
    #     except ValueError:
    #         # print('excep')
    #         return 0

    def bell(self, x, y, wid_x, wid_y):
        bell_amp = 1.0
        return bell_amp * (1.0/(1.0+(abs((x-self.mu_x)/wid_x)) ** self.flat)) * (1.0/(1.0+abs((y-self.mu_y)/wid_y) ** self.flat))

    def compute_shape(self, x, y):
        return self.power_gain * self.compute_branches(x=x, y=y)

    def emission_function(self, points, direction, spectrum, world, ray, primitive, to_local, to_world):
        # print(type(points[0][0]))
        return self.compute_shape(x=points[0], y=points[1])
        # spectrum.samples += 0.0 #self.compute_shape(x=points[0], y=points[1])
        # return spectrum


class GenericSinglePlasma(InhomogeneousVolumeEmitter):
    def __init__(self, contour_type, stacking_type, grid_side):

        self.INTEG_STEP = 0.005
        self.CONTOUR_NUM = 50
        self.Z_NUM = 30

        super().__init__(integrator=NumericalIntegrator(step=self.INTEG_STEP, min_samples=1))
        # super().__init__(integrator=MyTestNumIntegrator(step=self.INTEG_STEP, min_samples=2))

        self.contour_type = contour_type
        self.stacking = stacking_type
        self.grid_side = grid_side

        self.make_profile()
        self.create_KDTree()
        self.create_grid()

    def get_contour(self):
        contour_name = self.contour_type + '.npy'
        tck = np.load(contour_name)
        xi, yi = interpolate.splev(np.linspace(0, 1, self.CONTOUR_NUM), tck)
        contour_x = 2.0 * (xi - np.mean(xi)) / np.max(xi)
        contour_y = 2.0 * (yi - np.mean(yi)) / np.max(yi)
        return 0.03*contour_x, 0.03*contour_y
        # plt.figure()
        # # plt.plot(x_new, y_new, 'or')
        # plt.plot(xi, yi, 'ob')

    @staticmethod
    def invert_gaussian(z, sigma=1.0):
        return sqrt(-2 * sigma * sigma * log(z))

    @staticmethod
    def invert_quad(z):
        return sqrt(1.0 - z)

    @staticmethod
    def invert_bell(z):
        return 0.1 * (log(-1.0 + 1.0 / z) - np.log(0.010101)) / log(2)

    @staticmethod
    def invert_flat_top(z):
        return sqrt(sqrt(-log(z)))

    def make_profile(self):
        contour_x, contour_y = self.get_contour()

        if self.stacking == 'gaussian':
            stack_func = self.invert_gaussian
        elif self.stacking == 'quadratic':
            stack_func = self.invert_quad
        elif self.stacking == 'bell':
            stack_func = self.invert_bell
        elif self.stacking == 'flat_top':
            stack_func = self.invert_flat_top
        else:
            print("ERROR on stacking function attribution")

        x, y, z = [], [], []
        for z_temp in np.linspace(0.0000001, 0.9999999, num=self.Z_NUM, endpoint=True):
            amp_x = stack_func(z_temp)
            amp_y = amp_x
            for x_temp, y_temp in zip(amp_x * contour_x, amp_y * contour_y):
                x.append(x_temp)
                y.append(y_temp)
                z.append(z_temp)

        self.x = np.array(x)
        self.y = np.array(y)
        self.z = np.array(z)

        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')
        # ax.plot_trisurf(x, y, z, color='r')

    def create_KDTree(self):
        xy = np.dstack((self.x, self.y))[0]
        self.mytree = cKDTree(xy)

    def create_grid(self):
        xx = np.linspace(-0.1, 0.1, num=self.grid_side)
        yy = np.linspace(-0.1, 0.1, num=self.grid_side)
        xm, ym = np.meshgrid(xx, yy)
        xm_raveled = np.ravel(xm)
        ym_raveled = np.ravel(ym)
        xy = np.dstack((ym_raveled, xm_raveled))[0]
        self.z_grid = self.z[self.mytree.query(xy)[1]].reshape((self.grid_side, self.grid_side))

    def change_profile(self, center, width, power_gain):
        center = np.array(center[0])
        self.center_x = center[0]
        self.center_y = center[1]
        self.power_gain = power_gain
        width = 1.0/np.array(width[0])
        self.wid_x_up = width[0][0]
        self.wid_x_dn = width[0][1]
        self.wid_y_up = width[1][1]
        self.wid_y_dn = width[1][0]

    def compute_shape(self, x_array, y_array):
        x = x_array - self.center_x
        y = y_array - self.center_y

        upup_idx = np.intersect1d(np.where(x >= 0.0)[0], np.where(y > 0.0)[0])
        updn_idx = np.intersect1d(np.where(x > 0.0)[0], np.where(y <= 0.0)[0])
        dnup_idx = np.intersect1d(np.where(x < 0.0)[0], np.where(y >= 0.0)[0])
        dndn_idx = np.intersect1d(np.where(x <= 0.0)[0], np.where(y < 0.0)[0])

        x[upup_idx] *= self.wid_x_up
        y[upup_idx] *= self.wid_y_up
        x[updn_idx] *= self.wid_x_up
        y[updn_idx] *= self.wid_y_dn
        x[dnup_idx] *= self.wid_x_dn
        y[dnup_idx] *= self.wid_y_up
        x[dndn_idx] *= self.wid_x_dn
        y[dndn_idx] *= self.wid_y_dn

        idx_x = ((x / 0.2 + 0.5) * (self.grid_side)).astype(dtype=int)
        idx_y = ((y / 0.2 + 0.5) * (self.grid_side)).astype(dtype=int)

        x_over = np.where(idx_x >= self.grid_side)
        idx_x[x_over] = -100
        y_over = np.where(idx_y >= self.grid_side)
        idx_y[y_over] = -100

        x_neg = np.where(idx_x < 0.0)
        y_neg = np.where(idx_y < 0.0)

        idx_final = (idx_x, idx_y)
        out = self.power_gain * self.z_grid[idx_final]

        out[x_neg] = 0.0
        out[y_neg] = 0.0

        return out

    def old_compute_shape(self, x_array, y_array):
        # for sanity check
        z_out = []
        for x_val, y_val in zip(x_array, y_array):
            x = x_val - self.center_x
            y = y_val - self.center_y

            if x_val > self.center_x and y_val > self.center_y:
                x *= self.wid_x_up
                y *= self.wid_y_up
            elif x_val > self.center_x and y_val < self.center_y:
                x *= self.wid_x_up
                y *= self.wid_y_dn
            elif x_val < self.center_x and y_val > self.center_y:
                x *= self.wid_x_dn
                y *= self.wid_y_up
            else:
                x *= self.wid_x_dn
                y *= self.wid_y_dn

            idx_x = int((x / 0.2 + 0.5) * self.grid_side)
            idx_y = int((y / 0.2 + 0.5) * self.grid_side)

            if idx_x > 0 and idx_x < self.grid_side and idx_y > 0 and idx_y < self.grid_side:
                out_val = self.power_gain * self.z_grid[idx_x, idx_y]
            else:
                out_val = 0.0
            z_out.append(out_val)
        return np.array(z_out)

    def emission_function(self, point, direction, spectrum, world, ray, primitive, to_local, to_world):

        x = point.x - self.center_x
        y = point.y - self.center_y
        if point.x > self.center_x and point.y > self.center_y:
            x *= self.wid_x_up
            y *= self.wid_y_up
        elif point.x > self.center_x and point.y < self.center_y:
            x *= self.wid_x_up
            y *= self.wid_y_dn
        elif point.x < self.center_x and point.y > self.center_y:
            x *= self.wid_x_dn
            y *= self.wid_y_up
        else:
            x *= self.wid_x_dn
            y *= self.wid_y_dn

        idx_x = int((x/0.2 + 0.5)*self.grid_side)
        idx_y = int((y/0.2 + 0.5)*self.grid_side)

        if idx_x > 0 and idx_x < self.grid_side and idx_y > 0 and idx_y < self.grid_side:
            spectrum.samples += self.power_gain * self.z_grid[idx_x, idx_y]
        else:
            spectrum.samples += 0
        return spectrum

    
    # vectorization is slower!!! -> only helps for tiny integrations steps
    # def emission_function(self, points, direction, spectrum, world, ray, primitive, to_local, to_world):
    #     x = points[0] - self.center_x
    #     y = points[1] - self.center_y
    #
    #     upup_idx = np.intersect1d(np.where(x >= 0.0)[0], np.where(y > 0.0)[0])
    #     updn_idx = np.intersect1d(np.where(x > 0.0)[0], np.where(y <= 0.0)[0])
    #     dnup_idx = np.intersect1d(np.where(x < 0.0)[0], np.where(y >= 0.0)[0])
    #     dndn_idx = np.intersect1d(np.where(x <= 0.0)[0], np.where(y < 0.0)[0])
    #
    #     x[upup_idx] *= self.wid_x_up
    #     y[upup_idx] *= self.wid_y_up
    #     x[updn_idx] *= self.wid_x_up
    #     y[updn_idx] *= self.wid_y_dn
    #     x[dnup_idx] *= self.wid_x_dn
    #     y[dnup_idx] *= self.wid_y_up
    #     x[dndn_idx] *= self.wid_x_dn
    #     y[dndn_idx] *= self.wid_y_dn
    #
    #     # idx_x = ((x / 0.2 + 0.5) * (self.grid_side)).astype(dtype=int)
    #     # idx_y = ((y / 0.2 + 0.5) * (self.grid_side)).astype(dtype=int)
    #     #
    #     # x_over = np.where(abs(idx_x) >= self.grid_side)
    #     # idx_x[x_over] = -100
    #     # y_over = np.where(abs(idx_y) >= self.grid_side)
    #     # idx_y[y_over] = -100
    #     #
    #     # x_neg = np.where(idx_x < 0.0)
    #     # y_neg = np.where(idx_y < 0.0)
    #     #
    #     # idx_final = (idx_x, idx_y)
    #     # out = self.power_gain * self.z_grid[idx_final]
    #     #
    #     # out[x_neg] = 0.0
    #     # out[y_neg] = 0.0
    #
    #     # spectrum.samples += self.power_gain * self.z_grid[idx_x, idx_y]
    #
    #     return np.zeros(self.z_grid.shape)




# from mpl_toolkits.mplot3d import Axes3D
# import matplotlib.pyplot as plt
# import time
#
# start = time.time()
#
# # plasma = GenericSinglePlasma(type='test1', amp_x_init=0.1, amp_y_init=0.1, stacking='yooo', contour_num=50, z_num=20)
# # x, y, z = plasma.make_profile()
# # x, y, z = np.array(x), np.array(y), np.array(z)
# # order = np.argsort(x)
# # x = x[order]
# # y = y[order]
# # z = z[order]
#
# x, y, z = [], [], []
# for x_val in np.linspace(-0.1, 0.1, 20):
#     for y_val in np.linspace(-0.1, 0.1, 200):
#         x.append(x_val)
#         y.append(y_val)
#         z.append(-x_val*x_val+y_val*y_val)
#
# # print(x)
# # print()
# # print(y)
# x, y, z = np.array(x), np.array(y), np.array(z)
#
# permut = np.random.permutation(x.shape[0])
# print(permut)
# x = x[permut]
# y = y[permut]
# z = z[permut]
#
# print(x)
#
# f = interpolate.interp2d(x, y, z, kind='linear')
#
# x_lin = np.linspace(-0.1, 0.1, 50)
# y_lin = np.linspace(-0.1, 0.1, 50)
# xx, yy = np.meshgrid(x_lin, y_lin)
# zz = f(x_lin, y_lin)
#
# end = time.time()
# print(end - start)
#
#
# fig1 = plt.figure()
# ax1 = fig1.add_subplot(111, projection='3d')
# ax1.plot_trisurf(x, y, z, color='b')
#
# fig2 = plt.figure()
# ax2 = fig2.add_subplot(111, projection='3d')
# ax2.plot_surface(xx, yy, zz, color='r')
#
# plt.show()
#
