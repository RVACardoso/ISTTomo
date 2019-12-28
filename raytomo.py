import numpy as np
import matplotlib.pyplot as plt
from math import degrees, radians, asin

from raysect.optical import World, translate, rotate, Point3D, d65_white
from raysect.optical.material import Lambert, UniformSurfaceEmitter, AbsorbingSurface
from raysect.optical.library import RoughNickel
from raysect.optical.library.spectra.colours import red, green, blue, yellow, purple, orange, cyan, red_orange
from raysect.optical.observer import PinholeCamera, RGBPipeline2D, Pixel, PowerPipeline0D, TargettedPixel
from raysect.primitive import Box, Sphere
from raysect.primitive.mesh import import_stl

from vacuum import Vacuum

from raysect.core import SerialEngine

class RayTomo:

    def __init__(self):
        world = World()
        self.world = world
        self.top_px = []
        self.out_px = []
        self.top_power = []
        self.out_power = []

        self.vessel_width = 0.060
        self.vessel_in_rad = 0.100
        self.vessel_out_rad = 0.101

        self.lid_height = 0.0016
        # self.cam_out_radius = 0.0189
        self.cam_in_radius = 0.0184
        self.tube_height = 0.0401
        # self.lid_top = 0.00110       # como parece certo pela folha
        # self.lid_outer = 0.00155
        #
        # # self.lid_top = 0.00155     # como estava inicialmente
        # # self.lid_outer = 0.00110

        self.lid_top = 0.00155
        self.lid_outer = 0.00155

        # Top calib 2
        self.x_shift_top = -0.005
        self.y_shift_top = 0.09685
        self.x_shift_outer = -0.1091

        self.top_twist = 0.000 + 0.0003
        self.top_px_first_x = self.x_shift_top + 0.0001 + 0.000375 + 7 * 0.00095 + 0.0
        self.top_px_first_y = self.y_shift_top + 0.009 + self.top_twist - 0.0003
        self.top_px_z = self.vessel_width / 2

        self.out_twist = 0.0007 - 0.0003
        self.out_px_first_y = 0.0001 + 0.000375 + 7 * 0.00095 - 0.00015 + 0.0001
        self.out_px_first_x = self.x_shift_outer - 0.009 + self.out_twist + 0.00035 - 0.0005
        self.out_px_z = self.vessel_width / 2

        self.top_data = []
        self.out_data = []

        self.real_data = []

    def check_scene(self, max_iter=200):

        self.vessel.material = Lambert(blue)
        self.camera_outer.material = Lambert(yellow)
        self.camera_top.material = Lambert(yellow)
        self.source.material = Lambert(green)
        self.top_pinhole.material = Lambert(green)
        self.out_pinhole.material = Lambert(green)


        # cube walls
        bottom = Box(lower=Point3D(-0.99, -1.02, -0.99), upper=Point3D(0.99, -1.01, 0.99), parent=self.world,
                      material=Lambert(red))
        # top = Box(lower=Point3D(-0.99, 1.01, -0.99), upper=Point3D(0.99, 1.02, 0.99), parent=self.world,
        #           material=Lambert(red))
        left = Box(lower=Point3D(1.01, -0.99, -0.99), upper=Point3D(1.02, 0.99, 0.99), parent=self.world,
                    material=Lambert(yellow))
        # right = Box(lower=Point3D(-1.02, -0.99, -0.99), upper=Point3D(-1.01, 0.99, 0.99), parent=self.world,
        #             material=Lambert(purple))
        back = Box(lower=Point3D(-0.99, -0.99, 1.01), upper=Point3D(0.99, 0.99, 1.02), parent=self.world,
                   material=Lambert(orange))

        # various wall light sources
        light_front = Box(lower=Point3D(-1.5, -1.5, -10.1), upper=Point3D(1.5, 1.5, -10), parent=self.world,
                          material=UniformSurfaceEmitter(d65_white, 1.0))
        light_top = Box(lower=Point3D(-0.99, 1.01, -0.99), upper=Point3D(0.99, 1.02, 0.99), parent=self.world,
                        material=UniformSurfaceEmitter(d65_white, 1.0), transform=translate(0, 1.0, 0))

        light_bottom = Box(lower=Point3D(-0.99, -3.02, -0.99), upper=Point3D(0.99, -3.01, 0.99), parent=self.world,
                        material=UniformSurfaceEmitter(d65_white, 1.0), transform=translate(0, 1.0, 0))

        light_right = Box(lower=Point3D(-1.92, -0.99, -0.99), upper=Point3D(-1.91, 0.99, 0.99), parent=self.world,
                          material=UniformSurfaceEmitter(d65_white, 1.0))

        light_left = Box(lower=Point3D(1.91, -0.99, -0.99), upper=Point3D(1.92, 0.99, 0.99), parent=self.world,
                          material=UniformSurfaceEmitter(d65_white, 1.0))

        # Process the ray-traced spectra with the RGB pipeline.
        rgb = RGBPipeline2D()

        # camera
        pix = 1000
        camera = PinholeCamera((pix, pix), pipelines=[rgb], transform=translate(-0.01, 0.0, -0.25) * rotate(0, 0, 0))
        # camera = PinholeCamera((pix, pix), pipelines=[rgb], transform=translate(0.0, 0.0, 0.4) * rotate(180, 0, 0))
        # top view
        # camera = PinholeCamera((pix, pix), pipelines=[rgb], transform=translate(0.0, self.vessel_out_rad+0.15, self.vessel_width/2)*rotate(0, -90, 0))
        # prof
        camera = PinholeCamera((pix, pix), pipelines=[rgb], transform=translate(-0.13, 0.13, -0.2) * rotate(-25, -25.0, 0))

        # camera top side
        # camera = PinholeCamera((pix, pix), pipelines=[rgb], transform=translate(self.x_shift_top, self.top_px_first_y+0.0004, self.top_px_z-self.cam_in_radius+0.005)*rotate(0, 0, 0))
        # camera top down-up
        # camera = PinholeCamera((pix, pix), pipelines=[rgb], transform=translate(self.x_shift_top, self.top_px_first_y-0.01, self.vessel_width/2)*rotate(0, 90, 0))
        # camera top up-down
        # camera = PinholeCamera((pix, pix), pipelines=[rgb], transform=translate(self.x_shift_top-0.004, self.top_px_first_y+self.lid_top+self.tube_height-0.01, self.vessel_width/2)*rotate(0, -90, 0))

        # camera out side
        # camera = PinholeCamera((pix, pix), pipelines=[rgb], transform=translate(-self.vessel_out_rad-0.015, 0.000, self.vessel_width/2-self.cam_in_radius/2+0.0001))
        # camera out down-up
        # camera = PinholeCamera((pix, pix), pipelines=[rgb], transform=translate(self.out_px_first_x+0.005+0.005, 0.0, self.vessel_width/2)*rotate(90, 0, 0))
        # camera out up-down
        # camera = PinholeCamera((pix, pix), pipelines=[rgb], transform=translate(-self.vessel_out_rad-self.tube_height-0.01, 0.0, self.vessel_width/2-0.005)*rotate(-90, 0, 0))

        # camera - pixel sampling settings
        camera.fov = 60  # 45
        camera.pixel_samples = 10

        # camera - ray sampling settings
        camera.spectral_rays = 1
        camera.spectral_bins = 25
        camera.parent = self.world

        plt.ion()
        p = 1
        while not camera.render_complete:
            print("Rendering pass {}...".format(p))
            camera.observe()
            print()
            p += 1
            if p > max_iter:
                break

        plt.ioff()
        rgb.display()

    def add_isttok(self, pixel_samples=1000000):
        nickel_roughness = 0.23
        min_wl = 50
        max_wl = 51

        # add vessel and cameras
        self.vessel = import_stl("../isttok_3d/vessel5_stl.stl", scaling=1, mode='binary', parent=self.world,
                                 material=RoughNickel(nickel_roughness),  #AbsorbingSurface(),
                                 transform=translate(0, 0, 0) * rotate(0, 0, 0))

        self.camera_top = import_stl("../isttok_3d/camera_top3_stl.stl", scaling=1, mode='binary', parent=self.world,
                                     material=AbsorbingSurface(),
                                     transform=translate(self.x_shift_top,
                                                         self.y_shift_top + self.tube_height + self.lid_top,
                                                         self.vessel_width / 2.0) * rotate(0, -90, 0))

        # self.camera_outer = import_stl("../isttok_3d/camera_outer3_stl.stl", scaling=1, mode='binary', parent=self.world,
        #                                material=AbsorbingSurface(),
        #                                transform=translate(self.x_shift_outer - self.tube_height - self.lid_outer,
        #                                                    0.0,
        #                                                    self.vessel_width / 2.0) * rotate(-90, 0, 0))

        self.camera_outer = import_stl("../isttok_3d/camera_outer4_newpin_stl.stl", scaling=1, mode='binary', parent=self.world,
                                       material=AbsorbingSurface(),
                                       transform=translate(self.x_shift_outer - self.tube_height - self.lid_outer,
                                                           0.0,
                                                           self.vessel_width / 2.0) * rotate(-90, 0, 0))

        pinhole_sphere_radius = 0.0005 # 0.0005
        self.top_pinhole = Sphere(radius=pinhole_sphere_radius, parent=self.world,
                                  transform=translate(self.x_shift_top, self.y_shift_top, self.vessel_width/2),
                                  material=Vacuum())
        pinhole_sphere_radius = 0.00035  # 0.0005
        self.out_pinhole = Sphere(radius=pinhole_sphere_radius, parent=self.world,
                                  transform=translate(self.x_shift_outer, 0.0, self.vessel_width / 2),
                                  material=Vacuum())

        for i in range(16):
            self.top_power.append(PowerPipeline0D(accumulate=False))
            self.out_power.append(PowerPipeline0D(accumulate=False))

            top_px_x = self.top_px_first_x - i * 0.00095
            top_px_y = self.top_px_first_y - i * 2 * (self.top_twist / 15)
            top_angle = degrees(asin(2 * self.top_twist / 0.01425))

            out_px_y = self.out_px_first_y - i * 0.00095
            out_px_x = self.out_px_first_x - i * 2 * (self.out_twist / 15)
            out_angle = -degrees(asin(2 * self.out_twist / 0.01425))

            self.top_px.append(TargettedPixel(targets=[self.top_pinhole], targetted_path_prob=1.0, pipelines=[self.top_power[i]],
                                     x_width=0.00075, y_width=0.00405,
                                     min_wavelength=min_wl, max_wavelength=max_wl,
                                     spectral_bins=1, pixel_samples=pixel_samples, parent=self.world, quiet=True,
                                     ray_importance_sampling=True, ray_important_path_weight=0.05, ray_max_depth=50,
                                     transform=translate(top_px_x, top_px_y, self.top_px_z) *
                                               rotate(0, 0, top_angle) * rotate(0, -90, 0)))

            self.out_px.append(TargettedPixel(targets=[self.out_pinhole], targetted_path_prob=1.0, pipelines=[self.out_power[i]],
                                     x_width=0.00075, y_width=0.00405,
                                     min_wavelength=min_wl, max_wavelength=max_wl,
                                     spectral_bins=1, pixel_samples=pixel_samples, parent=self.world, quiet=True,
                                     ray_importance_sampling=True, ray_important_path_weight=0.05, ray_max_depth=50,
                                     transform=translate(out_px_x, out_px_y, self.out_px_z) *
                                               rotate(0, 0, out_angle) * rotate(-90, 0, 90)))

            # self.top_px.append(Box(Point3D(-0.000375, -0.002025, 0.0), Point3D(0.000375, 0.002025, 0.0005), parent=self.world, material=UniformSurfaceEmitter(green, scale=0.01),
            # transform=translate(top_px_x, top_px_y, self.top_px_z)*rotate(0, 0, top_angle)*rotate(0, -90, 0)))
            #
            # self.out_px.append(Box(Point3D(-0.000375, -0.002025, 0.0), Point3D(0.000375, 0.002025, 0.0005), parent=self.world, material=UniformSurfaceEmitter(green, scale=0.01),
            # transform=translate(out_px_x, out_px_y, self.out_px_z)*rotate(0, 0, out_angle)*rotate(-90, 0, 90)))

    def place_source(self):
        # abstract method
        raise NotImplementedError

    def simulate_top_rays(self, pixels=range(16)):
        top_data = []
        for i in pixels:
            self.top_px[i].observe()
            # print("{}, {} -> {}".format(self.top_power[i].value.mean,
            #                                                  self.top_power[i].value.error(), self.top_power[i].value.error()*100/self.top_power[i].value.mean))
            top_data.append(self.top_power[i].value.mean)

        self.top_data_sum = sum(top_data)
        self.top_data = top_data

    def simulate_out_rays(self, pixels=range(16)):
        out_data = []
        for i in pixels:
            self.out_px[i].observe()
            # print("{}, {} -> {}".format(self.out_power[i].value.mean,
            #                                                  self.out_power[i].value.error(), self.out_power[i].value.error()*100/self.out_power[i].value.mean))
            out_data.append(self.out_power[i].value.mean)

        self.out_data_sum = sum(out_data)
        self.out_data = out_data

    def plot_data(self, fig_name):

        import matplotlib
        font = {'family': 'normal',
                'size': 32}
        matplotlib.rc('font', **font)

        plt.figure()
        # plt.subplot(211)
        # angles = np.linspace(0, 360, 360)
        # vessel_x = self.vessel_in_rad * np.cos(angles)
        # vessel_y = self.vessel_in_rad * np.sin(angles)
        # plt.scatter(x=vessel_x, y=vessel_y)
        # try:
        #     plt.scatter(x=self.lamp_radius * np.cos(radians(self.lamp_angle)),
        #             y=self.lamp_radius * np.sin(radians(self.lamp_angle)),
        #             c="r", s=100)
        # except AttributeError:
        #     pass
        # plt.scatter(x=self.vessel_in_rad, y=0, c='g', s=50)
        # plt.scatter(x=0.005, y=self.vessel_in_rad - 0.0005, c='g', s=50)
        # plt.gca().set_aspect('equal', adjustable='box')

        cam_nr = np.arange(16)

        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.5)

        plt.subplot(211)
        plt.bar(cam_nr, self.top_data)
        plt.xticks(cam_nr, cam_nr)
        plt.ylabel("Power (W)")
        # plt.ylabel("Fraction of TOP total power (P = " + str(self.top_data_sum) + ")")
        plt.xlabel("Detector number")

        plt.subplot(212)
        # plt.barh(cam_nr, np.flip(self.out_data))
        plt.bar(cam_nr, self.out_data)
        plt.xticks(cam_nr, cam_nr, rotation='vertical')
        plt.yticks(rotation='vertical')
        plt.ylabel("Power (W)")
        # plt.ylabel("Fraction of OUTER total power (P = " + str(self.out_data_sum) + ")")
        plt.xlabel("Detector number")

        # plt.suptitle("SIMULATED DATA - " + fig_name)
        # plt.savefig("./withlegs_test_sim/"+fig_name+".png")

        # cam_nr = np.arange(16)
        #
        font_size = 15

        plt.rc('xtick', labelsize=font_size)
        plt.rc('ytick', labelsize=font_size)

        # plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        #
        # plt.bar(cam_nr, self.out_data)
        # plt.xticks(cam_nr, cam_nr)
        # plt.xlabel("Pixel Number", fontsize=font_size)
        # plt.ylabel("Power received (W)", fontsize=font_size)
        # plt.text(2.0, 0.85e-5, 'OUTER camera', horizontalalignment='center', fontsize=25)


    def create_save(self):
        # abstract method
        raise NotImplementedError

    def load_data_to_save(self):
        # abstract method
        raise NotImplementedError
