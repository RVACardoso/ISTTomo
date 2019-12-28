from reallamp import RealLamp
from realplasma import  RealPlasma
from raylamp import RayLamp
from rayplasma import RayPlasma
from plasmaemit import SinglePlasmaEmit, NpySinglePlasmaEmit, GenericSinglePlasma
from geomatrix import GeoMatrix

import matplotlib.pyplot as plt
import numpy as np

import time


# geomat = GeoMatrix(matrix_file="geo_test50-px25000.npy")
# # geomat_dir = GeoMatrix(matrix_file="geo_test50-px25000-noreflect.npy")
# # geomat_ref = GeoMatrix(matrix_data=(geomat.geo_matrix - geomat_dir.geo_matrix))
#
# geomat = GeoMatrix(matrix_file="geo_mat60-newpin-ref-px1000.npy")
#
# # for i in range(32):
# #     plt.subplot(6,6,i+1)
# #     geomat.plot_matrix(detector_nr=i)
# #
# # 4 and 30
# # geomat.plot_matrix(detector_nr=4)
# # plt.title("Viewing cone for top detector #5")
# geomat.plot_matrix(detector_nr=5)
# plt.title("Viewing cone for top detector #5")
# plt.show()

# contour_type = 'test1'
# stacking_type = 'gaussian'
# #
# sinth_plasma = GenericSinglePlasma(contour_type=contour_type, stacking_type=stacking_type, grid_side=60)
# # # loop across center, width and power gain
# center = [[0.05, 0.02]]
# width = [[(1.0, 2.0), (1.0, 1.0)]]
# power_gain = 5000
# sinth_plasma.change_profile(center=center, width=width, power_gain=power_gain)
#
# start = time.time()
# geomat.get_tomogram_data(pixel_div=5, sinth_plasma=sinth_plasma)
# geomat.get_matrix_sensors()
#
# middle = time.time()
# print(middle-start)
#
# import matplotlib
# font = {'family' : 'normal',
#         'size'   : 23}
# matplotlib.rc('font', **font)
#
# # geomat.plot_data(fig_name='geomat fig name')
# # geomat.plot_tomogram()
#
# # plt.rcParams.update({'font.size': 18})
# cam_nr = np.arange(16)
# plt.figure()
# plt.subplot(121)
# plt.bar(cam_nr, geomat.top_data)
# plt.xticks(cam_nr, cam_nr)
# plt.ylabel("Power collected (W)")
# plt.text(3.0, 0.000004, 'Top camera', horizontalalignment='center', fontsize=25)
# # plt.xlabel("Pixel number")
# plt.title("Ray-tracing simulation")
# plt.xlabel("Detector number")
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#
# start2 = time.time()
#
# ray_plasma = RayPlasma(coord='cartesian', pixel_nr=60, count=1)
# ray_plasma.add_isttok(pixel_samples=30000)
# ray_plasma.place_source()
# sinth_plasma1 = GenericSinglePlasma(contour_type=contour_type, stacking_type=stacking_type, grid_side=60)
# sinth_plasma1.change_profile(center=center, width=width, power_gain=power_gain)
# ray_plasma.source.material = sinth_plasma1
# ray_plasma.simulate_top_rays()
# ray_plasma.simulate_out_rays()
# # ray_plasma.plot_data(fig_name='ray fig name')
# # ray_plasma.get_tomogram_data(pixel_div=10, sinth_plasma=sinth_plasma1)
# # ray_plasma.plot_tomogram()
# print(time.time()-start2)
#
# plt.subplot(122)
# plt.bar(cam_nr, ray_plasma.top_data)
# plt.xticks(cam_nr, cam_nr)
# plt.ylabel("Power colected (W)")
# plt.text(3.0, 0.000004, 'Top camera', horizontalalignment='center', fontsize=25)
# plt.xlabel("Detector number")
# plt.title("Matrix product")
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
# plt.suptitle("Comparison between ray-tracing and matrix product")
#
# print(np.abs(geomat.top_data-ray_plasma.top_data)/np.array(ray_plasma.top_data))
#
# plt.show()



#######################################################################################################

# All working using GenericSinglePlasma
# # RayPlasma.csv_file_name = 'csv_test1.csv'
# ray_plasma = RayPlasma(coord='cartesian', pixel_nr=60, count=1)
# # ray_plasma.create_save()
# ray_plasma.add_isttok(pixel_samples=10000)
# ray_plasma.place_source()
# # ray_plasma.check_scene()
# contour_type = 'test1'
# stacking_type = 'gaussian'
# sinth_plasma = GenericSinglePlasma(contour_type=contour_type, stacking_type=stacking_type, grid_side=100)
# # # loop across center, width and power gain
# # center = [[-0.03, -0.01]]
# # width = [[(1.0, 1.0), (1.0, 1.0)]]
# center = [[0.03, -0.03]]                                ## compare LoS
# # center = [[-0.025, -0.02]]
# # center = [[-0.025, 0.02]]
# width = [[(0.75, 0.75), (0.75, 0.75)]]                  ## compare LoS
#
# import matplotlib
# font = {'family' : 'normal',
#         'size'   : 20}
# matplotlib.rc('font', **font)
# #
# power_gain = 5.1
# sinth_plasma.change_profile(center=center, width=width, power_gain=power_gain)
# ray_plasma.source.material = sinth_plasma
#
# # from raysect.optical.material import UnityVolumeEmitter
# # ray_plasma.source.material = UnityVolumeEmitter()
#
# ray_plasma.simulate_top_rays()
# ray_plasma.simulate_out_rays()
# # ray_plasma.get_tomogram_data(pixel_div=10, sinth_plasma=sinth_plasma)
# ray_plasma.plot_compare_los()
# # ray_plasma.load_data_to_save()
# # ray_plasma.plot_data(fig_name='fig name')
# # print(ray_plasma.top_data)
# # print(ray_plasma.out_data)
# # ray_plasma.plot_tomogram()
# plt.show()

# old_pin_out = [1.2895686881252291e-09, 1.6799814074055568e-09, 2.0854084940942626e-09, 2.5058108684345996e-09, 2.958729093841035e-09, 3.355633069934302e-09, 3.5340404593937113e-09, 3.623638433393569e-09, 3.578610936080435e-09, 3.3530634848887927e-09, 3.0339896560323847e-09, 2.682312993932888e-09, 2.2093599751347877e-09, 1.83887453188438e-09, 1.4987782944877209e-09, 1.1848801952232872e-09]
# new_pin_out = [4.927117874284028e-10, 6.658060155723725e-10, 8.158440717833898e-10, 9.825189441004487e-10, 1.1551734787359194e-09, 1.2549745299977083e-09, 1.398434611383289e-09, 1.406995752863905e-09, 1.3869940604973499e-09, 1.3121613883662494e-09, 1.1823101594293927e-09, 1.047703852384906e-09, 8.773084732187701e-10, 7.301019033171488e-10, 5.710056256097953e-10, 4.620470983649682e-10]
#
# print(list(np.array(old_pin_out)/np.array(new_pin_out)))

# plt.rcParams.update({'font.size': 20})

# cam_nr = np.arange(16)
# plt.figure()
# plt.subplot(221)
# plt.bar(cam_nr, ray_plasma.top_data)
# plt.xticks(cam_nr, cam_nr)
# plt.ylabel("Power collected (W)")
# plt.text(12.0, 0.000005, 'Top camera', horizontalalignment='center')
# # plt.xlabel("Pixel number")
# plt.title("Max. emissivity: 5 $kW/m^{-3}$")
# plt.ylim(bottom=0, top=6.5e-6)
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
# top5 = ray_plasma.top_data
#
# plt.subplot(223)
# plt.bar(cam_nr, ray_plasma.out_data)
# plt.xticks(cam_nr, cam_nr)
# plt.ylabel("Power collected (W)")
# plt.ylim(bottom=0, top=2.6e-6)
# plt.text(12.0, 0.000002, 'Outer camera', horizontalalignment='center')
# plt.xlabel("Detector number")
# out5 = ray_plasma.out_data
#
#
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#
# power_gain = 500000
# sinth_plasma.change_profile(center=center, width=width, power_gain=power_gain)
# ray_plasma.source.material = sinth_plasma
# ray_plasma.simulate_top_rays()
# ray_plasma.simulate_out_rays()
# ray_plasma.get_tomogram_data(pixel_div=10, sinth_plasma=sinth_plasma)
#
# plt.subplot(222)
# plt.bar(cam_nr, ray_plasma.top_data, color='r')
# plt.xticks(cam_nr, cam_nr)
# plt.ylabel("Power collected (W)")
# plt.text(12.0, 0.0005, 'Top camera', horizontalalignment='center')
# # plt.xlabel("Pixel number")
# plt.title("Max. emissivity: 500 $kW/m^{-3}$")
# plt.ylim(bottom=0, top=6.5e-4)
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#
# plt.subplot(224)
# plt.bar(cam_nr, ray_plasma.out_data, color='r')
# plt.xticks(cam_nr, cam_nr)
# plt.ylabel("Power collected (W)")
# plt.text(12.0, 0.0002, 'Outer camera', horizontalalignment='center')
# plt.xlabel("Detector number")
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
# plt.ylim(bottom=0, top=2.6e-4)
# plt.suptitle("Linearity of sensor measurements with profile emissivity")
#
# print(np.array(ray_plasma.top_data)/np.array(top5))
# print(np.array(ray_plasma.out_data)/np.array(out5))
# plt.show()

# All working using SinglePlasmaEmit and NpyPlasmaEmit!
# RayPlasma.csv_file_name = 'csv_test1.csv'
# ray_plasma = RayPlasma(coord='cartesian', pixel_nr=50, count=1)
# ray_plasma.create_save()
# ray_plasma.add_isttok(pixel_samples=100)
# ray_plasma.place_source()
# # ray_plasma.check_scene()
# shape = "gaussian"
# center = [[-0.025, 0.0]]
# width = [[(0.07, 0.03), (0.035, 0.025)]]
# flatlen = 0.5
# power_gain = 5100
# start = time.time()
# sinth_plasma = SinglePlasmaEmit(shape=shape, center=center, width=width, flatlen=flatlen, power_gain=power_gain)
# ray_plasma.source.material = sinth_plasma
# ray_plasma.simulate_top_rays()
# ray_plasma.simulate_out_rays()
# np_sinth_plasma = NpySinglePlasmaEmit(shape=shape, center=center, width=width, flatlen=flatlen, power_gain=power_gain)
# ray_plasma.get_tomogram_data(pixel_div=10, sinth_plasma=np_sinth_plasma)
# ray_plasma.load_data_to_save()
# end = time.time()
# print(end - start)
# ray_plasma.plot_data(fig_name='fig name')
# ray_plasma.plot_tomogram()
# plt.show()

# # RealLamp working!
# real_lamp = RealLamp()
# real_lamp.get_data(lamp_radius=0.05, lamp_angle=210.0)    # reflections on ECPD article : lamp_radius=0.05, lamp_angle=225.0
# # real_lamp.convert_to_watt()
# # real_lamp.top_pixels = real_lamp.top_pixels/sum(real_lamp.top_pixels)
# # real_lamp.out_pixels = real_lamp.out_pixels/sum(real_lamp.out_pixels)
# real_lamp.plot_data()
# real_top, real_out = real_lamp.provide_data()
# print(real_top)
# print()
# print(real_out)
# plt.show()

# plt.rcParams.update({'font.size': 25})
# import matplotlib
# font = {'family' : 'normal',
#         'size'   : 22}
# matplotlib.rc('font', **font)
#
# cam_nr = np.arange(16)
#
# plt.figure()
# plt.subplot(211)
# plt.bar(cam_nr, real_lamp.top_pixels)
# plt.xticks(cam_nr, cam_nr)
# plt.ylabel("Fraction of power")
# plt.text(2.0, 0.7, 'Top camera', horizontalalignment='center')
# # plt.xlabel("Pixel number")
# plt.xlabel("Detector number")
#
#
# plt.subplot(212)
# sum_top_pixels = sum(real_lamp.out_pixels)
# plt.bar(cam_nr, real_lamp.out_pixels)
# plt.xticks(cam_nr, cam_nr)
# plt.ylabel("Fraction of power")
# plt.text(2.0, 0.4, 'Outer camera', horizontalalignment='center')
# plt.xlabel("Detector number")
#
# plt.suptitle("Fraction of total collected power")
# plt.subplots_adjust(hspace=0.3)
# plt.show()

# # RealPlasma working!
real_plasma = RealPlasma()
real_plasma.get_data(shotnr=47258)  # 46100 , 47181, 47706
# real_plasma.get_instant(time=304000)
real_plasma.get_window_avg(start=300000, end=305000)
real_plasma.new_plot_data()
real_plasma.plot_orig_signals()
# real_top, real_out = real_plasma.provide_data()
# print(real_top)
# print()
# print(real_out)
plt.show()


#
# # print(real_plasma.top_signals.shape)
# # print(real_plasma.top_time)
# top_signals = np.array(real_plasma.top_signals)
# out_signals = np.array(real_plasma.out_signals)
# top_time = np.array(real_plasma.top_time[0])
# out_time = np.array(real_plasma.out_time[0])
#
# dc_index = np.where(top_time > 2000)[0][0]
# # print(dc_index)
# top_dc, outer_dc = [], []
# for i in range(16):
#     top_dc.append(np.mean(top_signals[i][:dc_index]))
#     outer_dc.append(np.mean(out_signals[i][:dc_index]))
# top_dc = np.array(top_dc)
# outer_dc = np.array(outer_dc)
#
# top_idx = ( (top_time>280000) & (top_time<310000) )
# # print(top_idx)
# # print(top_time[top_idx])  #.shape)
# # print(top_signals[:,top_idx].shape)
# # print(top_signals[:,top_idx])
# # print()
# out_idx = ( (out_time>280000) & (out_time<310000) )
# # print(out_idx)
# # print(out_time[out_idx])  #.shape)
# # print(out_signals[:, out_idx].shape)
# # print(out_signals[:, out_idx])
#
# # print(np.where(top_time > 380000 and top_time < 390000))
# # print()
#
# # print(top_signals.shape)
# # print(top_time.shape)
# top_time = top_time[top_idx]
# top_signals = top_signals[:, top_idx]
# out_signals = out_signals[:, out_idx]
#
# top_signals = top_signals.transpose() - top_dc
# out_signals = out_signals.transpose() - outer_dc
#
# # for i in range(16):
# #     top_signals[i] -= top_dc[i]
# #     out_signals[i] -= outer_dc[i]
#
#
# out_data = np.hstack((top_time.reshape((299, 1)), top_signals, out_signals))
# print(out_data.shape)
# print()
# # print(out_data)
# np.save("movie47258", out_data)




# RayLamp working!
# RayLamp.csv_file_name = 'csv_test1.csv'
# ray_lamp = RayLamp()
# ray_lamp.add_isttok(pixel_samples=500)
# # ray_lamp.place_source(radius=0.05, angle=105.0, legs=True)
# ray_lamp.place_source(radius=0.03, angle=45.0, legs=True)
# # ray_lamp.check_scene()
# ray_lamp.simulate_top_rays()
# ray_lamp.simulate_out_rays()
# # ray_lamp.plot_data(fig_name='fig name')
# # ray_lamp.create_csv()
# # ray_lamp.load_data_to_csv()
# top_loss_val = ray_lamp.calculate_top_loss()
# out_loss_val = ray_lamp.calculate_out_loss()
# print(top_loss_val)
# print(out_loss_val)
# ray_lamp.plot_compare()
# plt.show()


# RayLamp: stats for quality of simulation
#
# top_stat = []
# out_stat = []
#
# real_top, real_out = [], []
# sim_top, sim_out = [], []
#
# for angle in range(0,360,45):
#         print(angle)
#         ray_lamp = RayLamp()
#         ray_lamp.add_isttok(pixel_samples=500)
#         ray_lamp.place_source(radius=0.03, angle=angle, legs=True)
#         ray_lamp.simulate_top_rays()
#         ray_lamp.simulate_out_rays()
#         real_top_frac, top_data_frac = ray_lamp.calculate_top_loss()
#         sim_top.append(top_data_frac)
#         real_top.append(real_top_frac)
#
#         real_out_frac, out_data_frac = ray_lamp.calculate_out_loss()
#         sim_out.append(out_data_frac)
#         real_out.append(real_out_frac)
#
#
# for radius in [0.05, 0.07]:
#         for angle in range(0,360,15):
#                 print(angle)
#                 ray_lamp = RayLamp()
#                 ray_lamp.add_isttok(pixel_samples=500)
#                 ray_lamp.place_source(radius=radius, angle=angle, legs=True)
#                 ray_lamp.simulate_top_rays()
#                 ray_lamp.simulate_out_rays()
#                 real_top_frac, top_data_frac = ray_lamp.calculate_top_loss()
#                 sim_top.append(top_data_frac)
#                 real_top.append(real_top_frac)
#
#                 real_out_frac, out_data_frac = ray_lamp.calculate_out_loss()
#                 sim_out.append(out_data_frac)
#                 real_out.append(real_out_frac)
#
#
# # top_stat = np.array(top_stat)
# # out_stat = np.array(out_stat)
#
# sim_top = np.array(sim_top).flatten()
# real_top = np.array(real_top).flatten()
# sim_out = np.array(sim_out).flatten()
# real_out = np.array(real_out).flatten()
#
# top_data = np.vstack((sim_top, real_top))
# out_data = np.vstack((sim_out, real_out))
# # print(top_data.shape)
#
# top_coeff = np.corrcoef(top_data)
# print(top_coeff)
#
# print()
#
# out_coeff = np.corrcoef(out_data)
# print(out_coeff)



# print(sim_top)
# print(sim_top.shape)
# print(real_top.shape)
#
#
# print(top_stat)
# print(top_stat.shape)
# print()
#
#
# print(np.mean(top_stat))
# print(np.std(top_stat))
#
# print()
#
# print(out_stat)
# print(np.mean(out_stat))
# print(np.std(out_stat))

# RayPlasma working!
# # RayPlasma.csv_file_name = 'csv_test1.csv'
# ray_plasma = RayPlasma()
# ray_plasma.add_isttok(pixel_samples=1000000)
# ray_plasma.place_source(shape='lorentzian', count=1, center=[[0.0, 0.0]], width=[0.05])
# # ray_plasma.check_scene()
# ray_plasma.simulate_top_rays()
# ray_plasma.simulate_out_rays()
# ray_plasma.plot_data(fig_name='fig name')
# ray_plasma.get_tomogram_data()
# ray_plasma.plot_tomogram()
# # ray_plasma.create_csv()
# # ray_plasma.load_data_to_csv()
# plt.show()
