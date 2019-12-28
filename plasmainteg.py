# from raysect.optical cimport new_point3d
# from libc.math cimport floor
# cimport cython

from raysect.optical.material.emitter import VolumeIntegrator

from math import floor, sqrt
from raysect.core.math.point import Point3D
import numpy as np


class MyTestNumIntegrator(VolumeIntegrator):

    def __init__(self, step, min_samples=5):
        self._step = step
        self._min_samples = min_samples

    @property
    def step(self):
        return self._step

    @step.setter
    def step(self, value):
        if value <= 0:
            raise ValueError("Numerical integration step size can not be less than or equal to zero")
        self._step = value

    @property
    def min_samples(self):
        return self._min_samples

    @min_samples.setter
    def min_samples(self, value):
        if value < 2:
            raise ValueError("At least two samples are required to perform the numerical integration.")
        self._min_samples = value

    def integrate_old(self, spectrum, world, ray, primitive, material, start_point,end_point,
                             world_to_primitive, primitive_to_world):

        # convert start and end points to local space
        start = start_point.transform(world_to_primitive)
        end = end_point.transform(world_to_primitive)

        # obtain local space ray direction and integration length
        integration_direction = start.vector_to(end)
        length = sqrt(integration_direction.length)   # length = integration_direction.get_length()

        # nothing to contribute?
        if length == 0.0:
            return spectrum

        integration_direction = integration_direction.normalise()
        ray_direction = -integration_direction #.neg()

        # calculate number of complete intervals (samples - 1)
        intervals = max(self._min_samples - 1, int( floor(length / self._step)))

        # adjust (increase) step size to absorb any remainder and maintain equal interval spacing
        step = length / intervals

        # create working buffers
        emission = ray.new_spectrum()
        emission_previous = ray.new_spectrum()

        # sample point and sanity check as bounds checking is disabled
        emission_previous = material.emission_function(start, ray_direction, emission_previous, world, ray, primitive, world_to_primitive, primitive_to_world)
        self._check_dimensions(emission_previous, spectrum.bins)

        # numerical integration
        c = 0.5 * step
        for interval in range(1, intervals):

            # calculate location of new sample point
            t = interval * step
            sample_point = Point3D(
                start.x + t * integration_direction.x,
                start.y + t * integration_direction.y,
                start.z + t * integration_direction.z
            )

            # sample point and sanity check as bounds checking is disabled
            emission = material.emission_function(sample_point, ray_direction, emission, world, ray, primitive, world_to_primitive, primitive_to_world)
            self._check_dimensions(emission, spectrum.bins)

            # trapezium rule integration
            # print(spectrum.bins)
            for index in range(spectrum.bins):
                # print(type(spectrum))
                spectrum.samples[index] += c * (emission.samples[index] + emission_previous.samples[index])

            # swap buffers and clear the active buffer
            temp = emission_previous
            emission_previous = emission
            emission = temp
            emission.clear()

        return spectrum

    def integrate(self, spectrum, world, ray,primitive, material, start_point,end_point,
                             world_to_primitive, primitive_to_world):

        # convert start and end points to local space
        start = start_point.transform(world_to_primitive)
        end = end_point.transform(world_to_primitive)

        # obtain local space ray direction and integration length
        integration_direction = start.vector_to(end)
        length = sqrt(integration_direction.length)   # length = integration_direction.get_length()

        # nothing to contribute?
        if length == 0.0:
            return spectrum

        integration_direction = integration_direction.normalise()
        ray_direction = -integration_direction #.neg()

        # calculate number of complete intervals (samples - 1)
        intervals = max(self._min_samples - 1, int(floor(length / self._step)))

        # adjust (increase) step size to absorb any remainder and maintain equal interval spacing
        step = length / intervals

        # create working buffers
        emission = ray.new_spectrum()

        sample_points_x = [start.x + interval*step * integration_direction.x for interval in range(intervals)]
        sample_points_y = [start.y + interval*step * integration_direction.y for interval in range(intervals)]
        sample_points = np.array([sample_points_x, sample_points_y])

        emission = material.emission_function(sample_points, ray_direction, emission, world, ray, primitive,
                                              world_to_primitive, primitive_to_world)

        pair_sum = (emission[:-1] + np.array(emission, copy=True)[1:])
        spectrum.samples[0] += np.sum(0.5 * step * pair_sum)

        return spectrum

    def _check_dimensions(self, spectrum, bins):
        if spectrum.samples.ndim != 1 or spectrum.samples.shape[0] != bins:
            raise ValueError("Spectrum returned by emission function has the wrong number of samples.")