from raysect.optical.material.emitter import HomogeneousVolumeEmitter


class Vacuum(HomogeneousVolumeEmitter):

    def emission_function(self, direction, spectrum, world, ray, primitive, to_local, to_world):
        spectrum.samples += 0.0
        return spectrum