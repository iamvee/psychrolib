from .unit_system import SI, IP
from .generic_unit import GenericUnit


class Pressure(GenericUnit):
    PASCAL_TO_PSI_RATIO = 1
    PASCAL_TO_ATMOSPHERE_RATIO = 1

    @property
    def pascal(self):
        return self._base_value

    @property
    def atmosphere(self):
        return self._base_value

    @property
    def psi(self):
        return self._base_value

    @property
    def barr(self):
        return self._base_value


class Pascal(Pressure, SI):
    pass


class Psi(Pressure, IP):
    def __init__(self, value):
        super().__init__(value * self.PASCAL_TO_PSI_RATIO)


class Atmosphere(Pascal):
    def __init__(self, value):
        super().__init__(value * self.PASCAL_TO_ATMOSPHERE_RATIO)



