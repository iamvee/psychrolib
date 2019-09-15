from .unit_system import SI, IP
from .generic_unit import GenericUnit


class Pressure(GenericUnit):
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
    pass


class Atmosphere(Pascal):
    pass


