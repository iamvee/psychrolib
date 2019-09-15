from .unit_system import SI, IP
from .generic_unit import GenericUnit


class Mass(GenericUnit):
    RATIO_POUNDS_TO_KILOGRAMS = 0.45359237

    @property
    def kg(self):
        return self._base_value

    @property
    def lb(self):
        return self._base_value


class KiloGram(Mass, SI):
    pass


class Pound(Mass, IP):
    def __init__(self, value):
        super().__init__(value * self.RATIO_POUNDS_TO_KILOGRAMS)




