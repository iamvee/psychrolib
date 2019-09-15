from .unit_system import SI, IP
from .generic_unit import GenericUnit


class Temperature(GenericUnit):
    ZERO_CELSIUS_IN_KELVIN = 273.15
    ZERO_FAHRENHEIT_IN_RANKINE = 491.67
    FAHRENHEIT_TO_CELSIUS_RATIO = 1.8

    @property
    def kelvin(self):
        return self._base_value

    @property
    def celsius(self):
        return self.kelvin - self.ZERO_CELSIUS_IN_KELVIN

    @property
    def fahrenheit(self):
        return self.FAHRENHEIT_TO_CELSIUS_RATIO * self.celsius + 32

    @property
    def rankine(self):
        return self.FAHRENHEIT_TO_CELSIUS_RATIO * self.kelvin


class Kelvin(Temperature, SI):
    pass


class Celsius(Kelvin):
    def __init__(self, value):
        super().__init__(value + self.ZERO_CELSIUS_IN_KELVIN)


class Rankine(Temperature, IP):
    def __init__(self, value):
        super().__init__(value / self.FAHRENHEIT_TO_CELSIUS_RATIO)


class Fahrenheit(Rankine):
    def __init__(self, value):
        super().__init__(value + self.ZERO_FAHRENHEIT_IN_RANKINE)

