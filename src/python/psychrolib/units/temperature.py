from .unit_system import SI, IP
from .generic_unit import GenericUnit


class Temperature(GenericUnit):
    ZERO_CELSIUS_IN_KELVIN = 273.15
    ZERO_FAHRENHEIT_IN_RANKINE = 491.67
    FAHRENHEIT_TO_CELSIUS_RATIO = 1.8

    absolute_measure = True

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

    k = kelvin
    r = rankine
    c = celsius
    f = fahrenheit

    def __sub__(self, other):
        return DeltaTemperature(self.kelvin - other.kelvin)


class DeltaTemperature:
    FAHRENHEIT_TO_CELSIUS_RATIO = 1.8

    def __init__(self, value):
        self._celsius = value

    @property
    def fahrenheit(self):
        return self.FAHRENHEIT_TO_CELSIUS_RATIO * self._celsius

    @property
    def celsius(self):
        return self._celsius

    def __repr__(self):
        return f"Delta ... {self._celsius} C"

    c = celsius
    f = fahrenheit




class Kelvin(Temperature, SI):
    def __repr__(self):
        return f"{self.kelvin} K"

class Celsius(Kelvin):
    def __init__(self, value):
        super().__init__(value + self.ZERO_CELSIUS_IN_KELVIN)

    def __repr__(self):
        return f"{self.celsius} °C"

class Rankine(Temperature, IP):
    def __init__(self, value):
        super().__init__(value / self.FAHRENHEIT_TO_CELSIUS_RATIO)
    def __repr__(self):
        return f"{self.rankine} R"

class Fahrenheit(Rankine):
    def __init__(self, value):
        super().__init__(value + self.ZERO_FAHRENHEIT_IN_RANKINE)
    def __repr__(self):
        return f"{self.fahrenheit} °F"