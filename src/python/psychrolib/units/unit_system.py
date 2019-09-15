from enum import Enum, auto


class MeasurementUnit(Enum):
    SI = auto()
    IP = auto()


class SI:
    unit = MeasurementUnit.SI

class IP:
    unit = MeasurementUnit.IP
