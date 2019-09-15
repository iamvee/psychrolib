from enum import Enum, auto


class MeasurementUnit(Enum):
    SI = auto()
    IP = auto()


class MeasurementUnitClass:
    pass


class SI(MeasurementUnitClass):
    unit = MeasurementUnit.SI


class IP(MeasurementUnitClass):
    unit = MeasurementUnit.IP


SI_SYSTEM = SI()
IP_SYSTEM = IP()
