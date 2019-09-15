from enum import Enum, auto
#######################################################################################################
# Helper functions
#######################################################################################################

# Unit system to use.
class UnitSystem(Enum):
    """
    Private class not exposed used to set automatic enumeration values.
    """
    IP = auto()
    SI = auto()


IP = UnitSystem.IP
SI = UnitSystem.SI

PSYCHROLIB_UNITS = None

PSYCHROLIB_TOLERANCE = 1.0


# Tolerance of temperature calculations

def SetUnitSystem(Units: UnitSystem) -> None:
    """
    Set the system of units to use (SI or IP).

    Args:
        Units: string indicating the system of units chosen (SI or IP)

    Notes:
        This function *HAS TO BE CALLED* before the library can be used

    """
    global PSYCHROLIB_UNITS
    global PSYCHROLIB_TOLERANCE

    if not isinstance(Units, UnitSystem):
        raise ValueError("The system of units has to be either SI or IP.")

    PSYCHROLIB_UNITS = Units

    # Define tolerance on temperature calculations
    # The tolerance is the same in IP and SI
    if Units == IP:
        PSYCHROLIB_TOLERANCE = 0.001 * 9. / 5.
    else:
        PSYCHROLIB_TOLERANCE = 0.001


def GetUnitSystem() -> Optional[UnitSystem]:
    """
    Return system of units in use.

    """
    return PSYCHROLIB_UNITS


def isIP() -> bool:
    """
    Check whether the system in use is IP or SI.

    """
    if PSYCHROLIB_UNITS == IP:
        return True
    elif PSYCHROLIB_UNITS == SI:
        return False
    else:
        raise ValueError('The system of units has not been defined.')


#######################################################################################################
# Conversion between temperature units
#######################################################################################################

def GetTRankineFromTFahrenheit(TFahrenheit: float) -> float:
    """
    Utility function to convert temperature to degree Rankine (°R)
    given temperature in degree Fahrenheit (°F).

    Args:
        TRankine: Temperature in degree Fahrenheit (°F)

    Returns:
        Temperature in degree Rankine (°R)

    Notes:
        Exact conversion.

    """
    # Zero degree Fahrenheit (°F) expressed as degree Rankine (°R)
    ZERO_FAHRENHEIT_AS_RANKINE = 459.67

    TRankine = TFahrenheit + ZERO_FAHRENHEIT_AS_RANKINE
    return TRankine


def GetTKelvinFromTCelsius(TCelsius: float) -> float:
    """
    Utility function to convert temperature to Kelvin (K)
    given temperature in degree Celsius (°C).

    Args:
        TCelsius: Temperature in degree Celsius (°C)

    Returns:
        Temperature in Kelvin (K)

    Notes:
        Exact conversion.

    """
    # Zero degree Celsius (°C) expressed as Kelvin (K)
    ZERO_CELSIUS_AS_KELVIN = 273.15

    TKelvin = TCelsius + ZERO_CELSIUS_AS_KELVIN
    return TKelvin

