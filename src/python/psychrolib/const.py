from .units import Fahrenheit, Celsius

#######################################################################################################
# Global constants
#######################################################################################################

R_DA_IP = 53.350
"""float: Universal gas constant for dry air (IP version)

    Units:
        ft lb_Force lb_DryAir⁻¹ R⁻¹

    Reference:
        ASHRAE Handbook - Fundamentals (2017) ch. 1
"""

R_DA_SI = 287.042
"""float: Universal gas constant for dry air (SI version)

    Units:
        J kg_DryAir⁻¹ K⁻¹

    Reference:
        ASHRAE Handbook - Fundamentals (2017) ch. 1
"""

MAX_ITER_COUNT = 100
"""int: Maximum number of iterations before exiting while loops.

"""

MIN_HUM_RATIO = 1e-7
"""float: Minimum acceptable humidity ratio used/returned by any functions.
          Any value above 0 or below the MIN_HUM_RATIO will be reset to this value.

"""
class Water:
    FREEZING_POINT = Celsius(0.0)
    TRIPLE_POINT = Celsius(0.01)


FREEZING_POINT_WATER_IP = Fahrenheit(32.0)
"""float: Freezing point of water in Fahrenheit.

"""

FREEZING_POINT_WATER_SI = Celsius(0.0)
"""float: Freezing point of water in Celsius.

"""

TRIPLE_POINT_WATER_IP = Fahrenheit(32.018)
"""float: Triple point of water in Fahrenheit.

"""

TRIPLE_POINT_WATER_SI = Celsius(0.01)
"""float: Triple point of water in Celsius.

"""
