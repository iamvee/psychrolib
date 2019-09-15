from .psychrolib import *


def CalcPsychrometricsFromTDewPoint(TDryBulb: float, TDewPoint: float, Pressure: float) -> tuple:
    """
    Utility function to calculate humidity ratio, wet-bulb temperature, relative humidity,
    vapour pressure, moist air enthalpy, moist air volume, and degree of saturation of air given
    dry-bulb temperature, dew-point temperature, and pressure.

    Args:
        TDryBulb : Dry-bulb temperature in °F [IP] or °C [SI]
        TDewPoint : Dew-point temperature in °F [IP] or °C [SI]
        Pressure : Atmospheric pressure in Psi [IP] or Pa [SI]

    Returns:
        Humidity ratio in lb_H₂O lb_Air⁻¹ [IP] or kg_H₂O kg_Air⁻¹ [SI]
        Wet-bulb temperature in °F [IP] or °C [SI]
        Relative humidity in range [0, 1]
        Partial pressure of water vapor in moist air in Psi [IP] or Pa [SI]
        Moist air enthalpy in Btu lb⁻¹ [IP] or J kg⁻¹ [SI]
        Specific volume of moist air in ft³ lb⁻¹ [IP] or in m³ kg⁻¹ [SI]
        Degree of saturation [unitless]

    """
    HumRatio = GetHumRatioFromTDewPoint(TDewPoint, Pressure)
    TWetBulb = GetTWetBulbFromHumRatio(TDryBulb, HumRatio, Pressure)
    RelHum = GetRelHumFromHumRatio(TDryBulb, HumRatio, Pressure)
    VapPres = GetVapPresFromHumRatio(HumRatio, Pressure)
    MoistAirEnthalpy = GetMoistAirEnthalpy(TDryBulb, HumRatio)
    MoistAirVolume = GetMoistAirVolume(TDryBulb, HumRatio, Pressure)
    DegreeOfSaturation = GetDegreeOfSaturation(TDryBulb, HumRatio, Pressure)
    return HumRatio, TWetBulb, RelHum, VapPres, MoistAirEnthalpy, MoistAirVolume, DegreeOfSaturation


def CalcPsychrometricsFromRelHum(TDryBulb: float, RelHum: float, Pressure: float) -> tuple:
    """
    Utility function to calculate humidity ratio, wet-bulb temperature, dew-point temperature,
    vapour pressure, moist air enthalpy, moist air volume, and degree of saturation of air given
    dry-bulb temperature, relative humidity and pressure.

    Args:
        TDryBulb : Dry-bulb temperature in °F [IP] or °C [SI]
        RelHum : Relative humidity in range [0, 1]
        Pressure : Atmospheric pressure in Psi [IP] or Pa [SI]

    Returns:
        Humidity ratio in lb_H₂O lb_Air⁻¹ [IP] or kg_H₂O kg_Air⁻¹ [SI]
        Wet-bulb temperature in °F [IP] or °C [SI]
        Dew-point temperature in °F [IP] or °C [SI].
        Partial pressure of water vapor in moist air in Psi [IP] or Pa [SI]
        Moist air enthalpy in Btu lb⁻¹ [IP] or J kg⁻¹ [SI]
        Specific volume of moist air in ft³ lb⁻¹ [IP] or in m³ kg⁻¹ [SI]
        Degree of saturation [unitless]

    """
    HumRatio = GetHumRatioFromRelHum(TDryBulb, RelHum, Pressure)
    TWetBulb = GetTWetBulbFromHumRatio(TDryBulb, HumRatio, Pressure)
    TDewPoint = GetTDewPointFromHumRatio(TDryBulb, HumRatio, Pressure)
    VapPres = GetVapPresFromHumRatio(HumRatio, Pressure)
    MoistAirEnthalpy = GetMoistAirEnthalpy(TDryBulb, HumRatio)
    MoistAirVolume = GetMoistAirVolume(TDryBulb, HumRatio, Pressure)
    DegreeOfSaturation = GetDegreeOfSaturation(TDryBulb, HumRatio, Pressure)
    return HumRatio, TWetBulb, TDewPoint, VapPres, MoistAirEnthalpy, MoistAirVolume, DegreeOfSaturation


def CalcPsychrometricsFromTWetBulb(TDryBulb: float, TWetBulb: float, Pressure: float) -> tuple:
    """
    Utility function to calculate humidity ratio, dew-point temperature, relative humidity,
    vapour pressure, moist air enthalpy, moist air volume, and degree of saturation of air given
    dry-bulb temperature, wet-bulb temperature, and pressure.

    Args:
        TDryBulb : Dry-bulb temperature in °F [IP] or °C [SI]
        TWetBulb : Wet-bulb temperature in °F [IP] or °C [SI]
        Pressure : Atmospheric pressure in Psi [IP] or Pa [SI]

    Returns:
        Humidity ratio in lb_H₂O lb_Air⁻¹ [IP] or kg_H₂O kg_Air⁻¹ [SI]
        Dew-point temperature in °F [IP] or °C [SI]
        Relative humidity in range [0, 1]
        Partial pressure of water vapor in moist air in Psi [IP] or Pa [SI]
        Moist air enthalpy in Btu lb⁻¹ [IP] or J kg⁻¹ [SI]
        Specific volume of moist air in ft³ lb⁻¹ [IP] or in m³ kg⁻¹ [SI]
        Degree of saturation [unitless]

    """
    HumRatio = GetHumRatioFromTWetBulb(TDryBulb, TWetBulb, Pressure)
    TDewPoint = GetTDewPointFromHumRatio(TDryBulb, HumRatio, Pressure)
    RelHum = GetRelHumFromHumRatio(TDryBulb, HumRatio, Pressure)
    VapPres = GetVapPresFromHumRatio(HumRatio, Pressure)
    MoistAirEnthalpy = GetMoistAirEnthalpy(TDryBulb, HumRatio)
    MoistAirVolume = GetMoistAirVolume(TDryBulb, HumRatio, Pressure)
    DegreeOfSaturation = GetDegreeOfSaturation(TDryBulb, HumRatio, Pressure)
    return HumRatio, TDewPoint, RelHum, VapPres, MoistAirEnthalpy, MoistAirVolume, DegreeOfSaturation


class Psychrometrics:
    def __init__(self,
                 temperature_of_dry_bulb: Any,
                 temperature_of_wet_bulb: Any,
                 pressure_value: Any,
                 system: MeasurementUnitClass = SI_SYSTEM):
        if temperature_of_wet_bulb > temperature_of_dry_bulb:
            raise ValueError("Wet bulb temperature is above dry bulb temperature")

        self.unit_system = system

        if self.unit_system == SI_SYSTEM:
            self.temperature_of_dry_bulb = Celsius(temperature_of_dry_bulb)
            self.temperature_of_wet_bulb = Celsius(temperature_of_wet_bulb)
            self.pressure = Pascal(pressure_value)
        else:
            self.temperature_of_dry_bulb = Fahrenheit(temperature_of_dry_bulb)
            self.temperature_of_wet_bulb = Fahrenheit(temperature_of_wet_bulb)
            self.pressure = Psi(pressure_value)

        self._temperature_dry = Fahrenheit(self.temperature_of_dry_bulb)
        self._temperature_wet = Fahrenheit(self.temperature_of_wet_bulb)

        if self.temperature_of_wet_bulb.celsius < -100 or self.temperature_of_dry_bulb.celsius > 200:
            raise ValueError("Dry bulb temperature must be in range [-100, 200]°C")

        if self.unit_system == IP_SYSTEM:
            self.vapour_bounds = [Fahrenheit(-148), Fahrenheit(392)]
        else:
            self.vapour_bounds = [Celsius(-100), Celsius(200)]

    @property
    def report(self) -> tuple:
        return (self.humidity_ratio,
                self.temperature_of_dew_point,
                self.rel_humidity,
                self.vap_pressure,
                self.MoistAirEnthalpy,
                self.MoistAirVolume,
                self.DegreeOfSaturation)

    @property
    def diff_temperature_dry_wet(self) -> DeltaTemperature:
        return self.temperature_of_dry_bulb - self.temperature_of_wet_bulb

    @property
    def humidity_ratio(self) -> float:
        if isIP():
            if self.temperature_of_wet_bulb >= FREEZING_POINT_WATER_IP:
                HumRatio = \
                    ((
                             1093 - 0.556 * self.temperature_of_wet_bulb.fahrenheit) * self.sat_hum_ratio - 0.240 * self.diff_temperature_dry_wet.fahrenheit) / (
                            1093 + 0.444 * self.temperature_of_dry_bulb.fahrenheit - self.temperature_of_wet_bulb.fahrenheit)
            else:
                HumRatio = ((
                                    1220 - 0.04 * self.temperature_of_wet_bulb.fahrenheit) * self.sat_hum_ratio - 0.240 * self.diff_temperature_dry_wet.fahrenheit) / (
                                   1220 + 0.444 * self.temperature_of_dry_bulb.fahrenheit - 0.48 * self.temperature_of_wet_bulb.fahrenheit)
        else:
            if self.temperature_of_dry_bulb >= FREEZING_POINT_WATER_SI:
                HumRatio = ((
                                    2501. - 2.326 * self.temperature_of_wet_bulb.celsius) * self.sat_hum_ratio - 1.006 * self.diff_temperature_dry_wet.celsius) / (
                                   2501. + 1.86 * self.temperature_of_dry_bulb.celsius - 4.186 * self.temperature_of_wet_bulb.celsius)
            else:
                HumRatio = ((
                                    2830. - 0.24 * self.temperature_of_wet_bulb.celsius) * self.sat_hum_ratio - 1.006 * self.diff_temperature_dry_wet.celsius) / (
                                   2830. + 1.86 * self.temperature_of_dry_bulb.celsius - 2.1 * self.temperature_of_wet_bulb.celsius)
        # Validity check.
        result = max(HumRatio, MIN_HUM_RATIO)

        # TODO
        # if result < 0:
        #     raise ValueError("Humidity ratio cannot be negative")

        return result

    @property
    def bounded_humidity_ratio(self) -> float:
        return max(self.humidity_ratio, MIN_HUM_RATIO)

    @property
    def temperature_of_dew_point(self):
        if self.humidity_ratio < 0:
            raise ValueError("Humidity ratio cannot be negative")

        TDewPoint = GetTDewPointFromVapPres(self.temperature_of_dry_bulb, self.vap_pressure)
        return TDewPoint

    @staticmethod
    def _calc_saturation_vapour_pressure(tmp: Temperature) -> float:
        TDryBulb = tmp
        if isIP():
            t = tmp.rankine

            if (tmp.fahrenheit <= TRIPLE_POINT_WATER_IP.fahrenheit):
                LnPws = (-1.0214165E+04 / t - 4.8932428 - 5.3765794E-03 * t + 1.9202377E-07 * t ** 2 \
                         + 3.5575832E-10 * math.pow(t, 3) - 9.0344688E-14 * math.pow(t, 4) + 4.1635019 * math.log(t))
            else:
                LnPws = -1.0440397E+04 / t - 1.1294650E+01 - 2.7022355E-02 * t + 1.2890360E-05 * t ** 2 \
                        - 2.4780681E-09 * math.pow(t, 3) + 6.5459673 * math.log(t)
        else:
            t = tmp.kelvin

            c = CoefficientPI

            if (tmp.celsius <= TRIPLE_POINT_WATER_SI.celsius):
                LnPws = c.C1 / t + c.C2 + c.C3 * t + c.C4 * t ** 2 + c.C5 * math.pow(t, 3) + c.C6 * math.pow(t,
                                                                                                             4) + c.C7 * math.log(
                    t)
            else:
                LnPws = c.C8 / t + c.C9 + c.C10 * t + c.C11 * t ** 2 + c.c12 * math.pow(t, 3) + c.C13 * math.log(t)

        SatVapPres = math.exp(LnPws)
        return SatVapPres

    def x(self):
        if self.unit_system == IP_SYSTEM:
            BOUNDS = [-148, 392]
        else:
            BOUNDS = [-100, 200]

        # Validity check -- bounds outside which a solution cannot be found
        if self.vap_pressure < GetSatVapPres(BOUNDS[0]) or self.vap_pressure > GetSatVapPres(BOUNDS[1]):
            raise ValueError("Partial pressure of water vapor is outside range of validity of equations")

        # We use NR to approximate the solution.
        # First guess
        TDewPoint = TDryBulb  # Calculated value of dew point temperatures, solved for iteratively
        lnVP = math.log(VapPres)  # Partial pressure of water vapor in moist air

        index = 1

        while True:
            TDewPoint_iter = TDewPoint  # TDewPoint used in NR calculation
            lnVP_iter = math.log(GetSatVapPres(TDewPoint_iter))

            # Derivative of function, calculated analytically
            d_lnVP = dLnPws_(TDewPoint_iter)

            # New estimate, bounded by the search domain defined above
            TDewPoint = TDewPoint_iter - (lnVP_iter - lnVP) / d_lnVP
            TDewPoint = max(TDewPoint, BOUNDS[0])
            TDewPoint = min(TDewPoint, BOUNDS[1])

            if ((math.fabs(TDewPoint - TDewPoint_iter) <= PSYCHROLIB_TOLERANCE)):
                break

            if (index > MAX_ITER_COUNT):
                raise ValueError("Convergence not reached in GetTDewPointFromVapPres. Stopping.")

            index = index + 1

        TDewPoint = min(TDewPoint, TDryBulb)
        return TDewPoint

    @property
    def sat_vap_pressure(self):
        return self._calc_saturation_vapour_pressure(self.temperature_of_dry_bulb)

    @property
    def sat_hum_ratio(self) -> float:
        SatHumRatio = 0.621945 * self.sat_vap_pressure / (self.pressure - self.sat_vap_pressure)

        return max(SatHumRatio, MIN_HUM_RATIO)

    @property
    def vap_pressure(self):
        if self.vap_pressure < 0:
            raise ValueError("Partial pressure of water vapor in moist air cannot be negative")

        # Validity check -- bounds outside which a solution cannot be found
        if (self.vap_pressure < self._calc_saturation_vapour_pressure(self.vapour_bounds[0]) or
                self.vap_pressure > self._calc_saturation_vapour_pressure(self.vapour_bounds[1])
        ):
            raise ValueError("Partial pressure of water vapor is outside range of validity of equations")

        return GetVapPresFromHumRatio(self.humidity_ratio, self.pressure)

    @property
    def rel_humidity(self):
        return self.vap_pressure / self.sat_vap_pressure

    @property
    def MoistAirEnthalpy(self) -> float:
        if isIP():
            return (0.240 * self.temperature_of_dry_bulb.fahrenheit + self.bounded_humidity_ratio * (
                    1061 + 0.444 * self.temperature_of_dry_bulb.fahrenheit))
        else:
            return (1.006 * self.temperature_of_dry_bulb.celsius + self.bounded_humidity_ratio * (
                    2501. + 1.86 * self.temperature_of_dry_bulb.celsius)) * 1000

    @property
    def MoistAirVolume(self):
        if isIP():
            return R_DA_IP * GetTRankineFromTFahrenheit(self.temperature_of_dry_bulb) * (
                    1 + 1.607858 * self.bounded_humidity_ratio) / (
                           144 * self.pressure)
        else:
            return R_DA_SI * GetTKelvinFromTCelsius(self.temperature_of_dry_bulb) * (
                    1 + 1.607858 * self.bounded_humidity_ratio) / self.pressure

    @property
    def DegreeOfSaturation(self):
        pass


if __name__ == '__main__':
    pass
