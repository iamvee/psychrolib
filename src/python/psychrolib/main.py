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

        if self.unit_system == IP_SYSTEM:
            self.vapour_bounds = [Fahrenheit(-148), Fahrenheit(392)]
        else:
            self.vapour_bounds = [Celsius(-100), Celsius(200)]

        # cached data
        self._humidity_ratio = None

    # todo
    def validate(self):
        if self.temperature_of_wet_bulb.celsius < -100 or self.temperature_of_dry_bulb.celsius > 200:
            raise ValueError("Dry bulb temperature must be in range [-100, 200]°C")

    @property
    def report(self) -> tuple:
        return (self.humidity_ratio,
                self.temperature_of_dew_point,
                self.rel_humidity,
                self.vap_pressure,
                self.MoistAirEnthalpy,
                self.MoistAirVolume,
                self.DegreeOfSaturation)

    # fixme : Static Methods
    @staticmethod
    def _calc_saturation_vapour_pressure(tmp: Temperature) -> Pressure:
        t = tmp.kelvin
        t2, t3, t4 = t ** 2, t ** 3, t ** 4
        log_t = math.log(t)
        c = CoefficientSI

        if tmp.celsius <= TRIPLE_POINT_WATER_SI.celsius:
            ln_pws = c.C01 / t + c.C02 + c.C03 * t + c.C04 * t2 + c.C05 * t3 + c.C06 * t4 + c.C07 * log_t
        else:
            ln_pws = c.C08 / t + c.C09 + c.C10 * t + c.C11 * t2 + c.C12 * t3 + c.C13 * log_t

        res = math.exp(ln_pws)
        return Pascal(res)

    @staticmethod
    def _calc_vapour_pressure_from_humidity_ratio(humidity_ratio: HumidityRatio,
                                                  pressure_: Pressure,
                                                  bounded_humidity_ratio: HumidityRatio = None):
        if bounded_humidity_ratio is None:
            bounded_humidity_ratio = HumidityRatio(max(humidity_ratio.value, MIN_HUM_RATIO))
        return pressure_.pascal * bounded_humidity_ratio.value / (0.621945 + bounded_humidity_ratio.value)

    @staticmethod
    def _calculate_dev_point_from_vapour_pressure(temperature_of_dry_bulb:Temperature,
                                                  vapour_pressure: Pressure,
                                                  unit_system: UnitSystem,
                                                  vapour_bounds: list):

        temperature_of_dew_point = temperature_of_dry_bulb
        lnVP = math.log(vapour_pressure.pascal)  # Partial pressure of water vapor in moist air

        index = 1

        while True:
            temperature_of_dew_point_iter = temperature_of_dew_point
            lnVP_iter = math.log(GetSatVapPres(temperature_of_dew_point_iter))

            # Derivative of function, calculated analytically
            d_lnVP = dLnPws_(temperature_of_dew_point_iter)

            # New estimate, bounded by the search domain defined above
            TDewPoint = temperature_of_dew_point_iter - (lnVP_iter - lnVP) / d_lnVP
            TDewPoint = max(TDewPoint, vapour_bounds[0])
            TDewPoint = min(TDewPoint, vapour_bounds[1])

            if ((math.fabs(temperature_of_dew_point - temperature_of_dew_point_iter) <= PSYCHROLIB_TOLERANCE)):
                break

            if (index > MAX_ITER_COUNT):
                raise ValueError("Convergence not reached in GetTDewPointFromVapPres. Stopping.")

            index = index + 1

        temperature_of_dew_point = min(temperature_of_dew_point, temperature_of_dry_bulb)
        return temperature_of_dew_point

    # fixme : proerties
    @property
    def diff_temperature_dry_wet(self) -> DeltaTemperature:
        return self.temperature_of_dry_bulb - self.temperature_of_wet_bulb

    @property
    def humidity_ratio(self) -> HumidityRatio:
        if self._humidity_ratio is None:
            if self.temperature_of_dry_bulb >= Water.FREEZING_POINT:
                a = (2501. - 2.326 * self.temperature_of_wet_bulb.celsius) * self.sat_hum_ratio.value
                b = 1.006 * self.diff_temperature_dry_wet.celsius
                c = 2501. + 1.86 * self.temperature_of_dry_bulb.celsius - 4.186 * self.temperature_of_wet_bulb.celsius

            else:
                a = (2830. - 0.24 * self.temperature_of_wet_bulb.celsius) * self.sat_hum_ratio.value
                b = 1.006 * self.diff_temperature_dry_wet.celsius
                c = 2830. + 1.86 * self.temperature_of_dry_bulb.celsius - 2.1 * self.temperature_of_wet_bulb.celsius

            hr = (a - b) / c
            self._humidity_ratio = max(hr, MIN_HUM_RATIO.value)

        return self._humidity_ratio

    @property
    def bounded_humidity_ratio(self) -> HumidityRatio:
        val = max(self.humidity_ratio.value, MIN_HUM_RATIO.value)
        return HumidityRatio(value=val)

    @property
    def temperature_of_dew_point(self):
        if self.humidity_ratio < 0:
            raise ValueError("Humidity ratio cannot be negative")

        TDewPoint = GetTDewPointFromVapPres(self.temperature_of_dry_bulb, self.vap_pressure)
        return TDewPoint

    @property
    def sat_vap_pressure(self) -> Pressure:
        return self._calc_saturation_vapour_pressure(self.temperature_of_dry_bulb)

    @property
    def sat_hum_ratio(self) -> HumidityRatio:
        SatHumRatio = 0.621945 * self.sat_vap_pressure.pascal / (self.pressure.pascal - self.sat_vap_pressure.pascal)
        res = max(SatHumRatio, MIN_HUM_RATIO)
        return HumidityRatio(value=res)

    @property
    def vap_pressure(self) -> Pressure:
        if self.vap_pressure < 0:
            raise ValueError("Partial pressure of water vapor in moist air cannot be negative")

        # Validity check -- bounds outside which a solution cannot be found
        # fixme :
        if (self.vap_pressure < self._calc_saturation_vapour_pressure(self.vapour_bounds[0]) or
                self.vap_pressure > self._calc_saturation_vapour_pressure(self.vapour_bounds[1])):
            raise ValueError("Partial pressure of water vapor is outside range of validity of equations")

        return self._calc_vapour_pressure_from_humidity_ratio(self.humidity_ratio, self.pressure,
                                                              bounded_humidity_ratio=self.bounded_humidity_ratio)

    @property
    def rel_humidity(self):
        return self.vap_pressure.pascal / self.sat_vap_pressure.pascal

    @property
    def MoistAirEnthalpy(self) -> float:
        if isIP():
            return (0.240 * self.temperature_of_dry_bulb.fahrenheit + self.bounded_humidity_ratio.value * (
                    1061 + 0.444 * self.temperature_of_dry_bulb.fahrenheit))
        else:
            return (1.006 * self.temperature_of_dry_bulb.celsius + self.bounded_humidity_ratio.value * (
                    2501. + 1.86 * self.temperature_of_dry_bulb.celsius)) * 1000

    @property
    def MoistAirVolume(self):
        if isIP():
            return R_DA_IP * self.temperature_of_dry_bulb.rankine * (1 + 1.607858 * self.bounded_humidity_ratio.value) \
                   / (144 * self.pressure.psi)
        else:
            return R_DA_SI * self.temperature_of_dry_bulb.kelvin * (1 + 1.607858 * self.bounded_humidity_ratio.value) \
                   / self.pressure.pascal

    @property
    def DegreeOfSaturation(self):
        return self.bounded_humidity_ratio.value / self.sat_hum_ratio.value



if __name__ == '__main__':
    pass
