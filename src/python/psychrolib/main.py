from .psychrolib import *
from typing import Optional, Any


def whatever(TDryBulb, HumRatio, Pressure):
    VapPres = GetVapPresFromHumRatio(HumRatio, Pressure)
    MoistAirEnthalpy = GetMoistAirEnthalpy(TDryBulb, HumRatio)
    MoistAirVolume = GetMoistAirVolume(TDryBulb, HumRatio, Pressure)
    DegreeOfSaturation = GetDegreeOfSaturation(TDryBulb, HumRatio, Pressure)
    return VapPres, MoistAirEnthalpy, MoistAirVolume, DegreeOfSaturation


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

    VapPres, MoistAirEnthalpy, MoistAirVolume, DegreeOfSaturation = whatever(TDryBulb, HumRatio, Pressure)
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

    VapPres, MoistAirEnthalpy, MoistAirVolume, DegreeOfSaturation = whatever(TDryBulb, HumRatio, Pressure)

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

    VapPres, MoistAirEnthalpy, MoistAirVolume, DegreeOfSaturation = whatever(TDryBulb, HumRatio, Pressure)

    return HumRatio, TDewPoint, RelHum, VapPres, MoistAirEnthalpy, MoistAirVolume, DegreeOfSaturation


class Psychrometrics:
    def __init__(self,
                 temperature_of_dry_bulb: Any = None,
                 temperature_of_wet_bulb: Any = None,
                 pressure_value: Any = None,
                 temperature_of_dew_point: Any = None,
                 rel_humidity: Any = None,
                 system: MeasurementUnitClass = SI_SYSTEM):

        self._humidity_ratio = None
        self._temperature_of_dew_point = None
        self._temperature_of_wet_bulb = None
        self._sat_vap_pressure = None
        self._sat_vap_pressure__dry_bulb = None
        self._sat_hum_ratio = None
        self._sat_hum_ratio__dry_bulb = None
        self._vap_pressure = None
        self._rel_humidity = None
        self._moist_air_enthalpy = None
        self._moist_air_volume = None


        self.unit_system = system

        self._temperature_of_dry_bulb = Celsius(temperature_of_dry_bulb)

        if temperature_of_wet_bulb is not None:
            if temperature_of_wet_bulb > temperature_of_dry_bulb:
                raise ValueError("Wet bulb temperature is above dry bulb temperature")
            self._temperature_of_wet_bulb = Celsius(temperature_of_wet_bulb)

        if temperature_of_dew_point is not None:
            self._temperature_of_dew_point = Celsius(temperature_of_dew_point)




        if rel_humidity is not None:
            self._rel_humidity = rel_humidity

        self.pressure = Pascal(pressure_value)

        self.vapour_bounds = [Celsius(-100), Celsius(200)]

        # cached data

    # todo
    def validate(self):
        if self.temperature_of_wet_bulb.celsius < -100 or self.temperature_of_dry_bulb.celsius > 200:
            raise ValueError("Dry bulb temperature must be in range [-100, 200]°C")

        if self.humidity_ratio < 0:
            raise ValueError("Humidity ratio cannot be negative")

        if self.vap_pressure < 0:
            raise ValueError("Partial pressure of water vapor in moist air cannot be negative")

        if (self.vap_pressure < self._calc_saturation_vapour_pressure(self.vapour_bounds[0]) or
                self.vap_pressure > self._calc_saturation_vapour_pressure(self.vapour_bounds[1])):
            raise ValueError("Partial pressure of water vapor is outside range of validity of equations")

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
            bounded_humidity_ratio = HumidityRatio(max(humidity_ratio.value, MIN_HUM_RATIO.value))
        return pressure_.pascal * bounded_humidity_ratio.value / (0.621945 + bounded_humidity_ratio.value)

    @staticmethod
    def _calculate_dew_point_from_vapour_pressure(temperature_of_dry_bulb: Temperature,
                                                  vapour_pressure: Pressure,
                                                  vapour_bounds: list) -> Temperature:

        temperature_of_dew_point = temperature_of_dry_bulb
        lnVP = Pascal(math.log(vapour_pressure.pascal))  # Partial pressure of water vapor in moist air

        index = 1

        while True:
            temperature_of_dew_point_iter = temperature_of_dew_point
            sat_vap_pres_result = Psychrometrics._calc_saturation_vapour_pressure(temperature_of_dew_point_iter)
            lnVP_iter = Pascal(math.log(sat_vap_pres_result.pascal))

            # Derivative of function, calculated analytically
            d_lnVP = derivation_of_ln_saturation_vapour_pressure(temperature_of_dew_point_iter)

            # New estimate, bounded by the search domain defined above
            TDewPoint = temperature_of_dew_point_iter.celsius - (lnVP_iter.pascal - lnVP.pascal) / d_lnVP.pascal
            TDewPoint = max(TDewPoint, vapour_bounds[0].celsius)
            TDewPoint = min(TDewPoint, vapour_bounds[1].celsius)
            temperature_of_dew_point = Celsius(TDewPoint)
            temp_diff = DeltaTemperature(
                math.fabs(temperature_of_dew_point.kelvin - temperature_of_dew_point_iter.kelvin))
            if temp_diff.celsius <= PSYCHROLIB_TOLERANCE_TEMPERATURE.celsius:
                break

            if index > MAX_ITER_COUNT:
                raise ValueError("Convergence not reached in GetTDewPointFromVapPres. Stopping.")

            index = index + 1

        temperature_of_dew_point = min(temperature_of_dew_point.celsius, temperature_of_dry_bulb.celsius)
        return Celsius(temperature_of_dew_point)

    @staticmethod
    def _calculate_hum_ratio_from_temperature_of_dew_point(temperature_of_dew_point: Temperature) -> HumidityRatio:


    # fixme : properties
    @property
    def diff_temperature_dry_wet(self) -> DeltaTemperature:
        return self.temperature_of_dry_bulb - self.temperature_of_wet_bulb

    @property
    def humidity_ratio(self) -> HumidityRatio:

        if self._humidity_ratio is None:
            if self.temperature_of_dry_bulb.celsius >= Water.FREEZING_POINT.celsius:
                a = (2501. - 2.326 * self.temperature_of_wet_bulb.celsius) * self.sat_hum_ratio.value
                b = 1.006 * self.diff_temperature_dry_wet.celsius
                c = 2501. + 1.86 * self.temperature_of_dry_bulb.celsius - 4.186 * self.temperature_of_wet_bulb.celsius

            else:
                a = (2830. - 0.24 * self.temperature_of_wet_bulb.celsius) * self.sat_hum_ratio.value
                b = 1.006 * self.diff_temperature_dry_wet.celsius
                c = 2830. + 1.86 * self.temperature_of_dry_bulb.celsius - 2.1 * self.temperature_of_wet_bulb.celsius

            hr = (a - b) / c
            self._humidity_ratio = HumidityRatio(value=max(hr, MIN_HUM_RATIO.value))

        return self._humidity_ratio

    @property
    def bounded_humidity_ratio(self) -> HumidityRatio:
        val = max(self.humidity_ratio.value, MIN_HUM_RATIO.value)
        return HumidityRatio(value=val)

    @property
    def temperature_of_dew_point(self) -> Temperature:
        if self._temperature_of_dew_point is None:
            self._temperature_of_dew_point = self._calculate_dew_point_from_vapour_pressure(
                self.temperature_of_dry_bulb,
                self.vap_pressure,
                self.vapour_bounds)
        return self._temperature_of_dew_point

    @property
    def temperature_of_dry_bulb(self) -> Temperature:
        if self._temperature_of_dry_bulb is None:
            pass
        return self._temperature_of_dry_bulb


    # TODO : fix this region
    @property
    def temperature_of_wet_bulb(self) -> Temperature:
        if self._temperature_of_wet_bulb is None:
            BoundedHumRatio = self.bounded_humidity_ratio

            TDewPoint = GetTDewPointFromHumRatio(self.temperature_of_dry_bulb,
                                                 self.bounded_humidity_ratio,
                                                 self.pressure)

            # Initial guesses
            TWetBulbSup = TDryBulb
            TWetBulbInf = TDewPoint
            TWetBulb = (TWetBulbInf + TWetBulbSup) / 2

            index = 1
            # Bisection loop
            while ((TWetBulbSup - TWetBulbInf) > PSYCHROLIB_TOLERANCE):

                # Compute humidity ratio at temperature Tstar
                Wstar = GetHumRatioFromTWetBulb(TDryBulb, TWetBulb, Pressure)

                # Get new bounds
                if Wstar > BoundedHumRatio:
                    TWetBulbSup = TWetBulb
                else:
                    TWetBulbInf = TWetBulb

                # New guess of wet bulb temperature
                TWetBulb = (TWetBulbSup + TWetBulbInf) / 2

                if (index >= MAX_ITER_COUNT):
                    raise ValueError("Convergence not reached in GetTWetBulbFromHumRatio. Stopping.")

                index = index + 1
            return TWetBulb

        return self._temperature_of_wet_bulb

    @property
    def sat_vap_pressure(self) -> Pressure:
        if self._sat_vap_pressure is None:
            # self._sat_vap_pressure = self._calc_saturation_vapour_pressure(self.temperature_of_dry_bulb)
            self._sat_vap_pressure = self._calc_saturation_vapour_pressure(self.temperature_of_wet_bulb)
        return self._sat_vap_pressure

    @property
    def sat_vap_pressure__dry_bulb(self) -> Pressure:
        if self._sat_vap_pressure__dry_bulb is None:
            # self._sat_vap_pressure = self._calc_saturation_vapour_pressure(self.temperature_of_dry_bulb)
            self._sat_vap_pressure__dry_bulb = self._calc_saturation_vapour_pressure(self.temperature_of_dry_bulb)
        return self._sat_vap_pressure__dry_bulb

    @property
    def sat_hum_ratio(self) -> HumidityRatio:
        if self._sat_hum_ratio is None:
            SatHumRatio = 0.621945 * self.sat_vap_pressure.pascal / (
                    self.pressure.pascal - self.sat_vap_pressure.pascal)
            res = max(SatHumRatio, MIN_HUM_RATIO)
            self._sat_hum_ratio = HumidityRatio(value=res)

        return self._sat_hum_ratio

    @property
    def sat_hum_ratio__dry_bulb(self) -> HumidityRatio:
        if self._sat_hum_ratio__dry_bulb is None:
            SatHumRatio = 0.621945 * self.sat_vap_pressure__dry_bulb.pascal / (
                    self.pressure.pascal - self.sat_vap_pressure__dry_bulb.pascal)
            res = max(SatHumRatio, MIN_HUM_RATIO)
            self._sat_hum_ratio__dry_bulb = HumidityRatio(value=res)

        return self._sat_hum_ratio__dry_bulb

    @property
    def vap_pressure(self) -> Pressure:
        if self._vap_pressure is None:
            # Validity check -- bounds outside which a solution cannot be found
            # fixme :

            self._vap_pressure = self._calc_vapour_pressure_from_humidity_ratio(
                self.humidity_ratio, self.pressure, bounded_humidity_ratio=self.bounded_humidity_ratio)

        return Pascal(self._vap_pressure)

    @property
    def rel_humidity(self):
        if self._rel_humidity is None:
            self._rel_humidity = self.vap_pressure.pascal / self.sat_vap_pressure__dry_bulb.pascal

        return self._rel_humidity

    @property
    def MoistAirEnthalpy(self) -> float:
        if self._moist_air_enthalpy is None:
            if isIP():
                self._moist_air_enthalpy = (
                        0.240 * self.temperature_of_dry_bulb.fahrenheit + self.bounded_humidity_ratio.value *
                        (1061 + 0.444 * self.temperature_of_dry_bulb.fahrenheit))
            else:
                self._moist_air_enthalpy = (
                        1.006 * self.temperature_of_dry_bulb.celsius + self.bounded_humidity_ratio.value *
                        (2501.0 + 1.86 * self.temperature_of_dry_bulb.celsius))
                self._moist_air_enthalpy *= 1000

        return self._moist_air_enthalpy

    @property
    def MoistAirVolume(self):
        if self._moist_air_volume is None:
            if isIP():
                self._moist_air_volume = (
                        R_DA_IP * self.temperature_of_dry_bulb.rankine *
                        (1 + 1.607858 * self.bounded_humidity_ratio.value) / (144 * self.pressure.psi))
            else:
                self._moist_air_volume = (
                        R_DA_SI * self.temperature_of_dry_bulb.kelvin *
                        (1 + 1.607858 * self.bounded_humidity_ratio.value) / self.pressure.pascal)
        return self._moist_air_volume

    @property
    def DegreeOfSaturation(self):
        return self.bounded_humidity_ratio.value / self.sat_hum_ratio__dry_bulb.value


if __name__ == '__main__':
    pass
