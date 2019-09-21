import psychrolib
from psychrolib import main

args = 40, 20, 101325
x = main.Psychrometrics(temperature_of_dry_bulb=40,
                        temperature_of_wet_bulb=20,
                        pressure_value = 101325,
                        system = main.SI_SYSTEM)
main.SetUnitSystem(main.SI)

print(f"{'humidity ratio':>40} |", x.humidity_ratio)
print(f"{'t dewpoint':>40} |", x.temperature_of_dew_point)
print(f"{'rel hum':>40} |", x.rel_humidity)
print(f"{'vap pressure':>40} |", x.vap_pressure)
print(f"{'moist air enthalpy':>40} |", x.MoistAirEnthalpy)
print(f"{'moist air volume':>40} |", x.MoistAirVolume)
print(f"{'degree of saturation':>40} |", x.DegreeOfSaturation)

# print(main.CalcPsychrometricsFromTWetBulb(*args))

y = main.Psychrometrics(temperature_of_dry_bulb=40,
                        temperature_of_dew_point=x.temperature_of_dew_point.celsius,
                        pressure_value = 101325,
                        system = main.SI_SYSTEM)

print(f"{'humidity ratio':>40} |", y.humidity_ratio)
print(f"{'t dewpoint':>40} |", y.temperature_of_dew_point)
print(f"{'rel hum':>40} |", y.rel_humidity)
print(f"{'vap pressure':>40} |", y.vap_pressure)
print(f"{'moist air enthalpy':>40} |", y.MoistAirEnthalpy)
print(f"{'moist air volume':>40} |", y.MoistAirVolume)
print(f"{'degree of saturation':>40} |", y.DegreeOfSaturation)