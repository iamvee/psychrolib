import psychrolib
from psychrolib import main

args = 40, 20, 101325
x = main.Psychrometrics(*args)
main.SetUnitSystem(main.SI)

print(x.humidity_ratio)
print(main.CalcPsychrometricsFromTWetBulb(*args))