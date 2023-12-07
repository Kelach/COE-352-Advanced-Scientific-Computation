from math import pi, e, sin
userInput = {
    "globalNodes":11,
    "spatialDomain": [0, 1],
    "timeDomain": [0, 1],
    "BCs": [0, "NA", 0, "NA"],
    "timeStepSize": 1/551,
    "IC": lambda x: sin(pi*x),
    "f": lambda t, x:  (pi**2 - 1)*(e**-t) * sin(pi*x),
    "m": lambda x: (1 - x)*(1 + x)/4,
    "k": lambda x: 1,
    "quadraturePoints": 2,
    "timeDescretizationMethod": "backwardEuler"
}