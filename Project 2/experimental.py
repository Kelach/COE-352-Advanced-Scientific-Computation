import numpy as np
import matplotlib.pyplot as plt
from input import userInput
from FEA import FEA

def fowardEulerWithVaryingTimeSteps():
    """
    Runs FEA with different time step sizes using F. Euler time discretization
    """
    for i in range(551, 538,  -4):
        userInput["timeStepSize"] = 1/i
        userInput["timeDescretizationMethod"] = "forwardEuler"
        approxSolutionAtFinalTime = FEA(**userInput).tolist()
        plt.plot(np.linspace(userInput.get("spatialDomain")[0], userInput.get("spatialDomain")[1], userInput.get("globalNodes")), approxSolutionAtFinalTime, label=f"Δt = 1/{i}")
    
    plt.ylabel("Temperature")
    plt.xlabel("Position") 
    plt.legend(loc='best')
    plt.title(f"1D Heat Equation: Approximate Solution (t=1) with varying Δt")
    plt.show()

def backwardEulerWithVaryingTimeSteps():
    """
    Runs FEA with different time step sizes using B. Euler time discretization
    """
    for i in [551, 50, 10, 5, 1, 0.1,]:
        userInput["timeStepSize"] = 1/i
        userInput["timeDescretizationMethod"] = "backwardEuler"
        approxSolutionAtFinalTime = FEA(**userInput).tolist()
        plt.plot(np.linspace(userInput.get("spatialDomain")[0], userInput.get("spatialDomain")[1], userInput.get("globalNodes")), approxSolutionAtFinalTime, label=f"Δt = {1/i :.3f}")
    
    plt.ylabel("Temperature")
    plt.xlabel("Position") 
    plt.legend(loc='best')
    plt.title(f"1D Heat Equation: Approximate Solution (t=1) with varying Δt")
    plt.show()

def fowardEulerWithVaryingNodes():
    """
    Runs FEA with different number of nodes using F. Euler time discretization
    """
    for i in [11, 9, 5, 3]:
        userInput["globalNodes"] = i
        userInput["timeDescretizationMethod"] = "forwardEuler"
        approxSolutionAtFinalTime = FEA(**userInput).tolist()
        plt.plot(np.linspace(userInput.get("spatialDomain")[0], userInput.get("spatialDomain")[1], userInput.get("globalNodes")), approxSolutionAtFinalTime, label=f"# Nodes = {i}")
    
    plt.ylabel("Temperature")
    plt.xlabel("Position") 
    plt.legend(loc='best')
    plt.title(f"1D Heat Equation: Approximate Solution (t=1) with varying Nodes")
    plt.show()


if __name__ == "__main__":
    """
    Returns graphing plots to compare different time step sizes and number of nodes
    """
    # fowardEulerWithVaryingTimeSteps()
    # backwardEulerWithVaryingTimeSteps()
    # fowardEulerWithVaryingNodes()
    pass