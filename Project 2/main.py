from math import exp, pi, sin
from input import userInput
from FEA import FEA
import matplotlib.pyplot as plt
import numpy as np

def main():
    """
    Main function that runs the FEA
    """
    print(f"""
{"="*10} Running Continous Galerkin Finite Element Method on Heat Equation {"="*10}
Configurations:
    Number of Nodes: {userInput.get("globalNodes")}
    Delta Time: {userInput.get("timeStepSize")}
    Spatial Bounds: {userInput.get("spatialDomain")}
    Time Bounds: {userInput.get("timeDomain")}
    Boundary Conditions: {userInput.get("BCs")}
    Quadrature Points: {userInput.get("quadraturePoints")}
    Time Discretization Method: {userInput.get("timeDescretizationMethod")}

""")
    x_appox = np.linspace(userInput.get("spatialDomain")[0], userInput.get("spatialDomain")[1], userInput.get("globalNodes"))
    x_exact = np.linspace(userInput.get("spatialDomain")[0], userInput.get("spatialDomain")[1], 100)
    approxSolutionAtFinalTime = FEA(**userInput).tolist()
    exactSolutionAtFinalTime = exp(-userInput.get("timeDomain")[1]) * np.array([sin(pi*node) for node in x_exact])  
    plt.ylabel("Temperature")
    plt.xlabel("Position") 
    plt.plot(x_exact, exactSolutionAtFinalTime, 'r', label="Exact Solution", linewidth=3)
    plt.plot(x_appox, approxSolutionAtFinalTime, '--b', label="Approximate Solution", linewidth=2)
    plt.legend()
    plt.title("1D Heat Equation: Exact vs. Approximate Solution t=1")
    plt.show()

if __name__ == "__main__":
    main()