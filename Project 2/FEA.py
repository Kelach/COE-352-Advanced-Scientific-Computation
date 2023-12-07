from collections import defaultdict
from integrator import integrationWithQaudrature, getMappedFunction
from scipy.linalg import inv
import numpy as np

def getBoundaryConditions(N):
    """
    Function that gets the boundary conditions from the user
    """
    BCs = defaultdict(list)
    for i in range(2):
        for type in ["", "_prime"]:
            try:
                BC = eval(input(f"What is the boundary condition at the boundary x{i}{type}? (type 'NA' if nonw"))
                if i == 1:
                    BCs[N-1].append(BC)
                else:
                    BCs[0].append(BC)
            except NameError:
                pass
    return BCs


def FEA(**kwargs):
    """
    1-D Finite Element Analysis Function
    Requires the following parameters:
        - globalNodes: number of global nodes
        - spatialDomain: bounds of the spatial domain
        - timeDomain: bounds of the time domain
        - BCs: boundary conditions
        - timeStepSize: time step size
        - IC: initial condition
        - f: function that represents the external force
        - mappedBasis: basis + test function mapped to reference space
        - quadraturePoints: number of quadrature points
    """
    # externalForce = eval(input("What is the external force?"))
    globalNodesCount = kwargs.get("globalNodes", 11)
    spatialDomain = kwargs.get("spatialDomain", [0, 1])
    timeDomain = kwargs.get("timeDomain", [0, 1])
    timeStepSize = kwargs.get("timeStepSize", 0.1)
    explicitMethod = True if kwargs.get("timeDescretizationMethod", "forwardEuler") == "forwardEuler" else False
    
    # gets boundary conditions from kwargs 
    BCs = [kwargs.get("BCs", [0, 0])[0:2], kwargs.get("BCs", [0, 0])[2:4] ]
    IC = kwargs.get("IC", [0 for _ in range(globalNodesCount)])
    globalElements = globalNodesCount - 1
    stepSize = (spatialDomain[1] - spatialDomain[0]) / globalElements
    
    # defining the grid space
    globalNodes = [spatialDomain[0] + i * stepSize for i in range(globalNodesCount)]
    
    # defining connectivity map
    connectivityMap = [[i, i + 1] for i in range(globalElements)]

    # defining the stiffness matrix
    stiffnessMatrix = np.array([[0.0 for _ in range(globalNodesCount)] for _ in range(globalNodesCount)])
    massMatrix = np.array([[0.0 for _ in range(globalNodesCount)] for _ in range(globalNodesCount)])
    forceVector = np.array([0.0 for _ in range(globalNodesCount)])
    localStiffnessMatrix = [[0, 0] for _ in range(2)]
    localMassMatrix = [[0, 0] for _ in range(2)]
    localForceVector = [0, 0]
    quadraturePoints = kwargs.get("quadraturePoints", 2)
    f = kwargs["f"]     # function f to be integrated over the reference domain [-1, 1]]
    k = kwargs["k"]     # function k to be integrated over the reference domain [-1, 1]]
    m = kwargs["m"]     # function m to be integrated over the reference domain [-1, 1]]

    for element in range(globalElements):
        # Local Element Calculation
        for localNodeIndex in range(2):
            for i in range(2):
                basisFunction1 = lambda x: (1 - x) / 2 if localNodeIndex == 0 else (1 + x ) / 2
                basisFunction2 = lambda x: (1 - x) / 2 if i == 0 else (1 + x ) / 2
                basisFunctionPrime1 = lambda x: -1/stepSize if localNodeIndex == 0 else 1/stepSize 
                basisFunctionPrime2 = lambda x: -1/stepSize if i == 0 else 1/stepSize 
                functionToIntegrate =  lambda x: basisFunctionPrime1(x)*basisFunctionPrime2(x)*(stepSize/2) 
                localStiffnessMatrix[localNodeIndex][i] = integrationWithQaudrature(quadraturePoints, functionToIntegrate)
                functionToIntegrate = lambda x: (basisFunction1(x)*basisFunction2(x))*stepSize/2
                localMassMatrix[localNodeIndex][i] = integrationWithQaudrature(quadraturePoints, functionToIntegrate)
        
        # Finite Element Assembly
        for localNodeIndex in range(2):
            globalNodeIndex = connectivityMap[element][localNodeIndex]
            
            for i in range(2):
                globalNodeIndex2 = connectivityMap[element][i]
                stiffnessMatrix[globalNodeIndex][globalNodeIndex2] += localStiffnessMatrix[localNodeIndex][i]
                massMatrix[globalNodeIndex][globalNodeIndex2] += localMassMatrix[localNodeIndex][i]

                # if globalNodeIndex != globalNodeIndex2:
                #     stiffnessMatrix[globalNodeIndex][globalNodeIndex2] *= -1

    # Solving reduced system (assumes homogenous BC)
    displacementVector = np.array([ IC(xi) for xi in globalNodes[1:-1] ]) # initializes the vector with the initial condition at t=0
    reducedMassMatrix = massMatrix[1:-1, 1:-1]
    inverseReducedMassMatrix = inv(reducedMassMatrix)
    reducedStiffnessMatrix = stiffnessMatrix[1:-1, 1:-1]
    inverseReducedBMatrix = inv(reducedMassMatrix/timeStepSize + reducedStiffnessMatrix)
    currentTime = 0

    # Iterates over time domain by a the timeStepSize
    while currentTime < timeDomain[1]:
        # Building time-dependent f vector
        forceVector = np.array([0.0 for _ in range(globalNodesCount)])
        for globalElementIndex in range(globalElements):
            # Local Element Calculation
            for localNodeIndex in range(2):
                # defining mapped functions
                basisFunction = lambda eta: (1 - eta) / 2 if localNodeIndex == 0 else (1 + eta) / 2
                newTime = currentTime + timeStepSize if explicitMethod is False else currentTime # jumps to next time step if using backward euler
                # defines bounds of local element to be mapped to [-1, 1]
                a, b = globalNodes[globalElementIndex], globalNodes[globalElementIndex+1] 
                fMapped = getMappedFunction(lambda x: f(newTime, x), a, b, stepSize)
                # performing integration
                functionToIntegrate = lambda x: fMapped(x) * basisFunction(x)
                localForceVector[localNodeIndex] = integrationWithQaudrature(quadraturePoints, functionToIntegrate)
            
            # Finite Element Assembly
            for localNodeIndex in range(2):
                globalNodeIndex = connectivityMap[globalElementIndex][localNodeIndex]
                forceVector[globalNodeIndex] += localForceVector[localNodeIndex]
        
        # Updating solution
        if explicitMethod is True:
            # Forward Euler update equation
            displacementVector = (np.eye(len(reducedMassMatrix)) - timeStepSize*inverseReducedMassMatrix @ reducedStiffnessMatrix) @ displacementVector + (timeStepSize * inverseReducedMassMatrix @ forceVector[1:-1])
        else:
            # Backward Euler update equation
            displacementVector = inverseReducedBMatrix @ forceVector[1:-1] + inverseReducedBMatrix @ (reducedMassMatrix @ displacementVector/timeStepSize)
        
        currentTime += timeStepSize

    # applying essential boundary conditions
    displacementVector = np.append(displacementVector, float(BCs[1][0]))
    displacementVector = np.insert(displacementVector, 0, float(BCs[0][0]))

    return displacementVector