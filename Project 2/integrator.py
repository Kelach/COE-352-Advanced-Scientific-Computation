quadratureLookUp = {
    1: {
        "points": [0],
        "weights": [2]
    },
    2: {
        "points": [-0.5773502691896257, 0.5773502691896257],
        "weights": [1, 1]
    },
    3: {
        "points": [-0.7745966692414834, 0, 0.7745966692414834],
        "weights": [0.5555555555555556, 0.8888888888888888, 0.5555555555555556]
    },
    4: {
        "points": [-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526],
        "weights": [0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538]
    },
    5: {
        "points": [-0.906179845938664, -0.5384693101056831, 0, 0.5384693101056831, 0.906179845938664],
        "weights": [0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891]
    }
}
def integrationWithQaudrature(N: int, f: callable):
    """
    Function that performs integration with quadrature
    """
    # returns sum of all N calculations of weight * f(points) 
    return sum ([ quadratureLookUp[N]["weights"][i]*f(quadratureLookUp[N]["points"][i]) for i in range(N) ])

def getMappedFunction(f: callable, a : int, b : int, h: int):
    """
    Maps given function f in to reference domain [-1, 1]
    """
    # return lambda x: (h/2) * f( a + (b-a)*(x+1) / 2)
    return lambda eta: (h/2) * f (  (eta + 1)*(h/2) + a )