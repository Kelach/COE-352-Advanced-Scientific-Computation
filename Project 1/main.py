# Write MATLAB or PYTHON software that will calculate (1) the equilibrium displacements, (2) the internal stresses and (3) the elongations  of a spring/mass system using the 3 steps discussed in class and in the book. Solve the force balance first, then back-calculate the elongation and internal stress vectors using the displacements.  Your code should be created in a general way where the user will input the following: the number of springs/masses, the spring constants for each spring, the masses, and which boundary condition to apply.   Your code should allow the user to apply one or two fixed ends. 

# For this project, you will need to solve Ku=f.  You will calculate the condition number of K by  (1) calculating and printing the singular values (and eigenvalues) using YOUR SVD algorithm, and (2) calculating and printing a l2-condition number.  Do not use a software call to calculate a condition number.  USE YOUR SVD ROUTINE TO SOLVE THE KU=F SYSTEM!
def main():
    """
    Main function for Project 1. This function will calculate the equilibrium
    displacements, the internal stresses, and the elongations of a spring/mass
    system using thew 3 steps discussed in class and in the book. Solves the
    force balance first, then back-calculates the elongation and internal
    stress vectors using the displacements. This function is created in a
    general way where the user will input the following: the number of
    springs/masses, the spring constants for each spring, the masses, and
    which boundary condition to apply. This function allows the user to apply
    one or two fixed ends.
    """
    
    N = eval(input("Enter the number of springs/masses: "))
    spring_constants = eval(input("Enter the spring constants for each spring: "))
    masses = eval(input("Enter the masses: "))
    conditions = eval(input("Enter the boundary conditions to apply: "))



    

    # main function
    pass