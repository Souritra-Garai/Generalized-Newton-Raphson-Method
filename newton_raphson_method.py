import numpy as np
import numpy.linalg as la

# Implementation of Generalized Newton Raphson Method
# Based from Turns, Stephen R. "An Introduction to Combustion", pp - 710-712
# Appendix E, "Generalized Newton’s Method for the Solution of Nonlinear Equations"

debug = False

tolerance = 1e-7

e_tol = 1e-5

def partial_derivative(f, x_vector, j) :

    x_j = x_vector[j]

    epsilon = e_tol
    
    if x_j > 1 :
        
        epsilon = x_j * e_tol

    x_vector_2 = np.copy(x_vector)

    x_vector_2[j] += epsilon

    return ( f(x_vector_2) - f(x_vector) ) / epsilon

def stop_criterion(d_vector, x_vector) :

    return ( ( np.abs(x_vector) > tolerance )*( np.abs( d_vector / x_vector ) < tolerance ) + ( np.abs(x_vector) < tolerance )*( np.abs( d_vector ) < tolerance ) ).all()

def f_vector(f_array, x_vector, n) :

    vector = np.zeros_like(x_vector)

    for i in range(n) :

        vector[i] = f_array[i](x_vector)

    return vector

def norm(f_vector) :

    return np.abs(f_vector).sum()

def Jacobian(f_array, x_vector, n) :

    J = np.zeros((n,n))

    for i in range(n) :

        for j in range(n) :

            J[i][j] = partial_derivative(f_array[i], x_vector, j)

    return J
    
def guess(f_array, x_vector, n) :

    J = Jacobian(f_array, x_vector, n)
    F = -f_vector(f_array, x_vector, n)

    if debug :

        print('Jacobian : ')
        print(J)
        print('- f_vector : ')
        print(F)

    return la.solve(J, F), F

def find_roots(f_array_input, first_guess) :

    n = len(first_guess)

    f_array = np.copy(f_array_input)
    
    x_vector = np.copy(first_guess)

    if debug :

        iteration = 0
        print('iteration : ', iteration)
        print('x_vector : ', x_vector)

    d_vector, f_vector = guess(f_array, x_vector, n)

    prev_norm = curr_norm = norm(f_vector)
        
    while not stop_criterion(d_vector, x_vector) :

        if curr_norm > prev_norm :
            
            d_vector /= 5

        x_vector = x_vector + d_vector

        if debug :

            iteration += 1
            print('iteration : ', iteration)
            print('x_vector : ', x_vector)

        d_vector, f_vector = guess(f_array, x_vector, n)

        prev_norm = curr_norm

        curr_norm = norm(f_vector)

    return x_vector