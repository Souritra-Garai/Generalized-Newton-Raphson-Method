# Generalized-Newton-Raphson-Method

Required python packages : numpy, numpy.linalg

The file newton_raphson_method.py contains a implementation of Generalized Newton’s Method for the Solution of Nonlinear Equations.
It is designed to solve system of equations of the kind
$$\vec{f}(\vec{x}) = \vec{0} $$

For using Newton Raphson method to solve the above equation numerically,
Use find_roots function, with the arguments
- f_vector : array of functions
- first_guess : array of real number to be used as first guess

Based on
Turns, Stephen R. "An Introduction to Combustion", pp - 710-712
Appendix E, "Generalized Newton’s Method for the Solution of Nonlinear Equations"
