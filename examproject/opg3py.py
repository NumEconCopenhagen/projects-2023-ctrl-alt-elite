import numpy as np
from scipy.optimize import minimize

def griewank(x):
    return griewank_(x[0],x[1])
    
def griewank_(x1,x2):
    A = x1**2/4000 + x2**2/4000
    B = np.cos(x1/np.sqrt(1))*np.cos(x2/np.sqrt(2))
    return A-B+1

def refined_global_optimizer(bounds, tolerance, warmup_iterations, max_iterations):
    # Set the seed for the random uniform distribution
    np.random.seed(101)

    # Define current best solution and its corresponding initial guess
    x_ast = None
    x0_ast = None

    # List of effective initial guesses
    x0_list = []

    # Loop over the maximum number of iterations
    for k in range(max_iterations):
        # Generate a random initial guess for x with uniform distribution
        x_k = np.random.uniform(bounds[0], bounds[1], size=2)

        # If the number of iterations is greater than the number of warm-up iterations,
        if k >= warmup_iterations:
            # Calculate chi and update the effective initial guess x^k0
            chi_k = 0.50 * (2 / (1 + np.exp((k - warmup_iterations) / 100)))
            x_k0 = chi_k * x_k + (1 - chi_k) * x_ast
        else:
            # Otherwise, set the effective initial guess x^k0 to the random initial guess x^k
            x_k0 = x_k

        # Append the effective initial guess x^k0 to the list of effective initial guesses
        x0_list.append(x_k0)

        # Minimize the Griewank function using the BFGS method
        res = minimize(griewank, x_k0, method='BFGS', tol=tolerance)

        # If the current best solution is 0 or the minimum value of the Griewank function is less than the current best solution,
        if k == 0 or res.fun < x_ast:
            x_ast = res.fun
            x0_ast = res.x

        # If the current best solution is less than the tolerance,
        if x_ast < tolerance:
            break
    # Return the optimal solution and the list of effective initial guesses
    return x0_ast, x0_list