import numpy as np
from scipy import optimize
import sympy as sm
from sympy import Symbol
from sympy.solvers import solve
sm.init_printing(use_unicode=True) # for pretty printing
from IPython.display import display
import matplotlib.pyplot as plt # baseline module
import ipywidgets as widgets
from types import SimpleNamespace

class ConsModel():
    def __init__(self,do_print=True):

        if do_print: print('initializing the model')

        self.par = SimpleNamespace()
        self.val = SimpleNamespace()
        
        if do_print: 
            print('calling.setup()')
        self.setup()
    
    def setup(self):
        """define the baseline parameters of the Solow model"""
        par = self.par
        val = self.val

        par.L = sm.symbols('L', positive=True)
        par.C = sm.symbols('C')
        par.G = sm.symbols('G', positive=True)

        par.alpha = sm.symbols('alpha')
        par.kappa = sm.symbols('kappa', positive=True)
        par.v = sm.symbols('v', positive=True)
        par.w = sm.symbols('w', positive=True)
        par.tau = sm.symbols('tau')
        par.wtilde = sm.symbols('wtilde')

    def solve_analytical(self):
        par = self.par
        
        # Define the utility function:
        V_func = sm.log(par.C**par.alpha * par.G**(1-par.alpha))-par.v*(par.L**2)/2

        # Define the consumption function:
        #par.w_tilde = (1-par.tau) * par.w
        C_func = par.kappa+par.wtilde*par.L

        #Substitute the consumption function into the utility function
        V_sub = V_func.subs(par.C, C_func)

        # Calculate the first derivative of the utility function with respect to L
        dV_dL = sm.diff(V_sub, par.L)
        eq = sm.Eq(dV_dL,0)

        # Solve the first-order condition
        optimal_L = sm.solve(eq, par.L)

        # Print the optimal labor input
        #print("Optimal Labor Input (L):", optimal_L)
        return optimal_L
    
    def plot1(self):

        par = self.par

        # Set the range of w values
        w_values = np.linspace(0.00001,6,100)

        # Create an empty list to store results
        optimal_L_values =[]

        # Iterate over w values
        for w in w_values:
            par.wtilde = (1 - par.tau) * w
            optimal_L = self.solve_analytical()
            optimal_L_values.append(optimal_L[0])

        # Plot the result
        plt.plot(w_values, optimal_L_values)
        plt.xlabel('w')
        plt.ylabel('L')
        plt.title('Optimal Labor Input as a Function of w')
        plt.grid(True)
        plt.show()
 

