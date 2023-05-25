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

        par.L = sm.symbols('L')
        par.C = sm.symbols('C')
        par.G = sm.symbols('G')

        par.alpha = sm.symbols('alpha')
        par.kappa = sm.symbols('kappa')
        par.v = sm.symbols('v')
        par.w = sm.symbols('w')
        par.tau = sm.symbols('tau')
        par.w = sm.symbols('w')

    def solve_analytical(self):
        par = self.par
        # Define the utility function:
        V_func = sm.log(par.C**par.alpha * par.G**(1-par.alpha))-par.v*(par.L**2)/2

        # Define the consumption function:
        #w_tilde = (1-par.tau) * par.w
        C_func = par.kappa+par.w*par.L

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

  #val.alpha = 0.5
        #val.kappa = 1
        #val.v = 1/(2*16**2)
        #val.w = 1
        #val.tau = 0.3