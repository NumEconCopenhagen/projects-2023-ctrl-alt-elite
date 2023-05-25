from typing import Any
import numpy as np
from scipy import linalg
from scipy import optimize
import sympy as sm
from types import SimpleNamespace
import matplotlib.pyplot as plt


class ProfitClass():
    def __init__(self,do_print=True):
        """create the model yay"""

        if do_print: print('initializing the model')

        self.par = SimpleNamespace()
        self.val = SimpleNamespace()
        self.sim = SimpleNamespace()
        
        if do_print: print('calling.setup()')
        self.setup()
    
    def setup(self):
        """define parameters of the model"""
        val = self.val
        par = self.par
        sim = self.sim

        #model parameters
        par.eta = sm.symbols('eta')
        par.kappa = sm.symbols('kappa')
        par.w = sm.symbols('w')
        par.ell_t = sm.symbols('ell_t')
        par.Pi = sm.symbols('Pi')


        val.eta = 0.5
        val.kappa = 1.0
        val.w = 1.0
        val.ell_t = sm.symbols('ell_t') 
        
    def solve_numerical_kappa1(self):
        par = self.par
        val = self.val
    
        # Choose a specific value for kappa
        val.kappa = 1.0
        
        Pi = val.kappa*val.ell_t**(1-val.eta) -val.w*val.ell_t
        
        # Calculate the derivative of the profit equation with respect to ell_t
        d_Pi = sm.diff(Pi, val.ell_t)
        
        # Find the value of ell_t that maximizes profits
        optimal_ell = sm.solve(d_Pi, val.ell_t)[0]
        
        # Substitute the parameter values into the optimal_ell expression
        optimal_ell_val = optimal_ell.subs([(val.kappa, val.kappa), (val.eta, val.eta), (val.w, val.w)])
        
        optimal_ell_val

        print("Optimal ell_t value for kappa = 1.0:", optimal_ell_val)

    
    def solve_numerical_kappa2(self):
        par = self.par
        val = self.val
        
        # Choose a specific value for kappa
        val.kappa = 2.0
        
        Pi = val.kappa*val.ell_t**(1-val.eta) -val.w*val.ell_t
        
        # Calculate the derivative of the profit equation with respect to ell_t
        d_Pi = sm.diff(Pi, val.ell_t)
        
        # Find the value of ell_t that maximizes profits
        optimal_ell = sm.solve(d_Pi, val.ell_t)[0]
        
        # Substitute the parameter values into the optimal_ell expression
        optimal_ell_val = optimal_ell.subs([(val.kappa, val.kappa), (val.eta, val.eta), (val.w, val.w)])
        
        optimal_ell_val

        print("Optimal ell_t value for kappa = 2.0:", optimal_ell_val)
