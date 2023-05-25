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
        par.rho = sm.symbols('rho')
        par.iota = sm.symbols('iota')
        par.sigma_epsilon = sm.symbols('sigma_epsilon')
        par.R = sm.symbols('R')
        par.delta = sm.symbols('delta')


        val.eta = 0.5
        val.kappa = 1.0
        val.w = 1.0
        val.ell_t = sm.symbols('ell_t') 
        val.rho = 0.90
        val.iota = 0.01
        val.sigma_epsilon = 0.10   
        val.R = (1 + 0.01) ** (1 / 12) 
        val.delta = 0.05
        
        sim.T = 120

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

    def calculate_H(self):
        
        par = self.par
        val = self.val
        sim = self.sim

        T = sim.T

        # Initialize variables
        epsilon = np.random.normal(-0.5 * val.sigma_epsilon ** 2, val.sigma_epsilon, size=T)
        kappa = np.zeros(T)
        kappa[0] = 1.0
        ell = np.zeros(T)
        ell[0] = ((1 - val.eta) * kappa[0] / val.w) ** (1 / val.eta)
        value_function = np.zeros(T)

        # Simulate the shock series and compute the value function
        for t in range(1, T):
         kappa[t] = np.exp(val.rho * np.log(kappa[t - 1]) + epsilon[t])
         ell[t] = ((1 - val.eta) * kappa[t] / val.w) ** (1 / val.eta)
         value_function[t] = val.R ** -t * (kappa[t] * ell[t] ** (1 - val.eta) - val.w * ell[t] - (ell[t] != ell[t - 1]) * val.iota)

        # Calculate the expected value of the salon
        K = 1  # Number of shock series
        expected_value = np.mean([np.sum(value_function) for _ in range(K)])

        print("Expected value of the salon (H):", expected_value)

    def calculate_H_K(self, K):
        par = self.par
        val = self.val
        sim = self.sim

        T = sim.T
        results = []

        for _ in range(K):
            # Initialize variables
            epsilon = np.random.normal(-0.5 * val.sigma_epsilon ** 2, val.sigma_epsilon, size=T)
            kappa = np.zeros(T)
            kappa[0] = 1.0
            ell = np.zeros(T)
            ell[0] = ((1 - val.eta) * kappa[0] / val.w) ** (1 / val.eta)
            value_function = np.zeros(T)

            # Simulate the shock series and compute the value function
            for t in range(1, T):
                kappa[t] = np.exp(val.rho * np.log(kappa[t - 1]) + epsilon[t])
                ell[t] = ((1 - val.eta) * kappa[t] / val.w) ** (1 / val.eta)
                value_function[t] = val.R ** -t * (kappa[t] * ell[t] ** (1 - val.eta) - val.w * ell[t] - (ell[t] != ell[t - 1]) * val.iota)

            # Calculate the expected value of the salon
            expected_value = np.sum(value_function)
            results.append(expected_value)

        average_H = np.mean(results)
        return average_H

    
    def H_Delta(self, K):
        par = self.par
        val = self.val
        sim = self.sim

        # Initialize the sum of salon values
        sum_h = 0

        # Perform Monte Carlo simulation
        for k in range(K):
            # Generate shock series
            np.random.seed(k)  # Set the seed for reproducibility
            epsilon = np.random.normal(loc=-0.5 * val.sigma_epsilon**2, scale=val.sigma_epsilon, size=sim.T)

            # Initialize variables
            kappa_prev = 1
            ell_prev = 0
            h = 0

            # Calculate the value of the salon for each period
            for t in range(sim.T):
                kappa = np.exp(val.rho * np.log(kappa_prev) + epsilon[t])
                ell_ast = ((1 - val.eta) * kappa / val.w)**(1 / val.eta)

                if t == 0 or abs(ell_prev - ell_ast) > val.delta:
                    ell = ell_ast
                else:
                    ell = ell_prev

                profit = kappa * ell**(1 - val.eta) - val.w * ell - (ell != ell_prev) * val.iota
                h += val.R**(-t) * profit

                # Update variables for the next period
                kappa_prev = kappa
                ell_prev = ell

            sum_h += h

           # Calculate the expected value of the salon
            H = sum_h / K
            return H