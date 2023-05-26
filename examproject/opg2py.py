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

    #Question 1

    def solve_numerical_kappa1(self):
        par = self.par
        val = self.val
    
        # value for kappa
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
        
        # value for kappa
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

    #Question 2

    def Calculate_H(self):
        
        par = self.par
        val = self.val
        sim = self.sim

        # Specify the number of simulations
        K = 1000

        # Initialize variables
        kappa_prev = 1.0
        ell_prev = 0
        ex_post_value = 0

        # Simulation loop for random shock series
        for k in range(K):
            ex_post_value_k = 0
            kappa_t = np.zeros(120)
            
            # Simulation loop for each month
            for t in range(120):
                # Generate random shock
                epsilon_t = np.random.normal(-0.5 * val.sigma_epsilon**2, val.sigma_epsilon)
                    
                # Compute demand shock for current month
                kappa_t = np.exp(val.rho * np.log(kappa_prev) + epsilon_t)
                    
                # Calculate optimal number of hairdressers
                ell_t = ((1 - 0.5) * kappa_t / 1.0) ** (1 / 0.5)
                    
                # Compute ex post value of the salon for current month
                ex_post_value_k += val.R**(-t) * (kappa_t * ell_t**(1-0.5) - 1.0 * ell_t - (ell_t != ell_prev) * val.iota)
                    
                # Update previous number of hairdressers and demand-shock
                ell_prev = ell_t
                kappa_prev = kappa_t
                
            # Add ex post value of the salon for current simulation to total ex post value
            ex_post_value += ex_post_value_k

        # Calculate expected value of the salon
        H = ex_post_value / K

        # Print the result
        print("Expected value of the salon (H):", H)

    #Question 3

    def calculate_H_Delta(self):
        par = self.par
        val = self.val
        sim = self.sim

        # Initialize variables
        kappa_prev = 1.0
        ell_prev = 0
        ex_post_value = 0

        # Specify the number of simulations
        K = 1000

        # Simulation loop for random shock series
        for k in range(K):
            ex_post_value_k = 0
            kappa_t = np.zeros(120)
            
            # Simulation loop for each month
            for t in range(120):
                # Generate random shock
                epsilon_t = np.random.normal(-0.5 * val.sigma_epsilon**2, val.sigma_epsilon)
                
                # Compute demand shock for current month
                kappa_t = np.exp(val.rho * np.log(kappa_prev) + epsilon_t)
                
                # Calculate optimal number of hairdressers
                ell_star = ((1 - 0.5) * kappa_t / 1.0) ** (1 / 0.5)
                
                # Adjust number of hairdressers based on policy
                if abs(ell_prev - ell_star) > val.delta:
                    ell_t = ell_star
                else:
                    ell_t = ell_prev
                
                # Compute ex post value of the salon for current month
                ex_post_value_k += val.R**(-t) * (kappa_t * ell_t**(1-0.5) - 1.0 * ell_t - (ell_t != ell_prev) * val.iota)
                
                # Update previous number of hairdressers and demand-shock
                ell_prev = ell_t
                kappa_prev = kappa_t
            
            # Add ex post value of the salon for current simulation to total ex post value
            ex_post_value += ex_post_value_k

        # Calculate expected value of the salon
        H = ex_post_value / K

        # Print the result
        print("Expected value of the salon (H):", H)

   # Question 4

    def Maximize_H_Delta(self):
        par = self.par
        val = self.val
        sim = self.sim

        # Specify the number of simulations
        K = 1000

        # Specify the range of values for Delta
        delta_values = np.linspace(0, 0.3, 50) 

        # Create list for storing H values
        H_values = []

        # Initialize variables
        optimal_delta = None
        max_H = float('-inf')

        # Perform grid search
        for delta in delta_values:
            # Initialize variables
            kappa_prev = 1.0
            ell_prev = 0
            ex_post_value = 0
            
            # Simulation loop
            for k in range(K):
                ex_post_value_k = 0
                kappa_t = np.zeros(120)
                
                for t in range(120):
                    # Generate random shock
                    epsilon_t = np.random.normal(-0.5 * val.sigma_epsilon**2, val.sigma_epsilon)
                    
                    # Compute demand shock for current month
                    kappa_t = np.exp(val.rho * np.log(kappa_prev) + epsilon_t)
                    
                    # Calculate optimal number of hairdressers
                    ell_star = ((1 - 0.5) * kappa_t / 1.0) ** (1 / 0.5)
                    
                    # Adjust number of hairdressers based on policy
                    if abs(ell_prev - ell_star) > delta:
                        ell_t = ell_star
                    else:
                        ell_t = ell_prev
                    
                    # Compute ex post value of the salon for current month
                    ex_post_value_k += val.R**(-t) * (kappa_t * ell_t**(1-0.5) - 1.0 * ell_t - (ell_t != ell_prev) * val.iota)
                    
                    # Update previous number of hairdressers and demand-shock
                    ell_prev = ell_t
                    kappa_prev = kappa_t
                
                # Add ex post value of the salon for current simulation to total ex post value
                ex_post_value += ex_post_value_k
            
            # Calculate expected value of the salon
            H = ex_post_value / K

            # Add H to list of H values
            H_values.append(H)
            
            # Update optimal delta if current H is higher
            if H > max_H:
                max_H = H
                optimal_delta = delta

        # Plot the H values
        plt.plot(delta_values, H_values)
        plt.xlabel('Delta')
        plt.ylabel('Expected Value of the Salon (H)')
        plt.title('Optimal Delta for Maximizing Profitability')
        plt.axvline(x=optimal_delta, color='r', linestyle='--', label='Optimal Delta')
        plt.legend()
        plt.grid(True)
        plt.show()

        # Print the result
        print("Optimal Delta:", optimal_delta)
        print("Maximum Expected value of the salon (H):", max_H)

    #Question 5