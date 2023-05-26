import numpy as np
from scipy import optimize
from scipy.optimize import *
from scipy.optimize import root
import sympy as sm
from sympy import Symbol
from sympy.solvers import solve
sm.init_printing(use_unicode=True) # for pretty printing
from IPython.display import display
import matplotlib.pyplot as plt # baseline module
import ipywidgets as widgets
from types import SimpleNamespace
import math

class ConsModel():
    def __init__(self,do_print=False):

        if do_print: print('initializing the model')

        self.par = SimpleNamespace()
        self.val = SimpleNamespace()
        
        if do_print: 
            print('calling.setup()')
        self.setup()
    
    def setup(self):
        """define the baseline parameters of the model"""
        par = self.par

        par.L = sm.symbols('L', nonnegative=True)
        par.C = sm.symbols('C', nonnegative=True)
        par.G = sm.symbols('G', positive=True)

        par.alpha = sm.symbols('alpha', nonnegative=True)
        par.kappa = sm.symbols('kappa', positive=True)
        par.v = sm.symbols('v', positive=True)
        par.w = sm.symbols('w', positive=True)
        par.tau = sm.symbols('tau', nonnegative=True)
        par.wtilde = sm.symbols('wtilde', positive=True)

        par.rho = sm.symbols('rho')
        par.sigma = sm.symbols('sigma')
        par.eps = sm.symbols('epsilon')

        opt_tau = None

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
        #eq = sm.Eq(dV_dL,0)

        # Solve the first-order condition
        optimal_L = sm.solve(dV_dL, par.L)

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

    def plot2(self):

        par = self.par
        par.alpha = 0.5
        par.kappa = 1
        par.v = 1/(2*(16**2))
        par.w = 1

        # We set the values of the tax close to 0 and 1 as the model cannot solve for tau = 0 and tau = 1
        tau_values = np.linspace(0.001,0.999,100)

        # Create an empty list to store results
        optimal_L_values =[]
        G_func_values =[]
        V_func_values =[]
        C_func_values =[]

        # Iterate over tau values
        for tau in tau_values:
            par.tau = tau
            par.wtilde = (1 - par.tau) * par.w
            
            optimal_L = self.solve_analytical()
            optimal_L_values.append(optimal_L[0])

            G_func = par.tau * par.w * optimal_L[0]
            G_func_values.append(G_func)

            C_func = par.kappa+par.wtilde*optimal_L[0]
            C_func_values.append(C_func)

            V_func = sm.log(C_func**par.alpha * G_func **(1 - par.alpha)) - par.v * (optimal_L[0]**2) / 2
            V_func_values.append(V_func)

        # Plot optimal labor input as a function of tau
        plt.subplot(311)
        plt.plot(tau_values, optimal_L_values)
        plt.xlim(0, 1)
        plt.xlabel('tau')
        plt.ylabel('L')
        plt.title('Optimal labor input as a function of tau')
        plt.grid(True)

        # Plot G_func as a function of tau
        plt.subplot(312)
        plt.plot(tau_values, G_func_values)
        plt.xlim(0, 1)
        plt.xlabel('tau')
        plt.ylabel('G')
        plt.title('Government spenditure as a function of tau')
        plt.grid(True)

        # Plot V_func as a function of tau
        plt.subplot(313)
        plt.plot(tau_values, V_func_values)
        plt.xlim(0, 1)
        plt.xlabel('tau')
        plt.ylabel('V')
        plt.title('Utility as a function of tau')
        plt.grid(True)

        # Adjust the spacing between subplots
        plt.subplots_adjust(wspace=0.5, hspace=0.5)

        # Set the figure size
        fig = plt.gcf()
        fig.set_size_inches(8, 12)

        # Display the plots
        plt.tight_layout()
        plt.show()

    
    def optimal_labour(self): #replace wtilde
        par = self.par
        optimal_L = self.solve_analytical()[0]
        optimal_L = optimal_L.subs(par.wtilde, (1 - par.tau) * par.w)
        return optimal_L
    
    def optimal_tax(self):
        global opt_tau  # Use the global opt_tau parameter

        par = self.par

        par.alpha = 0.5
        par.kappa = 1
        par.v = 1/(2*(16**2))
        par.w = 1
        par.tau = par.tau

        L_opt = self.optimal_labour()
        C_func = par.kappa + (1 - par.tau) * par.w*par.L
        Cl = C_func.subs(par.L, L_opt)

        G_func = par.tau * par.w * par.L
        Gl = G_func.subs(par.L, L_opt)

        V_func = sm.log(Cl**par.alpha * Gl**(1 - par.alpha)) - par.v * (par.L**2) / 2
        Vl = V_func.subs(par.L, L_opt)

        diff_Vl = sm.diff(Vl, par.tau)

        opt_tau= sm.solve(diff_Vl, par.tau)

        #return opt_tau
        print('The optimal tax is:')

        return opt_tau

    def plot3(self):

        par = self.par
        par.alpha = 0.5
        par.kappa = 1
        par.v = 1/(2*(16**2))
        par.w = 1

        # We set the values of the tax close to 0 and 1 as the model cannot solve for tau = 0 and tau = 1
        tau_values = np.linspace(0.001,0.999,100)

        # Create an empty list to store results
        optimal_L_values =[]
        G_func_values =[]
        V_func_values =[]
        C_func_values =[]

        # Iterate over tau values
        for tau in tau_values:
            par.tau = tau
            par.wtilde = (1 - par.tau) * par.w
            
            optimal_L = self.solve_analytical()
            optimal_L_values.append(optimal_L[0])

            G_func = par.tau * par.w * optimal_L[0]
            G_func_values.append(G_func)

            C_func = par.kappa+par.wtilde*optimal_L[0]
            C_func_values.append(C_func)

            V_func = sm.log(C_func**par.alpha * G_func **(1 - par.alpha)) - par.v * (optimal_L[0]**2) / 2
            V_func_values.append(V_func)

        # Plot optimal labor input as a function of tau
        plt.subplot(311)
        plt.plot(tau_values, optimal_L_values)
        plt.axvline(x=opt_tau, color='red', linestyle='--')
        plt.xlim(0, 1)
        plt.xlabel('tau')
        plt.ylabel('L')
        plt.title('Optimal labor input as a function of tau')
        plt.grid(True)

        # Plot G_func as a function of tau
        plt.subplot(312)
        plt.plot(tau_values, G_func_values)
        plt.axvline(x=opt_tau, color='red', linestyle='--')
        plt.xlim(0, 1)
        plt.xlabel('tau')
        plt.ylabel('G')
        plt.title('Government spenditure as a function of tau')
        plt.grid(True)

        # Plot V_func as a function of tau
        plt.subplot(313)
        plt.plot(tau_values, V_func_values)
        plt.axvline(x=opt_tau, color='red', linestyle='--')
        plt.xlim(0, 1)
        plt.xlabel('tau')
        plt.ylabel('V')
        plt.title('Utility as a function of tau')
        plt.grid(True)

        # Adjust the spacing between subplots
        plt.subplots_adjust(wspace=0.5, hspace=0.5)

        # Set the figure size
        fig = plt.gcf()
        fig.set_size_inches(8, 12)

        # Display the plots
        plt.tight_layout()
        plt.show()

    def solve_ces(self):
        par = self.par

        par.alpha = 0.5
        par.kappa = 1
        par.v = 1/(2*(16**2))
        par.w = 1
        par.eps = 1
        par.tau = 0.514531123095038

        scenarios = [
        {"sigma": 1.001, "rho": 1.001},
        {"sigma": 1.5, "rho": 1.5}]

        results = []
    
        for scenario in scenarios:
            par.sigma = scenario["sigma"]
            par.rho = scenario["rho"]

            # Define the utility function
            def utility_func(L):
                C = par.kappa + (1 - par.tau) * par.w * L
                G = par.tau * par.w * L
                return -(((((par.alpha * C ** ((par.sigma - 1) / par.sigma)) + (1 - par.alpha) * G ** ((par.sigma - 1) / par.sigma)) ** (par.sigma / (par.sigma - 1))) ** (1 - par.rho) - 1) / (1 - par.rho) - par.v * (L ** (1 + par.eps) / (1 + par.eps)))

            # Solve for the optimal value of L
            result = minimize(utility_func, x0=10, method='L-BFGS-B', bounds=[(0,24)])  # Use an appropriate initial guess for L
            optimal_L = result.x[0]

            # Calculate the corresponding value of G
            optimal_G = par.tau * par.w * optimal_L

            results.append((optimal_L, optimal_G))

        return results
    
    def tax_ces(self):
        par = self.par

        par.alpha = 0.5
        par.kappa = 1
        par.v = 1 / (2 * (16 ** 2))
        par.w = 1
        par.eps = 1

        scenarios = [
        {"sigma": 1.001, "rho": 1.001},
        {"sigma": 1.5, "rho": 1.5}]

        results = []
    
        for scenario in scenarios:
            par.sigma = scenario["sigma"]
            par.rho = scenario["rho"]

            def V_func(tau, L):
                C = par.kappa + (1 - tau) * par.w * L
                G = tau * par.w * L
                return ((((par.alpha * C ** ((par.sigma - 1) / par.sigma)) + (1 - par.alpha) * G ** ((par.sigma - 1) / par.sigma)) ** (par.sigma / (par.sigma - 1))) ** (1 - par.rho) - 1) / (1 - par.rho) - par.v * (L ** (1 + par.eps) / (1 + par.eps))

            def optimize_tau(tau):
                result_L = minimize_scalar(lambda L: -V_func(tau, L), bounds=(0, 24), method='bounded')
                optimal_L = result_L.x
                return -result_L.fun, optimal_L

            result_tau = minimize_scalar(lambda tau: -optimize_tau(tau)[0], bounds=(0, 1), method='bounded')
            optimal_tau = result_tau.x
            results.append(optimal_tau)

        return results


    
