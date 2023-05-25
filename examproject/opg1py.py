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

        par.L = sm.symbols('L', nonnegative=True)
        par.C = sm.symbols('C', nonnegative=True)
        par.G = sm.symbols('G', positive=True)

        par.alpha = sm.symbols('alpha', nonnegative=True)
        par.kappa = sm.symbols('kappa', positive=True)
        par.v = sm.symbols('v', positive=True)
        par.w = sm.symbols('w', positive=True)
        par.tau = sm.symbols('tau', nonnegative=True)
        par.wtilde = sm.symbols('wtilde', positive=True)

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

        # We set the values of the tax close to 0 and 1 as the model cannot solve for tau = 0 and tau = 1
        tau_values = np.linspace(0.001,0.999,100)

        # Create an empty list to store results
        optimal_L_values =[]
        G_func_values =[]
        V_func_values =[]
        C_func_values =[]

        par.alpha = 0.5
        par.kappa = 1
        par.v = 1/(2*(16**2))
        par.w = 1

        # Iterate over tau values
        for tau in tau_values:
            par.tau = tau
            par.wtilde = (1 - par.tau) * par.w
            
            optimal_L = self.solve_analytical()
            optimal_L_values.append(optimal_L[0])

            G_func = par.tau * par.w * optimal_L[0] * par.wtilde
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
 

