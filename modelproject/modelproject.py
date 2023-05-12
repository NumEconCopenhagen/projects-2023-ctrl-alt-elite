from pyexpat import model
from scipy import optimize
import numpy as np
import scipy.optimize as opt
from scipy import linalg
import sympy as sm
from sympy.solvers import solve
from sympy import Symbol
import matplotlib.pyplot as plt
from types import SimpleNamespace
from tabulate import tabulate
import ipywidgets as widgets
from scipy.optimize import root_scalar
from ipywidgets import interact, interactive, fixed, interact_manual

class OilSolowModelClass():
    def __init__(self,do_print=True):
        """create the model"""

        if do_print: print('initializing the model')

        self.par = SimpleNamespace()
        self.val = SimpleNamespace()
        self.sim = SimpleNamespace()
        
        if do_print: print('calling.setup()')
        self.setup()
    
    def setup(self):
        """define the baseline parameters of the Solow model"""

        val = self.val
        par = self.par
        sim = self.sim

        par.k = sm.symbols('k')
        par.alpha = sm.symbols('alpha')
        par.delta = sm.symbols('delta')
        par.beta = sm.symbols('beta')
        par.etha = sm.symbols('etha')
        par.phi = sm.symbols('phi')
        par.s = sm.symbols('s')
        par.se = sm.symbols('s_E')
        par.n = sm.symbols('n')
        par.g = sm.symbols('g')
        par.d = sm.symbols('D')
        par.dT = sm.symbols('dT')
        par.R = sm.symbols('R')

        # model parameters
    
        val.se = 0.005     # oil udvingings rate
        val.s = 0.3     # savings rate
        val.n = 0.01    # population growth rate
        val.delta = 0.05    # capital depreciation rate
        val.alpha = 0.2     # capital share of output
        val.g = 0.027    # technological progress rate
        val.beta =  0.6   # labor share of output
        val.etha = 0.2    # energy share of output
        val.phi =0.5   # climate damage parameter
        val.phi2 = 0
        val.k_0 = 1         # initial capital stock
        val.L_0 = 1         # initial labor force
        val.A_0 = 1         # initial total factor productivity
        val.R_0 = 1         # initial oil reserves
        val.e = 1
        val.Y_0 = 1
       

    
        #simulation parameters
        par.simT =100  # 100 perioder
        sim.D = np.zeros(par.simT)  # climate economic damages 
        sim.E = np.zeros(par.simT) # Oil usage
        sim.R = np.zeros(par.simT) # Oil stock left
        sim.k = np.zeros(par.simT) # Physical capital
        sim.L = np.zeros(par.simT) # Workforce
        sim.A = np.zeros(par.simT) # total factor productivit
        sim.Y = np.zeros(par.simT) # Output - GDP
        sim.Phi = np.zeros(par.simT) # climate damage parameter
        sim.g = np.zeros(par.simT) # technological progress rate
        sim.z = np.zeros(par.simT)

    def solve_steady_state(self):
        val = self.val
        par = self.par

        # steady state values
        par.k_star = (val.s * val.A_0 ** (1 - val.beta) * (val.L_0 * val.R_0) ** val.beta / (val.delta + val.g)) ** (1 / (1 - val.alpha - val.beta))
        par.L_star = val.L_0 * (1 + val.n) ** (par.simT - 1)
        par.A_star = val.A_0 * (1 + val.g) ** (par.simT - 1)
        par.R

    # How does climate damages develop over time and how much og the production is lost do to this?
    def solve_climate_damage(self):
        par = self.par
        sim = self.sim
        val = self.val

        #looping over each period t

        sim.R[0] = 1
        for t in range(par.simT-1):
            sim.E[t] = sim.R[t] * val.se
            sim.R[t+1] = sim.R[t] * (1 - val.se)
            sim.D[t] = 1 - ((sim.R[0] * (1 - val.se) ** t) / sim.R[t]) ** val.phi
            # calculate output Y

    def calculate_D_t(self, t=50, se=0.005, phi=1):
        R_0 = 1 # initial value of R
        R_t = R_0 * (1 - se) ** t
        D_t = 1 - ((R_0 * (1 - se) ** t) / R_0) ** phi
        return D_t
    
    # find the size of phi
    def find(self):
       val = self.val
       # growth rate equation with phi isolated
       phi = (val.beta + val.etha) * ((val.beta / (val.beta + val.etha)) * val.g - (val.etha / (val.beta + val.etha)) * val.n - (val.etha / (val.beta + val.etha)) * val.se) / (-1 * val.se)
       return abs(phi)
       sum = find()
       print("phi =",sum)

    # find balanced growth phi=0.5
    def balance_climate(self):
        val = self.val

        gy =  (( val.beta)/ (val.beta + val.etha) )* val.g - ((val.etha)/(val.beta + val.etha)) * val.n - ((val.etha)/(val.beta + val.etha)) * val.se - ((val.phi)/(val.beta + val.etha)) * val.se
        return abs(gy)
        sum = balance_climate()
        print ("Phi= 0.5",balance_climate)

    def balance_no_climate(self):
        val = self.val

        gy2 =  (( val.beta)/ (val.beta + val.etha) )* val.g - ((val.etha)/(val.beta + val.etha)) * val.n - ((val.etha)/(val.beta + val.etha)) * val.se - ((val.phi2)/(val.beta + val.etha)) * val.se
        return abs(gy2)
        sum = balance_no_climate()
        print ("Phi= 0",balance_no_climate)

    def solve_output(self): 
        par = self.par
        sim = self.sim
        val = self.val

        # Initial values
        val.k_0 = 1.0
        val.e_0 = 1.0
        val.A_0 = 1.0
        val.Y_0 = sim.Y[t] = (1 - sim.D[t]) * sim.k[t] ** val.alpha * sim.A[t] ** val.beta * val.e ** val.etha 

        # Number of periods to simulate
        T = 100

        #         Arrays to store the results
        sim.k = np.zeros(T+1)
        sim.A = np.zeros(T+1)
        sim.Y = np.zeros(T+1)

        sim.k[0] = val.k_0
        sim.L[0] = val.L_0
        sim.A[0] = val.A_0
        sim.Y[0] = val.Y_0

        # Simulate the model
        for t in range(T):
            sim.Y[t] = sim.A[t] * (sim.k[t]**val.alpha) * (sim.L[t]**(1-val.alpha))
            sim.k[t+1] = val.s * sim.Y[t] * (1-val.delta) + (1-val.delta) * sim.k[t]
            sim.L[t+1] = (1+val.n) * sim.L[t]
            sim.A[t+1] = (1+val.g) * sim.A[t]

        # Plot the output over time
        fig, ax = plt.subplots()
        ax.plot(range(T+1), sim.Y)
        ax.set_xlabel('Time')
        ax.set_ylabel('Output')
        ax.set_title('Output over Time')
        plt.show()


      
