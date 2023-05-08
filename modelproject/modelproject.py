from scipy import optimize
import numpy as np
from scipy import optimize
from scipy import linalg
import sympy as sm
from sympy.solvers import solve
from sympy import Symbol
import matplotlib.pyplot as plt
from types import SimpleNamespace
from tabulate import tabulate
import ipywidgets as widgets
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
        par.s = sm.symbols('s')
        par.n = sm.symbols('n')
        par.g = sm.symbols('g')
        par.d = sm.symbols('D')
        par.dT = sm.symbols('dT')

        # model parameters
    
        val.se = 0.005     # oil udvingings rate
        val.s = 0.3     # savings rate
        val.n = 0.01    # population growth rate
        val.delta = 0.05    # capital depreciation rate
        val.alpha = 0.2     # capital share of output
        val.g = 0.027    # technological progress rate
        val.beta =  0.6 
        val.etha = 0
        val.phi =0
    
        #simulation parameters
        par.simT =100  # 100 perioder
        sim.D = np.zeros(par.simT)  # climate economic damages 
        sim.E = np.zeros(par.simT) # Oil usage
        sim.R = np.zeros(par.simT) # Oil stock left
        sim.k = np.zeros(par.simT) # Physical capital
        sim.L = np.zeros(par.simT) # Workforce
        sim.A = np.zeros(par.simT) # total factor productivit
        sim.Y = np.zeros(par.simT) # Output - GDP