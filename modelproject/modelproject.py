from multiprocessing import Value
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
    
        val.s_E = 0.005     # oil udvingings rate
        val.s = 0.3     # savings rate
        val.n = 0.01    # population growth rate
        val.delta = 0.05    # capital depreciation rate
        val.alpha = 0.2     # capital share of output
        val.g = 0.027    # technological progress rate
        val.beta =  0.6   # labor share of output
        val.etha = 0.2    # energy share of output
        val.phi = 0.5  # climate damage parameter
        val.phi2 = 0
        val.k_0 = 1         # initial capital stock
        val.L_0 = 1         # initial labor force
        val.A_0 = 1         # initial total factor productivity
        val.R = 1         # initial oil reserves
        val.e = 1
        val.Y_0 = 1
        val.D = 0
       

    
        #simulation parameters
        par.simT =100  # 100 perioder
        sim.D = np.zeros(par.simT)  # climate economic damages 
        sim.E = np.zeros(par.simT) # Oil usage
        sim.R = np.zeros(par.simT) # Oil stock left
        sim.K = np.zeros(par.simT) # Physical capital
        sim.L = np.zeros(par.simT) # Workforce
        sim.A = np.zeros(par.simT) # total factor productivit
        sim.Y = np.zeros(par.simT) # Output - GDP
        sim.Phi = np.zeros(par.simT) # climate damage parameter
        sim.g = np.zeros(par.simT) # technological progress rate
        sim.z = np.zeros(par.simT)
        sim.fracY = np.zeros(par.simT)
        sim.fracYD = np.zeros(par.simT)
        sim.fracYDgrowth = np.zeros(par.simT)
        sim.fracY_ext = np.zeros(par.simT)

    def simulate (self, Phi_param):
         #function that makes a simulation of all scenarios 

         par = self.par
         val = self.val
         sim = self.sim

         # looping over each period t
         for t in range(par.simT):
            if t == 0: 
                #Setting values and equations common for all simulations in period 0. 
                K_lag = 1
                L_lag = 1
                A_lag = 1
                R_lag = 1 
                D_lag = 1- (val.R /R_lag)**val.phi
                E_lag = val.s_E*R_lag
                Y_lag = (1 - D_lag) * K_lag** val.alpha *(A_lag*L_lag)**val.beta*E_lag**val.etha


                L = sim.L[t] = L_lag
                A = sim.A[t] = A_lag

                # setting the values for period 0
                if phi_param == 0:
                    # Setting the equations for period 0
                    D_lag = 1- (val.R /R_lag)**phi_param
                    



       





    def simulate1(self):
        par = self.par
        val = self.val
        sim = self.sim
        
     # simulating without climate change
     
        val.phi = 0
        for t in range(par.simT):
            if t == 0: 
             # setting the values for period 0
             A_lag = 1
             K_lag = 1 
             L_lag = 1 
             R_lag = 1 
             D_lag = 1- (val.R /R_lag)**val.phi
             E_lag = val.s_E*R_lag
             Y_lag = (1 - D_lag) * K_lag** val.alpha *(A_lag*L_lag)**val.beta*E_lag**val.etha

             #the model equations for period 0

             L=sim.L[t] = L_lag
             A=sim.A[t] = A_lag
             K= sim.K[t] = K_lag
             R = sim.R[t] =R_lag
             D = sim.D[t] = 1- (val.R /R_lag)**val.phi
             E =sim.E[t] = val.s_E* R_lag
             Y = sim.Y[t] = (1-sim.D[t]) * sim.K[t]** val.alpha *(sim.A[t]*sim.L[t])**val.beta*sim.E[t]**val.etha

            else:
            # setting the lagged values from period t=1 to t=100
             A_lag = sim.A[t-1] = (1+val.g)*A_lag
             K_lag = sim.K[t-1] = val.s*sim.Y[t-1] +(1-val.delta)*K_lag
             L_lag = sim.L[t-1] = (1+val.n)*L_lag
             R_lag = sim.R[t-1] = R_lag - sim.E[t-1]
             D_lag = sim.D[t-1] = 1- (sim.R[t-1] /R_lag)**val.phi
             E_lag = sim.E[t-1] = val.s_E* sim.R[t-1]
             Y_lag = sim.Y[t-1] = (1-sim.D[t-1]) * sim.K[t-1]**val.alpha * (sim.A[t-1]*sim.L[t-1])**val.beta * sim.E[t-1]**val.etha

            # the model equations for period t = 1 to t = 100
             A = sim.A[t] = (1+val.g)*A_lag
             K = sim.K[t] = val.s*Y_lag +(1-val.delta)*K_lag
             L = sim.L[t] = (1+val.n)*L_lag
             R = sim.R[t] = R_lag - E_lag
             D = sim.D[t] = 1- (val.R /R_lag)**val.phi
             E = sim.E[t] = val.s_E* R_lag
             Y = sim.Y[t] = 1-sim.D[t] * sim.K[t]** val.alpha *(sim.A[t]*sim.L[t])**val.beta*sim.E[t]**val.etha
    
        # calculating the relative growth in GDP
        sim.fracY[t] = (sim.Y[t]/sim.L[t])/(sim.Y[0]/sim.L[0])
    
    def simulate2(self):
        par = self.par
        val = self.val
        sim = self.sim
        
     # simulating without climate change
     
        val.phi = 0.5
        for t in range(par.simT):
            if t == 0: 
             # setting the values for period 0
             A_lag = 1
             K_lag = 1 
             L_lag = 1 
             R_lag = 1 
             D_lag = 1- (val.R /R_lag)**val.phi
             E_lag = val.s_E*R_lag
             Y_lag = (1 - D_lag) * K_lag** val.alpha *(A_lag*L_lag)**val.beta*E_lag**val.etha

             #the model equations for period 0

             L=sim.L[t] = L_lag
             A=sim.A[t] = A_lag
             K= sim.K[t] = K_lag
             R = sim.R[t] =R_lag
             D = sim.D[t] = 1- (val.R /R_lag)**val.phi
             E =sim.E[t] = val.s_E* R_lag
             Y = sim.Y[t] = (1-sim.D[t]) * sim.K[t]** val.alpha *(sim.A[t]*sim.L[t])**val.beta*sim.E[t]**val.etha

            else:
            # setting the lagged values from period t=1 to t=100
             A_lag = sim.A[t-1] = (1+val.g)*A_lag
             K_lag = sim.K[t-1] = val.s*sim.Y[t-1] +(1-val.delta)*K_lag
             L_lag = sim.L[t-1] = (1+val.n)*L_lag
             R_lag = sim.R[t-1] = R_lag - sim.E[t-1]
             D_lag = sim.D[t-1] = 1- (sim.R[t-1] /R_lag)**val.phi
             E_lag = sim.E[t-1] = val.s_E* sim.R[t-1]
             Y_lag = sim.Y[t-1] = (1-sim.D[t-1]) * sim.K[t-1]**val.alpha * (sim.A[t-1]*sim.L[t-1])**val.beta * sim.E[t-1]**val.etha


            # the model equations for period t = 1 to t = 100
             A = sim.A[t] = (1+val.g)*A_lag
             K = sim.K[t] = val.s*Y_lag +(1-val.delta)*K_lag
             L = sim.L[t] = (1+val.n)*L_lag
             R = sim.R[t] = R_lag - E_lag
             D = sim.D[t] = 1- (val.R /R_lag)**val.phi
             E = sim.E[t] = val.s_E* R_lag
             Y = sim.Y[t] = 1-sim.D[t] * sim.K[t]** val.alpha *(sim.A[t]*sim.L[t])**val.beta*sim.E[t]**val.etha
    
        # calculating the relative growth in GDP with climate change  
        sim.fracYD[t] = (sim.Y[t]/sim.L[t])/(sim.Y[0]/sim.L[0])
    
    def simulate3(self):
        par = self.par
        val = self.val
        sim = self.sim
        
     # simulating without climate change
     
        val.phi = 1
        for t in range(par.simT):
            if t == 0: 
             # setting the values for period 0
             A_lag = 1
             K_lag = 1 
             L_lag = 1 
             R_lag = 1 
             D_lag = 1- (val.R /R_lag)**val.phi
             E_lag = val.s_E*R_lag
             Y_lag = (1 - D_lag) * K_lag** val.alpha *(A_lag*L_lag)**val.beta*E_lag**val.etha

             #the model equations for period 0

             L=sim.L[t] = L_lag
             A=sim.A[t] = A_lag
             K= sim.K[t] = K_lag
             R = sim.R[t] =R_lag
             D = sim.D[t] = 1- (val.R /R_lag)**val.phi
             E =sim.E[t] = val.s_E* R_lag
             Y = sim.Y[t] = (1-sim.D[t]) * sim.K[t]** val.alpha *(sim.A[t]*sim.L[t])**val.beta*sim.E[t]**val.etha

            else:
            # setting the lagged values from period t=1 to t=100
             A_lag = sim.A[t-1] = (1+val.g)*A_lag
             K_lag = sim.K[t-1] = val.s*sim.Y[t-1] +(1-val.delta)*K_lag
             L_lag = sim.L[t-1] = (1+val.n)*L_lag
             R_lag = sim.R[t-1] = R_lag - sim.E[t-1]
             D_lag = sim.D[t-1] = 1- (sim.R[t-1] /R_lag)**val.phi
             E_lag = sim.E[t-1] = val.s_E* sim.R[t-1]
             Y_lag = sim.Y[t-1] = (1-sim.D[t-1]) * sim.K[t-1]**val.alpha * (sim.A[t-1]*sim.L[t-1])**val.beta * sim.E[t-1]**val.etha


            # the model equations for period t = 1 to t = 100
             A = sim.A[t] = (1+val.g)*A_lag
             K = sim.K[t] = val.s*Y_lag +(1-val.delta)*K_lag
             L = sim.L[t] = (1+val.n)*L_lag
             R = sim.R[t] = R_lag - E_lag
             D = sim.D[t] = 1- (val.R /R_lag)**val.phi
             E = sim.E[t] = val.s_E* R_lag
             Y = sim.Y[t] = 1-sim.D[t] * sim.K[t]** val.alpha *(sim.A[t]*sim.L[t])**val.beta*sim.E[t]**val.etha
    
          # calculating the relative growth in GDP with growing climate change 
        sim.fracYDgrowth[t] = (sim.Y[t]/sim.L[t])/(sim.Y[0]/sim.L[0])


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

        gy =  (( val.beta)/ (val.beta + val.etha) )* val.g - ((val.etha)/(val.beta + val.etha)) * val.n - ((val.etha)/(val.beta + val.etha)) * val.s_E - ((val.phi)/(val.beta + val.etha)) * val.s_E
        return abs(gy)
        sum = balance_climate()
        print ("Phi= 0.5",balance_climate)

    def balance_no_climate(self):
        val = self.val

        gy2 =  (( val.beta)/ (val.beta + val.etha) )* val.g - ((val.etha)/(val.beta + val.etha)) * val.n - ((val.etha)/(val.beta + val.etha)) * val.s_E - ((val.phi2)/(val.beta + val.etha)) * val.s_E
        return abs(gy2)
        sum = balance_no_climate()
        print ("Phi= 0",balance_no_climate)




