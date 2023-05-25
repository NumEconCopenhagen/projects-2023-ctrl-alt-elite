
# Imports
from types import SimpleNamespace

import numpy as np
from scipy import optimize

import pandas as pd 
import matplotlib.pyplot as plt

# Define class
class HouseholdSpecializationModelClass:

    def __init__(self):
        """ setup model """

        # a. create namespaces
        par = self.par = SimpleNamespace()
        sol = self.sol = SimpleNamespace()

        # b. preferences
        par.rho = 2.0
        par.nu = 0.001
        par.epsilon = 1.0
        par.omega = 0.5 

        # c. household production
        par.alpha = 0.5
        par.sigma = 1.0

        # d. wages
        par.wM = 1.0
        par.wF = 1.0
        par.wF_vec = np.linspace(0.8,1.2,5)

        # e. targets
        par.beta0_target = 0.4
        par.beta1_target = -0.1

        # f. solution
        sol.LM_vec = np.zeros(par.wF_vec.size)
        sol.HM_vec = np.zeros(par.wF_vec.size)
        sol.LF_vec = np.zeros(par.wF_vec.size)
        sol.HF_vec = np.zeros(par.wF_vec.size)

        sol.beta0 = np.nan
        sol.beta1 = np.nan

    def calc_utility(self,LM,HM,LF,HF):
        """ calculate utility """

        par = self.par
        sol = self.sol

        # a. consumption of market goods
        C = par.wM*LM + par.wF*LF

        # b. home production
        if par.sigma == 1:
            H = HM**(1-par.alpha)*HF**par.alpha
        elif par.sigma == 0:
            H = min(HM,HF)
        else:
            H = ((1-par.alpha)*HM**((par.sigma-1)/par.sigma) + par.alpha*HF**((par.sigma-1)/par.sigma))**(par.sigma/(par.sigma-1))

        # c. total consumption utility
        Q = C**par.omega*H**(1-par.omega)
        utility = np.fmax(Q,1e-8)**(1-par.rho)/(1-par.rho)

        # d. disutility of work
        epsilon_ = 1+1/par.epsilon
        TM = LM+HM
        TF = LF+HF
        disutility = par.nu*(TM**epsilon_/epsilon_+TF**epsilon_/epsilon_)
        
        return utility - disutility

    def solve_discrete(self,do_print=False):
        """ solve model discretely """
        
        par = self.par
        sol = self.sol
        opt = SimpleNamespace()
        
        # a. all possible choices
        x = np.linspace(0,24,49)
        LM,HM,LF,HF = np.meshgrid(x,x,x,x) # all combinations
    
        LM = LM.ravel() # vector
        HM = HM.ravel()
        LF = LF.ravel()
        HF = HF.ravel()

        # b. calculate utility
        u = self.calc_utility(LM,HM,LF,HF)
    
        # c. set to minus infinity if constraint is broken
        I = (LM+HM > 24) | (LF+HF > 24) # | is "or"
        u[I] = -np.inf
    
        # d. find maximizing argument
        j = np.argmax(u)
        
        opt.LM = LM[j]
        opt.HM = HM[j]
        opt.LF = LF[j]
        opt.HF = HF[j]

        # e. print
        if do_print:
            for k,v in opt.__dict__.items():
                print(f'{k} = {v:6.4f}')

        return opt
    

    # Question 3  
    def solve_cont(self, do_print=False):
        """ solve model continously """
        par = self.par
        sol = self.sol
        opt = SimpleNamespace()
        
        def obj(x): #objective function to maximize utility
            LM, HM, LF, HF = x
            return -self.calc_utility(LM, HM, LF, HF)

        # constraints and bounds
        def constraints(x):
            LM, HM, LF, HF = x
            return [24 - LM-HM, 24 -LF-HF]

        constraints = ({'type':'ineq', 'fun' :constraints})

        bounds = ((0,24),(0,24),(0,24),(0,24))

        # call solver
        x0=[6,6,6,6]
        result = optimize.minimize(obj, x0, method='nelder-mead', bounds=bounds, constraints=constraints)

        opt.LM = result.x[0]
        opt.HM = result.x[1]
        opt.LF = result.x[2]
        opt.HF = result.x[3]

        return opt

    #Question 4
    def solve_wF_vec(self,discrete=False):
        """ solve model for vector of female wages """
        par = self.par
        sol = self.sol

        for i, wF in enumerate(par.wF_vec):

            par.wF = wF

            if discrete:
                result = self.solve_discrete()
            else:
                result = self.solve_cont()
            
            sol.LM_vec[i] = result.LM
            sol.HM_vec[i] = result.HM
            sol.LF_vec[i] = result.LF
            sol.HF_vec[i] = result.HF

        return sol

    def run_regression(self):
        """ run regression """

        par = self.par
        sol = self.sol

        x = np.log(par.wF_vec)
        y = np.log(sol.HF_vec/sol.HM_vec)
        A = np.vstack([np.ones(x.size),x]).T
        sol.beta0,sol.beta1 = np.linalg.lstsq(A,y,rcond=None)[0]

        return sol.beta0, sol.beta1

    
    def estimate(self,alpha=0.5,sigma=0.5):
        """ estimate alpha and sigma """

        par = self.par
        sol = self.sol

        # Set the desired alpha and sigma
        par.alpha = alpha
        par.sigma = sigma

       # Solve the model for the vector of female wages in continuous time
        result = self.solve_wF_vec(discrete=False)

        # Perform the regression
        x = np.log(par.wF_vec)
        y = np.log(result.HF_vec / result.HM_vec)
        A = np.vstack([np.ones(x.size), x]).T
        sol.beta0, sol.beta1 = np.linalg.lstsq(A, y, rcond=None)[0]

        # Plot the results
        fig, ax = plt.subplots()
        ax.scatter(x, y, color='blue', label='Data')
        ax.plot(x, A.dot([sol.beta0, sol.beta1]), color='red', label='Regression Line')
        ax.set_xlabel('log(wF)')
        ax.set_ylabel('log(HF/HM)')
        ax.set_title('Regression: log(HF/HM) vs. log(wF)')
        ax.legend()
        

        
        print( sol.beta0, sol.beta1)
        plt.show()
        

    def solve_wF_alpha(self, discrete=False):
     """ Solve model for vector of female wages """
     par = self.par
     sol = self.sol

     for i, wF in enumerate(par.wF_vec):
         par.wF = wF

          # Set alpha to 0.5
         par.alpha = 0.5

         if discrete:
             result = self.solve_discrete()
         else:
             result = self.solve_cont()

         sol.LM_vec[i] = result.LM
         sol.HM_vec[i] = result.HM
         sol.LF_vec[i] = result.LF
         sol.HF_vec[i] = result.HF

         return sol

    def run_regression_alpha(self):
     """ Run regression """
     par = self.par
     sol = self.sol

     x = np.log(par.wF_vec)
     y = np.log(sol.HF_vec / sol.HM_vec)
     A = np.vstack([np.ones(x.size), x]).T
     sol.beta0, sol.beta1 = np.linalg.lstsq(A, y, rcond=None)[0]

     return sol.beta0, sol.beta1

        
