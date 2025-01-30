import numpy as np
from scipy.optimize import fsolve
import sys, os
from typing import Callable
import matplotlib.pyplot as plt

import dyMES.rfunctions as rf
import dyMES.default_models as dm


class model:
    def __init__(self, initial_state:float, params:dict = {}, transition_function:Callable = None, num_groups:str = None) -> None:
        """
        Args:
            initial_state: Value of state at t=0
            params: Dict storing all parameters used in transition function
            transition_function: Transition function as described in paper, look at run_model.ipynb for example usage
            num_groups: key to variable in "params" which represents the number of groups in the system
        """

        self.states = [initial_state]
        self.lambdas = [[-0.000001,0]] #Default lambdas used as a starting point to find true values
        self.params = params
        self.time = [0]
        self.num_groups = num_groups

        if transition_function is not None:
            self.func = transition_function

        #Set default model to the original ecological toy model
        else:
            print("Using Default Transition Function")
            print("Steady State at N = 100")
            print("Parameters:", dm.eco_params)
            self.func = dm.eco_transition_function
            self.params = dm.eco_params
            self.num_groups = 'S'
            
        if self.num_groups == None:
            self.num_groups = "No Groups"
            self.params[self.num_groups] = 1
        


        #Initialize lambdas, lambda 2 should be 0, so self.derivatives is set to 0
        self.derivatives = [0]
        self.lambdas = [self.brute_force_update(init=True)]
        self.lambdas[0][1] = 0 
        

        #Calculate initial derivatives, should be 0 at steady state
        Z = rf.R_mean(self.func, self.states[-1], self.params, self.lambdas[-1])
        f_mean = rf.R_mean(self.func, self.states[-1], self.params, self.lambdas[-1], lambda n: self.func(n, self.states[-1], self.params))/Z 
        self.derivatives[-1] = f_mean * self.params[self.num_groups]
        



    def update(self, time:float, dt=0.1, error_lim=float("inf")) -> None:
        """Updates model for "time" duration using timesteps of size "dt"

        Args:
            time: Length of time to perform updates for
            dt: Timestep to use, dt=0.1 seems to be a good value to start off with
            error_lim: Error limit to accept before printing out a warning.
        """

        num_timesteps = int(time/dt)


        #Main Update Loop
        for timestep in range(num_timesteps):
            
            self.time.append(self.time[-1] + dt)
            self.states.append(self.states[-1] + dt * self.derivatives[-1])

            #Update derivatives
            Z = rf.R_mean(self.func, self.states[-1], self.params, self.lambdas[-1])
          
            f_mean = rf.R_mean(self.func, self.states[-1], self.params, self.lambdas[-1], lambda n: self.func(n, self.states[-1], self.params))/Z 
           
            self.derivatives.append(f_mean * self.params[self.num_groups])

     
            new_lambda = self.brute_force_update()
            self.check_constraints(new_lambda, error_lim = error_lim)
                
    
            
            self.lambdas.append(new_lambda)


            
     

    
    def brute_force_update(self, init : bool = False) -> list:
        """Returns lambda calculated with brute force method

        Args:
            single_lambda: set to true when finding initial lambdas, indicates
                that lambda2 should be 0
        Returns:
            list containing new lambdas
            

        """
        def constraints(lambdas: np.array, init_lambdas = False) -> list:
            """Returns error between actual <n>, <nf> and the value calculated by the given lambdas

            Args:
                lambdas: A 1d array of the two lambdas

            Returns:
                2 element list with the difference in <n> and <nf>
            
            """
            if(init_lambdas):
                lambdas[1] = 0
            Z = rf.R_mean(self.func, self.states[-1], self.params, lambdas) #Calculate normalization factor

            n_mean = rf.Rn_mean(self.func, self.states[-1], self.params, lambdas)/Z
            f_mean = rf.R_mean(self.func, self.states[-1], self.params, lambdas, mean_func=lambda n: self.func(n, self.states[-1], self.params))/Z

            if init_lambdas:
                return [n_mean * self.params[self.num_groups] - self.states[-1], 0]
            
            return [n_mean * self.params[self.num_groups] - self.states[-1], f_mean * self.params[self.num_groups] - self.derivatives[-1]]
        

        if init:
           
            new_lambdas = fsolve(constraints, self.lambdas[-1], args=True)
            return new_lambdas
        
        starting_point = [self.lambdas[-1][0] - 0.0000000001, self.lambdas[-1][1] - 0.001]
        new_lambdas = fsolve(constraints, starting_point)
       
     

        return new_lambdas



    def find_steady_state_params(self, param_key : str) -> float:
        """Finds value for param_key such that <f> = 0

        Args: 
            param_key: key of parameter in self.params to be tuned

        Returns:
            value of parameter such that the system is in steady state.

        """
        
        params_copy = dict(self.params)

        def get_derivatives(par_val):

            params_copy[param_key] = par_val
        
            lambdas = self.lambdas[0] #Lambda 1 is set to 0 at initial iteration

            Z = rf.R_mean(self.func, self.states[-1], params_copy, lambdas) #Calculate normalization factor

            
            f_mean = rf.R_mean(self.func, self.states[-1], params_copy, lambdas, mean_func=lambda n: self.func(n, self.states[-1], params_copy))/Z
            
            return f_mean #Want this to be 0
     
        return fsolve(get_derivatives, 0.02)[0]
    
    
    def update_param(self, param_key : str, new_val : float):
        """ Update parameters according to method outlined in paper

        Args:
            param_key: Key in self.params corresponding to value to be changed
            new_val: new value
        """
        
        self.params[param_key] = new_val

        Z = rf.R_mean(self.func, self.states[-1], self.params, self.lambdas[-1]) #Calculate normalization factor 
        f_mean = rf.R_mean(self.func, self.states[-1], self.params, self.lambdas[-1], mean_func=lambda n: self.func(n, self.states[-1], self.params))/Z
        self.derivatives[-1] = f_mean * self.params[self.num_groups]



    def check_constraints(self, lambdas: list, error_lim: float):
        """ Asserts that calculated value of <n> using current lambdas is close to self.state[-1]
        
        Args:
            lambdas: list containing lambdas
            error_lim: error tolerance 
        """
        Z = rf.R_mean(self.func, self.states[-1], self.params, lambdas) #Calculate normalization factor

        n_mean = rf.Rn_mean(self.func, self.states[-1], self.params, lambdas)/Z
        f_mean = rf.R_mean(self.func, self.states[-1], self.params, lambdas, mean_func=lambda n: self.func(n, self.states[-1], self.params))/Z

        iter_num = len(self.states)
        assert (n_mean * self.params[self.num_groups] - self.states[-1]) ** 2 < error_lim**2, "Constraints not satisfied at iteration " + str(iter_num)
        

    def graph(self) -> None:
        """ Displays graph of state vs time, derivatives vs time, lambdas vs time
        """


        fig, axs = plt.subplots(2, 2,constrained_layout=True)
       
        axs[0, 0].plot(self.time,self.states)
        axs[0, 0].set_title('State')
        axs[0, 1].plot(self.time,self.derivatives)
        axs[0, 1].set_title('Derivative')
        axs[1, 0].plot(self.time, np.array(self.lambdas)[:,0])
        axs[1, 0].set_title('Lambda 1')
        axs[1, 1].plot(self.time, np.array(self.lambdas)[:,1])
        axs[1, 1].set_title('Lambda 2')

        for ax in axs.flat:
            ax.set(xlabel='Time')

        plt.show()
