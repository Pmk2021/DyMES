import numpy as np
from typing import Callable

def R_mean(transition_function : Callable, N : float, m_eff: float, params : dict, lambdas : list, mean_func = lambda n: 1) -> float:
    n_array = np.arange(0,m_eff + 1)
    f_array = mean_func(n_array)
    return np.sum(R(transition_function, N, m_eff, params, lambdas) * f_array)

def Rn_mean(transition_function : Callable, N : float, m_eff: float, params : dict, lambdas : list, mean_func = lambda n: 1) -> float:
    n_array = np.arange(0,m_eff + 1)

    return np.sum(R(transition_function, N, m_eff, params, lambdas) * n_array)

def Rn2_mean(transition_function : Callable, N : float, m_eff: float, params : dict, lambdas : list, mean_func = lambda n: 1) -> float:
    n_array = np.arange(0,m_eff + 1)**2
    return np.sum(R(transition_function, N, m_eff, params, lambdas) * n_array)

def R(transition_function : Callable, N : float, m_eff: float, params : dict, lambdas : list, mean_func = lambda n: 1) -> np.array:
    n_array = np.arange(0,m_eff + 1)

    return np.exp(lambdas[0] * n_array + lambdas[1] * transition_function(n_array, N, params))
    