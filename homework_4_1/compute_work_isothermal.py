# gas law calculations for isothermal work
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid

# function ideal gas
def ideal_gas (Vi, R, T, n):

    #generate grid and blank matrix 
    work_iso = []
    N = 100
    dv_value = np.linspace(Vi, 3*Vi, N) 

    #intergal using multiple volumes to find work at each point then store
    for dv in dv_value:
        V = np.linspace(Vi, dv, N)
        P = (n*R*T)/V
        work = -trapezoid(P, V)
        work_iso.append(work)
    return work_iso, dv_value
