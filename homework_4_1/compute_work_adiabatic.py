import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid

#function adiotic pressure
def adiabatic_gas(Vi, R, gama, n, T):

    #grid and space to store values
    work_adi = []
    N = 100
    dv_value = np.linspace(Vi, 3*Vi, N) 
    Pi = n*R*T/Vi

    #intergal using multiple volumes to find work at each point then store
    for dv in dv_value:
        V = np.linspace(Vi, dv, N)
        P = Pi*(Vi/V)**gama
        work = -trapezoid(P, V)
        work_adi.append(work)

    return work_adi, dv_value

