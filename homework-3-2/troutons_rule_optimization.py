import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats
from scipy.optimize import minimize
from lecture_07_regression import ols_slope, ols_intercept, ols

df = pd.read_csv("C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\trouton.csv")
x = df['T_B (K)']
y = df['H_v (kcal/mol)']*4184 #convert to j/mol

def objective(parm):
    a, b = parm
    res = y - (a*x + b)
    return np.sum(res ** 2)

result = minimize(
    fun= objective, # function to minimize
    x0 = [103,0], # initial guess
    method="Nelder-Mead", # optimization method
    tol=1e-6 # tolerance for termination
)

ols_slope(x,y)
ols_intercept(x,y)
ols(x,y)

a, b = result["x"]
y_min = (a*x) + b

#find the slope using the functions and generate a line
slope, intercept = ols(x,y)
intercept = intercept/1000
mse = np.mean(np.square(result["x"][0] - slope))
print(f"ols slope {slope:.3f} and intercept {intercept:.3f}")
print(f"\nminimize slope {a:.3f} and intercept {b/1000:.3f}")
print(f"\ndifference in slope values between the two methods {slope - a}")
print("\n\nThere is very little difference between the two values either method is acceptable in this model")

#should just import function from previous code but ran out of time with p-chem test (for grpah too)
perfect_liquids_x, perfect_liquids_y = [], []
imperfect_liquids_x, imperfect_liquids_y = [], []
quantum_x, quantum_y = [], []
metal_x, metal_y = [], []
for i, element in enumerate(df['Class']):
    x_i = x.iloc[i]
    y_i = y.iloc[i]
    if element == 'Perfect liquids':
        perfect_liquids_x.append(x_i)
        perfect_liquids_y.append(y_i)
    elif element == 'Imperfect liquids':
        imperfect_liquids_x.append(x_i)
        imperfect_liquids_y.append(y_i)
    elif element == 'Liquids subject to quantum effects':
        quantum_x.append(x_i)
        quantum_y.append(y_i)
    else:
        metal_x.append(x_i)
        metal_y.append(y_i)

#generate plot using different classes from the matrix values stored earlier
plt.scatter(perfect_liquids_x, perfect_liquids_y, color = 'green', label = 'perfect liquids')
plt.scatter(imperfect_liquids_x, imperfect_liquids_y,  color = 'blue', label = 'imperfect liquids')
plt.scatter(quantum_x, quantum_y, color = 'teal', label = 'elemets subject to quantum effects')
plt.scatter(metal_x, metal_y, color = 'grey', label = 'metals')
plt.plot(x, y_min, color = 'purple', label = 'min (Hv = a * TB + b)' )
plt.text(
    x.max() * 0.7, y.min() * 1.05,
    f"a = {a:.2f} J/mol·K\n"
    f"b = {b/1000:.2f} kJ/mol",
    fontsize=10,
    bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5')
    )
plt.title('Troutons Rule')
plt.ylabel('H_v jol/mol')
plt.xlabel('T_B boiling point in kelvin')
plt.grid(True)
plt.legend()
plt.savefig('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\homework3_2\\Trouton_rule_optimization', dpi=300)

plt.show()
