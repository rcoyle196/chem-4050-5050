import numpy as np
import scipy as sci
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import t
from lecture_07_regression import ols_slope, ols_intercept, ols

#linnear regression
#Hv = a · TB + b (slope of a)
#Trouton’s Rule, Regression, and Uncertainty Analysis
#ln P = -(d_H/R)(1/T)+C
#d_s = 4.5R + RlnT

#read the csv file and get teh column name
df = pd.read_csv("C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\trouton.csv")
list(df.columns)
df['H_v (kcal/mol)'].sort_values
#print(df)
#print(df.head())
#print(df.dtypes)

# set x and y equal to the correct values and convert units
x = df['T_B (K)']
y = df['H_v (kcal/mol)']*4184 #convert to j/mol

#def ols_slope(x, y):
#    x_mean = np.mean(x)
#    y_mean = np.mean(y)
#    sum_n = np.sum((x-x_mean) * (y-y_mean))
#    sum_d = np.sum((x-x_mean) **2)
#    return sum_n/sum_d

#def ols_intercept(x, y):
#    x_mean = np.mean(x)
#    y_mean = np.mean(y)
#    slope = ols_slope(x, y)
#    
#    return y_mean - slope * x_mean
#def ols(x, y):
#    slope = ols_slope(x, y)
#    intercept = ols_intercept(x, y)
#    return slope, intercept

# using the imported functions
ols_slope(x,y)
ols_intercept(x,y)
ols(x,y)

#find the slope using the functions and generate a line
slope, intercept = ols(x,y)
intercept = intercept/1000
print(f"intercept = {intercept:.4f}\nslope(a) = {slope:.4f}")
line = slope * x + intercept

# set matrix to store values for graphs
perfect_liquids_x, perfect_liquids_y = [], []
imperfect_liquids_x, imperfect_liquids_y = [], []
quantum_x, quantum_y = [], []
metal_x, metal_y = [], []

#run loop to get correspoinding x and y values from the data frams and appended them to the empty matrix
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
plt.scatter(perfect_liquids_x, perfect_liquids_y, color = 'purple', label = 'perfect liquids')
plt.scatter(imperfect_liquids_x, imperfect_liquids_y,  color = 'blue', label = 'imperfect liquids')
plt.scatter(quantum_x, quantum_y, color = 'teal', label = 'elemets subject to quantum effects')
plt.scatter(metal_x, metal_y, color = 'grey', label = 'metals')
plt.plot(x, line, color = 'red', label = 'Hv = a · TB + b' )
plt.title('Troutons Rule')
plt.xlabel('H_v jol/mol')
plt.ylabel('T_B in kelvin')
plt.grid(True)
plt.legend()
plt.show()

# show our graph vs with the troutons value expectations imposed on the same graph
exp_value = 88 * x
plt.scatter(perfect_liquids_x, perfect_liquids_y, color = 'purple', label = 'perfect liquids')
plt.scatter(imperfect_liquids_x, imperfect_liquids_y,  color = 'blue', label = 'imperfect liquids')
plt.scatter(quantum_x, quantum_y, color = 'teal', label = 'elemets subject to quantum effects')
plt.scatter(metal_x, metal_y, color = 'grey', label = 'metals')
plt.plot(x, line, color = 'red', label = 'Hv = a · TB + b' )
plt.plot(x, exp_value, '--', color = 'black', label = 'Trouts rule H_v * 88')
plt.title('Troutons Rule')
plt.ylabel('enthalpy of vaporization H_v jol/mo')
plt.xlabel('boiling point T_B in kelvin')
plt.grid(True)
plt.legend()
plt.show()

#find the mean squared root between our expected value and actual/ our expected and troutons expected
mse_t = np.mean(np.square(exp_value - line))
mse = np.mean(np.square(y - line))
print(f"\n\nMSE Troutons rule vs our model {mse_t:.3f}\nMSE actual vs predicted {mse:.3f}")


# find the coffidence intervals and for a fun excersise lets find the confidence interval between our expected values and the troutons expected values
res = y - slope
def sse(res):
    return np.sum(res **2)

def var(res):
    return sse(res)/(len(res)-2)

def se_slope(x, res):
    s_squared = var(res)
    x_mean = np.mean(x)
    return np.sqrt(s_squared / np.sum((x - x_mean) ** 2))

def se_intercept(x, res):
    s_squared = var(res)
    x_mean = np.mean(x)
    n = len(x)
    return np.sqrt(s_squared * (1/n + (x_mean**2) / np.sum((x - x_mean) ** 2)))


def conf_int_slope(x, res, confidence_level):
    se = se_slope(x, res)
    num = len(x)
    deg_f = num - 2
    alpha = 1 - confidence_level
    critical_t_value = t.ppf(1 - alpha/2, deg_f)
    return critical_t_value * se

def conf_int_intercept(x, res, confidence_level):
    se = se_intercept(x, res)
    num = len(x)
    deg_f = num - 2
    alpha = 1 - confidence_level
    critical_t_value = t.ppf(1 - alpha/2, deg_f)
    return critical_t_value * se

print(f"actual vs expected\nslope {slope:.3f}  +/- {conf_int_slope(x, res, 0.95):.3f}")
print(f"intercept {intercept:.3f}   +/- {conf_int_intercept(x, res, 0.95):.3f}")

#Troutons exp values vs our expected values
res_Trouton = exp_value - line
print(f"\n\nexpected vs troutons rule\nslope {slope:.3f} +/- {conf_int_slope(x , res_Trouton, 0.95):.3f}")
print(f"intercept {intercept:.3f}   +/- {conf_int_intercept(x, res_Trouton, 0.95):.3f}\n\n")

print("our model is not great but it does model the trouton rule slope nicely based on confidence interval and eye test\nneither model seems to do a great job especially around the metals and our slope is stepper\nour model will mapp the metals a little better then troutons rule")

#print(sse(res))
#print(var(res))
#print(se_slope(x, res))
#print(se_intercept(x, res))


#final plot brining it all together
plt.scatter(perfect_liquids_x, perfect_liquids_y, color = 'purple', label = 'perfect liquids')
plt.scatter(imperfect_liquids_x, imperfect_liquids_y,  color = 'blue', label = 'imperfect liquids')
plt.scatter(quantum_x, quantum_y, color = 'teal', label = 'elemets subject to quantum effects')
plt.scatter(metal_x, metal_y, color = 'grey', label = 'metals')
plt.plot(x, line, color = 'red', label = 'Hv = a · TB + b' )
plt.plot(x, exp_value, '--', color = 'black', label = 'Trouts rule H_v * 88')
#add text to my plot
plt.text(
    x.max() * 0.7, y.min() * 1.05,
    f"a = {slope:.2f} J/mol·K ± {conf_int_slope(x, res, 0.95):.2f}\n"
    f"b = {intercept:.2f} kJ/mol ± {conf_int_intercept(x, res, 0.95):.2f}",
    fontsize=10,
    bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5')
)
plt.title('Troutons Rule')
plt.ylabel('enthalpy of vaporization H_v jol/mo')
plt.xlabel('boiling point T_B in kelvin')
plt.grid(True)
plt.legend()
plt.savefig('C:\\Users\\rcoyl\\OneDrive\\Documents\\git_hub_wexler\\chem-4050-5050\\homework_3_1\\Troutons_plot', dpi=300)
plt.show()