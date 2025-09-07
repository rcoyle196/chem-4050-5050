import pandas as pd
import matplotlib.pyplot as plt
import numpy as py

df = pd.read_csv(r"C:\Users\rcoyl\OneDrive\Documents\git_hub_wexler\--1\volume_pressure_data.csv")
print (df.head()) #test to see if its working

#stats of data
mean_v = df['Volume'].mean()
mean_p = df['Pressure'].mean()
median_v = df['Volume'].median()
median_p = df['Pressure'].median()
std_V = df['Volume'].std()
std_P = df['Pressure'].std()
print()
print (f"\t volume \tpressure \n mean \t {mean_v} \t\t{mean_p:.2f}\n median  {median_v} \t\t{median_p:.2f}\n std \t {std_V:.2f} \t\t{std_P:.2f}")

#print(type (df['Volume'])) #making sure type is correct
plt.plot(df['Volume'], df['Pressure'], color = 'red', linestyle = '-') #use df['col'] as x and y input no need to set varriables at this time
plt.xlabel('volume')
plt.ylabel('Pressure')
plt.title('pressure/volume')
plt.grid(True)
plt.show()

#quadratic
x=df['Volume'] 
y=df['Pressure']
quad = py.polyfit(x, y, 2) # setting up polynomial
a, b, c = quad #setting coeffients 
y_fit = a*x**2 + b*x + c #setting up an equation
plt.scatter(df['Volume'], df['Pressure'], color = 'red', label = 'csv #s', marker = 'o')
plt.plot(df['Volume'], y_fit, color = 'blue', label = 'quadratic fit')

# 2 differnt plots that are imposed on one graph
plt.xlabel('volume')
plt.ylabel('Pressure')
plt.title('pressure/volume')
plt.grid(True)
plt.legend()
plt.show() # can set the plots to a varriable to only show one at a time
print()
print(f"y = {a:.2f}x**2 {b:.2f}x {c:.2f}")

print()

#linear fit would be
inv = py.polyfit(x, y, 1)
m, k = inv
y_lin_fit = x*m + k

#graph adding using graph model above
plt.scatter(x, y, color = 'red', label = 'csv #s', marker = 'o')
plt.plot(x, y_fit, color = 'blue', label = 'quadratic fit')
plt.plot(x, y_lin_fit, color = 'green', label = 'inverse fit')
plt.xlabel('volume')
plt.ylabel('Pressure')
plt.title('pressure/volume')
plt.grid(True)
plt.legend()
plt.show() 
print(f"y = {m:.2f}x + {b:.2f}")

#root mean squared
rms_q = py.sqrt(py.mean((y-y_fit)**2))
print(f"quadratic root mean squared \n{rms_q:.2f}")
rms_l = py.sqrt(py.mean((y-y_lin_fit)**2))
print(f"root mean squared inverse \n{rms_l:.2f}")

#quadratic seems to be a better fit, and it does not seem to fit the typical model of p*v-c = 0
#the quadratic is a better fit but it still seems to be off with a rsm above 3 (maybe a exponential would be better)