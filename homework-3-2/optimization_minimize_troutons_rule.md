# Ols method vs Minimization 

## Minimization code 
    result = minimize (
        fun= funcion #place function here with set varriable
        x0 = [varriables] #initial guess
        method = "Nelder-Mead" #method
        tol=1e-6
    )

## OLS code
    def ols_slope(x, y):
        x_mean = np.mean(x)
        y_mean = np.mean(y)
        sum_n = np.sum((x-x_mean) * (y-y_mean))
        sum_d = np.sum((x-x_mean) **2)
        return sum_n/sum_d
    def ols_intercept(x, y):
        x_mean = np.mean(x)
        y_mean = np.mean(y)
        slope = ols_slope(x, y)    
        return y_mean - slope * x_mean
    def ols(x, y):
        slope = ols_slope(x, y)
        intercept = ols_intercept(x, y)
        return slope, intercept

|Minimize function|OLS method|
|-----------------|----------|
|slower: can be thrown off by intial guess: method matters: | only can be linear: only can handel limited data set
|-----------------|----------|
flexible: can fit most data sets: costom parameters and bounds: works with non anylitical models| faster: stable no adjusted paramaters set by user 
     

