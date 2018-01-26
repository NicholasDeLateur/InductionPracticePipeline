import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import least_squares
from scipy.integrate import odeint


# Make it work for Python 2+3 and with Unicode
import io
try:
    to_unicode = unicode
except NameError:
    to_unicode = str

# Using some induction data, get k_loss

'''Import some induction time course data'''

params = []

path = "degradation_timecourse.xlsx"
timecourse_pL21433_100nM = pd.read_excel(path)
timecourse_pL21433_100nM_minutes = timecourse_pL21433_100nM.iloc[:,0]
timecourse_pL21433_100nM_GFP = timecourse_pL21433_100nM.iloc[:,1]
timecourse_pL21433_100nM

'''Write the model''' 

# We treat this like constitutive model; we'll figure out what k_express is when broken down as a function of AHL and LuxR later
# All we want is to get k_loss which is independent of k_express, it's just the shape of the induction curve
def constant_expression(y, t, *args):
    # Unpack the arguments
    k_express = args[0] # A summation of transcription and translation
    k_loss = args[1] # A summation of dilution and degredation
    
    # Assign current values
    GFP = y
    
    # Calculate change in values
    dGFP = k_express - GFP*k_loss
    
    # Apply change in values to the current values in the next time step
    return dGFP


'''ODE parameter estimate the model to get k_loss''' 
 
# Timesteps
n_steps = 100000 #fineness of timesteps

sim1_t = np.linspace(min(timecourse_pL21433_100nM_minutes), max(timecourse_pL21433_100nM_minutes), n_steps) 

def residuals_timecourse_pL21433_100nM(p):
    p = tuple(p)
    sim_P = odeint(constant_expression, init_P, timecourse_pL21433_100nM_minutes, args = p)[:,0]
    res = sim_P - timecourse_pL21433_100nM_GFP
    return res

init_P = [0] # Initial condition of GFP
initial_guess = (100, 0.1)
low_bounds = [0, 0]
up_bounds = [1000, 1]
fitted_params = least_squares(residuals_timecourse_pL21433_100nM, initial_guess, bounds=(low_bounds, up_bounds)).x
k_loss = fitted_params[1]
vals = {}
vals["k_loss"] = k_loss
vals["part"] = "cds"
vals["id"] = "GFP"
params.append(vals)
vals = {}
vals["k_loss"] = k_loss
vals["part"] = "cds"
vals["id"] = "LuxR"
params.append(vals)


plt.plot(sim1_t, odeint(constant_expression, init_P, sim1_t, args = tuple(initial_guess))[:,0], c='r', label='GFP - Initial Param Simlulation')
plt.scatter(timecourse_pL21433_100nM_minutes, timecourse_pL21433_100nM_GFP, c='b', label='GFP - Experimental')
plt.plot(sim1_t, odeint(constant_expression, init_P, sim1_t, args = tuple(fitted_params))[:,0], c='g', label='GFP - Fitted Param Simlulation')
plt.legend(loc = 'best')
plt.xlabel('Time (minutes)')
plt.ylabel('GFP (GeoMean MEFL)')
plt.grid()
plt.yscale('log')
plt.show()


# Now that we know GFP (and other protein) k_loss, we can derive k_express for pCons


'''Import the steady state data for pLacIq promoters'''

path = "pLacIq_SS.xlsx"
pLacIq_SS = pd.read_excel(path)
pLacIq_SS


'''Write the model'''

def calculate_kexpress(SS, k_loss):
    return (SS * k_loss)

'''Solve for k_express'''
for x in range(3,-1,-1):
    P018U015_SS = pLacIq_SS.iloc[-1,1+x]
    P018U015_rate = calculate_kexpress(P018U015_SS, k_loss)

    vals = {}
    vals["k_express"] = P018U015_rate
    vals["part"] = "promoter"
    vals["type"] = "constitutive"
    vals["id"] = "P018U01" + str(5+x)
    params.append(vals)

'''Let us quickly check if our time courses and our derived parameters match up for P018U015'''

sim2_t = np.linspace(min(pLacIq_SS.minutes), max(pLacIq_SS.minutes), n_steps) 
init_P = [P018U015_SS]
y = odeint(constant_expression, init_P, sim2_t, args = (P018U015_rate, k_loss))[:,0]

plt.scatter(pLacIq_SS.minutes, pLacIq_SS.P018U015, c='b', label='GFP - Experimental')
plt.plot(sim2_t, y, c='g', label='GFP - Fitted Param Simlulation')
plt.legend(loc = 'best')
plt.xlabel('Time (minutes)')
plt.ylabel('GFP (GeoMean MEFL)')
plt.grid()
plt.yscale('log')
plt.ylim([10,10000])
plt.show()

# To simulate pLux we need to know max, basal, n, and Kd (ignoring different LuxR give different basal/max right now)

'''Import the steady state data for induced expression for a given LuxR amount under a constitutive promoter'''

'''Write the model'''

def hill_function(x, basal, maximal, Kd, n):
    return basal + maximal * (x**n / (Kd + x**n))

def log_equation_to_fit(x, basal, maximal, Kd, n):
    return np.log10(hill_function(x, basal, maximal, Kd, n))
# Unfortunately this is a necessary hack because the curve_fit function doesn't allow you to input function(function)

list_of_params_to_fit = ['basal_rate', 'max_rate', 'K_d', 'n']

def report_paramaters(fit_param_names, fit_param_values, fit_param_stdevs):
    vals = {}
    for each in range(len(fit_param_names)):
        print(fit_param_names[each], 'is', fit_param_values[each], 'with a standard deviation of', fit_param_stdevs[each])
        vals[fit_param_names[each]] = fit_param_values[each]
    vals["part"] = "promoter"
    vals["type"] = "activated"
    return vals

for x in range(7,-1,-1):
    # This is the data for a LuxR under control of P018U015
    print("paramaters for pL2f14" + str(x+33) + "_SS.xlsx")
    path = "pL2f14" + str(x+33) + "_SS.xlsx"
    pL2f1433_SS = pd.read_excel(path)
    pL2f1433_SS


    '''Curve fit to find the Basal, Max, Kd, and n'''

    init_guess = [100, 100, 100, 2]
    low_bounds = [0, 0, 0, 0]
    up_bounds = [1000000, 1000000, 1000000, 10]

    fit_params, covar = curve_fit(log_equation_to_fit, pL2f1433_SS.iloc[:,0], np.log10(pL2f1433_SS.iloc[:,1]), p0 = init_guess, bounds=(low_bounds, up_bounds))
    # Fitting the log of the data to the log of the equation removes bias towards the high points by the residuals
    # Treats the data more equally; allow accurate basal parameter estimation with minimal loss of accuracy of the maximal expression

    std_dev_error_of_fit_params = np.sqrt(np.diag(covar))
    vals = report_paramaters(list_of_params_to_fit, fit_params, std_dev_error_of_fit_params)
    vals["id"] = "pL2f14" + str(x+33)
    params.append(vals)


with io.open('params.json', 'w', encoding='utf8') as outfile:
    str_ = json.dumps(params,
                      indent=4, sort_keys=True,
                      separators=(',', ': '), ensure_ascii=False)
    outfile.write(to_unicode(str_))


pL2f1433_basal_SS = fit_params[0]
pL2f1433_basal_rate = pL2f1433_basal_SS * k_loss
pL2f1433_max_SS = fit_params[1]
pL2f1433_max_rate = pL2f1433_max_SS * k_loss
pL2f1433_Kd = fit_params[2]
pL2f1433_n = fit_params[3]

sim3_x = np.linspace(min(pL2f1433_SS.iloc[:,0]), max(pL2f1433_SS.iloc[:,0]), n_steps)
plt.scatter(pL2f1433_SS.iloc[:,0], pL2f1433_SS.iloc[:,1], c='b', label='Data')
plt.plot(sim3_x, hill_function(sim3_x, *fit_params), c='r', label='fit')
plt.xlabel(pL2f1433_SS.columns[0])
plt.ylabel(pL2f1433_SS.columns[1])
plt.legend(loc = 'best')
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.show()

'''Import the data for time course of induced expression'''

# Specifically we want the time series data of pL2f1433 over a variety of AHL
path = "pL2f1433_tc.xlsx"
pL2f1433_tc = pd.read_excel(path)
pL2f1433_tc

'''Write the model'''

# define your modules
def P018U015():
    return P018U015_rate # MEFL minute**-1
    
def P007U015(LuxR, AHL):
    P007U015_production_rate = pL2f1433_basal_rate + (pL2f1433_max_rate * (AHL**pL2f1433_n / (pL2f1433_Kd + AHL**pL2f1433_n))) 
    return P007U015_production_rate # MEFL minute**-1

def pL2f1433(y, t):
    # Unpack your current amount of each species
    LuxR, GFP, AHL = y
    
    # Determine the change in each species
    dLuxR = P018U015() - k_loss*LuxR
    dGFP = P007U015(LuxR, AHL) - k_loss*GFP
    dAHL = 0 # for now we're assuming AHL was added exogenously and never degrades
    
    # Return the change in each species; make sure same order as your init values
    # scipy.odeint will take these values and apply them to the current value of each species in the next time step for you
    return [dLuxR, dGFP, dAHL]

'''Simulate the induced gene expression'''
# Initial Conditions
# LuxR, GFP, AHL
init_P = [P018U015_SS, pL2f1433_basal_SS, 1000]
  
sim4_t = np.linspace(min(pL2f1433_tc['minutes']), max(pL2f1433_tc['minutes']), n_steps)

# Example simulation at AHL = 1000
sim_P = odeint(pL2f1433, [P018U015_SS, pL2f1433_basal_SS, 1000], sim4_t)

'''compare to time course data'''

# Unpack the multi dimensional data
t_exp = pL2f1433_tc.loc[pL2f1433_tc['lasAHL (nM)'] == 0]
t_exp = t_exp.iloc[:,1]
y_1 = pL2f1433_tc.loc[pL2f1433_tc['lasAHL (nM)'] == 0]
y_1 = y_1.iloc[:,2]
y_2 = pL2f1433_tc.loc[pL2f1433_tc['lasAHL (nM)'] == 0.01]
y_2 = y_2.iloc[:,2]
y_3 = pL2f1433_tc.loc[pL2f1433_tc['lasAHL (nM)'] == 0.1]
y_3 = y_3.iloc[:,2]
y_4 = pL2f1433_tc.loc[pL2f1433_tc['lasAHL (nM)'] == 1]
y_4 = y_4.iloc[:,2]
y_5 = pL2f1433_tc.loc[pL2f1433_tc['lasAHL (nM)'] == 10]
y_5 = y_5.iloc[:,2]
y_6 = pL2f1433_tc.loc[pL2f1433_tc['lasAHL (nM)'] == 100]
y_6 = y_6.iloc[:,2]

plt.scatter(t_exp, y_1, c='g', label='0')
plt.plot(sim4_t, odeint(pL2f1433, [P018U015_SS, pL2f1433_basal_SS, 0], sim4_t)[:,1], c='g', label = 'sim 0')
plt.scatter(t_exp, y_2, c='r', label='0.01')
plt.plot(sim4_t, odeint(pL2f1433, [P018U015_SS, pL2f1433_basal_SS, 0.01], sim4_t)[:,1], c='r', label = 'sim 0.01')
plt.scatter(t_exp, y_3, c='b', label='0.1')
plt.plot(sim4_t, odeint(pL2f1433, [P018U015_SS, pL2f1433_basal_SS, 0.1], sim4_t)[:,1], c='b', label = 'sim 0.1')
plt.scatter(t_exp, y_4, c='y', label='1')
plt.plot(sim4_t, odeint(pL2f1433, [P018U015_SS, pL2f1433_basal_SS, 1], sim4_t)[:,1], c='y', label = 'sim 1')
plt.scatter(t_exp, y_5, c='plum', label='10')
plt.plot(sim4_t, odeint(pL2f1433, [P018U015_SS, pL2f1433_basal_SS, 10], sim4_t)[:,1], c='plum', label = 'sim 10')
plt.scatter(t_exp, y_6, c='black', label='100')
plt.plot(sim4_t, odeint(pL2f1433, [P018U015_SS, pL2f1433_basal_SS, 100], sim4_t)[:,1], c='black', label = 'sim 100')
plt.legend(loc = 'best')
plt.xlabel('Time (minutes)')
plt.ylabel('GFP (GeoMean MEFL)')
plt.grid()
plt.yscale('log')
plt.xlim([-10, 600])
plt.ylim([20, 4000])
plt.show()
