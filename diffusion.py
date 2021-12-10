"""
Created in 2021
@author: Marie Romedenne, Rishi Pillai
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from numpy import *
import pylab as p
import os
import glob
import pandas as pd
from scipy.optimize import curve_fit
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from numpy import *
import pylab as p
import os
import glob
import pandas as pd
from scipy.optimize import curve_fit
import csv
import numpy as np
from scipy import integrate
import argparse

# Main Code #
def column(matrix, i):
    return [row[i] for row in matrix]

def compute_g(u,D,dt):
    global stepsdone

    corr_var=2

    """ given a u (x , t ) in array , compute g (x,t)= D * d ^2 u / dx ^2
    using central differences with spacing h ,
    and return g (x , t ). """
    stepsdone += 1

    d2u_dx2=np.zeros(u.shape, float)

    for i in range (1,len(u)-1):

        d2u_dx2[i] = (u[i+1]-2*u[i]+u[i-1])/(b/N)**2
        if corr_var==1:
            flux = -((1 / density)) * (n1 * (f1) ** (n1)) * ((stepsdone * dt + 1e-6) ** (n1 - 1)) - (f2) / (density)
        else:
            flux = -((1 / density)) * (2 ** (-0.5)) *(A2 ** (0.5)) * ((stepsdone * dt + 1e-6) ** (-0.5)) * Cr_act * 0.43
        u[0] = -1 / 3 * (-1 / D * flux * 2 * (b / N) + u[2] - 4 * u[1])

    i=len(u)-1
    u[i] = u[i-1]

    return D*d2u_dx2

def func(z, n1, f1, f2):
    return ((f1*z)**(n1)+(f2*z))

def advance_time(u,g,dt):
    global counter

    time_passed = stepsdone * dt / 3600  # in hours
    """ Given the array u , the rate of change array g ,
    and a timestep dt , compute the solution for u
    after t , using simple Euler method . """
    
    counter += 1
    if counter >2:
        sys.exit

    u = u + dt *g

    u_list.append(u)
    return u

def surfaceconcentration(u):

    """ Computes surface concentration in wt.% """
    return u[2]

def depthofattack(u,x):

    """ Computes surface concentration in wt.% """
    for w in range(0, len(x)):
        if u[w] < 0:
            attack=x[w]
    return attack



def do_steps(i, nsteps=10):

    global u

    for i in range(nsteps):
        time_passed = stepsdone * dt / 3600  # in hours

        if time_passed <= final_time:

            g = compute_g(u,D,dt)
            u = advance_time(u,g,dt)
            time_passed_list.append(time_passed)

        
            with open(simu_output_file+'-result-u-with-time-step.csv', 'w+') as f1, open(simu_output_file+'-x-inmicrons.csv', 'w', newline='') as f2, open(simu_output_file+'-time_passed_list_inhrs.csv', 'w', newline='') as f3:
                my_writer1 = csv.writer(f1, delimiter=';')
                my_writer2 = csv.writer(f2, delimiter=';')
                my_writer3 = csv.writer(f3, delimiter=';')

                my_writer1.writerow(["concentration profiles in weight fraction"])
                my_writer1.writerows(u_list)

                my_writer2.writerow(["x in microns"])
                my_writer2.writerow(x * 10000)

                my_writer3.writerow(["time_passed_inhrs"])
                my_writer3.writerow(time_passed_list)


        else:
            sys.exit()

    print (" stepsdone=%5d, time_passed =%8gh, surfaceconcentration(u)=%8g, depthofattack(u,x)=%8g" %
            (stepsdone,time_passed,surfaceconcentration(u),depthofattack(u,x)))
    l.set_ydata(u) # update data in plot
    return l,

""" Return corest with kp as function of nature of oxide """
def dX_dt(X, t=0):
        return array([b1 / a1 * A1 / X[1], #dWm/dt
                    b1 / a1 * (b1 / a1 * A1 / X[1] - B1), #dWr/dt
                    1 / a1 * b1 / a1 * A1 / X[1] - b1 / a1 * B1, #dW/dt
                    -b1 / a1 * B1]) #dWs/dt
def run():

    # global variables used across multiple functions
    global simu_output_file
    global n1
    global u
    global f1
    global f2
    global u_list
    global time_passed_list
    global x
    global l
    global Metal_loss
		
    """Integrating the ODE using scipy.integrate"""

    t = linspace(0, final_time, 1000)  # time in h and number of points it is an array
    X0 = array([0.0000001, 0.0000001, 0.0000001, 0.0000001])  # initials conditions
    X, infodict = integrate.odeint(dX_dt, X0, t, full_output=True)

    time_todataframe = pd.DataFrame(t,columns = ['time in hrs'])
    time_todataframe.to_csv(simu_output_file+"-time_hrs.csv")

    """a function that helps print rows of a matrix here X"""
    Wm=column(X, 0)
    Wr=column(X, 1)
    W=column(X, 2)
    Ws=column(X, 3)

    file = open(simu_output_file, 'w+', newline='')

    # writing the data into the file
    with file:
        write = csv.writer(file)
        write.writerow(["Wm - mg.cm-2","Wr- mg.cm-2","W- mg.cm-2","Ws-mg.cm-2"])
        write.writerows(X)

    """plotting"""
    Wm, Wr, W, Ws = X.T
    Metal_loss=Wm/density*10000
    print(Metal_loss)
    
    time_sec=t*3600
    Wm_data=np.array(Wm)

    xdata = time_sec
    ydata = Wm_data

    popt, pcov = curve_fit(func, xdata, ydata, bounds=(0, [0.5, 0.0001, 0.0001]))
    
    print("this is n1, f1 and f2 from the fit", popt)
    n1=popt[0] #unit in nothing
    f1=popt[1] #unit in mg^(1/n1).cm^(-2/n1).s^-1
    f2=popt[2] #unit in mg.cm-2.s-1
    x=np.linspace(a,b,N)

    u_list=[]
    time_passed_list = []
    u = np.zeros(x.shape, float)
    
    for w in range (0,len(x)):
        u[w] = cr_ini/100

    u_list.append(u)
    u_list.append(u)

    # First set up the figure, the axis, and the plot element we want to animate
    fig1 = plt.figure ()  # setup animation
    ax = plt.axes(xlim=(a, b), ylim=(0, cr_ini/100+0.1))
    l ,= plt.plot (x ,u, 'b-o') # plot initial u (x , t )

    # then compute solution and animate

    line_ani = animation.FuncAnimation(fig1, do_steps, range(10000))
    plt.show()
                   
if __name__=='__main__':
    
    stepsdone=0
    counter=0

    parser = argparse.ArgumentParser(description='Process simulation (diffusion)')
    parser.add_argument('input_file', help='simulation input file name')
    parser.add_argument('output_file', help='simulation output file name')

    args = parser.parse_args()    

    # Read simulation input from file
    simu_input_file = args.input_file
    simu_output_file = args.output_file
    simu_input = open(simu_input_file, 'r')
    simu_values = simu_input.readlines()
    simu_input.close()
    
    final_time              = int(simu_values[0])    # hours
    Temp                    = int(simu_values[1])    # Temperature in C entered globally by user
    kp0                     = float(simu_values[2])  # Pre-exponential:Parabolic rate constant for Cr2O3 in mg2cm-4h-1
    Q_kp					= float(simu_values[3]) # Activation energy:Parabolic rate constant for Cr2O3 in
    kp0_ms                  = float(simu_values[4])  # Pre-exponential:Parabolic rate constant for Cr dissolution in KCl-MgCl2 mg2cm-4h-1
    Q_kp_ms					= float(simu_values[5]) # Activation energy:Parabolic rate constant for Cr dissolution in KCl-MgCl2 kJ/mol
    a1                      = float(simu_values[6])  # 2.17, stoichiometric constant for Cr2O3
    b1                      = float(simu_values[7])  # 3.17, , stoichiometric constant for Cr2O3
    a                       = float(simu_values[8])  
    b                       = float(simu_values[9])  # Half component thickness in cm
    N                       = int(simu_values[10])   # number of nodes
    cr_ini                  = float(simu_values[11]) # initial Cr concentration in wt%
    D0                      = float(simu_values[12]) # Pre-exponential:Diffusion coefficient of Cr in m2/s
    Q_D                     = float(simu_values[13]) # Pre-exponential:Activation energy Diffusion coefficient of Cr in kJ/mol
    dt                      = float(simu_values[14]) # seconds
    density					= float(simu_values[15]) # alloy density in kg/m3
    const_Cract             = float(simu_values[16]) # constant for Cr-activity
    slope_Cract             = float(simu_values[17]) # slope for Cr-activity
    
    """calculates parabolic constant in weight of metal / Cr mg2cm-4h-1"""
    A1 = a1 ** 2 * (kp0*exp(-Q_kp*1e3/(8.314*(Temp+273.15)))) / 2
    A2 = kp0_ms*exp(-Q_kp_ms*1e3/(8.314*(Temp+273.15)))
    Cr_act=slope_Cract*(1/(Temp+273.15))+const_Cract
    print(Cr_act)
    B1=0
    D=D0*exp(-Q_D*1e3/(8.314*(Temp+273.15)))

    run()
    print("* Calculation completed.")

 


