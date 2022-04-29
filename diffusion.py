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

def compute_g(u,D,dt, corr_var):
    global stepsdone

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

def surfaceconcentration(u,x):
    """ Computes surface concentration in wt.% """
    surface=u[0] # is this ok?
    for w in range(0, len(x)):
        if u[w] < 0:
            surface=u[w+2]
    return surface

def depthofattack(u,x):

    """ Computes surface concentration in wt.% """
    for w in range(0, len(x)):
        if u[w] < 0.01*cr_ini:
            attack=1e4*x[w]
    return attack



def do_steps(i, fig1, corr_var, nsteps=10):

    global u
    global time_at_below
    global attack
    time_passed = stepsdone * dt / 3600  # in hours
    if time_passed <= final_time:
        
        for i in range(nsteps):
            g = compute_g(u,D,dt, corr_var)
            u = advance_time(u,g,dt)
            time_passed_list.append(time_passed)

            with open(simu_output_file+'-result-u-with-time-step.csv', 'w+') as f1, open(simu_output_file+'-x-inmicrons.csv', 'w', newline='') as f2, open(simu_output_file+'-time_passed_list_inhrs.csv', 'w', newline='') as f3:
                my_writer1 = csv.writer(f1, delimiter=';')
                my_writer2 = csv.writer(f2, delimiter=';')
                my_writer3 = csv.writer(f3, delimiter=';')

                #my_writer1.writerow(["concentration profiles in weight fraction"])
                my_writer1.writerows(u_list)


                my_writer2.writerow(["x in microns"])
                my_writer2.writerow(x * 10000)

                my_writer3.writerow(["time_passed_inhrs"])
                my_writer3.writerow(time_passed_list)
            
            if surfaceconcentration(u,x)<0.05:
                print("corrosion lifetime criteria met at ", time_passed)

            if stepsdone%50==0:
                print (" stepsdone=%5d, time_passed =%8gh, surfaceconcentration(u,x)=%8g, depthofattack(u,x)=%8g" %
                        (stepsdone,time_passed,surfaceconcentration(u,x),depthofattack(u,x)))
            
            time_at_below = time_passed
            attack = depthofattack(u,x)
    else:
        plt.close(fig1)
    
    l.set_ydata(u) # update data in plot
    return l,

""" Return corest with kp as function of nature of oxide """
def dX_dt(X, t=0):
        return array([b1 / a1 * A1 / X[1], #dWm/dt
                    b1 / a1 * (b1 / a1 * A1 / X[1] - B1), #dWr/dt
                    1 / a1 * b1 / a1 * A1 / X[1] - b1 / a1 * B1, #dW/dt
                    -b1 / a1 * B1]) #dWs/dt
import os

def prepend_line(file_name, line):
    """ Insert given string as a new line at the beginning of a file """
    # define name of temporary dummy file
    dummy_file = file_name + '.bak'
    # open original file in read mode and dummy file in write mode
    with open(file_name, 'r') as read_obj, open(dummy_file, 'w') as write_obj:
        # Write given line to the dummy file
        write_obj.write(line + '\n')
        # Read lines from original file one by one and append them to the dummy file
        for line in read_obj:
            write_obj.write(line)
    # remove original file
    os.remove(file_name)
    # Rename dummy file as the original file
    os.rename(dummy_file, file_name)

def run(corr_var):

    global A1, A2, Cr_act, b1, D, dt, B1
    A1 = a1 ** 2 * (kp0*exp(-Q_kp*1e3/(8.314*(Temp+273.15)))) / 2
    A2 = kp0_ms*exp(-Q_kp_ms*1e3/(8.314*(Temp+273.15)))
    Cr_act=slope_Cract*(1/(Temp+273.15))+const_Cract
    print(Cr_act)
    B1=0
    
    D=D0*exp(-Q_D*1e3/(8.314*(Temp+273.15)))
    dt = final_time*3600/500
    
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

    """plotting"""
    Wm, Wr, W, Ws = X.T
    Metal_loss=Wm/density*10000
    
    file = open(simu_output_file, 'w+', newline='')

    # writing the data into the file
    with file:
        write = csv.writer(file)
        write.writerow(["Wm - mg.cm-2","Wr- mg.cm-2","W- mg.cm-2","Ws-mg.cm-2"])
        write.writerows(X)

        masschangevalues = [
            ["Wm - mg.cm-2","Wr- mg.cm-2","W- mg.cm-2","Ws-mg.cm-2","Metal_Loss"],
            [Wm, Wr, W, Ws, Metal_loss]
        ]

    with open(simu_output_file+'-Masschange_values.csv', 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        # write the data
        writer.writerows(masschangevalues)


    
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

    global time_at_below
    time_at_below = 999999
    global attack
    attack = -1
    # then compute solution and animate
    for i in range(0,10000):
        do_steps(i, fig1, corr_var, nsteps=10)

    plt.show()
    
    df = pd.read_csv(simu_output_file+'-result-u-with-time-step.csv', delimiter=';')
    last_row = df.iloc[-1].values.tolist()
    j = 0
    for i in range(0, len(last_row)):
        if float(last_row[i])>0:
            j = i
            break
    df = df.iloc[:,j+1:]
    df.to_csv(simu_output_file+'-result-u-with-time-step.csv',index=False)
    
    prepend_line(simu_output_file+'-result-u-with-time-step.csv','concentration profiles in weight fraction')
    return time_at_below, attack

if __name__=='__main__':
    
    stepsdone=0
    counter=0

    parser = argparse.ArgumentParser(description='Process simulation (diffusion)')
    parser.add_argument('input_file', help='simulation input file name')
    parser.add_argument('output_file', help='simulation output file name')
    parser.add_argument('--corr_var', help='set corr_var 1 or 2', default=2)

    args = parser.parse_args()    

    # Read simulation input from file
    simu_input_file = args.input_file
    simu_output_file = args.output_file
    simu_input = open(simu_input_file, 'r')
    simu_values = simu_input.readlines()
    simu_input.close()
    corr_var = int(args.corr_var)
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
    density					= float(simu_values[15]) # alloy density in kg/m3
    const_Cract             = float(simu_values[16]) # constant for Cr-activity
    slope_Cract             = float(simu_values[17]) # slope for Cr-activity
    print("corr_var is set to =", corr_var)
    """calculates parabolic constant in weight of metal / Cr mg2cm-4h-1"""
    run(corr_var)
    print("* Calculation completed.")

 


