from tornado.web import Application, RequestHandler
from tornado.ioloop import IOLoop
import json
import pandas as pd
import os
import random
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
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from numpy import *
import pylab as p
import os
import glob
import pandas as pd
from scipy.optimize import curve_fit
from scipy.integrate import solve_ivp
import csv
import numpy as np
from scipy import integrate
import argparse

global stepsdone
global counter
stepdone = 0
counter = 0

# Main Code #
def column(matrix, i):
    return [row[i] for row in matrix]

def compute_g(u,D,dt, b, density, n1, f1, f2, A2, N, Cr_act, corr_var):
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

def advance_time(u,g,dt,u_list):
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

def dY_dt(t, Y, C3_e, eps0, sig0, C3_Q, R, T, h, H_star, C4, K_GP_0, con_p, d, sig_p, Kp, stress):                
    return array([(C3_e*eps0*exp(-C3_Q/(R*T))/(1-Y[2]))*sinh(Y[5]*(1-Y[1])/(sig0*(1-Y[3])*(1-Y[4]))), # dY[0]/dt corresponds to Y[0]=y0 here or =y1 in MatLab
                        (h/Y[5])*(1-Y[1]/H_star)*(C3_e*eps0*exp(-C3_Q/(R*T))/(1-Y[2]))*sinh(Y[5]*(1-Y[1])/(sig0*(1-Y[3])*(1-Y[4]))), #dY[1]/dt
                        C4*((1-Y[2])**2)*(C3_e*eps0*exp(-C3_Q/(R*T))/(1-Y[2]))*sinh(Y[5]*(1-Y[1])/(sig0*(1-Y[3])*(1-Y[4]))), #dY[2]/dt dislocation multiplication
                        (exp(K_GP_0*con_p+sig_p*log(Y[5])))/(d**2*(Y[3]+1e-4)), #dY[3]/dt
                        (Kp/3)*(1-Y[4])**4, #dY[4]/dt
                        stress*(C3_e*eps0*exp(-C3_Q/(R*T))/(1-Y[2]))*sinh(Y[5]*(1-Y[1])/(sig0*(1-Y[3])*(1-Y[4])))]) #dY[5]/dt

def timetocreep1(time,creepstrain):

    """ Computes surface concentration in wt.% """
    idx=(np.abs(np.asarray(creepstrain)-0.01).argmin())
    return time[idx]

def timetocreep2(time,creepstrain):

    """ Computes surface concentration in wt.% """
    idx=(np.abs(np.asarray(creepstrain)-0.02).argmin())
    return time[idx]

def do_steps(i, fig1, density, f1, f2, b, n1, A2, N, Cr_act, simu_output_file, time_passed_list, D, dt, final_time, u, u_list, x, l, corr_var, nsteps=10):
    global time_at_below
    global attack
    time_passed = stepsdone * dt / 3600  # in hours
    if time_passed <= final_time:
        for i in range(nsteps):      
            time_passed = stepsdone * dt / 3600  # in hours
            
            g = compute_g(u, D, dt, b, density, n1, f1, f2, A2, N, Cr_act, corr_var)
            u = advance_time(u, g, dt, u_list)
            time_passed_list.append(time_passed)
            
            if time_passed <=final_time:                
                with open(simu_output_file+'-result-u-with-time-step.csv', 'w+') as file1, open(simu_output_file+'-x-inmicrons.csv', 'w', newline='') as file2, open(simu_output_file+'-time_passed_list_inhrs.csv', 'w', newline='') as file3:
                    my_writer1 = csv.writer(file1, delimiter=';')
                    my_writer2 = csv.writer(file2, delimiter=';')
                    my_writer3 = csv.writer(file3, delimiter=';')

                    my_writer1.writerow(["concentration profiles in weight fraction"])
                    my_writer1.writerows(u_list)

                    my_writer2.writerow(["x in microns"])
                    my_writer2.writerow(x * 10000)

                    my_writer3.writerow(["time_passed_inhrs"])
                    my_writer3.writerow(time_passed_list)
                
                if surfaceconcentration(u)<0.05:
                    print("corrosion lifetime criteria met at ", time_passed)
                    if time_passed< time_at_below:
                        time_at_below = time_passed
                    
                print (" stepsdone=%5d, time_passed =%8gh, surfaceconcentration(u)=%8g, depthofattack(u,x)=%8g" %
                (stepsdone,time_passed,surfaceconcentration(u),depthofattack(u,x)))
                attack = depthofattack(u,x)
    else:
        plt.close(fig1)

    l.set_ydata(u) # update data in plot
    return l,

""" Return corest with kp as function of nature of oxide """
def dX_dt(X, t, a1, b1, A1, B1):
    
    return array([b1 / a1 * A1 / X[1], #dWm/dt
                b1 / a1 * (b1 / a1 * A1 / X[1] - B1), #dWr/dt
                1 / a1 * b1 / a1 * A1 / X[1] - b1 / a1 * B1, #dW/dt
                -b1 / a1 * B1]) #dWs/dt

def execute_creep(simu_output_file, stress, C3_e, eps0, sig0, C3_Q, R, T, h, H_star, C4, K_GP_0, con_p, d, sig_p, Kp):

    """Integrating the ODE using scipy.integrate"""
    sol = solve_ivp(lambda t, Y: dY_dt(t, Y, C3_e, eps0, sig0, C3_Q, R, T, h, H_star, C4, K_GP_0, con_p, d, sig_p, Kp, stress), [0,15863], [0.0, 0.0, 0.0, 0.0, 0.0, stress], method='Radau', rtol=1e-5, atol=1e-5,first_step=1e-3, dense_output=bool)
    creep=sol.y
    time=sol.t
    creepstrain=sol.y[0]
    solution=np.vstack([time,creep]) #stacks time and creep
    creep_solution=np.transpose(solution) #transpose array to save into CSV file

    file = open(simu_output_file+'-result_time_yi.csv', 'w+', newline='')

    #writing the data into the file
    with file:
        write = csv.writer(file)
        write.writerow(["time hrs","y1","y2","y3","y4","y5","y6"])
        #write.writerows(map(lambda x: [x], creep))
        #write.writerows(map(lambda x: [x], time))
        write.writerows(creep_solution)

    #- time to 1% creep strain
	#- time to 2% creep strain
	#- time to creep rupture

    return timetocreep1(time,creepstrain), timetocreep2(time,creepstrain)

def execute_corr(simu_output_file, final_time, Temp, kp0, Q_kp, kp0_ms, Q_kp_ms, a1, b1, a, b, N, cr_ini, D0, Q_D, dt, density, const_Cract, slope_Cract, A1, A2, Cr_act, B1, D, corr_var):
    global stepsdone
    global counter    
    
    stepsdone = 0
    counter = 0
    
    """Integrating the ODE using scipy.integrate"""

    t = linspace(0, 30000, 1000)  # time in h and number of points it is an array
    X0 = array([0.0000001, 0.0000001, 0.0000001, 0.0000001])  # initials conditions
    X, infodict = integrate.odeint(dX_dt, X0, t, args=(a1, b1, A1, B1), full_output=True)

    time_todataframe = pd.DataFrame(t, columns = ['time in hrs'])
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
    global time_at_below
    time_at_below = 999999
    global attack
    attack = -1
    #line_ani = animation.FuncAnimation(fig1, do_steps, range(10000), repeat = False, fargs=(fig1, density, f1, f2, b, n1, A2, N, Cr_act, simu_output_file, time_passed_list, D, dt, final_time, u, u_list, x, l, corr_var))
    #
    for i in range(0,10000):
        do_steps(i, fig1, density, f1, f2, b, n1, A2, N, Cr_act, simu_output_file, time_passed_list, D, dt, final_time, u, u_list, x, l, corr_var, nsteps=10)
    #plt.show()
    print("Answer:", time_at_below)
    print("Depth of Attack:", attack)
    return time_at_below, attack

# For ignoring cross origin issues
class NewHandler(RequestHandler):
     def set_default_headers(self):
          self.set_header("Access-Control-Allow-Origin", "*")
          self.set_header("Access-Control-Allow-Methods", "POST,GET,PUT,OPTIONS")
          self.set_header("Access-Control-Allow-Headers", "Content-Type, Depth, User-Agent, X-File-Size,X-Requested-With, X-Requested-By,If-Modified-Since,X-File-Name,Cache-Control")

class RequestCorrosionRun(RequestHandler):
    def post(self):

        request = json.loads(self.request.body)
        request = json.loads(request)
        thickness = request['thickness']
        temperature = request['temperature']
        endtime = request['endtime']
        material = request['material']
        output_file = request['output_file']
        environment = request['environment']
    
        if environment == "sCO2":
            corr_var = 1
        else:
            corr_var = 2
        result, time_at_below, attack = run_corrosion(output_file, thickness, temperature, endtime, material, environment, corr_var)

        output = {}
        output['output'] = result
        output['time_at_below'] = time_at_below
        output['attack'] = attack
        self.write(output)

class RequestCreepRun(RequestHandler):
    def post(self):
        request = json.loads(self.request.body)
        request = json.loads(request)
        stress = request['stress']
        temperature = request['temperature']
        material = request['material']
        output_file = request['output_file']
        thickness = request['thickness']

        simu_output_file, timetocreep1, timetocreep2 = run_creep(output_file, thickness, temperature, material, stress)
        output = {}
        output['output'] = simu_output_file
        output['timetocreep1'] =timetocreep1
        output['timetocreep2'] =timetocreep2

        self.write(output)

def run_creep(output_file, thickness, temperature, material, stress):
    if material == "740H":
        simu_input_file = "../data/CreepProp_282.txt"
    elif material == "282":
        simu_input_file = "../data/CreepProp_282.txt"
    elif material == "625":
        simu_input_file = "../data/CreepProp_282.txt"
    else:
        print(material, "error!")

    if output_file.strip()=="":
        output_file = "default"

    simu_output_file = "output/"+output_file
    simu_input = open(simu_input_file, 'r')
    simu_values = simu_input.readlines()
    simu_input.close()
    test_pass = False

    print("*stress =", stress)
    print("*material =", material)
    print("*temperature =", temperature)
    print("*thickness =", thickness)

    try:
        temp                    = int(temperature)        
        d                       = float(thickness)  
        stress                  = float(stress)  

        test_pass = True
    
    except:
        temp                    = int(simu_values[0])    
        d                       = float(simu_values[4]) 
        stress                       = float(simu_values[3]) 

    Ts                      = int(simu_values[1]) # GP solvus temperature in K from TC
    R                       = float(simu_values[2]) # gas constant 8.314 J.k-1.mol-1
    a_p                     = float(simu_values[5]) # pre-exp coarsening rate constant h-1
    Q_p                     = float(simu_values[6]) # Activation energy coarsening kJ/mol
    a_sp                    = float(simu_values[7]) # pre-exp depletion of GP microns2 h-1
    Q_sp                    = float(simu_values[8]) # Activation energy depletion of GP kJ/mol
    sig_p                   = float(simu_values[9]) # Stress dependence parameter 1 for depletion of GP
    con_p                   = float(simu_values[10]) # Stress dependence parameter 2 for depletion of GP
    a_ys                    = float(simu_values[11]) #  constant for temperature dependence of yield strength
    b_ys                    = float(simu_values[12]) #  constant for temperature dependence of yield strength
    c_ys                    = float(simu_values[13]) #  constant for temperature dependence of yield strength
    d_ys                    = float(simu_values[14]) #  constant for temperature dependence of yield strength
    e_ys                    = float(simu_values[15]) #  constant for temperature dependence of yield strength
    a_E                     = float(simu_values[16])  # constant for temperature dependence of Young's modulus
    b_E                    = float(simu_values[17])  # constant for temperature dependence of Young's modulus
    c_E                    = float(simu_values[18])  # constant for temperature dependence of Young's modulus
    d_E                    = float(simu_values[19])  # constant for temperature dependence of Young's modulus
    e_E                    = float(simu_values[20])  #  constant for temperature dependence of Young's modulus                     =
    f_E                    = float(simu_values[21])  #  constant for temperature dependence of Young's modulus
    a_s                    = float(simu_values[22])  #  constant for temperature dependence of strengthening phase vol. fraction
    b_s                    = float(simu_values[23])  #  constant for temperature dependence of strengthening phase vol. fraction
    c_s                    = float(simu_values[24])  #  constant for temperature dependence of strengthening phase vol. fraction
    d_s                    = float(simu_values[25])  # constant for temperature dependence of strengthening phase vol. fraction
    par1                   = float(simu_values[26])  # parameter1 related to alloy
    par2                   = float(simu_values[27])  # parameter2 related to alloy
    par3                   = float(simu_values[28])  # parameter3 related to alloy
    par4                   = float(simu_values[29])  # parameter4 related to alloy
    par5                   = float(simu_values[30])  # parameter4 related to alloy
    C1                     = float(simu_values[31]) # Back stress constant 1 MPa
    C2                     = float(simu_values[32]) # Back stress constant 2 kJ/mol
    C3_Q                    = float(simu_values[33]) # Activation energy for intrinsic strain rate kJ/mol
    C4                      = float(simu_values[34]) # Dislocation damage constant

    """calculates values from input file"""

    T =temp+273.15 #Temperature in Kelvin
    E=1e3*(f_E*temp**5+e_E*temp**4+d_E*temp**3+c_E*temp**2+b_E*temp+a_E) #Young Modulus in MPa
    ys = e_ys* (temp) ** 4 + d_ys * temp ** 3 + c_ys * temp ** 2 + b_ys * temp + a_ys  # Yield strength MPa
    vol_gp=d_s*temp**3+c_s*temp**2+b_s*temp+a_s # vol. fraction of gamma prime

    Kp=a_p*exp(-Q_p*1e3/(R*T)) #coarsening constant in h-1

    K_GP_0 = a_sp * exp(-Q_sp*1e3 / R / T)  # depth of GP depletion in 740 in microns

    sig0=C1*(1-exp(-C2/(R*Ts))*(Ts/T-1))
    eps0=2*sqrt(vol_gp)*(1-vol_gp)*(sqrt(pi/4)-sqrt(vol_gp))
    h=E*vol_gp
    H_star=2*vol_gp/(1+vol_gp)

    C3_e = exp(par1 + par2 * log(stress / ys) + par3 * sinh(stress / T) + par4 * log(1 / T) + par5 * (1 / T))

    timetocreep1, timetocreep2 = execute_creep(simu_output_file, stress, C3_e, eps0, sig0, C3_Q, R, T, h, H_star, C4, K_GP_0, con_p, d, sig_p, Kp)
    return simu_output_file, timetocreep1, timetocreep2
    
def run_corrosion(output_file, thickness, temperature, endtime, material, environment, corr_var):
    if material == "740H":
        simu_input_file = "../data/CorrProp_740H.txt"
    elif material == "282":
        simu_input_file = "../data/CorrProp_282.txt"
    elif material == "625":
        simu_input_file = "../data/CorrProp_625.txt"
    else:
        print(material, "error!")

    if output_file.strip()=="":
        output_file = "default"
    simu_output_file = "output/"+output_file
    simu_input = open(simu_input_file, 'r')
    simu_values = simu_input.readlines()
    simu_input.close()
    test_pass = False

    print("*endtime =", endtime)
    print("*temperature =", temperature)
    print("*thickness =", thickness)

    try:
        final_time              = int(endtime)    # hours
        Temp                    = float(temperature)  # Temperature in C entered globally by user      
        b                       = float(thickness)  # Half component thickness in cm
        test_pass = True
    
    except:
        print("** Warning: Please enter correct values for component thickness, temperature, and end time")
        final_time              = int(simu_values[0])    # hours
        Temp                    = int(simu_values[1])    # Temperature in C entered globally by user
        b                       = float(simu_values[9])  # Half component thickness in cm


    print("*final_time =", final_time)
    print("*Temp =", Temp)
    print("*thickness =", b)
    

    if test_pass!=True:
        final_time              = int(simu_values[0])    # hours
        Temp                    = int(simu_values[1])    # Temperature in C entered globally by user
        b                       = float(simu_values[9])  # Half component thickness in cm    
    
    kp0                     = float(simu_values[2])  # Pre-exponential:Parabolic rate constant for Cr2O3 in mg2cm-4h-1
    Q_kp					= float(simu_values[3]) # Activation energy:Parabolic rate constant for Cr2O3 in
    kp0_ms                  = float(simu_values[4])  # Pre-exponential:Parabolic rate constant for Cr dissolution in KCl-MgCl2 mg2cm-4h-1
    Q_kp_ms					= float(simu_values[5]) # Activation energy:Parabolic rate constant for Cr dissolution in KCl-MgCl2 kJ/mol
    a1                      = float(simu_values[6])  # 2.17, stoichiometric constant for Cr2O3
    b1                      = float(simu_values[7])  # 3.17, , stoichiometric constant for Cr2O3
    a                       = float(simu_values[8])  
    N                       = int(simu_values[10])   # number of nodes
    cr_ini                  = float(simu_values[11]) # initial Cr concentration in wt%
    D0                      = float(simu_values[12]) # Pre-exponential:Diffusion coefficient of Cr in m2/s
    Q_D                     = float(simu_values[13]) # Pre-exponential:Activation energy Diffusion coefficient of Cr in kJ/mol
    dt                      = float(simu_values[14]) # seconds
    density					= float(simu_values[15]) # alloy density in kg/m3
    const_Cract             = float(simu_values[16]) # constant for Cr-activity
    slope_Cract             = float(simu_values[17]) # slope for Cr-activity
    
    print("**a1=",a1)
    print("**b1=",b1)

    """calculates parabolic constant in weight of metal / Cr mg2cm-4h-1"""
    A1 = a1 ** 2 * (kp0*exp(-Q_kp*1e3/(8.314*(Temp+273.15)))) / 2
    A2 = kp0_ms*exp(-Q_kp_ms*1e3/(8.314*(Temp+273.15)))
    Cr_act=slope_Cract*(1/(Temp+273.15))+const_Cract
    print(Cr_act)
    B1=0
    D=D0*exp(-Q_D*1e3/(8.314*(Temp+273.15)))

    time_at_below, attack = execute_corr(simu_output_file, final_time, Temp, kp0, Q_kp, kp0_ms, Q_kp_ms, a1, b1, a, b, N, cr_ini, D0, Q_D, dt, density, const_Cract, slope_Cract, A1, A2, Cr_act, B1, D, corr_var)
    
    print("* Calculation completed.")
    print(b, Temp, final_time, material, environment, attack)
    return simu_output_file, time_at_below, attack

def make_app():
    urls = [
        ("/corrosion-run/", RequestCorrosionRun),
        ("/creep-run/", RequestCreepRun),
    ]
    return Application(urls, debug=True)

if __name__ == '__main__':

    print("* Ready")
    
    app = make_app()
    app.listen(8888)
    IOLoop.instance().start()
    
'''
# POST test
curl -d '{"input":"test"}' -H "Content-Type: application/json" -X POST http://localhost:8888/test/
  '''
