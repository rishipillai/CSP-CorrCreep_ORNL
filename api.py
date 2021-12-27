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
import creep
import diffusion

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
        simu_input_file = "data/CreepProp_740H.txt"
    elif material == "282":
        simu_input_file = "data/CreepProp_282.txt"
    elif material == "625":
        simu_input_file = "data/CreepProp_625.txt"
    else:
        print(material, "error!")

    if output_file.strip()=="":
        output_file = "default"

    creep.simu_output_file = "output/"+output_file
    simu_input = open(simu_input_file, 'r')
    simu_values = simu_input.readlines()
    simu_input.close()
    test_pass = False

    print("*stress =", stress)
    print("*material =", material)
    print("*temperature =", temperature)
    print("*thickness =", thickness)

    try:
        creep.temp                    = int(temperature)        
        creep.d                       = float(thickness)  
        creep.stress                  = float(stress)  

        test_pass = True
    
    except:
        creep.temp                    = int(simu_values[0])    
        creep.d                       = float(simu_values[4]) 
        creep.stress                       = float(simu_values[3]) 

    creep.Ts                      = int(simu_values[1]) # GP solvus temperature in K from TC
    creep.R                       = float(simu_values[2]) # gas constant 8.314 J.k-1.mol-1
    creep.a_p                     = float(simu_values[5]) # pre-exp coarsening rate constant h-1
    creep.Q_p                     = float(simu_values[6]) # Activation energy coarsening kJ/mol
    creep.a_sp                    = float(simu_values[7]) # pre-exp depletion of GP microns2 h-1
    creep.Q_sp                    = float(simu_values[8]) # Activation energy depletion of GP kJ/mol
    creep.a_ys                    = float(simu_values[9]) #  constant for temperature dependence of yield strength
    creep.b_ys                    = float(simu_values[10]) #  constant for temperature dependence of yield strength
    creep.c_ys                    = float(simu_values[11]) #  constant for temperature dependence of yield strength
    creep.d_ys                    = float(simu_values[12]) #  constant for temperature dependence of yield strength
    creep.e_ys                    = float(simu_values[13]) #  constant for temperature dependence of yield strength
    creep.a_E                     = float(simu_values[14])  # constant for temperature dependence of Young's modulus
    creep.b_E                    = float(simu_values[15])  # constant for temperature dependence of Young's modulus
    creep.c_E                    = float(simu_values[16])  # constant for temperature dependence of Young's modulus
    creep.d_E                    = float(simu_values[17])  # constant for temperature dependence of Young's modulus
    creep.e_E                    = float(simu_values[18])  #  constant for temperature dependence of Young's modulus                     =
    creep.f_E                    = float(simu_values[19])  #  constant for temperature dependence of Young's modulus
    creep.a_s                    = float(simu_values[20])  #  constant for temperature dependence of strengthening phase vol. fraction
    creep.b_s                    = float(simu_values[21])  #  constant for temperature dependence of strengthening phase vol. fraction
    creep.c_s                    = float(simu_values[22])  #  constant for temperature dependence of strengthening phase vol. fraction
    creep.d_s                    = float(simu_values[23])  # constant for temperature dependence of strengthening phase vol. fraction
    creep.par1                   = float(simu_values[24])  # parameter1 related to alloy
    creep.par2                   = float(simu_values[25])  # parameter2 related to alloy
    creep.par3                   = float(simu_values[26])  # parameter3 related to alloy
    creep.par4                   = float(simu_values[27])  # parameter4 related to alloy
    creep.par5                   = float(simu_values[28])  # parameter4 related to alloy

    creep.C1                     = float(simu_values[29]) # Back stress constant 1 MPa
    creep.C2                     = float(simu_values[30]) # Back stress constant 2 kJ/mol
    creep.C3_Q                    = float(simu_values[31]) # Activation energy for intrinsic strain rate kJ/mol
    creep.C4                      = float(simu_values[32]) # Dislocation damage constant

    creep.LMP_1                   = float(simu_values[33]) # Stress dependence of LMP - parameter 1
    creep.LMP_2                   = float(simu_values[34]) # Stress dependence of LMP - parameter 2
    creep.LMP_3                   = float(simu_values[35])  # Stress dependence of LMP - parameter 3
    creep.LMP_4                   = float(simu_values[36])  # Stress dependence of LMP - parameter 4

    creep.LMP_GP1                   = float(simu_values[37]) # GP dependence of LMP - parameter 1
    creep.LMP_GP2                   = float(simu_values[38]) # GP dependence of LMP - parameter 2
    creep.LMP_GP3                   = float(simu_values[39]) # GP dependence of LMP - parameter 3
    creep.LMP_GP4                   = float(simu_values[40]) # GP dependence of LMP - parameter 4
    creep.LMP_GP5                   = float(simu_values[41])  # GP dependence of LMP - pa
        
    """calculates values from input file"""

    timetocreep1, timetocreep2 = creep.run()

    return creep.simu_output_file, timetocreep1, timetocreep2
    
def run_corrosion(output_file, thickness, temperature, endtime, material, environment, corr_var):
    if material == "740H":
        simu_input_file = "data/CorrProp_740H.txt"
    elif material == "282":
        simu_input_file = "data/CorrProp_282.txt"
    elif material == "625":
        simu_input_file = "data/CorrProp_625.txt"
    else:
        print(material, "error!")

    if output_file.strip()=="":
        output_file = "default"
    diffusion.stepsdone=0
    diffusion.counter=0
    diffusion.simu_output_file = "output/"+output_file
    simu_input = open(simu_input_file, 'r')
    simu_values = simu_input.readlines()
    simu_input.close()
    test_pass = False

    print("*endtime =", endtime)
    print("*temperature =", temperature)
    print("*thickness =", thickness)

    try:
        diffusion.final_time              = int(endtime)    # hours
        diffusion.Temp                    = float(temperature)  # Temperature in C entered globally by user      
        diffusion.b                       = float(thickness)  # Half component thickness in cm
        diffusion.test_pass = True
    
    except:
        print("** Warning: Please enter correct values for component thickness, temperature, and end time")
        diffusion.final_time              = int(simu_values[0])    # hours
        diffusion.Temp                    = int(simu_values[1])    # Temperature in C entered globally by user
        diffusion.b                       = float(simu_values[9])  # Half component thickness in cm


    print("*final_time =", diffusion.final_time)
    print("*Temp =", diffusion.Temp)
    print("*thickness =", diffusion.b)
    

    if test_pass!=True:
        diffusion.final_time              = int(simu_values[0])    # hours
        diffusion.Temp                    = int(simu_values[1])    # Temperature in C entered globally by user
        diffusion.b                       = float(simu_values[9])  # Half component thickness in cm    
    
    diffusion.kp0                     = float(simu_values[2])  # Pre-exponential:Parabolic rate constant for Cr2O3 in mg2cm-4h-1
    diffusion.Q_kp					= float(simu_values[3]) # Activation energy:Parabolic rate constant for Cr2O3 in
    diffusion.kp0_ms                  = float(simu_values[4])  # Pre-exponential:Parabolic rate constant for Cr dissolution in KCl-MgCl2 mg2cm-4h-1
    diffusion.Q_kp_ms					= float(simu_values[5]) # Activation energy:Parabolic rate constant for Cr dissolution in KCl-MgCl2 kJ/mol
    diffusion.a1                      = float(simu_values[6])  # 2.17, stoichiometric constant for Cr2O3
    diffusion.b1                      = float(simu_values[7])  # 3.17, , stoichiometric constant for Cr2O3
    diffusion.a                       = float(simu_values[8])  
    diffusion.N                       = int(simu_values[10])   # number of nodes
    diffusion.cr_ini                  = float(simu_values[11]) # initial Cr concentration in wt%
    diffusion.D0                      = float(simu_values[12]) # Pre-exponential:Diffusion coefficient of Cr in m2/s
    diffusion.Q_D                     = float(simu_values[13]) # Pre-exponential:Activation energy Diffusion coefficient of Cr in kJ/mol
    #dt                      = float(simu_values[14]) # seconds
    diffusion.density					= float(simu_values[15]) # alloy density in kg/m3
    diffusion.const_Cract             = float(simu_values[16]) # constant for Cr-activity
    diffusion.slope_Cract             = float(simu_values[17]) # slope for Cr-activity

    time_at_below, attack = diffusion.run()
    
    print("* Calculation completed.")
    print(diffusion.b, diffusion.Temp, diffusion.final_time, material, environment, attack)
    return diffusion.simu_output_file, time_at_below, attack

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
