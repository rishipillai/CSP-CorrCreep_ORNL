"""
Created in 2021
@author: 6mr
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
from scipy.integrate import solve_ivp
import csv
import argparse

""" Returns creep"""


# Y is an array containing Y[0], Y[1], Y[2], Y[3], Y[4] and Y[5] which corresponds respectively to y1, y2, y3, y4, y5 and y6 in the Matlab code
# The following function will calculate the derivative against time, t, for all yi in the array
def dY_dt(t, Y):
    return array(
        [(C3_e * eps0 * exp(-C3_Q / (R * T)) / (1 - Y[2])) * sinh(Y[5] * (1 - Y[1]) / (sig0 * (1 - Y[3]) * (1 - Y[4]))),
         # dY[0]/dt corresponds to Y[0]=y0 here or =y1 in MatLab
         (h / Y[5]) * (1 - Y[1] / H_star) * (C3_e * eps0 * exp(-C3_Q / (R * T)) / (1 - Y[2])) * sinh(
             Y[5] * (1 - Y[1]) / (sig0 * (1 - Y[3]) * (1 - Y[4]))),  # dY[1]/dt
         C4 * ((1 - Y[2]) ** 2) * (C3_e * eps0 * exp(-C3_Q / (R * T)) / (1 - Y[2])) * sinh(
             Y[5] * (1 - Y[1]) / (sig0 * (1 - Y[3]) * (1 - Y[4]))),  # dY[2]/dt dislocation multiplication
         (K_GP) / (d ** 2 * (Y[3] + 1e-4)),  # dY[3]/dt
         (Kp / 3) * (1 - Y[4]) ** 4,  # dY[4]/dt
         stress * (C3_e * eps0 * exp(-C3_Q / (R * T)) / (1 - Y[2])) * sinh(
             Y[5] * (1 - Y[1]) / (sig0 * (1 - Y[3]) * (1 - Y[4])))])  # dY[5]/dt


def lmpstress(stress, temp):
    LMP = LMP_1 + LMP_2 * log10(stress) + LMP_3 * (log10(stress)) ** 2 + LMP_4 * (log10(stress)) ** 3
    rupturetime = 10 ** ((LMP / (temp + 273.15)) - 20)
    return rupturetime


def lmpgp(vol_gp, temp):
    LMP_GP = LMP_GP1 + LMP_GP2 * log10(vol_gp) + LMP_GP3 * log10(vol_gp) ** 2 + LMP_GP4 * log10(
        vol_gp) ** 3 + LMP_GP5 * log10(vol_gp) ** 4
    rupturetime_GP = 10 ** ((LMP_GP / (temp + 273.15)) - 20)
    return rupturetime_GP


def timetocreep1(time, creepstrain):
    """ Computes surface concentration in wt.% """
    idx = (np.abs(np.asarray(creepstrain) - 0.01).argmin())
    return time[idx]


def timetocreep2(time, creepstrain):
    """ Computes surface concentration in wt.% """
    idx = (np.abs(np.asarray(creepstrain) - 0.02).argmin())
    return time[idx]


def run():
    """Integrating the ODE using scipy.integrate"""
    sol = solve_ivp(dY_dt, [0, 15863], [0.0, 0.0, 0.0, 0.0, 0.0, stress], method='Radau', rtol=1e-5, atol=1e-5,
                    first_step=1e-3, dense_output=bool)
    creep = sol.y
    time = sol.t
    creepstrain = sol.y[0]

    solution = np.vstack([time, creep])  # stacks time and creep
    creep_solution = np.transpose(solution)  # transpose array to save into CSV file

    file = open(args.output_file + '-result_time_yi.csv', 'w+', newline='')

    # writing the data into the file
    with file:
        write = csv.writer(file)
        write.writerow(["time hrs", "y1", "y2", "y3", "y4", "y5", "y6"])
        # write.writerows(map(lambda x: [x], creep))
        # write.writerows(map(lambda x: [x], time))
        write.writerows(creep_solution)

        print(timetocreep1(time, creepstrain), timetocreep2(time, creepstrain))

        lifetimevalues = [
            ['Time to 1% creep (h)', 'Time to 2% creep (h)', 'Rupture time (h)', 'Rupture time based on gp (h)'],
            [timetocreep1(time, creepstrain), timetocreep2(time, creepstrain), lmpstress(stress, temp),
             lmpgp(vol_gp, temp)]
        ]

    with open('lifetime_values.csv', 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        # write the data
        writer.writerows(lifetimevalues)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process simulation (creep)')
    parser.add_argument('input_file', help='simulation input file name')
    parser.add_argument('output_file', help='simulation output file name')

    args = parser.parse_args()

    """Read simulation input from file"""
    simu_input_file = args.input_file
    simu_output_file = args.output_file
    simu_input = open(simu_input_file, 'r')
    simu_values = simu_input.readlines()
    simu_input.close()

    temp = int(simu_values[0])  # Temperature Celcius
    Ts = int(simu_values[1])  # GP solvus temperature in K from TC
    R = float(simu_values[2])  # gas constant 8.314 J.k-1.mol-1
    stress = float(simu_values[3])  # Stress in MPa
    d = float(simu_values[4])  # component thickness in microns
    a_p = float(simu_values[5])  # pre-exp coarsening rate constant h-1
    Q_p = float(simu_values[6])  # Activation energy coarsening kJ/mol
    a_sp = float(simu_values[7])  # pre-exp depletion of GP microns2 h-1
    Q_sp = float(simu_values[8])  # Activation energy depletion of GP kJ/mol
    a_ys = float(simu_values[9])  # constant for temperature dependence of yield strength
    b_ys = float(simu_values[10])  # constant for temperature dependence of yield strength
    c_ys = float(simu_values[11])  # constant for temperature dependence of yield strength
    d_ys = float(simu_values[12])  # constant for temperature dependence of yield strength
    e_ys = float(simu_values[13])  # constant for temperature dependence of yield strength
    a_E = float(simu_values[14])  # constant for temperature dependence of Young's modulus
    b_E = float(simu_values[15])  # constant for temperature dependence of Young's modulus
    c_E = float(simu_values[16])  # constant for temperature dependence of Young's modulus
    d_E = float(simu_values[17])  # constant for temperature dependence of Young's modulus
    e_E = float(simu_values[18])  # constant for temperature dependence of Young's modulus                     =
    f_E = float(simu_values[19])  # constant for temperature dependence of Young's modulus
    a_s = float(simu_values[20])  # constant for temperature dependence of strengthening phase vol. fraction
    b_s = float(simu_values[21])  # constant for temperature dependence of strengthening phase vol. fraction
    c_s = float(simu_values[22])  # constant for temperature dependence of strengthening phase vol. fraction
    d_s = float(simu_values[23])  # constant for temperature dependence of strengthening phase vol. fraction
    par1 = float(simu_values[24])  # parameter1 related to alloy
    par2 = float(simu_values[25])  # parameter2 related to alloy
    par3 = float(simu_values[26])  # parameter3 related to alloy
    par4 = float(simu_values[27])  # parameter4 related to alloy
    par5 = float(simu_values[28])  # parameter4 related to alloy

    C1 = float(simu_values[29])  # Back stress constant 1 MPa
    C2 = float(simu_values[30])  # Back stress constant 2 kJ/mol

    C3_Q = float(simu_values[31])  # Activation energy for intrinsic strain rate kJ/mol
    C4 = float(simu_values[32])  # Dislocation damage constant

    LMP_1 = float(simu_values[33])  # Stress dependence of LMP - parameter 1
    LMP_2 = float(simu_values[34])  # Stress dependence of LMP - parameter 2
    LMP_3 = float(simu_values[35])  # Stress dependence of LMP - parameter 3
    LMP_4 = float(simu_values[36])  # Stress dependence of LMP - parameter 4

    LMP_GP1 = float(simu_values[37])  # GP dependence of LMP - parameter 1
    LMP_GP2 = float(simu_values[38])  # GP dependence of LMP - parameter 2
    LMP_GP3 = float(simu_values[39])  # GP dependence of LMP - parameter 3
    LMP_GP4 = float(simu_values[40])  # GP dependence of LMP - parameter 4
    LMP_GP5 = float(simu_values[41])  # GP dependence of LMP - pa

    """calculates values from input file"""

    T = temp + 273.15  # Temperature in Kelvin
    E = 1e3 * (f_E * temp ** 5 + e_E * temp ** 4 + d_E * temp ** 3 + c_E * temp ** 2 + b_E * temp + a_E)  # Young Modulus in MPa
    ys = e_ys * (temp) ** 4 + d_ys * temp ** 3 + c_ys * temp ** 2 + b_ys * temp + a_ys  # Yield strength MPa
    vol_gp = d_s * temp ** 3 + c_s * temp ** 2 + b_s * temp + a_s  # vol. fraction of gamma prime

    Kp = a_p * exp(-Q_p * 1e3 / (R * T))  # coarsening constant in h-1

    K_GP = exp(a_sp) * exp(-Q_sp * 1e3 / R / T)  # depth of GP depletion in 740 in microns

    sig0 = C1 * (1 - exp(-C2 / (R * Ts) * (Ts / T - 1)))
    eps0 = 2 * sqrt(vol_gp) * (1 - vol_gp) * (sqrt(pi / 4) - sqrt(vol_gp))
    h = E * vol_gp
    H_star = 2 * vol_gp / (1 + vol_gp)

    C3_e = exp(par1 + par2 * (1 / T) + par3 * sinh(stress / T) + par4 * log(1 / T) + par5 * log(stress / ys))

    run()
    print("* Calculation completed.")
