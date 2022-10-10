"""
Functions to help process TCP Data AT DIFFERENT SLIDING SPEEDS
Functions Present:
- Getting intact columns from processed experiments
- Plotting fitted exponential graphs of decay of intact molecules
- Get friction values from fc.ave.dump file
- Processing fc_ave.dump file to get average shear/normal stress and average mu

"""

import pandas as pd
import matplotlib.pyplot as plt
import scipy
from scipy import optimize
import numpy as np
import statistics
import os
import os.path
import sys
import glob
import sys
from numpy import mean as m

Speeds = ["1ms", "20ms", "50ms"]
Pressures = ['3GPa', '4GPa', '5GPa']
speed = "1ms"
pressures = ['3GPa', '4GPa', '5GPa']

def get_intact_columns_constant_temperature(Temperature, Pressures):

    """
    speed: Speed directory that you ran the experiments at
    :pressures: Different pressures that you ran the experiments at with constant speed
    :return: Dataframe with concatenated and renamed columns from each of the experiments at the input pressures and speed
    """

    #Makes a dataframe from the intact molecuels csv of the first speed/temperature/pressure experiment in your defined list
    All_Dataframe = pd.read_csv('F:/PhD/TCPDecompositionExperiments/Completed/Fe3O4/{}/1GPa/'
                                'processed/Fe2O3-200-iso-octane_IntactMols_0.3'.format(Temperature), sep='\t')

    #Take the correct columns from the above datafram using iloc
    Big_Dataframe_NotProcessed = All_Dataframe.iloc[:, [0, 2]]

    #Rename the column to the pressure of the first experiment
    Big_Dataframe = Big_Dataframe_NotProcessed.rename(columns={'Intact_molecules_noomit' : '1GPa'})

    #Quick for loop to get the rest of the columns from the remaining

    Pressures = Pressures
    for P in Pressures:
        Dataframe = pd.read_csv('F:/PhD/TCPDecompositionExperiments/Completed/Fe3O4/{}/{}/'
                                'processed/Fe2O3-200-iso-octane_IntactMols_0.3'.format(Temperature, P), sep='\t')
        UnnamedBig_DataframeP = Dataframe.iloc[:, [0, 2]]
        #print(Big_DataframeP)
        Big_DataframeP = UnnamedBig_DataframeP.rename(columns= {'Timestep' : 'Timestep_{}'.format(P),
                                                                'Intact_molecules_noomit' : '{}'.format(P)})
        Big_Dataframe = pd.concat([Big_Dataframe, Big_DataframeP], axis =1)

    # Using .dropna from pandas library to get rid of rows with NaNa (which mean the simulation ran out of time
    Big_Dataframe = Big_Dataframe.dropna(axis=0)

    return Big_Dataframe

def get_intact_columns_constant_pressure(Pressure, Temperatures):

    """
    speed: Speed directory that you ran the experiments at
    :pressures: Different pressures that you ran the experiments at with constant speed
    :return: Dataframe with concatenated and renamed columns from each of the experiments at the input temperatures
    """

    # Makes a dataframe from the intact molecuels csv of the first speed/temperature/pressure experiment in your defined list
    All_Dataframe = pd.read_csv('F:/PhD/TCPDecompositionExperiments/Completed/Fe3O4/400K/{}/'
                                'processed/Fe2O3-200-iso-octane_IntactMols_0.3'.format(Pressure), sep='\t')

    # Take the correct columns from the above datafram using iloc
    Big_Dataframe_NotProcessed = All_Dataframe.iloc[:, [0, 2]]

    # Rename the column to the pressure of the first experiment
    Big_Dataframe = Big_Dataframe_NotProcessed.rename(columns={'Intact_molecules_noomit': '400K'})

    # Quick for loop to get the rest of the columns from the remaining

    for T in Temperatures:
        Dataframe = pd.read_csv('F:/PhD/TCPDecompositionExperiments/Completed/Fe3O4/{}/{}/'
                                'processed/Fe2O3-200-iso-octane_IntactMols_0.3'.format(T, Pressure), sep='\t')
        UnnamedBig_DataframeP = Dataframe.iloc[:, [0, 2]]
        # print(Big_DataframeP)
        Big_DataframeP = UnnamedBig_DataframeP.rename(columns={'Timestep': 'Timestep_{}'.format(T),
                                                               'Intact_molecules_noomit': '{}'.format(T)})
        Big_Dataframe = pd.concat([Big_Dataframe, Big_DataframeP], axis=1)

    # Using .dropna from pandas library to get rid of rows with NaNa (which mean the simulation ran out of time
    Big_Dataframe = Big_Dataframe.dropna(axis=0)

    return Big_Dataframe

def get_fitted_plots_constant_temperature(BigDataframe, speed, Temperature):
    MoleculesIntact = BigDataframe

    Timestep_List = MoleculesIntact['Timestep'].to_list()
    OneGPa = MoleculesIntact['1GPa'].to_list()
    TwoGPa = MoleculesIntact['2GPa'].to_list()
    ThreeGPa = MoleculesIntact['3GPa'].to_list()
    FourGPa = MoleculesIntact['4GPa'].to_list()
    FiveGPa = MoleculesIntact['5GPa'].to_list()

    # Setting up time variable for plots later on
    Time = []
    for x in Timestep_List:
        Time.append(((x - 400000) / int(4 * 10 ** 6)))

    def monoExp(x, m, t):
        return m * np.exp(-t * x)

    # perform the fit
    p0 = (48, .1)  # start with values near those we expect
    sigma = np.ones(len(Time))
    sigma[[0]] = 48
    params1GPa, cv = scipy.optimize.curve_fit(monoExp, Time, OneGPa, p0, sigma=sigma)
    m, t = params1GPa
    #print(params1GPa)

    def monoExpOne(x, a, b):
        return a * np.exp(-b * x)

    # perform the fit
    p01 = (48, .1)  # start with values near those we expect
    sigma = np.ones(len(Time))
    sigma[0] = 48
    params2GPa, cv = scipy.optimize.curve_fit(monoExpOne, Time, TwoGPa, p01, sigma=sigma)
    a, b = params2GPa
    #print(params2GPa)

    def monoExpTwo(x, c, d):
        return c * np.exp(-d * x)

    # perform the fit
    p02 = (48, .1)  # start with values near those we expect
    sigma = np.ones(len(Time))
    sigma[0] = 48
    params3GPa, cv = scipy.optimize.curve_fit(monoExpTwo, Time, ThreeGPa, p02, sigma=sigma)
    c, d = params3GPa
    #print(params3GPa)

    def monoExpThree(x, e, f):
        return e * np.exp(-f * x)

    # perform the fit
    p03 = (48, .1)  # start with values near those we expect
    sigma = np.ones(len(Time))
    sigma[0] = 48
    params4GPa, cv = scipy.optimize.curve_fit(monoExpThree, Time, FourGPa, p03, sigma=sigma)
    e, f = params4GPa
    #print(params4GPa)

    def monoExpFour(x, g, h):
        return g * np.exp(-h * x)

    # perform the fit
    p03 = (100, 2)  # start with values near those we expect
    sigma = np.ones(len(Time))
    sigma[0] = 48
    params5GPa, cv = scipy.optimize.curve_fit(monoExpFour, Time, FiveGPa, p03, sigma=sigma)
    g, h = params5GPa
    #print(params5GPa)

    fitted_functionOneGPa = []
    fitted_functionTwoGPa = []
    fitted_functionThreeGPa = []
    fitted_functionFourGPa = []
    fitted_functionFiveGPa = []

    for x in Time:
        fitted_functionOneGPa.append(m * np.exp(-t * x))

    for x in Time:
        fitted_functionTwoGPa.append(a * np.exp(-b * x))

    for x in Time:
        fitted_functionThreeGPa.append(c * np.exp(-d * x))

    for x in Time:
        fitted_functionFourGPa.append(e * np.exp(-f * x))

    for x in Time:
        fitted_functionFiveGPa.append(g * np.exp(-h * x))

    # plot the results
    plt.plot(Time, OneGPa, '--', label="1GPa", color='blue')
    plt.plot(Time, fitted_functionOneGPa, color='blue')
    plt.plot(Time, TwoGPa, '--', label="2GPa", color='red')
    plt.plot(Time, fitted_functionTwoGPa, color='red')
    plt.plot(Time, ThreeGPa, '--', label="3GPa", color='orange')
    plt.plot(Time, fitted_functionThreeGPa, color='orange')
    plt.plot(Time, FourGPa, '--', label="4GPa", color='green')
    plt.plot(Time, fitted_functionFourGPa, color='green')
    plt.plot(Time, FiveGPa, '--', label="5GPa", color='purple')
    plt.plot(Time, fitted_functionFiveGPa, color='purple')
    plt.title("TCP Graphs With Exponential Fits {} at {}".format(Temperature, speed))
    plt.ylabel('Intact Molecules')
    plt.legend()
    plt.xlabel('Time (ns)')
    plt.savefig('F:/PhD/TCPDecompositionExperiments/Completed/Fe3O4/'
               'IntactMoleculeComparison_Pressure_{}.png'.format(Temperature))
    plt.close()
    #plt.show()

    return Timestep_List, fitted_functionOneGPa, fitted_functionTwoGPa, fitted_functionThreeGPa, fitted_functionFourGPa, fitted_functionFiveGPa

def get_fitted_plots_constant_pressure(BigDataframe, speed, Pressure):
    MoleculesIntact = BigDataframe

    Timestep_List = MoleculesIntact['Timestep'].to_list()
    OneGPa = MoleculesIntact['300K'].to_list()
    TwoGPa = MoleculesIntact['400K'].to_list()
    ThreeGPa = MoleculesIntact['500K'].to_list()
    FourGPa = MoleculesIntact['600K'].to_list()
    FiveGPa = MoleculesIntact['700K'].to_list()

    # Setting up time variable for plots later on
    Time = []
    for x in Timestep_List:
        Time.append(((x - 400000) / int(4 * 10 ** 6)))

    def monoExp(x, m, t):
        return m * np.exp(-t * x)

    # perform the fit
    p0 = (48, .1)  # start with values near those we expect
    sigma = np.ones(len(Time))
    sigma[[0, -300]] = .2
    params1GPa, cv = scipy.optimize.curve_fit(monoExp, Time, OneGPa, p0, sigma=sigma)
    m, t = params1GPa
    #print(params1GPa)

    def monoExpOne(x, a, b):
        return a * np.exp(-b * x)

    # perform the fit
    p01 = (48, .1)  # start with values near those we expect
    sigma = np.ones(len(Time))
    sigma[[0, -300]] = .2
    params2GPa, cv = scipy.optimize.curve_fit(monoExpOne, Time, TwoGPa, p01, sigma=sigma)
    a, b = params2GPa
    #print(params2GPa)

    def monoExpTwo(x, c, d):
        return c * np.exp(-d * x)

    # perform the fit
    p02 = (48, .1)  # start with values near those we expect
    sigma = np.ones(len(Time))
    sigma[[0, -300]] = .4
    params3GPa, cv = scipy.optimize.curve_fit(monoExpTwo, Time, ThreeGPa, p02, sigma=sigma)
    c, d = params3GPa
    #print(params3GPa)

    def monoExpThree(x, e, f):
        return e * np.exp(-f * x)

    # perform the fit
    p03 = (48, .1)  # start with values near those we expect
    sigma = np.ones(len(Time))
    sigma[[0]] = .8
    params4GPa, cv = scipy.optimize.curve_fit(monoExpThree, Time, FourGPa, p03, sigma=sigma)
    e, f = params4GPa
    #print(params4GPa)

    def monoExpFour(x, g, h):
        return g * np.exp(-h * x)

    # perform the fit
    p03 = (100, 2)  # start with values near those we expect
    sigma = np.ones(len(Time))
    sigma[[0]] = .3
    params5GPa, cv = scipy.optimize.curve_fit(monoExpFour, Time, FiveGPa, p03, sigma=sigma)
    g, h = params5GPa
    #print(params5GPa)

    fitted_function_300K = []
    fitted_function_400K = []
    fitted_function_500K = []
    fitted_function_600K = []
    fitted_function_700K = []

    for x in Time:
        fitted_function_300K.append(m * np.exp(-t * x))

    for x in Time:
        fitted_function_400K.append(a * np.exp(-b * x))

    for x in Time:
        fitted_function_500K.append(c * np.exp(-d * x))

    for x in Time:
        fitted_function_600K.append(e * np.exp(-f * x))

    for x in Time:
        fitted_function_700K.append(g * np.exp(-h * x))

    # plot the results
    plt.plot(Time, OneGPa, '--', label="300K", color='blue')
    plt.plot(Time, fitted_function_300K, color='blue')
    plt.plot(Time, TwoGPa, '--', label="400K", color='red')
    plt.plot(Time, fitted_function_400K, color='red')
    plt.plot(Time, ThreeGPa, '--', label="500K", color='orange')
    plt.plot(Time, fitted_function_500K, color='orange')
    plt.plot(Time, FourGPa, '--', label="600K", color='green')
    plt.plot(Time, fitted_function_600K, color='green')
    plt.plot(Time, FiveGPa, '--', label="700K", color='purple')
    plt.plot(Time, fitted_function_700K, color='purple')
    plt.title("TCP Graphs With Exponential Fits {} at {}".format(Pressure, speed))
    plt.ylabel('Intact Molecules')
    plt.legend()
    plt.xlabel('Time (ns)')
    plt.savefig('F:/PhD/TCPDecompositionExperiments/Completed/Fe3O4/'
               'IntactMoleculeComparison_Pressure_{}.png'.format(Pressure))
    plt.close()
    #plt.show()

    return Timestep_List , fitted_function_300K, fitted_function_400K, fitted_function_500K, fitted_function_600K, fitted_function_700K

def get_average_shear_normal_stress_and_average_mu_constant_temperature(Temperature, Pressures):
    Friction_Coefficient_Dataframe_Unnamed = pd.read_csv('F:/PhD/TCPDecompositionExperiments/Completed/Fe3O4/{}/1GPa/'
                                'fc_ave.dump'.format(Temperature), sep=' ')
    Friction_Coefficient_Dataframe = Friction_Coefficient_Dataframe_Unnamed.rename(columns={'v_s_bot' : 'Shear Stress 1GPa', 'v_p_bot' : 'Normal Stress 1GPa'})

    for P in Pressures:
        Dataframe = pd.read_csv('F:/PhD/TCPDecompositionExperiments/Completed/Fe3O4/{}/{}/'
                                'fc_ave.dump'.format(Temperature, P), sep=' ')
        Big_DataframeP = Dataframe.rename(columns= {'Timestep': 'Timestep {}'.format(P),
                                                        'v_s_bot': 'Shear Stress {}'.format(P),
                                                        'v_p_bot': 'Normal Stress {}'.format(P)})

        Friction_Coefficient_Dataframe = pd.concat([Friction_Coefficient_Dataframe, Big_DataframeP], axis =1)
        Friction_Coefficient_Dataframe = Friction_Coefficient_Dataframe.dropna()


    #print(Friction_Coefficient_Dataframe)
    Mu_Final_Dataframe = Friction_Coefficient_Dataframe.iloc[:, [0, 1, 2, 4, 5, 7, 8, 10, 11, 13, 14]]
    print(Mu_Final_Dataframe)


    #print(Mu_Final_Dataframe)

    ShearStressMeans = Mu_Final_Dataframe[['Shear Stress 1GPa', 'Shear Stress 2GPa', 'Shear Stress 3GPa', 'Shear Stress 4GPa', 'Shear Stress 5GPa']].mean()
    Average_Shear_Stress_Dictionary = ShearStressMeans.to_dict()
    #print(ShearStressMeans)
    NormalStressMeans = Mu_Final_Dataframe[['Normal Stress 1GPa', 'Normal Stress 2GPa', 'Normal Stress 3GPa', 'Normal Stress 4GPa', 'Normal Stress 5GPa']].mean()
    NormalStressMeans = NormalStressMeans.to_dict()
    #print(NormalStressMeans)


    Average_Mu_Dictionary = {}

    Normal_Stress = NormalStressMeans.get('Normal Stress 1GPa')
    Shear_Stress = ShearStressMeans.get('Shear Stress 1GPa')
    Average_Mu = Shear_Stress / Normal_Stress
    Average_Mu_Dictionary.update({'Average Mu 1GPa': Average_Mu})

    for P in Pressures:

        Normal_Stress = NormalStressMeans.get('Normal Stress {}'.format(P))
        Shear_Stress = ShearStressMeans.get('Shear Stress {}'.format(P))
        Average_Mu = Shear_Stress / Normal_Stress
        Average_Mu_Dictionary.update({'Average Mu {}'.format(P): Average_Mu})

    Average_Shear_Stress_List = list(Average_Shear_Stress_Dictionary.values())
    Average_Mu_List = list(Average_Mu_Dictionary.values())
    Average_Shear_Stress_List = [x / 10000 for x in Average_Shear_Stress_List]

    return Average_Shear_Stress_List, Average_Mu_List, NormalStressMeans

def get_average_shear_normal_stress_and_average_mu_constant_pressure(Pressure, Temperatures):
    Friction_Coefficient_Dataframe_Unnamed = pd.read_csv('F:/PhD/TCPDecompositionExperiments/Completed/Fe3O4/400K/{}/'
                                'fc_ave.dump'.format(Pressure), sep=' ')
    Friction_Coefficient_Dataframe = Friction_Coefficient_Dataframe_Unnamed.rename(columns={'v_s_bot' : 'Shear Stress 400K', 'v_p_bot' : 'Normal Stress 400K'})

    for T in Temperatures:
        Dataframe = pd.read_csv('F:/PhD/TCPDecompositionExperiments/Completed/Fe3O4/{}/{}/'
                                'fc_ave.dump'.format(T, Pressure), sep=' ')
        Big_DataframeP = Dataframe.rename(columns= {'Timestep': 'Timestep {}'.format(T),
                                                        'v_s_bot': 'Shear Stress {}'.format(T),
                                                        'v_p_bot': 'Normal Stress {}'.format(T)})

        Friction_Coefficient_Dataframe = pd.concat([Friction_Coefficient_Dataframe, Big_DataframeP], axis =1)
        Friction_Coefficient_Dataframe = Friction_Coefficient_Dataframe.dropna()


    #print(Friction_Coefficient_Dataframe)
    Mu_Final_Dataframe = Friction_Coefficient_Dataframe.iloc[:, [0, 1, 2, 4, 5, 7, 8, 10, 11, 13, 14]]

    ShearStressMeans = Mu_Final_Dataframe[['Shear Stress 400K', 'Shear Stress 500K', 'Shear Stress 600K']].mean()
    Average_Shear_Stress_Dictionary = ShearStressMeans.to_dict()
    #print(ShearStressMeans)
    NormalStressMeans = Mu_Final_Dataframe[['Normal Stress 400K', 'Normal Stress 500K', 'Normal Stress 600K']].mean()
    NormalStressMeans = NormalStressMeans.to_dict()
    #print(NormalStressMeans)

    Average_Mu_Dictionary = {}

    Normal_Stress = NormalStressMeans.get('Normal Stress 400K')
    #print(Normal_Stress)
    Shear_Stress = ShearStressMeans.get('Shear Stress 400K')
    #print(Shear_Stress)
    Average_Mu = Shear_Stress / Normal_Stress
    Average_Mu_Dictionary.update({'Average Mu 400K': Average_Mu})

    for T in Temperatures:

        Normal_Stress = NormalStressMeans.get('Normal Stress {}'.format(T))
        Shear_Stress = ShearStressMeans.get('Shear Stress {}'.format(T))
        Average_Mu = Shear_Stress / Normal_Stress
        Average_Mu_Dictionary.update({'Average Mu {}'.format(T): Average_Mu})

    Average_Shear_Stress_List = list(Average_Shear_Stress_Dictionary.values())
    Average_Mu_List = list(Average_Mu_Dictionary.values())
    Average_Shear_Stress_List = [x / 10000 for x in Average_Shear_Stress_List]

    return Average_Shear_Stress_List, Average_Mu_List, NormalStressMeans

def plot_shear_stress_vs_normal_stress(Average_Shear_Stress_List_1, Average_Shear_Stress_List_2, Average_Shear_Stress_List_3, Average_Shear_Stress_List_4, Average_Shear_Stress_List_5, Temp1, Temp2, Temp3, Temp4, Temp5):
    x = np.array([1, 2, 3, 4, 5])
    a, b = np.polyfit(x, Average_Shear_Stress_List_1, 1)
    c, d = np.polyfit(x, Average_Shear_Stress_List_2, 1)
    e, f = np.polyfit(x, Average_Shear_Stress_List_3, 1)
    g, h = np.polyfit(x, Average_Shear_Stress_List_4, 1)
    i, j = np.polyfit(x, Average_Shear_Stress_List_5, 1)

    fig1, ax2 = plt.subplots()
    ax2.set_title('Shear Stress vs Normal Stress at Different Temperatures')
    ax2.set_xlabel('Normal Stress (GPa)')
    ax2.set_ylabel('Shear Stress (GPa)')
    ax2.scatter(x, Average_Shear_Stress_List_1)
    ax2.scatter(x, Average_Shear_Stress_List_2)
    ax2.scatter(x, Average_Shear_Stress_List_3)
    ax2.scatter(x, Average_Shear_Stress_List_4)
    ax2.scatter(x, Average_Shear_Stress_List_5)
    ax2.plot(x, a * x + b, label=Temp1)
    ax2.plot(x, c * x + d, label=Temp2)
    ax2.plot(x, e * x + f, label=Temp3)
    ax2.plot(x, g * x + h, label=Temp4)
    ax2.plot(x, i * x + j, label=Temp5)
    ax2.legend()
    plt.show()

def plot_variation_in_mu(Average_Mu_List_1, Average_Mu_List_2, Average_Mu_List_3, Average_Mu_List_4, Average_Mu_List_5, temp1, temp2, temp3, temp4, temp5):
    x = np.array([1, 2, 3, 4, 5])
    a, b = np.polyfit(x, Average_Mu_List_1, 1)
    c, d = np.polyfit(x, Average_Mu_List_2, 1)
    e, f = np.polyfit(x, Average_Mu_List_3, 1)
    g, h = np.polyfit(x, Average_Mu_List_4, 1)
    i, j = np.polyfit(x, Average_Mu_List_5, 1)

    fig1, ax1 = plt.subplots()
    ax1.set_title('Normal Stress vs Mu at Different Speeds')
    ax1.set_xlabel('Normal Stress (GPa)')
    ax1.set_ylabel('Shear Stress (GPa)')
    ax1.scatter(x, Average_Mu_List_1)
    ax1.scatter(x, Average_Mu_List_2)
    ax1.scatter(x, Average_Mu_List_3)
    ax1.scatter(x, Average_Mu_List_4)
    ax1.scatter(x, Average_Mu_List_5)
    ax1.plot(x, a * x + b, label=temp1)
    ax1.plot(x, c * x + d, label=temp2)
    ax1.plot(x, e * x + f, label=temp3)
    ax1.plot(x, g * x + h, label=temp4)
    ax1.plot(x, i * x + j, label=temp5)
    ax1.legend()
    plt.show()


def get_MATLABFIT_dissociation_rates(TimestepList, Coefficient, Cutoff=None):
    Time = []
    for x in TimestepList:
        Time.append(((x - 400000) / int(4 * 10 ** 6)))
    Extrapolation_Constant = 999 / len(Time)
    FittedCurveValues = []
    for x in Time:
        Value = 48 * np.exp(Coefficient*(x))
        FittedCurveValues.append(Value)
    UnextrapolatedRate = (np.log(FittedCurveValues[-1] / 48)) / Time[-1]
    NanosecondRate = Extrapolation_Constant * UnextrapolatedRate * -1
    LnRate = np.log(NanosecondRate)
    Value = list(range(0, len(Time)))

    if Cutoff != None:
        TimeCutoff = []
        [TimeCutoff.append(Time[x]) for x in Value if Time[x] <= Cutoff]
        fitted_function = []
        for x in TimeCutoff:
            fitted_function.append(48 * np.exp(Coefficient * x))

        ShortIntactMoleculeList = []
        [ShortIntactMoleculeList.append(fitted_function[x]) for x in Value[:len(TimeCutoff)]]
        UnextrapolatedRate = (np.log((ShortIntactMoleculeList[-1] / 48))) / TimeCutoff[-1]

        Extrapolation_Constant = 999 / len(TimeCutoff)
        NanosecondRate = Extrapolation_Constant #* (UnextrapolatedRate * -1)
        LnRate = np.log(NanosecondRate)
    else:
        Extrapolation_Constant = 999 / len(Time)
        FittedCurveValues = []
        for x in Time:
            Value = 48 * np.exp(Coefficient*(x))
            FittedCurveValues.append(Value)
        #print(Time[-1])
        UnextrapolatedRate = (np.log((FittedCurveValues[-1] / 48))) / Time[-1]
        NanosecondRate = (UnextrapolatedRate * -1) #* Extrapolation_Constant
        LnRate = np.log(NanosecondRate)
    return NanosecondRate, LnRate

def plot_variation_in_shear_stress_constanttemp(temperature, pressures):
    Friction_Coefficient_Dataframe_Unnamed = pd.read_csv('F:/PhD/TCPDecompositionExperiments/Completed/Fe3O4/{}/1GPa/'
                                'fc_ave.dump'.format(temperature), sep=' ')
    Friction_Coefficient_Dataframe = Friction_Coefficient_Dataframe_Unnamed.rename(columns={'v_s_bot' : 'Shear Stress 1GPa', 'v_p_bot' : 'Normal Stress 1GPa'})

    for P in Pressures:
        Dataframe = pd.read_csv('F:/PhD/TCPDecompositionExperiments/Completed/Fe3O4/{}/{}/'
                                'fc_ave.dump'.format(temperature, P), sep=' ')
        Big_DataframeP = Dataframe.rename(columns= {'TimeStep': 'Timestep {}'.format(P),
                                                        'v_s_bot': 'Shear Stress {}'.format(P),
                                                        'v_p_bot': 'Normal Stress {}'.format(P)})

        Friction_Coefficient_Dataframe = pd.concat([Friction_Coefficient_Dataframe, Big_DataframeP], axis =1)
        Friction_Coefficient_Dataframe = Friction_Coefficient_Dataframe.dropna()

    Timestep = Friction_Coefficient_Dataframe.TimeStep.tolist()

    Shear_Stress_1GPa = Friction_Coefficient_Dataframe['Shear Stress 1GPa'].tolist()
    Shear_Stress_2GPa = Friction_Coefficient_Dataframe['Shear Stress 2GPa'].tolist()
    Shear_Stress_3GPa = Friction_Coefficient_Dataframe['Shear Stress 3GPa'].tolist()
    Shear_Stress_4GPa = Friction_Coefficient_Dataframe['Shear Stress 4GPa'].tolist()
    Shear_Stress_5GPa = Friction_Coefficient_Dataframe['Shear Stress 5GPa'].tolist()

    Shear_Stress_1GPa = [x / 10000 for x in Shear_Stress_1GPa]
    Shear_Stress_2GPa = [x / 10000 for x in Shear_Stress_2GPa]
    Shear_Stress_3GPa = [x / 10000 for x in Shear_Stress_3GPa]
    Shear_Stress_4GPa = [x / 10000 for x in Shear_Stress_4GPa]
    Shear_Stress_5GPa = [x / 10000 for x in Shear_Stress_5GPa]

    Time = []
    for x in Timestep:
        Time.append(((x - 400000) / int(4 * 10 ** 6)))

    fig3, ax3 = plt.subplots()
    ax3.set_title('Variation in Shear Stress at {}'.format(temperature))
    ax3.set_xlabel('Time (ns)')
    ax3.set_ylabel('Shear Stress (GPa)')
    ax3.plot(Time, Shear_Stress_1GPa, label='1GPa')
    ax3.plot(Time, Shear_Stress_2GPa, label='2GPa')
    ax3.plot(Time, Shear_Stress_3GPa, label='3GPa')
    ax3.plot(Time, Shear_Stress_4GPa, label='4GPa')
    ax3.plot(Time, Shear_Stress_5GPa, label='5GPa')
    ax3.legend()
    plt.show()