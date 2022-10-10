import pandas as pd
import matplotlib.pyplot as plt
import scipy
from scipy import optimize
import numpy as np
import sys
from numpy import mean as m
import os
import os.path
import sys
import statistics
import glob
from HelperFunctions import get_intact_columns_constant_pressure
from HelperFunctions import get_intact_columns_constant_temperature
from HelperFunctions import plot_shear_stress_vs_normal_stress
from HelperFunctions import get_fitted_plots_constant_temperature
from HelperFunctions import get_fitted_plots_constant_pressure
from HelperFunctions import get_average_shear_normal_stress_and_average_mu_constant_temperature
from HelperFunctions import get_average_shear_normal_stress_and_average_mu_constant_pressure
from mpl_toolkits import mplot3d as Axes3D
from HelperFunctions import plot_shear_stress_vs_normal_stress, plot_variation_in_mu
from HelperFunctions import get_MATLABFIT_dissociation_rates
#import Dissociation_Rates
from scipy.stats.distributions import t

Temperatures = ["500K", "600K"]
Pressures = ['2GPa', '3GPa', '4GPa', '5GPa']

Big_Dataframe_400K = get_intact_columns_constant_temperature("400K", Pressures)
Big_Dataframe_500K = get_intact_columns_constant_temperature("500K", Pressures)
Big_Dataframe_600K = get_intact_columns_constant_temperature("600K", Pressures)

Big_Dataframe_1GPa = get_intact_columns_constant_pressure('1GPa', Temperatures)
Big_Dataframe_2GPa = get_intact_columns_constant_pressure('2GPa', Temperatures)
Big_Dataframe_3GPa = get_intact_columns_constant_pressure('3GPa', Temperatures)
Big_Dataframe_4GPa = get_intact_columns_constant_pressure('4GPa', Temperatures)
Big_Dataframe_5GPa = get_intact_columns_constant_pressure('5GPa', Temperatures)

Big_Dataframe_1GPa.to_csv('F:/PhD/TCPDecompositionExperiments/Completed/Fe3O4/OneGPaComparison.csv')
Big_Dataframe_2GPa.to_csv('F:/PhD/TCPDecompositionExperiments/Completed/Fe3O4/TwoGPaComparison.csv')
Big_Dataframe_3GPa.to_csv('F:/PhD/TCPDecompositionExperiments/Completed/Fe3O4/ThreeGPaComparison.csv')
Big_Dataframe_4GPa.to_csv('F:/PhD/TCPDecompositionExperiments/Completed/Fe3O4/FourGPaComparison.csv')
Big_Dataframe_5GPa.to_csv('F:/PhD/TCPDecompositionExperiments/Completed/Fe3O4/FiveGPaComparison.csv')

#Big_Dataframe_300K.to_csv('F:/PhD/TCPDecompositionExperiments/Completed/Fe3O4/300KComparison.csv')
Big_Dataframe_400K.to_csv('F:/PhD/TCPDecompositionExperiments/Completed/Fe3O4/400KComparison.csv')
Big_Dataframe_500K.to_csv('F:/PhD/TCPDecompositionExperiments/Completed/Fe3O4/500KComparison.csv')
Big_Dataframe_600K.to_csv('F:/PhD/TCPDecompositionExperiments/Completed/Fe3O4/600KComparison.csv')
Big_Dataframe_600K.to_csv('F:/PhD/TCPDecompositionExperiments/Completed/Fe3O4/600KComparison.csv')
#Big_Dataframe_700K.to_csv('F:/PhD/TCPDecompositionExperiments/Completed/Fe3O4/700KComparison.csv')

############## Get Timestep lists for use with MATLAB dissociation rate function ######################

#Timestep_300K = Big_Dataframe_300K['Timestep']
Timestep_400K = Big_Dataframe_400K['Timestep']
Timestep_500K = Big_Dataframe_500K['Timestep']
Timestep_600K = Big_Dataframe_600K['Timestep']
#Timestep_700K = Big_Dataframe_700K['Timestep']
Timestep_1GPa = Big_Dataframe_1GPa['Timestep']
Timestep_2GPa = Big_Dataframe_2GPa['Timestep']
Timestep_3GPa = Big_Dataframe_3GPa['Timestep']
Timestep_4GPa = Big_Dataframe_4GPa['Timestep']
Timestep_5GPa = Big_Dataframe_5GPa['Timestep']
