# Import libraries
# You should try an import the bare minimum of modules
import sys # access system routines
import os
import glob
import re

import math
import scipy
import numpy as np
import matplotlib.pyplot as plt

# add path to our file
sys.path.append('c:/Users/robertsheehan/Programming/Python/Common/')
sys.path.append('c:/Users/robertsheehan/Programming/Python/Plotting/')

import Common
import Plotting

MOD_NAME_STR = "Plots" # use this in exception handling messages
       
def sg_smoothing_plot(filename):
    # make a plot of the data from the Savitzky-Golay filter tests
    # R. Sheehan 29 - 3 - 2022

    FUNC_NAME = ".sg_smoothing_plot()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
            
        if glob.glob(filename):
            # import the dataset
            hv_data = []; labels = []; marks = []; 
            hv_data_1 = []; labels_1 = []; marks_1 = []; 
            data = np.loadtxt(filename, delimiter = ',', unpack = True)
            hv_data.append([data[0], data[1]]); marks.append(Plotting.labs_lins[1]); labels.append('original'); 
            hv_data.append([data[0], data[2]]); marks.append(Plotting.labs_pts[2]); labels.append('noisy'); 
            hv_data.append([data[0], data[3]]); marks.append(Plotting.labs_lins[0]); labels.append('smoothed');    
            
            # plot the original data with the smoothed function                
            args = Plotting.plot_arg_multiple()
            
            #args.loud = True
            args.crv_lab_list = labels
            args.mrk_list = marks
            args.x_label = 'X'
            args.y_label = 'Y'
            args.fig_name = filename.replace('.txt','')
            args.plt_title = filename.replace('.txt','')
            
            Plotting.plot_multiple_curves(hv_data, args)
        else:
            ERR_STATEMENT = ERR_STATEMENT + "\nFile: " + filename + " not found"
            raise Exception
    except Exception as e:
        print(ERR_STATEMENT)
        print(e)
        
def sg_smoothing_derivative_plot(filename):
    # make a plot of the data from the Savitzky-Golay filter tests
    # R. Sheehan 29 - 3 - 2022

    FUNC_NAME = ".sg_smoothing_plot()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
            
        if glob.glob(filename):
            # import the dataset
            hv_data = []; labels = []; marks = []; 
            hv_data_1 = []; labels_1 = []; marks_1 = []; 
            data = np.loadtxt(filename, delimiter = ',', unpack = True)
            hv_data.append([data[0], data[4]]); marks.append(Plotting.labs_lins[1]); labels.append('original'); 
            hv_data.append([data[0], data[5]]); marks.append(Plotting.labs_pts[2]); labels.append('noisy'); 
            hv_data.append([data[0], data[3]]); marks.append(Plotting.labs_dashed[0]); labels.append('smoothed');    
            
            # plot the original data with the smoothed function                
            args = Plotting.plot_arg_multiple()
            
            args.loud = False
            args.crv_lab_list = labels
            args.mrk_list = marks
            args.x_label = 'X'
            args.y_label = 'dY / dX'
            args.plt_range = [0, 100, -0.2, 0.2]
            args.fig_name = filename.replace('.txt','')
            args.plt_title = filename.replace('.txt','')
            
            Plotting.plot_multiple_curves(hv_data, args)
        else:
            ERR_STATEMENT = ERR_STATEMENT + "\nFile: " + filename + " not found"
            raise Exception
    except Exception as e:
        print(ERR_STATEMENT)
        print(e)

def main():
    pass

if __name__ == '__main__':
    main()

    pwd = os.getcwd() # get current working directory
    
    print(pwd)
    
    # SG Filtered Gauss Data
    plot_Gauss_data = False
    if plot_Gauss_data:
        os.chdir('fgauss_test')
        
        m = 2
        ld = 0
        nl = nr = 4
        while m < 7:
            filename = "SG_Filter_fgauss_m_%(v1)d_ld_%(v2)d_nl_%(v3)d_nr_%(v3)d.txt"%{"v1":m, "v2":ld, "v3":nl}
            sg_smoothing_plot(filename)
            m = m + 2
        
        m = 4
        ld = 0
        nl = 4
        while nl < 65:
            nr = nl
            filename = "SG_Filter_fgauss_m_%(v1)d_ld_%(v2)d_nl_%(v3)d_nr_%(v3)d.txt"%{"v1":m, "v2":ld, "v3":nl}
            sg_smoothing_plot(filename)
            nl = nl * 2  

        os.chdir(pwd)
        
    # SG Filtered Ring Down Data
    plot_Ring_Down_data = True
    if plot_Ring_Down_data:
        m = 4
        ld = 1
        nl = 10
        
        os.chdir('ringdown_test') if ld == 0 else os.chdir('ringdown_der_test')
        
        while m < 9:
            filename = "SG_Filter_ring_down_m_%(v1)d_ld_%(v2)d_nl_%(v3)d_nr_%(v3)d.txt"%{"v1":m, "v2":ld, "v3":nl}
            sg_smoothing_plot(filename) if ld == 0  else sg_smoothing_derivative_plot(filename)
            m = m + 2
            
        m = 10
        nl = 8
        while nl < 129:
            filename = "SG_Filter_ring_down_m_%(v1)d_ld_%(v2)d_nl_%(v3)d_nr_%(v3)d.txt"%{"v1":m, "v2":ld, "v3":nl}
            sg_smoothing_plot(filename) if ld == 0  else sg_smoothing_derivative_plot(filename)
            nl = nl * 2
            
        os.chdir(pwd)
    
