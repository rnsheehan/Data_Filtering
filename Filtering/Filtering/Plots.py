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
            hv_data.append([data[0], data[1]]); marks.append(Plotting.labs_pts[0]); labels.append('data'); 
            hv_data.append([data[0], data[2]]); marks.append(Plotting.labs_lins[1]); labels.append('smoothing'); 
            #hv_data.append([data[0], data[3]]); marks.append(Plotting.labs_pts[2]); labels.append('derivative'); 
            
            # plot the original data with the smoothed function                
            args = Plotting.plot_arg_multiple()
            
            args.loud = True
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

def main():
    pass

if __name__ == '__main__':
    main()

    pwd = os.getcwd() # get current working directory
    
    print(pwd)
    
    #filename = "SG_Filter_fgauss_m_2_ld_0_nl_4_nr_4.txt"
    #filename = "SG_Filter_fgauss_m_4_ld_0_nl_4_nr_4.txt"
    #filename = "SG_Filter_fgauss_m_4_ld_0_nl_8_nr_8.txt"
    #filename = "SG_Filter_fgauss_m_4_ld_0_nl_16_nr_16.txt"
    
    #filename = "SG_Filter_ring_down_m_4_ld_0_nl_8_nr_8.txt"
    #filename = "SG_Filter_ring_down_m_4_ld_1_nl_8_nr_8.txt"
    #filename = "SG_Filter_ring_down_m_6_ld_1_nl_10_nr_10.txt"
    #filename = "SG_Filter_ring_down_m_6_ld_1_nl_20_nr_20.txt"
    filename = "SG_Filter_ring_down_m_6_ld_0_nl_40_nr_40.txt"
    #filename = "SG_Filter_ring_down_m_8_ld_0_nl_40_nr_40.txt"
    
    sg_smoothing_plot(filename)
    # data = np.loadtxt(filename, delimiter = ',', unpack = True)
    # print(data.shape)
    # print(data[0])
    
