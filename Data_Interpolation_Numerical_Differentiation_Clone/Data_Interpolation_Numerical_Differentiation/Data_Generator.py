# Generate data for use in beam propagation calculations
# position data should be randomly spaced if required
# R. Sheehan 17 - 9 - 2014

# Import libraries
import math
import numpy as np
import random
import matplotlib.pyplot as plt
import os
import glob
import re

# function definitions
# should really learn how to do this so it links to a file somewhere

def count_lines(thedata, thepath):
    # count the number of lines in a file that has been opened and read into memory
    # thedata is the stream that contains the data from an open file
    # how do you know if thedata contains data?
    # assume that you only call count_lines inside another function for reading data
    # thepath is the name of the file containing the data
    # R. Sheehan 26 - 4 - 2014

    nlines=0
    for lines in thedata:
        nlines = nlines + 1
##    print "There are %(nlines)d lines in %(path)s"%{"nlines":nlines,"path":thepath}
    return nlines

def open_file(thepath):
    # given a filename thepath, open it for reading and return a file object that can be accessed
    # it may be an idea to return a dictionary or something
    # R. Sheehan 26 - 4 - 2014

    thefile = file(thepath,"r") # open file for reading

    # check that the files are available
    if thefile.closed:
        print "%(path)s could not be opened"%{"path":thepath}
        return 0
    else:
        print "%(path)s is open"%{"path":thepath}
        return thefile

def read_data(thepath):
    # read data from an open file
    # R. Sheehan 26 - 4 - 2014

    datapts = 0 # old C++ habits die hard

    thefile = file(thepath,"r") # open file for reading

    # check that the files are available
    if thefile.closed:
        print "%(path)s could not be opened"%{"path":thepath}
        datapts = -1
    else:
        print "%(path)s is open"%{"path":thepath}

        thedata = thefile.readlines() # read the data from the file

        nlines = count_lines(thedata, thepath) # count the number of data points

        datapts = np.zeros([nlines]) # create an array of zeros of length nlines

        i=0
        for lines in thedata:
            datapts[i] = lines
            i = i + 1

        del thedata # clear thedata stream

    del thefile

    return datapts
	
def write_data(thepath, thedata):
    # write a vector containing data to a file
    # R. Sheehan 7 - 8 - 2014

    thefile = file(thepath,"w") # create a file for writing

    # check that the file is available
    if thefile.closed:
        print "%(path)s could not be opened"%{"path":thepath}
        datapts = -1
    else:
        # Write the data to the files
        print "%(path)s is open"%{"path":thepath}
        for data in thedata:
            thefile.write("%(num)0.9f\n"%{"num":data})

    # delete the file objects
    del thefile
    
    return 0

def func(pos):
    # function that returns a value at pos
    if abs(pos) < 1.0e-12:
        return 1.0
    else:
        arg = (pos)**2
        return math.cos(arg)
##        return math.cos(pos)

def dfunc(pos):
    # derivative of the function func
    if abs(pos) < 1.0e-12:
        return 0.0
    else:
        arg = (pos)**2
        return -2.0*(pos)*math.sin(arg)
##        return -1.0*math.sin(pos)

def output_interpolation_data(xl, xu, npoints, make_plot):
    # generate interpolation data over the range (xl, xu)
    # plot the data, output to files

##    X = np.linspace(xl, xu, npoints) # same dx for all
    X = np.zeros([npoints])
    for i in range(0,npoints, 1):
        X[i] = random.uniform(xl, xu)
    X = np.sort(X)

    #print X

    # generate function data
    Y = np.zeros([npoints])
    DY = np.zeros([npoints])
    for i in range(0,npoints,1):
        Y[i] = func(X[i])
        DY[i] = dfunc(X[i])

    maxval = max(np.max(Y), np.max(DY))
    minval = min(np.min(Y), np.min(DY) )

    if make_plot:
        # plot the data for the lolz
        fig = plt.figure()

        ax = fig.add_subplot(111)

        ax.plot(X, Y, 'r+-', label = 'f(x)')
        ax.plot(X, DY, 'gx-', label = 'f\'(x)')
        ax.legend(loc='upper center')

        plt.xlabel('x', fontsize = 20)
        plt.ylabel('y', fontsize = 20)
        plt.title('Function and its derivative')
        plt.axis( [xl, xu, minval, maxval] )

        plt.savefig('Interpolation_Data.png')

        plt.show()

        plt.clf()
        plt.cla()
        plt.close()

    # write the data to a files  
    xfile = "posdata.txt"
    write_data(xfile, X)
    
    yfile = "valdata.txt"
    write_data(yfile, Y)

    yfile = "derivdata.txt"
    write_data(yfile, DY)

    del X, Y, DY

    return 0

def gauss_func(pos, peak_loc, peak_width):
    # evaluate gaussian distribution function

    arg = (pos - peak_loc) / peak_width
    argsqr = arg**2

    return math.exp(-1.0*argsqr)

def nasty_func(pos):

    t1 = math.cos(3.0*pos)
    t2 = (math.sin((pos)**3))**2

    peaks = [0.2, 0.4, 0.6, 0.8]
    widths = [0.007, 0.01, 0.02, 0.04]

    retval = t1*t2

    for i in range(0,len(peaks),1):
        retval = retval+4.0*gauss_func(pos, peaks[i], widths[i])
    

    return retval

def output_filter_data(xl, xu, npoints, make_plot, add_noise = 0):
    # generate filtering data over the range (xl, xu)
    # plot the data, output to files

    X = np.linspace(xl, xu, npoints) # same dx for all
##    X = np.zeros([npoints])
##    for i in range(0,npoints, 1):
##        X[i] = random.uniform(xl, xu)
##    X = np.sort(X)

    #print X

    # parameters for gaussian noise
    mu = 0.0 # expected value
    sigma = 2.0 # variance

    # generate function data
    Y = np.zeros([npoints])
    for i in range(0,npoints,1):
        Y[i] = nasty_func(X[i])
        if add_noise:
            Y[i] = Y[i] + random.gauss(mu, sigma)

    maxval = np.max(Y)
    minval = np.min(Y)

    if make_plot:
        # plot the data for the lolz
        fig = plt.figure()

        ax = fig.add_subplot(111)

        ax.plot(X, Y, 'r+-', label = 'f(x)')
        
        ax.legend(loc='upper center')

        plt.xlabel('x', fontsize = 20)
        plt.ylabel('y', fontsize = 20)
        plt.title('Function')
        plt.axis( [xl, xu, minval, maxval] )

        plt.savefig('Filter_Data.png')

        plt.show()

        plt.clf()
        plt.cla()
        plt.close()

    # write the data to a files  
    xfile = "fil_posdata.txt"
    write_data(xfile, X)
    
    yfile = "fil_valdata.txt"
    write_data(yfile, Y)

    del X, Y

    return 0

def main():
    pass

if __name__ == '__main__':
    main()

    pwd = os.getcwd() # get current working directory
    os.chdir(pwd)

    print pwd

# Firstly, I want to generate data to test an interpolation
# and numerical differentiation code
# so I want to make sure that I can compute function and its derivative exactly
# then compare with numerical results

##    output_interpolation_data(0, 5, 150, 1)


    output_filter_data(0, 2.5, 1024, 1, 1)
    
