import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime, os, glob, json, string, inspect, emcee
from scipy.stats import median_abs_deviation
from matplotlib.patches import Circle
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.io import ascii, fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.table import Table
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)

def filter_times(lc):
    '''
    Take an input lc array and return
    a time filtered array, input/output arrrays:
        0 -- Time (MJD)
        1 -- Rate (counts/s/cm^2)
        2 -- error on rate (counts/s/cm^2)
    '''
    
    ti = 59340 # May 6 2021
    tf = 59511 # Nov 01 2021
    
    # Get index
    index = np.where((lc[:,0] > ti) & (lc[:,0] < tf))[0]
        
    #Return filtered array
    return lc[index,:]

def filter_errors(lc):
    '''
    Take the input lc and filter any data points
    that are greater than 2 standard deviations
    from the median value
    '''
    
    index = np.where(lc[:,2] < np.median(lc[:,2]) + 3.0 * np.std(lc[:,2], ddof = 1))[0]
    
    # Return good values 
    return lc[index,:]
    
def match_times(maxi,bat):
    '''
    Take the two input filtered lcs
    bat & maxi, and match them by times
    returning a single array that has index:
    0 -- Time (MJD)
    1 -- 4-10keV rate
    2 -- 4-10keV error
    3 -- 15-50keV rate
    4 -- 15-50keV error
    '''
    
    maxi_times = np.floor(maxi[:,0]) # This will remove the 0.5 that maxi uses
    bat_times  = bat[:,0]

    # Find the matching indexes
    index = np.in1d(bat_times,maxi_times)
    bat = bat[index,:]
    bat_times  = bat[:,0]
    index = np.in1d(maxi_times, bat_times)
    maxi = maxi[index,:]
    maxi[:,0] = np.floor(maxi[:,0])
        
    return np.array([maxi[:,0] + 0.5, maxi[:,1], maxi[:,2], bat[:,1], bat[:,2]])

def weighted_average(y,dy):
    '''
    Calculate the (variance)-weighted average  
    '''
    avg = np.average(y, weights = dy ** (-2))
    err = (np.sum(dy ** (-2))) ** (-0.5)
    
    return avg, err

def get_bin_indexes(arr):
    '''
    Here we are going to solve for bins ensuring that both atleast one of maxi or bat have a 3-sigma detection, input arr has the following indexes:
        0 -- Time (MJD)
        1 -- 4-10keV rate
        2 -- 4-10keV error
        3 -- 15-50keV rate
        4 -- 15-50keV error
    '''    
    
    # Binning arrays 
    left_index  = []
    right_index = []
    
    # Iterate through the arrays and group days until Swift-BAT has a 3-sigma detection within the bin
    k = 0
    while k < len(arr[0]):
        
        # If both data points are 3-sigma detections append the data
        if arr[3][k]/arr[4][k] >= 3.0 and arr[1][k] > 0.0:
            left_index.append(k)
            right_index.append(k+1)
            k += 1
            
        else:
            i = 0 
            bat_avg = arr[3][k]
            bat_avg_err = arr[4][k]
            maxi_avg = arr[1][k]
            maxi_avg_err = arr[2][k]
            while bat_avg/bat_avg_err < 3 or maxi_avg < 0.0:
                i += 1
                maxi_avg, maxi_avg_err = weighted_average(arr[1][k:k+(i+1)],arr[2][k:k+(i+1)])
                bat_avg, bat_avg_err = weighted_average(arr[3][k:k+(i+1)],arr[4][k:k+(i+1)])
                                
                if k + (i + 1) > len(arr[0]): # Force a break at the edge of the data set
                    break
                    
            if k + (i + 1) > len(arr[0]):
                left_index.append(int(k))
                right_index.append(len(arr[0]))
            else:
                left_index.append(int(k))
                right_index.append(int(k+i+1))
            k = k + (i + 1)
        
    return left_index, right_index

def group_days(arr,left_index,right_index):
    '''
    Here we are going to group days into independant time bins window
    ensuring that both atleast one of maxi or bat have a 3-sigma detection, input arr has the following indexes:
        0 -- Time (MJD)
        1 -- 4-10keV rate
        2 -- 4-10keV error
        3 -- 15-50keV rate
        4 -- 15-50keV error
    '''
    
    # Grouped arrays
    time_ctr   = []
    n_bins     = []
    maxi_rate  = []
    maxi_err   = []
    bat_rate   = []
    bat_err    = []
    key        = [] # 0 -- detection, 1 -- lower lim (maxi - no-det, bat - det)

    # Bin the data 
    for li, ri in zip(left_index, right_index):
        
        # Calculate the averages
        maxi_avg, maxi_avg_err = weighted_average(arr[1][li:ri],arr[2][li:ri])
        bat_avg, bat_avg_err = weighted_average(arr[3][li:ri],arr[4][li:ri])
        
        # Append the values
        time_ctr.append(np.mean(arr[0][li:ri]))
        n_bins.append(ri-li)
        bat_rate.append(bat_avg)
        bat_err.append(bat_avg_err)
                        
        # Append the keys that categorize the data

        # Both instruments are detected at 3 sigma
        if maxi_avg/maxi_avg_err >= 3: 
            key.append(0)
            maxi_rate.append(maxi_avg)
            maxi_err.append(maxi_avg_err)
                        
        # Only BAT has a 3 sigma detection, set conservative lower limit by using 3 x MAXI rms as signal  
        elif maxi_avg/maxi_avg_err < 3:
            key.append(1)
            maxi_rate.append(3.0 * maxi_avg_err)
            maxi_err.append(maxi_avg_err)
    
    return np.array([time_ctr, n_bins, key, maxi_rate, maxi_err, bat_rate, bat_err])


def logprob(p,x,y,x_err,y_err):
    '''
    log probability for Bayesian fit
    '''
    theta = p[0]
    if np.abs(theta-np.pi/4) > np.pi/4:
        return -np.inf
    Delta = (np.cos(theta)*y - np.sin(theta)*x)**2
    Sigma = (np.sin(theta))**2*x_err**2+(np.cos(theta))**2*y_err**2
    lp = -0.5*np.nansum(Delta/Sigma)-0.5*np.nansum(np.log(Sigma))

    return lp

def confidence_interval(y):
    '''
    Confidence interval, to get error from 
    MCMC samples
    '''
    median=np.median(y)
    pct15=np.percentile(y,15)
    pct85=np.percentile(y,85)
    list1=np.array([median,median-pct15,pct85-median])
    return list1


def MCMC_HR_fit(x,y,x_err,y_err,nWalkers=10,nBurn=1000,nSample=10000):
    '''
    Fit the HR ratio following the MCMC prescription
    from the WATCHDOG pipeline (Tetarenko et al. 2016)
    '''
    
    ndim = 1
    p0 = np.zeros((nWalkers,ndim))
    p0[:,0] = np.pi/4+np.random.randn(nWalkers)*0.1

    sampler = emcee.EnsembleSampler(nWalkers, ndim, logprob, 
                                    args=[x,y,x_err,y_err])
    pos, prob, state = sampler.run_mcmc(p0, nBurn)
    sampler.reset()
    sampler.run_mcmc(pos, nSample)
    samples=np.sort(np.tan(sampler.flatchain),axis=None) 
    
    CI = confidence_interval(samples)
    return CI

def plot_binned_lc(group_array, matched_array, l, r):

    '''
    Here we are going to plot the binned data and the matched data
    '''
    fig,ax = plt.subplots(figsize=(12,8))
    ax.errorbar(matched_array[0],matched_array[-2], yerr = matched_array[-1], fmt='.', label='Swift-BAT Unbinned')
    ax.errorbar(group_array[0], group_array[-2], yerr = group_array[-1], fmt='.', label='Swift-BAT binned')
    ax.errorbar(matched_array[0],matched_array[-4], yerr = matched_array[-3], fmt='.', label='MAXI/GSC Unbinned')
    ax.errorbar(group_array[0], group_array[-4], yerr = group_array[-3], fmt='.', label='MAXI/GSC binned')
    for li, ri in zip(l, r):
        if ri-li != 1:
            ax.axvline(matched_array[0][li],c = 'k' ,ls = ':')
            ax.axvline(matched_array[0][ri-1],c = 'r' ,ls = ':')
    ax.set_ylabel('Count Rate\n(Crab)')
    ax.set_xlabel('Observing Data (MJD)')
    ax.legend()
    plt.savefig('../plots/HR_lc_Binning.png')
    plt.clf()
    plt.close()
    
    
def main():
    # Load in the the "vanilla" light curves
    print('Loading in Watchdog [crab-calibrated] light curves')
    bat = np.loadtxt('../SAXJ1810_watchdog/BATObs/SAXJ1810.8m2609/SAXJ1810.8m2609_BAT_lc_15_50_crab.txt')
    maxi = np.loadtxt('../SAXJ1810_watchdog/MAXIObs/SAXJ1810.8m2609/SAXJ1810.8m2609_MAXI_lc_4_10_crab.txt')

    # Filter the relevant time range
    print('Filtering times and anomously large errors')
    bat_filtered = filter_times(bat)
    bat_filtered = filter_errors(bat_filtered)
    maxi_filtered = filter_times(maxi)
    maxi_filtered = filter_errors(maxi_filtered)

    # Match and group the array
    print('Matching the times of MAXI/GSC and Swift-BAT data')
    matched_array = match_times(maxi_filtered, bat_filtered)

    # Get the binning indexes:
    print('Solving for binning indexes and applying manual changes')
    l, r = get_bin_indexes(matched_array)

    # Apply some manual editing to prevent averaging across gaps when MAXI/GSC drops off
    # This could be automated BUT I am lazy
    r[32] += 1
    l[33] += 1

    # Group the arrays
    print('Grouping lightcurves to ensure a 3-sigma Swift-BAT detection')
    group_array = group_days(matched_array, l, r)

    # Plot the binned data overtop the unbinned data
    print('Plotting lightcurves')
    plot_binned_lc(group_array, matched_array, l, r)

    # Calculate the HR ratios for each time bin

    print('Calculating HR ratios')
    HR = []
    HR_neg = []
    HR_pos = []
    for k in range(len(group_array[2])):
        output = MCMC_HR_fit(group_array[-4,k],group_array[-2,k],group_array[-3,k],group_array[-1,k])
        HR.append(output[0])
        HR_neg.append(output[1])
        HR_pos.append(output[2])

    # Save the HR_fits as a text file
    np.savetxt('../results/HR_fits.txt', np.array([group_array[0], HR, HR_neg, HR_pos, group_array[2]]).T)

if __name__=="__main__":
    main()
