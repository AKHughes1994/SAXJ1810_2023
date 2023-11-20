import numpy as np

# Functions from the Hardness analysis that filters/Cleans the light curves for MAXI/GSC and Swift-BAT
def filter_times(lc):
    '''
    Take an input lc array and return
    a time filtered array, input/output arrrays:
        0 -- Time (MJD)
        1 -- Rate (counts/s/cm^2)
        2 -- error on rate (counts/s/cm^2)
    '''
    
    ti = 59340 # May 6 2021
    tf = 59519 # Nov 01 2021
    
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

def main():
    # Load in instrumental light curves (not in Crab units) and filter times and large errors

    # Swift-BAT light curve [15-50keV] -- columns: Time, Rate, Error
    bat  = np.loadtxt('../SAXJ1810_watchdog/BATObs/SAXJ1810.8m2609/SAXJ1810.8m2609_BAT_lc_15_50.txt')
    bat  = filter_errors(filter_times(bat))

    # MAXI/GSC light curve [4-10keV]
    maxi = np.loadtxt('../SAXJ1810_watchdog/MAXIObs/SAXJ1810.8m2609/SAXJ1810.8m2609_MAXI_lc_4_10.txt')
    maxi  = filter_errors(filter_times(maxi))

    # Save filtered lcs
    np.savetxt('../results/filtered_bat.txt',  bat)
    np.savetxt('../results/filtered_maxi.txt', maxi)

if __name__ == "__main__":
    main()
