import time, json, glob, os, argparse
import numpy as np
from astropy.io import fits
from xspec import *


def main():
    # Initialize XSpec
    Xset.abund = "wilm"    # Abundance of elements
    Xset.xsect = "vern"    # Cross-section
    Fit.query = "yes"      # default response to XSpec queries
    Xset.chatter = 0
    Fit.nIterations = 100  # fit iterations for each attemp
    Plot.xAxis = "KeV"     # X-axis for plotting set to energy instead of channel
    Fit.statTest = "chi"
    Fit.statMethod = "chi"
    Plot.device = '/null'
    Plot.yLog = True
    Plot.xAxis ='keV'
    Xset.parallel.error = 10

    # Array of the best fit value
    chi       = []
    dof       = []
    times     = []

    gamma     = []
    gamma_neg = []
    gamma_pos = []
    flux      = []
    flux_neg  = []
    flux_pos  = []

    bb_flux      = []
    bb_flux_neg  = []
    bb_flux_pos  = []
    bb_kt        = []
    bb_kt_neg    = []
    bb_kt_pos    = []

    disk_flux      = []
    disk_flux_neg  = []
    disk_flux_pos  = []
    disk_kt        = []
    disk_kt_neg    = []
    disk_kt_pos    = []

    # Load in the Spectra
    spectra_array = sorted(glob.glob('../spectra/*final.pi'))
    spectra_array = [spectra_array[2], spectra_array[6], spectra_array[11]] # These are the three epochs of interest


    # Load in the best parameters
    with open('../results/fit_dict_full_Emin0.5.json','r') as jfile:
        full_fit = json.load(jfile)
    nh    = full_fit['nh'][0] * 1e-22
    print(nh)

    complex_model = 'else' # As long as its not '1' or '2' this will perform the 3-component fitting 

    for i, spectrum in enumerate(spectra_array[:]):
        hdul = fits.open(spectrum)
        times.append(hdul[0].header['DATE-OBS'])  
        hdul.close()

        # Load in the spectrum  -- Trim bad energy values
        AllData(spectrum)
        AllData.ignore('*:10.0-**')
        AllData.ignore('*:**-0.5')
        AllData.notice('*:0.5-10.0')
        AllData.ignore('bad')

        Emin = 0.5
        Emax = 10.0

        if complex_model == '1':
            # Initialize the model
            Model('tbabs * (pegpwrlw + cflux * bbody)')
            AllModels(1).setPars({1:'{} -1'.format(nh), #TBABS
                2:'1.7,,0.0,0.0,5.0,5.0',3:Emin,4:Emax,5:'100', #pewgpwr_1
                6:Emin, 7:Emax, 8:'-10.', 9:'1.0,,0.05,0.05,5.0,5.0', 10:'1.0 -1' #thermal
            })

            # Perform the fit
            AllModels.systematic = 0.03
            for delt in [1e-2,1e-3,1e-4,1e-5]:
                Fit.delta = delt
                Fit.perform()

            Fit.error("2,5,8,9")

        elif complex_model == '2':
            # Initialize the model
            Model('tbabs * (cflux * bbody + cflux * diskbb)')
            AllModels(1).setPars({1:'{} -1'.format(nh), #TBABS
               2:Emin, 3:Emax, 4:'-10.', 5:'1.0,,0.05,0.05,5.0,5.0', 6:'1.0 -1', #thermal BB
               7:Emin, 8:Emax, 9:'-10.', 10:'0.5,,0.05,0.05,5.0,5.0', 11:'1.0 -1' #thermal DiskBB
            })

            # Perform the fit
            AllModels.systematic = 0.03
            for delt in [1e-2,1e-3,1e-4,1e-5]:
                Fit.delta = delt
                Fit.perform()

            Fit.error("4,5,9,10")

        else:
            # Initialize the model
            Model('tbabs * (pegpwrlw + cflux * bbody + cflux * diskbb)')
            AllModels(1).setPars({1:'{} -1'.format(nh), #TBABS
               2:'1.7,,0.0,0.0,5.0,5.0',3:Emin,4:Emax,5:'100', #pewgpwr_1
               6:Emin, 7:Emax, 8:'-10.', 9:'1.0,,0.05,0.05,5.0,5.0', 10:'1.0 -1', #thermal BB
               11:Emin, 12:Emax, 13:'-10.', 14:'0.5,,0.05,0.05,5.0,5.0', 15:'1.0 -1' #thermal DiskBB
            })

            # Perform the fit
            AllModels.systematic = 0.03
            for delt in [1e-2,1e-3,1e-4,1e-5]:
                Fit.delta = delt
                Fit.perform()

            Fit.error("2,5,8,9,13,14")
        print(spectrum)
        print('Chi/dof = %.2f/%.2f ~ %.2f' %(Fit.statistic, Fit.dof, Fit.statistic/Fit.dof))
    
        k=1
        # Append the values of interest
        gamma.append(AllModels(k)(2).values[0])
        gamma_neg.append(AllModels(k)(2).values[0] - AllModels(k)(2).error[0])
        gamma_pos.append(AllModels(k)(2).error[1] - AllModels(k)(2).values[0])

        flux.append(AllModels(k)(5).values[0] * 1e-12)
        flux_neg.append(AllModels(k)(5).values[0] * 1e-12 - AllModels(k)(5).error[0] * 1e-12)
        flux_pos.append(AllModels(k)(5).error[1] * 1e-12 - AllModels(k)(5).values[0] * 1e-12) 

        bb_kt.append(AllModels(k)(9).values[0])
        bb_kt_neg.append(AllModels(k)(9).values[0] - AllModels(k)(9).error[0])
        bb_kt_pos.append(AllModels(k)(9).error[1] - AllModels(k)(9).values[0])

        bb_flux.append(10 ** AllModels(k)(8).values[0])
        bb_flux_neg.append(np.abs(AllModels(k)(8).values[0] - AllModels(k)(8).error[0]) * np.abs(np.log10(10) * 10 ** AllModels(k)(8).values[0]))
        bb_flux_pos.append(np.abs(AllModels(k)(8).error[1] - AllModels(k)(8).values[0]) * np.abs(np.log10(10) * 10 ** AllModels(k)(8).values[0])) 

        disk_kt.append(AllModels(k)(14).values[0])
        disk_kt_neg.append(AllModels(k)(14).values[0] - AllModels(k)(14).error[0])
        disk_kt_pos.append(AllModels(k)(14).error[1] - AllModels(k)(14).values[0])

        disk_flux.append(10 ** AllModels(k)(13).values[0])
        disk_flux_neg.append(np.abs(AllModels(k)(13).values[0] - AllModels(k)(13).error[0]) * np.abs(np.log10(10) * 10 ** AllModels(k)(13).values[0]))
        disk_flux_pos.append(np.abs(AllModels(k)(13).error[1] - AllModels(k)(13).values[0]) * np.abs(np.log10(10) * 10 ** AllModels(k)(13).values[0])) 

        chi.append(Fit.statistic)
        dof.append(Fit.dof)

    out_dict = {'times':times, 'gamma':gamma, 
    'gamma_neg':gamma_neg, 'gamma_pos':gamma_pos,'flux':flux, 'flux_neg':flux_neg, 'flux_pos':flux_pos,
    'bb_kt':bb_kt, 'bb_kt_neg':bb_kt_neg, 'bb_kt_pos':bb_kt_pos,'bb_flux':bb_flux, 'bb_flux_neg':bb_flux_neg, 'bb_flux_pos':bb_flux_pos,
    'disk_kt':disk_kt, 'disk_kt_neg':disk_kt_neg, 'disk_kt_pos':disk_kt_pos,'disk_flux':disk_flux, 'disk_flux_neg':disk_flux_neg, 'disk_flux_pos':disk_flux_pos,
    'chi_tot':Fit.statistic, 'dof_tot':Fit.dof, 'chi':chi, 'dof':dof, 'Emin':AllModels(1)(3).values[0], 'Emax':AllModels(1)(4).values[0]}

    with open('../results/fit_dict_3comp.json','w') as jfile:
        json.dump(out_dict, jfile)

if __name__ == '__main__':
    main()
