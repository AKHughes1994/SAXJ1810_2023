import time, json, glob, os, argparse
import numpy as np
from astropy.io import fits
from xspec import *

def main(): 

    # Load in the Parameters (it's just the minimum energy)
    parser = argparse.ArgumentParser()
    parser.add_argument('--emin', \
                    help='minimum fit energy [keV]', \
                    dest='emin', default=0.5)
    parser.add_argument('--emax', \
                    help='minimum fit energy [keV]', \
                    dest='emax', default=10.0)

    args = parser.parse_args()
    emin = float(args.emin)
    emax = float(args.emax)

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

    x = []

    # Load in the Spectra
    spectra_array = sorted(glob.glob('../spectra/*final.pi'))[:-2]
    times = []

    with open('../results/fit_dict_individual_Emin{}.json'.format(emin),'r') as jfile:
        fit_bb = json.load(jfile)

    spectra_string = ''
    for i, spectrum in enumerate(spectra_array[:]):
        hdul = fits.open(spectrum)
        times.append(hdul[0].header['DATE-OBS'])
        spectra_string += f' {i+1}:{i+1} {spectrum}'

    # Load in the spectrum  -- Trim bad energy values
    AllData(spectra_string)
    AllData.ignore('*:10.0-**')
    AllData.ignore('*:**-0.5')
    AllData.notice('*:0.5-10.0')
    AllData.ignore('bad')
    Model('tbabs * (pegpwrlw + cflux * bbody)')

    # Initialize Arrays to store spectra values
    chi       = []
    dof       = []
    gamma     = []
    gamma_neg = []
    gamma_pos = []
    flux      = []
    flux_neg  = []
    flux_pos  = []
    nh        = []
    nh_neg    = []
    nh_pos    = []
    bb_flux      = []
    bb_flux_neg  = []
    bb_flux_pos  = []
    kt           = []
    kt_neg       = []
    kt_pos       = []

    # Initialize the model
    for i in range(len(spectra_array[:])):
        if i + 1 == 1:
            AllModels(i+1).setPars({1:'0.4', #TBABS
              2:'{},,0.0,0.0,5.0,5.0'.format(fit_bb["gamma"][i]),3:emin,4:emax,5:'{}'.format(fit_bb["flux"][i]/1e-12),
              6:emin, 7:emax, 8:'{}'.format(np.log10(fit_bb["bb_flux"][i])), 9:'{},,0.05,0.05,5.0,5.0'.format(fit_bb["kt"][i]), 10:'1.0 -1'})
        else:
            AllModels(i+1).setPars({1:'0.4', #TBABS
                2:'{},,0.0,0.0,5.0,5.0'.format(fit_bb["gamma"][i]),3:emin,4:emax,5:'{}'.format(fit_bb["flux"][i]/1e-12),
                6:emin, 7:emax, 8:'{}'.format(np.log10(fit_bb["bb_flux"][i])), 9:'{},,0.05,0.05,5.0,5.0'.format(fit_bb["kt"][i]), 10:'1.0 -1'})
            AllModels(i+1)(1).link = AllModels(1)(1)

    # Perform the fit
    AllModels.systematic = 0.03
    for delt in [1e-2,1e-3,1e-4,1e-5]:
        print('Fit with delt =',delt)
        Fit.delta = delt
        Fit.perform()

    for i in range(len(spectra_array[:])):
        if AllModels(i+1)(8).values[0] < -20: 
            AllModels(i+1).setPars({8:'-10.3', 9:'0.5'})
            x.append('yes')
        else:
            x.append('no')

    for delt in [1e-2,1e-3,1e-4,1e-5]:
        print('Fit with delt =',delt)
        Fit.delta = delt
        Fit.perform()

    # Record the chi-square value here becasue parallelizing error records the chi 0.0 for each fit
    for i in range(len(spectra_array)):
        k = i + 1
        chi.append(AllData(k).statistic)

    # Calculate 90% confidence intervals for the fits
    print('Calculating Error')
    n = len(spectra_array)
    error_range = '1'
    for i in range(len(spectra_array[:])):
        error_range += ',%s,%s,%s,%s' %(2 + 10*i, 5 + 10*i, 8 + 10*i, 9 + 10*i)
    Fit.error(error_range)

    # Append the spectra parameters
    for i in range(len(spectra_array)):
        k = i + 1

        nh.append(AllModels(k)(1).values[0] * 1e22)
        nh_neg.append(AllModels(k)(1).values[0] * 1e22 - AllModels(k)(1).error[0] * 1e22)
        nh_pos.append(AllModels(k)(1).error[1] * 1e22 - AllModels(k)(1).values[0] * 1e22) 

        gamma.append(AllModels(k)(2).values[0])
        gamma_neg.append(AllModels(k)(2).values[0] - AllModels(k)(2).error[0])
        gamma_pos.append(AllModels(k)(2).error[1] - AllModels(k)(2).values[0])

        flux.append(AllModels(k)(5).values[0] * 1e-12)
        flux_neg.append(AllModels(k)(5).values[0] * 1e-12 - AllModels(k)(5).error[0] * 1e-12)
        flux_pos.append(AllModels(k)(5).error[1] * 1e-12 - AllModels(k)(5).values[0] * 1e-12) 

        if AllModels(k)(8).values[0] > -20.: 
            kt.append(AllModels(k)(9).values[0])
            kt_neg.append(AllModels(k)(9).values[0] - AllModels(k)(9).error[0])
            kt_pos.append(AllModels(k)(9).error[1] - AllModels(k)(9).values[0])

            bb_flux.append(10 ** AllModels(k)(8).values[0])
            bb_flux_neg.append(np.abs(AllModels(k)(8).values[0] - AllModels(k)(8).error[0]) * np.abs(np.log10(10) * 10 ** AllModels(k)(8).values[0]))
            bb_flux_pos.append(np.abs(AllModels(k)(8).error[1] - AllModels(k)(8).values[0]) * np.abs(np.log10(10) * 10 ** AllModels(k)(8).values[0])) 
    
        else:
            kt.append(np.nan)
            kt_neg.append(np.nan)
            kt_pos.append(np.nan)

            bb_flux.append(np.nan)
            bb_flux_neg.append(np.nan)
            bb_flux_pos.append(np.nan)

    # Plot the spectra
    for i, spectrum in enumerate(spectra_array):
        Plot.commands = ()
        Plot.addCommand("rescale y1 0.01 11.0")
        Plot.addCommand("rescale x1 0.3 10.5")
        Plot.addCommand("label top {}".format(spectrum.split('output/')[-1]))
        Plot.addCommand(f"hardcopy ../images/full_Spectra{i+1}.ps/ps")
        Plot(f"{i+1}","data", "chi")
        dof.append(len(Plot.x())-4)

    # Combine the ps files into a single file
    ps_command = ' > ../images/full_JointPSFile_Emin{}.ps'.format(emin)
    for i in np.array(range(len(spectra_array)))[::-1]:
        ps_command = f'../images/full_Spectra{i+1}.ps '+ ps_command
    ps_command = 'psmerge ' + ps_command
    os.system(ps_command)
    os.system('rm -r ../images/full_Spectra*')

    # Save an output dictionary
    out_dict = {'times':times,
    'nh':nh, 'nh_neg':nh_neg, 'nh_pos':nh_pos,
    'gamma':gamma, 'gamma_neg':gamma_neg, 'gamma_pos':gamma_pos,
    'flux':flux, 'flux_neg':flux_neg, 'flux_pos':flux_pos,
    'kt':kt, 'kt_neg':kt_neg, 'kt_pos':kt_pos,
    'bb_flux':bb_flux, 'bb_flux_neg':bb_flux_neg, 'bb_flux_pos':bb_flux_pos,
    'chi_tot':Fit.statistic, 'dof_tot':Fit.dof, 'chi':chi, 'dof':dof, 'Emin':AllModels(1)(3).values[0], 'Emax':AllModels(1)(4).values[0]}

    with open('../results/fit_dict_full_Emin{}.json'.format(emin),'w') as jfile:
        json.dump(out_dict, jfile)

if __name__ == '__main__':
    main()
