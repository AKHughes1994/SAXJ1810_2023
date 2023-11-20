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

    # Initialize PyXspec parameters
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

    # Load in the Spectra
    spectra_array = sorted(glob.glob('../spectra/*final.pi'))[:-2]

    # Extract the values
    times = []
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

        # Initialize the model
        Model('tbabs * (pegpwrlw + cflux * bbody)')
        AllModels(1).setPars({1:'0.39', #TBABS
               2:'1.7,,0.0,0.0,5.0,5.0', 3:emin, 4:emax,5:'100', #pewgpwr_1
               6:emin, 7:emax, 8:'-10.', 9:'1.0,,0.05,0.05,5.0,5.0', 10:'1.0 -1' #thermal
        })

        # Perform the fit
        AllModels.systematic = 0.03
        for delt in [1e-2,1e-3,1e-4,1e-5]:
            Fit.delta = delt
            Fit.perform()

        print('\n {}'.format(i),spectrum) 
        Fit.error("1,2,5,8,9")
        print('Chi/dof = {:.2f}/{:.2f} ~ {:.2f}, nh = {:.2f}, Gamma = {:.2f}'.format(Fit.statistic, Fit.dof, Fit.statistic/Fit.dof, AllModels(1)(1).values[0], AllModels(1)(2).values[0]))

        k = 1
        chi.append(Fit.statistic)
        dof.append(Fit.dof)

        nh.append(AllModels(k)(1).values[0] * 1e22)
        nh_neg.append(AllModels(k)(1).values[0] * 1e22 - AllModels(k)(1).error[0] * 1e22)
        nh_pos.append(AllModels(k)(1).error[1] * 1e22 - AllModels(k)(1).values[0] * 1e22) 

        gamma.append(AllModels(k)(2).values[0])
        gamma_neg.append(AllModels(k)(2).values[0] - AllModels(k)(2).error[0])
        gamma_pos.append(AllModels(k)(2).error[1] - AllModels(k)(2).values[0])

        flux.append(AllModels(k)(5).values[0] * 1e-12)
        flux_neg.append(AllModels(k)(5).values[0] * 1e-12 - AllModels(k)(5).error[0] * 1e-12)
        flux_pos.append(AllModels(k)(5).error[1] * 1e-12 - AllModels(k)(5).values[0] * 1e-12) 

        kt.append(AllModels(k)(9).values[0])
        kt_neg.append(AllModels(k)(9).values[0] - AllModels(k)(9).error[0])
        kt_pos.append(AllModels(k)(9).error[1] - AllModels(k)(9).values[0])

        bb_flux.append(10 ** AllModels(k)(8).values[0])
        bb_flux_neg.append(np.abs(AllModels(k)(8).values[0] - AllModels(k)(8).error[0]) * np.abs(np.log10(10) * 10 ** AllModels(k)(8).values[0]))
        bb_flux_pos.append(np.abs(AllModels(k)(8).error[1] - AllModels(k)(8).values[0]) * np.abs(np.log10(10) * 10 ** AllModels(k)(8).values[0])) 

        # Plot the spectra + fits
        Plot.commands = ()
        Plot.addCommand("rescale y1 0.01 11.0")
        Plot.addCommand("rescale x1 0.3 11.0")
        Plot.addCommand("label top {}".format(spectrum.split('output/')[-1]))
        Plot.addCommand(f"hardcopy ../images/individual_Spectra{i+1}.ps/ps")
        Plot("data", "chi")

    # Combine the plotted ps files into a single file and remove the single epoch files
    ps_command = ' > ../images/individual_JointPSFile_Emin{}.ps'.format(emin)
    for i in np.array(range(len(spectra_array)))[::-1]:
        ps_command = f'../images/individual_Spectra{i+1}.ps '+ ps_command
    ps_command = 'psmerge ' + ps_command
    os.system(ps_command)
    os.system('rm -r ../images/individual_Spectra*')

    # Save an output dictionary
    out_dict = {'times':times,
    'nh':nh, 'nh_neg':nh_neg, 'nh_pos':nh_pos,
    'gamma':gamma, 'gamma_neg':gamma_neg, 'gamma_pos':gamma_pos,
    'flux':flux, 'flux_neg':flux_neg, 'flux_pos':flux_pos,
    'kt':kt, 'kt_neg':kt_neg, 'kt_pos':kt_pos,
    'bb_flux':bb_flux, 'bb_flux_neg':bb_flux_neg, 'bb_flux_pos':bb_flux_pos, 'chi':chi, 'dof':dof}

    with open('../results/fit_dict_individual_Emin{}.json'.format(emin),'w') as jfile:
        json.dump(out_dict, jfile)

if __name__ == '__main__':
    main()
