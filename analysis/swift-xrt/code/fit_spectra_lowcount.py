import time, json, glob, os, argparse
import numpy as np
from astropy.io import fits
from xspec import *

def main():
    
    # Load in the Parameters (it's just the minimum energy)
    parser = argparse.ArgumentParser()
    parser.add_argument('--emin', \
                    help='minimum fit energy [keV]', \
                    dest='emin',  default=0.5)
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
    Fit.statMethod = "cstat"
    Plot.device = '/null'
    Plot.yLog = True
    Plot.xAxis ='keV'  
    Xset.parallel.error = 10

    # Load in the Spectra
    spectra_array = sorted(glob.glob('../spectra/*final.pi'))[-2:]
    times = []
    
    # Load in the best parameters
    with open('../results/fit_dict_full_Emin{}.json'.format(emin),'r') as jfile:
        full_fit = json.load(jfile)

    gamma = np.average(full_fit['gamma'], weights = np.amax((full_fit['gamma_neg'],full_fit['gamma_pos']),axis=0) ** (-2))
    nh    = full_fit['nh'][0] * 1e-22

    # Extract the values
    chi       = []
    dof       = []
    flux      = []
    flux_neg  = []
    flux_pos  = []

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
        Model('tbabs * (pegpwrlw)')
        AllModels(1).setPars({1:'{} -1'.format(nh), #TBABS
               2:'{} -1'.format(gamma), 3:emin, 4:emax, 5:'1', #pewgpwr_1
        })

        # Perform the fit
        AllModels.systematic = 0.03
        for delt in [1e-2,1e-3,1e-4,1e-5]:
            Fit.delta = delt
            Fit.perform()

        print('\n %i' %i,spectrum) 
        Fit.error("5")
        print('Chi/dof = %.2f/%.2f ~ %.2f' %(Fit.testStatistic, Fit.dof, Fit.testStatistic/Fit.dof))

        k = 1
        chi.append(Fit.testStatistic)
        dof.append(Fit.dof)

        flux.append(AllModels(k)(5).values[0] * 1e-12)
        flux_neg.append(AllModels(k)(5).values[0] * 1e-12 - AllModels(k)(5).error[0] * 1e-12)
        flux_pos.append(AllModels(k)(5).error[1] * 1e-12 - AllModels(k)(5).values[0] * 1e-12) 

        # Plot the spectra + fits
        Plot.commands = ()
        #Plot.addCommand("rescale y1 0.01 11.0")
        Plot.addCommand("rescale x1 0.3 11.0")
        Plot.addCommand("label top %s" %spectrum.split('output/')[-1])
        Plot.addCommand(f"hardcopy ../images/lowcount_Spectra{i+1}.ps/ps")
        Plot("data", "chi")

    # Combine the plotted ps files into a single file and remove the single epoch files
    ps_command = ' > ../images/lowcount_JointPSFile.ps'
    for i in np.array(range(len(spectra_array)))[::-1]:
        ps_command = f'../images/lowcount_Spectra{i+1}.ps '+ ps_command
    ps_command = 'psmerge ' + ps_command
    os.system(ps_command)
    os.system('rm -r ../images/lowcount_Spectra*')

    # Save an output dictionary
    out_dict = {'times':times,'flux':flux, 'flux_neg':flux_neg, 'flux_pos':flux_pos, 'chi':chi, 'dof':dof}

    with open('../results/fit_dict_lowcount_Emin{}.json'.format(emin),'w') as jfile:
        json.dump(out_dict, jfile)

if __name__ == '__main__':
    main()
