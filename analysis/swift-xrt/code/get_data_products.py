from swifttools.xrt_prods import XRTProductRequest
import time, json, glob, os, argparse
import numpy as np

#Load in the arguments
parser = argparse.ArgumentParser(description='Read in the segments')
parser.add_argument('--segments', nargs='+', help='Pass the segments separated by a space')
parser.add_argument('--centroid',    action='store_true', help='Centroiding')
parser.add_argument('--no-centroid', action='store_false', dest='centroid', help='No Centroiding')
parser.add_argument('--targetID', default='00032459')
args = parser.parse_args()
segments = np.array(args.segments).astype(int)
centroid = args.centroid
targetID = str(args.targetID)

#Get current directory
working_directory = os.getcwd().split('code')[0]

#Make directory for pipeline output if one doesn't exist
if os.path.exists(working_directory + 'pipeline_output/') == False:
    print('# Making directory for Pipeline Output"')
    os.system('mkdir {}pipeline_output/'.format(working_directory))

#Format the segment list  to feed into the Pipeline commands
ObsID=''
for seg in segments:
    ObsID = ObsID + ',%s%03d' %(targetID,seg)
ObsID = ObsID[1:] #Remove leading comma

#Get the data products from the Swifttools pipeline
print('# Running the Swift XRT pipeline #')
myRequest = XRTProductRequest('hughes1@ualberta.ca', silent=False)

#Set the global Parameters
if centroid == True:
    myRequest.setGlobalPars(name='SAX J1810.8-2609',
                        targ=targetID,
                        SinceT0=False,
                        RA=272.68549,
                        Dec=-26.15030,
                        centroid=centroid,
                        useSXPS=False,
                        poserr=1)
else:
    myRequest.setGlobalPars(name='SAX J1810.8-2609',
                        targ=targetID,
                        SinceT0=False,
                        RA=272.68549,
                        Dec=-26.15030,
                        centroid=centroid,
                        useSXPS=False)


#Define the parameters to extract a light curve -- times are in seconds
myRequest.addSpectrum(whichData='user',
                        useObs=ObsID,
                        hasRedshift=False,
                        timeslice='obsid')

# Check if valid
if(myRequest.isValid()[0]==True):
    if not myRequest.submit():
        print (f"I couldn't submit error:{myRequest.submitError}")

    while not myRequest.complete:
        time.sleep(10)
    # Retrieve the products
    dictSpec = myRequest.retrieveSpectralFits()
else:
    print("BAD REQUEST:", myRequest.isValid()[1])
    exit()

print('# Extracting the pipeline output tar files #')
observations = [key for key in dictSpec.keys() if 'Obs' in key]
for observation in observations:
    pipeline_output_file = dictSpec[observation]['DataFile']
    os.system('curl {} -o {}spectra/{}.tar.gz'.format(pipeline_output_file, working_directory, observation))

print('# De-Tarring the files #')
for Obs in ObsID.split(','):
    os.system(f'tar -xf {working_directory}spectra/Obs_{Obs}.tar.gz -C {working_directory}spectra/')
