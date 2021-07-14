import os
import numpy as np
import re

tiltAngles = []
directions = []
filenames = []

print('List of *.csv files in the current directory: ')
with os.scandir(path='.') as directoryScan:
    for entry in directoryScan:
        if entry.name.endswith('.csv') and entry.is_file():
            filenames.append(entry.name) # Store filenames in list
            print(entry.name)
            angle_substring = re.search('__(.+?)deg__', entry.name)
            if angle_substring and (float(angle_substring.group(1)) not in tiltAngles):
                tiltAngles.append(float(angle_substring.group(1))) # Store tilt angles in list
            direction_substring = re.search('deg__(.+?)MHz__', entry.name)
#                 direction_substring = re.search('m__(.+?)MHz__', entry.name)
            if direction_substring and (int(direction_substring.group(1)) not in directions):
                directions.append(int(direction_substring.group(1))) # Store TX frequencies in list
                    
print('Grid made of ' + str(len(directions)) + 'x' + str(len(tiltAngles)) + ' points (directions x angles).')
print('Directions (TX frequencies): ' + str(directions))
print('Tilt angles: ' + str(tiltAngles))

FFTpeaks = np.ndarray((len(directions), len(tiltAngles)))

FFT_FREQ_BINS = 2**20
SAMPLING_FREQUENCY = 100e3 # According to "hrc-ps.py" script
print('FFT resolution: ' + str(SAMPLING_FREQUENCY/FFT_FREQ_BINS) + ' Hz')

directionIndex = 0
tiltAngleIndex = 0

IFI = False
IFQ = False

for filename in filenames:
    if (str(directions[directionIndex]) in filename) and (str(tiltAngles[tiltAngleIndex]) in filename):
        if ('ChA' in filename):
            rawSamples = np.genfromtxt(filename, delimiter = ',')
            IFI_mV = rawSamples[:,0]
            IFI = True
        elif ('ChB' in filename):
            rawSamples = np.genfromtxt(filename, delimiter = ',')
            IFQ_mV = rawSamples[:,0]
            IFQ = True
        else:
            print('Something wrong.')
        if IFI == True and IFQ == True:
            complexSignal_mV = np.array(len(IFI_mV))
            complexSignal_mV = np.add(IFI_mV, 1j*IFQ_mV)
            timeAxis_s = rawSamples[:,1]
            totalSamples = timeAxis_s.size
            # FFT computation
            FFT = np.fft.fft(complexSignal_mV, n = FFT_FREQ_BINS) # FFT of real signal
            FFT_mV = np.abs(1/(totalSamples)*FFT) # FFT magnitude
            FFT_dBV = 20*np.log10(FFT_mV/1000)
            freqAxis = np.fft.fftfreq(FFT_FREQ_BINS) # freqBins+1
            freqAxis_Hz = freqAxis * SAMPLING_FREQUENCY
            FFTpeaks[directionIndex, tiltAngleIndex] = np.amax(FFT_dBV)
            IFI = False
            IFQ = False
            directionIndex += 1
        if directionIndex >= len(directions):
            if tiltAngleIndex < (len(tiltAngles)-1):
                tiltAngleIndex += 1
            directionIndex = 0