import os
import numpy as np

print('List of *.csv files in the current directory: ')
with os.scandir(path='.') as directoryScan:
    for entry in directoryScan:
            if entry.name.endswith('.csv') and entry.is_file():
                print(entry.name)

scanDirections = input('Enter number of scanned directions: ')
scanAngles = input('Enter number of scanned angles: ')

FFTpeaks = np.ndarray((scanAngles, scanDirections))

anotherAngle = True

directionIndex = 0
angleIndex = 0

while anotherAngle == True:
    tiltAngle = input('Enter tilt angle to analyse: ')
    tiltAngle += 'deg'
    
    filenames = []
    
    with os.scandir(path='.') as directoryScan:
        for entry in directoryScan:
            if entry.name.endswith('.csv') and entry.is_file() and (tiltAngle in entry.name):
                filenames.append(entry.name)
                print(entry.name)
    
    currentCycle = 0
    
    for entry in filenames:
        rawSamples = np.genfromtxt(entry, delimiter = ',')
        voltageAxis_mV = rawSamples[:,0]
        timeAxis_s = rawSamples[:,1]
        totalSamples = timeAxis_s.size
        FFT_FREQ_BINS = 2**18 # Check "hrc-ps.py" script
        SAMPLING_FREQUENCY = 100e3 # Check "hrc-ps.py" script
        # FFT computation
        FFT = np.fft.rfft(voltageAxis_mV, n = FFT_FREQ_BINS) # FFT of real signal
        FFT_mV = np.abs(2/(totalSamples)*FFT) # FFT magnitude
        FFT_dBV = 20*np.log10(FFT_mV/1000)
        freqAxis = np.fft.rfftfreq(FFT_FREQ_BINS) # freqBins/2+1
        freqAxis_Hz = freqAxis * SAMPLING_FREQUENCY
        FFTpeaks.append(freqAxis_Hz[FFT_dBV.argmax()])
        print("{0:.0f} Hz".format(FFTpeaks[currentCycle]))
        currentCycle += 1

    whatNext = input('Another tilt angle? [y/n]:  ')
    continueTakingUserInput = True
    while continueTakingUserInput == True:
        if whatNext == 'y':
            anotherAngle = True
            continueTakingUserInput = False
        if whatNext == 'n':
            anotherAngle = False
            continueTakingUserInput = False
        if whatNext != 'y' and whatNext != 'n':
            print('Please type correctly [y/n].')
            continueTakingUserInput = True
    

