import os
import numpy as np
import re
import matplotlib.pyplot as plt

"""
HYPOTHESIS:
- Scan angle between +-15 DEG
- Scan angle linearly varying with TX frequency
"""

MAX_SCAN_ANGLE = np.deg2rad(15)

tiltAngles = []
directions = []
antennaHeights = []
filenames = []
swathsStep = []

print('List of *.csv files in the current directory: ')
with os.scandir(path='.') as directoryScan:
    for entry in directoryScan:
        if entry.name.endswith('.csv') and entry.is_file():
            filenames.append(entry.name) # Store filenames in list
            print(entry.name)
            angle_substring = re.search('__(.+?)deg__', entry.name)
            if angle_substring and (float(angle_substring.group(1)) not in tiltAngles):
                tiltAngles.append(float(angle_substring.group(1))) # Store tilt angles in list
                antennaHeight_substring = re.search('deg__(.+?)m__', entry.name)
                antennaHeights.append(float(antennaHeight_substring.group(1))) # Store antenna height
            direction_substring = re.search('m__(.+?)MHz__', entry.name)
            if direction_substring and (int(direction_substring.group(1)) not in directions):
                directions.append(int(direction_substring.group(1))) # Store TX frequencies in list

tiltAngles = np.asarray(tiltAngles)
directions = np.asarray(directions)
directions_DEG = np.zeros(len(directions))
for index in range(len(directions)):
    directions_DEG[index] = 15 - ((directions[index] - 23500) / 1000) * 30
antennaHeights = np.asarray(antennaHeights)

maxSwaths = np.array(2 * antennaHeights[:] / np.tan(np.deg2rad(tiltAngles[:])) * np.tan(MAX_SCAN_ANGLE))
swathsStep = np.array(maxSwaths[:] / (directions.size - 1))

TARGET_POSITION = np.array([0, 1.25]) # m,m [horiz., vert.]

print('Grid made of ' + str(len(directions)) + 'x' + str(len(tiltAngles)) + ' points (directions x angles):')
grid_horizontal = np.zeros((len(tiltAngles),len(directions_DEG)))
grid_vertical = np.zeros((len(tiltAngles),len(directions_DEG)))

for column in range(len(directions_DEG)):
    for row in range(len(tiltAngles)):
        grid_horizontal[row, column] = antennaHeights[row] / np.tan(np.deg2rad(tiltAngles[row]))
        grid_vertical[row, column] = - (maxSwaths[row] / 2) + column * swathsStep[row]
        
print('Directions (TX frequencies): ' + str(directions))
print('Tilt angles: ' + str(tiltAngles))
print('Swaths: ' + str(np.around(maxSwaths, decimals = 2)))
print('Swaths step: ' + str(np.around(swathsStep, decimals = 2)))

FFTpeaks = np.ndarray((len(directions), len(tiltAngles)))

FFT_FREQ_BINS = 2**20
SAMPLING_FREQUENCY = 100e3 # According to "hrc-ps.py" script
print('FFT resolution: ' + str(SAMPLING_FREQUENCY/FFT_FREQ_BINS) + ' Hz')

directionIndex = 0
tiltAngleIndex = 0

IFI = False
IFQ = False

TARGET_FREQUENCY = 34.4 / 2 * 28 # Hz
ARGMAX_RANGE = 5000 # bins
argmax_startBin = int(round((FFT_FREQ_BINS / (SAMPLING_FREQUENCY) * TARGET_FREQUENCY) - ARGMAX_RANGE / 2))
argmax_endBin = int(round((FFT_FREQ_BINS / (SAMPLING_FREQUENCY) * TARGET_FREQUENCY) + ARGMAX_RANGE / 2))

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
            FFT = np.fft.fft(complexSignal_mV, n = FFT_FREQ_BINS) # FFT of complex signal
            FFT_mV = np.abs(1/(totalSamples)*FFT) # FFT magnitude
            FFT_dBV = 20*np.log10(FFT_mV/1000)
            freqAxis = np.fft.fftfreq(FFT_FREQ_BINS) # freqBins+1
            freqAxis_Hz = freqAxis * SAMPLING_FREQUENCY
#             Plot FFT
#             plt.plot(np.fft.fftshift(freqAxis_Hz), np.fft.fftshift(FFT_dBV))
#             plt.ylabel('Spectrum magnitude (dBV)')
#             plt.xlabel('Frequency (Hz)')
#             plt.grid(True)
#             plt.show()
            # FFTpeaks update
            FFTpeaks[directionIndex, tiltAngleIndex] = np.amax(FFT_dBV[argmax_startBin:argmax_endBin])
            print('{0:.1f}'.format(FFTpeaks[directionIndex, tiltAngleIndex]) + ' dBV', end = ', ')
            print(filename)
            IFI = False
            IFQ = False
            directionIndex += 1
        if directionIndex >= len(directions):
            if (tiltAngleIndex + 1) < len(tiltAngles):
                tiltAngleIndex += 1
            directionIndex = 0
print('')

# Plot grid
figure, (axis1, axis2) = plt.subplots(1,2)
figure.suptitle('Grid scan analysis')
for row in range(len(tiltAngles)):
    for column in range(len(directions_DEG)):
        axis1.plot(grid_vertical[row,column], grid_horizontal[row,column],'bo')
axis1.plot(TARGET_POSITION[0], TARGET_POSITION[1],'go', markersize = 20)
axis1.set_title('Acquisition points and target position')
axis1.set_ylabel('Distance from radar [m]')
axis1.set_xlabel('Swath position [m]')

# Plot FFT peaks
FFTpeaks_rotated = np.rot90(FFTpeaks, k=3) # 3 rotation 90° counter-clockwise
im = axis2.imshow(FFTpeaks_rotated)
axis2.set_xlabel('Beam direction [DEG]')
axis2.set_ylabel('Tilt angle of the radar [DEG]')

# We want to show all ticks...
axis2.set_xticks(np.arange(len(directions_DEG)))
axis2.set_yticks(np.arange(len(tiltAngles)))
# ... and label them with the respective list entries
axis2.set_xticklabels(np.flip(directions_DEG))
axis2.set_yticklabels(tiltAngles)

# Rotate the tick labels and set their alignment.
plt.setp(axis2.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")

# Loop over data dimensions and create text annotations.
for i in range(len(tiltAngles)):
    for j in range(len(directions_DEG)):
        text = axis2.text(j, i, np.around(FFTpeaks_rotated[i, j],decimals=1),
                       ha="center", va="center", color="w")

axis2.set_title("Magnitude of received signal @ target frequency [dBV]")
figure.tight_layout()
plt.show()