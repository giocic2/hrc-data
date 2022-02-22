import easygui
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

# ANALYSIS SETTINGS
SAMPLING_FREQUENCY = 100e3 # According to "hrc-ps.py" script
FFT_RESOL = 1 # Hz
SMOOTHING_WINDOW = 10 # Hz
BANDWIDTH_THRESHOLD = 6 # dB
ZERO_FORCING = True # Enable forcing FFT to zero, everywhere except between FREQUENCY_MIN and FREQUENCY_MAX
FREQUENCY_MIN = 50 # Hz
FREQUENCY_MAX = 1_000 # Hz

# FFT bins and resolution
freqBins_FFT = int(2**np.ceil(np.log2(abs(SAMPLING_FREQUENCY/2/FFT_RESOL))))
print('FFT resolution: ' + str(SAMPLING_FREQUENCY/freqBins_FFT) + ' Hz')
print('FFT bins: ' + str(freqBins_FFT))
smoothingBins = int(round(SMOOTHING_WINDOW / (SAMPLING_FREQUENCY / freqBins_FFT)))
print('Size of smoothing window (moving average): ' + str(smoothingBins) + ' bins')
minBin = int(np.round(FREQUENCY_MIN / (SAMPLING_FREQUENCY/freqBins_FFT)))
FREQUENCY_MIN = minBin * SAMPLING_FREQUENCY/freqBins_FFT
print("Minimum frequency of interest: {:.1f} Hz".format(FREQUENCY_MIN))
maxBin = int(np.round(FREQUENCY_MAX / (SAMPLING_FREQUENCY/freqBins_FFT)))
FREQUENCY_MAX = maxBin * SAMPLING_FREQUENCY/freqBins_FFT
print("Maximum frequency of interest: {:.1f} Hz".format(FREQUENCY_MAX))

filename = None
while filename == None:
    filename = easygui.fileopenbox(title = "Choose *.csv file to analyse...", default = "*.csv")
print(filename)

rawSamples = np.genfromtxt(filename, delimiter = ',')

voltageAxis_mV = rawSamples[:,0]
timeAxis_s = rawSamples[:,1]
totalSamples = timeAxis_s.size

# Time-domain plot
plt.plot(timeAxis_s, voltageAxis_mV)
plt.ylabel('Voltage (mV)')
plt.xlabel('Time (s)')
plt.grid(True)
plt.show()

# FFT computation
voltageAxis_mV_win = voltageAxis_mV * np.hamming(totalSamples)
FFT = np.fft.rfft(voltageAxis_mV_win, n = freqBins_FFT ) # FFT of real signal
FFT_mV = np.abs(2/(totalSamples)*FFT) # FFT magnitude
if ZERO_FORCING == True:
    FFT_mV[0:minBin] = 0
    FFT_mV[maxBin:-1] = 0
FFT_max = np.amax(FFT_mV)
FFT_dBV = 20*np.log10(FFT_mV/1000)
freqAxis = np.fft.rfftfreq(freqBins_FFT ) # freqBins/2+1
freqAxis_Hz = freqAxis * SAMPLING_FREQUENCY
peakFreq = freqAxis_Hz[FFT_mV.argmax()]
# Plot FFT
plt.plot(freqAxis_Hz, FFT_dBV)
plt.ylabel('Spectrum magnitude (dBV)')
plt.xlabel('Frequency (Hz)')
if ZERO_FORCING == True:
    plt.xlim(FREQUENCY_MIN, FREQUENCY_MAX)
plt.grid(True)
plt.show()

# FFT computation: normalized and smoothed 
FFT_norm = FFT_mV / FFT_max
FFT_norm_dB = 20*np.log10(FFT_norm)
FFT_norm_smooth = np.convolve(FFT_norm, np.ones(smoothingBins), 'same') / smoothingBins
FFT_norm_dB_smooth = 20*np.log10(FFT_norm_smooth)
# Doppler centroid bandwidth
FFT_norm_dB_max = np.amax(FFT_norm_dB)
FFT_norm_dB_smooth_max = np.amax(FFT_norm_dB_smooth)
peakFreq = freqAxis_Hz[FFT_norm_dB.argmax()] # If two identical maxima, only the first occurrence is shown (negative frequency)
freqIndex = 0
stopIndex = 0
start_detected = False
startBand = 0
stopBand = 0
centroidDetected = False
while centroidDetected == False:
    for element in FFT_norm_dB_smooth:
        if element >= (FFT_norm_dB_smooth_max - BANDWIDTH_THRESHOLD):
            if start_detected == False:
                startBand = max(freqAxis_Hz[freqIndex],FREQUENCY_MIN)
                start_detected = True
            stopIndex = max(stopIndex,freqIndex)
            stopBand = freqAxis_Hz[stopIndex]
        freqIndex += 1
        if freqIndex >= (freqBins_FFT/2+1):
            centroidDetected = True
            break

print('Amplitude of FFT peak (norm.smooth.): {:.1f}'.format(FFT_norm_dB_smooth_max) + ' dB')
print('Bandwidth threshold (norm.smooth.): {:.1f}'.format(FFT_norm_dB_smooth_max - BANDWIDTH_THRESHOLD) + ' dB')
print('Bandwidth: {:.1f}'.format(stopBand - startBand) + ' Hz')
print('Bandwidth starts at {:.1f}'.format(startBand) + ' Hz')
print('Bandwidth stops at {:.1f}'.format(stopBand) + ' Hz')
print('Center of Doppler centroid: {:.1f}'.format((stopBand + startBand)/2) + ' Hz')

# Plot FFT: normalized and smoothed
plt.plot(freqAxis_Hz, FFT_norm_dB)
plt.plot(freqAxis_Hz, FFT_norm_dB_smooth)
plt.ylabel('Spectrum magnitude (dB)')
plt.xlabel('Frequency (Hz)')
if ZERO_FORCING == True:
    plt.xlim(FREQUENCY_MIN, FREQUENCY_MAX)
plt.grid(True)
plt.show()

# Spectrogram computation
f, t, Sxx = signal.spectrogram(voltageAxis_mV, fs = SAMPLING_FREQUENCY, noverlap=128 ,nperseg = 1024, nfft = 2**14, scaling = 'spectrum', detrend='constant')
plt.pcolormesh(t, f, Sxx, shading='gouraud')
plt.ylabel('Frequency [Hz]')
plt.xlabel('Time [sec]')
plt.axis([0, 1, 0, 1000])
plt.show()