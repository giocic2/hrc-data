import easygui
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.fft import fftshift

# ANALYSIS SETTINGS
SAMPLING_FREQUENCY = 2e3 # According to "hrc-ps.py" script
FFT_RESOL = 1 # Hz
ACQUISITION_TIME = 10 # s
SMOOTHING_WINDOW = 10 # Hz
BANDWIDTH_THRESHOLD = 6 # dB
ZERO_FORCING = True # Enable forcing FFT to zero, everywhere except between FREQUENCY_MIN and FREQUENCY_MAX
FREQUENCY_MIN = -500 # Hz
FREQUENCY_MAX = 500 # Hz
SPECTROGRAM = True
OFFSET_COMPENSATION = False

# FFT bins and resolution
freqBins_FFT = int(2**np.ceil(np.log2(abs(SAMPLING_FREQUENCY/2/FFT_RESOL))))
print('FFT resolution: ' + str(SAMPLING_FREQUENCY / freqBins_FFT) + ' Hz')
print('FFT bins: ' + str(freqBins_FFT))
smoothingBins = int(round(SMOOTHING_WINDOW / (SAMPLING_FREQUENCY / freqBins_FFT)))
print('Size of smoothing window (moving average): ' + str(smoothingBins) + ' bins')
minBin = int(freqBins_FFT/2 + np.round(FREQUENCY_MIN / (SAMPLING_FREQUENCY/freqBins_FFT)))
FREQUENCY_MIN = -SAMPLING_FREQUENCY/2 + minBin * SAMPLING_FREQUENCY/freqBins_FFT
print("Minimum frequency of interest: {:.1f} Hz".format(FREQUENCY_MIN))
maxBin = int(freqBins_FFT/2 + np.round(FREQUENCY_MAX / (SAMPLING_FREQUENCY/freqBins_FFT)))
FREQUENCY_MAX = -SAMPLING_FREQUENCY/2 + maxBin * SAMPLING_FREQUENCY/freqBins_FFT
print("Maximum frequency of interest: {:.1f} Hz".format(FREQUENCY_MAX))

filename_IFI = None
filename_IFQ = None
while filename_IFI == None:
    filename_IFI = easygui.fileopenbox(title = "Choose IFI *.csv file...", default = "*.csv")
print(filename_IFI)
while filename_IFQ == None:
    filename_IFQ = easygui.fileopenbox(title = "Choose IFQ *.csv file...", default = "*.csv")
print(filename_IFQ)

rawSamples_IFI = np.genfromtxt(filename_IFI, delimiter = ',')
rawSamples_IFQ = np.genfromtxt(filename_IFQ, delimiter = ',')

voltageAxis_IFI_mV = rawSamples_IFI[:,0]
voltageAxis_IFQ_mV = rawSamples_IFQ[:,0]
if OFFSET_COMPENSATION == True:
    voltageAxis_IFI_mV = voltageAxis_IFI_mV - np.mean(voltageAxis_IFI_mV)
    voltageAxis_IFQ_mV = voltageAxis_IFQ_mV - np.mean(voltageAxis_IFQ_mV)
timeAxis_s = rawSamples_IFI[:,1]
totalSamples = timeAxis_s.size

# Time-domain plot
plt.plot(timeAxis_s, voltageAxis_IFI_mV)
plt.plot(timeAxis_s, voltageAxis_IFQ_mV)
plt.ylabel('voltage (mV)')
plt.xlabel('time (s)')
plt.grid(True)
plt.legend(['IFI','IFQ'])
plt.show()

# FFT computation
complexSignal_mV = np.add(np.asarray(voltageAxis_IFI_mV), 1j*np.asarray(voltageAxis_IFQ_mV))
complexSignal_mV_win = complexSignal_mV * np.hamming(totalSamples)
FFT = np.fft.fftshift(np.fft.fft(complexSignal_mV_win, n = freqBins_FFT)) # FFT of complex signal
FFT_mV = np.abs(1/(totalSamples)*FFT) # FFT magnitude
if ZERO_FORCING == True:
    FFT_mV[0:minBin] = 0
    FFT_mV[maxBin:-1] = 0
FFT_max = np.amax(FFT_mV)
FFT_dBV = 20*np.log10(FFT_mV/1000)
freqAxis = np.fft.fftshift(np.fft.fftfreq(freqBins_FFT)) # freqBins+1
freqAxis_Hz = freqAxis * SAMPLING_FREQUENCY
# Plot FFT
plt.plot(freqAxis_Hz, FFT_dBV)
plt.ylabel('spectrum magnitude (dBV)')
plt.xlabel('frequency (Hz)')
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
        if freqIndex >= (freqBins_FFT+1):
            centroidDetected = True
            break

print('Amplitude of FFT peak: {:.1f}'.format(np.amax(FFT_dBV)) + ' dBV')
print('Amplitude of this FFT peak (norm.smooth.): {:.1f}'.format(FFT_norm_dB_smooth_max) + ' dB')
print('Bandwidth threshold (norm.smooth.): {:.1f}'.format(FFT_norm_dB_smooth_max - BANDWIDTH_THRESHOLD) + ' dB')
print('Bandwidth: {:.1f}'.format(stopBand - startBand) + ' Hz')
print('Bandwidth starts at {:.1f}'.format(startBand) + ' Hz')
print('Bandwidth stops at {:.1f}'.format(stopBand) + ' Hz')
print('Center of Doppler centroid: {:.1f}'.format((stopBand + startBand)/2) + ' Hz')

# Plot FFT: normalized and smoothed
plt.plot(freqAxis_Hz, FFT_norm_dB)
plt.plot(freqAxis_Hz, FFT_norm_dB_smooth)
plt.ylabel('spectrum magnitude (dB)')
plt.xlabel('frequency (Hz)')
if ZERO_FORCING == True:
    plt.xlim(FREQUENCY_MIN, FREQUENCY_MAX)
plt.grid(True)
plt.show()

# Spectrogram computation
if SPECTROGRAM == True:
    f, t, Sxx = signal.spectrogram(complexSignal_mV, fs = SAMPLING_FREQUENCY, noverlap=0, nperseg = 128, nfft = 2**10, scaling = 'spectrum', return_onesided=False, detrend=False)
    plt.pcolormesh(t, fftshift(f), fftshift(Sxx, axes=0), shading='gouraud')
    plt.ylabel('frequency (Hz)')
    plt.xlabel('time (s)')
    plt.axis([0, ACQUISITION_TIME, FREQUENCY_MIN, FREQUENCY_MAX])
    plt.show()