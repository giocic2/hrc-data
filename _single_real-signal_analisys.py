import easygui
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

# ANALYSIS SETTINGS
SAMPLING_FREQUENCY = 100e3 # According to "hrc-ps.py" script
FFT_RESOL = 1 # Hz
SMOOTHING_WINDOW = 10 # Hz
FREQUENCY_MIN = -50_000 # Hz (lower limit for doppler centroid estimation)
BANDWIDTH_THRESHOLD = 6 # dB

# FFT bins and resolution
freqBins_FFT = int(2**np.ceil(np.log2(abs(SAMPLING_FREQUENCY/2/FFT_RESOL))))
print('FFT resolution: ' + str(SAMPLING_FREQUENCY/freqBins_FFT) + ' Hz')
print('FFT bins: ' + str(freqBins_FFT))
smoothingBins = int(round(SMOOTHING_WINDOW / (SAMPLING_FREQUENCY / freqBins_FFT)))
print('Size of smoothing window (moving average): ' + str(smoothingBins) + ' bins')

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
FFT = np.fft.rfft(voltageAxis_mV, n = freqBins_FFT ) # FFT of real signal
FFT_mV = np.abs(2/(totalSamples)*FFT) # FFT magnitude
FFT_max = np.amax(FFT_mV)
FFT_dBV = 20*np.log10(FFT_mV/1000)
freqAxis = np.fft.rfftfreq(freqBins_FFT ) # freqBins/2+1
freqAxis_Hz = freqAxis * SAMPLING_FREQUENCY
peakFreq = freqAxis_Hz[FFT_mV.argmax()]
# Plot FFT
plt.plot(freqAxis_Hz, FFT_dBV)
plt.ylabel('Spectrum magnitude (dBV)')
plt.xlabel('Frequency (Hz)')
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
        else:
            start_detected = False
        freqIndex += 1
        if startBand < peakFreq and stopBand > peakFreq:
            centroidDetected = True
            break

print('Detected Doppler frequency: {:.1f}'.format(peakFreq) + ' Hz')
print('Amplitude of this FFT peak: {:.1f}'.format(20*np.log10(FFT_max/1000)) + ' dBV')
print('Bandwidth threshold: {:.1f}'.format(BANDWIDTH_THRESHOLD) + ' dB')
print('Bandwidth: {:.1f}'.format(stopBand - startBand) + ' Hz')
print('Bandwidth starts at {:.1f}'.format(startBand) + ' Hz')
print('Bandwidth stops at {:.1f}'.format(stopBand) + ' Hz')

# Plot FFT: normalized and smoothed
plt.plot(freqAxis_Hz, FFT_norm_dB)
plt.plot(freqAxis_Hz, FFT_norm_dB_smooth)
plt.ylabel('Spectrum magnitude (dB)')
plt.xlabel('Frequency (Hz)')
plt.grid(True)
plt.show()

# Spectrogram computation
f, t, Sxx = signal.spectrogram(voltageAxis_mV, fs = SAMPLING_FREQUENCY, noverlap=128 ,nperseg = 1024, nfft = 2**14, scaling = 'spectrum')
plt.pcolormesh(t, f, Sxx, shading='gouraud')
plt.ylabel('Frequency [Hz]')
plt.xlabel('Time [sec]')
plt.axis([0, 1, 0, 1000])
plt.show()