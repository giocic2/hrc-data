import easygui
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

# ANALYSIS SETTINGS
SAMPLING_FREQUENCY = 100e3 # According to "hrc-ps.py" script
FFT_RESOL = 1 # Hz
SMOOTHING_WINDOW = 10 # Hz

# FFT bins and resolution
freqBins_FFT = int(2**np.ceil(np.log2(abs(SAMPLING_FREQUENCY/2/FFT_RESOL))))
print('FFT resolution: ' + str(SAMPLING_FREQUENCY / freqBins_FFT) + ' Hz')
print('FFT bins: ' + str(freqBins_FFT))
smoothingBins = int(round(SMOOTHING_WINDOW / (SAMPLING_FREQUENCY / freqBins_FFT)))
print('Size of smoothing window (moving average): ' + str(smoothingBins) + ' bins')

filename_IFI = None
filename_IFQ = None
while filename_IFI == None:
    filename_IFI = easygui.fileopenbox(title = "Choose IFI *.csv file...", default = "*.csv")
print(filename_IFI)
while filename_IFQ == None:
    filename_IFQ = easygui.fileopenbox(title = "Choose IFQ *.csv file...", default = "*.csv")
print(filename_IFQ)

rawSamples_IFI = np.genfromtxt(filename_IFI, delimiter = ',')
rawSamples_IFQ = np.genfromtxt(filename_IFI, delimiter = ',')

voltageAxis_IFI_mV = rawSamples_IFI[:,0]
voltageAxis_IFQ_mV = rawSamples_IFQ[:,0]
timeAxis_s = rawSamples_IFI[:,1]
totalSamples = timeAxis_s.size

# FFT computation
complexSignal_mV = np.add(np.asarray(voltageAxis_IFI_mV), 1j*np.asarray(voltageAxis_IFQ_mV))
FFT = np.fft.fft(complexSignal_mV, n = freqBins_FFT) # FFT of complex signal
FFT_mV = np.abs(1/(totalSamples)*FFT) # FFT magnitude
FFT_max = np.amax(FFT_mV)
FFT_dBV = 20*np.log10(FFT_mV/1000)
freqAxis = np.fft.fftfreq(freqBins_FFT) # freqBins+1
freqAxis_Hz = freqAxis * SAMPLING_FREQUENCY
peakFreq = freqAxis_Hz[FFT_mV.argmax()]
# Plot FFT
plt.plot(np.fft.fftshift(freqAxis_Hz), np.fft.fftshift(FFT_dBV))
plt.ylabel('Spectrum magnitude (dBV)')
plt.xlabel('Frequency (Hz)')
plt.grid(True)
plt.show()

# FFT computation: normalized and smoothed 
FFT_norm = FFT_mV / FFT_max
FFT_norm_dB = 20*np.log10(FFT_norm)
FFT_norm_smooth = np.convolve(FFT_norm, np.ones(smoothingBins), 'same') / smoothingBins
FFT_norm_dB_smooth = 20*np.log10(FFT_norm_smooth)

print('Detected Doppler frequency: {:.1f}'.format(peakFreq) + ' Hz')
print('Amplitude of this FFT peak: {:.1f}'.format(20*np.log10(FFT_max/1000)) + ' dBV')

# Plot FFT: normalized and smoothed
plt.plot(np.fft.fftshift(freqAxis_Hz), np.fft.fftshift(FFT_norm_dB))
plt.plot(np.fft.fftshift(freqAxis_Hz), np.fft.fftshift(FFT_norm_dB_smooth))
plt.ylabel('Spectrum magnitude (dB)')
plt.xlabel('Frequency (Hz)')
plt.grid(True)
plt.show()

# Spectrogram computation
f, t, Sxx = signal.spectrogram(complexSignal_mV, fs = SAMPLING_FREQUENCY, nperseg = 2**10, nfft = 2**11, scaling = 'spectrum', return_onesided=False)
plt.pcolormesh(t, f, Sxx, shading='gouraud')
plt.ylabel('Frequency [Hz]')
plt.xlabel('Time [sec]')
##plt.axis([0, 2, 0, 1000])
plt.show()