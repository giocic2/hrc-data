import easygui
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

# ANALYSIS SETTINGS
SAMPLING_FREQUENCY = 100e3 # According to "hrc-ps.py" script
FFT_RESOL = 1 # Hz

# FFT bins and resolution
freqBins_FFT = int(2**np.ceil(np.log2(abs(SAMPLING_FREQUENCY/2/FFT_RESOL))))
print('FFT resolution: ' + str(SAMPLING_FREQUENCY/freqBins_FFT) + ' Hz')
print('FFT bins: ' + str(freqBins_FFT))

filename = None
while filename == None:
    filename = easygui.fileopenbox(title = "Choose *.csv file to analyse...", default = "*.csv")
print(filename)

rawSamples = np.genfromtxt(filename, delimiter = ',')

voltageAxis_mV = rawSamples[:,0]
timeAxis_s = rawSamples[:,1]
totalSamples = timeAxis_s.size

plt.plot(timeAxis_s, voltageAxis_mV)
plt.ylabel('Voltage (mV)')
plt.xlabel('Time (s)')
plt.grid(True)
plt.show()

# FFT computation
FFT = np.fft.rfft(voltageAxis_mV, n = freqBins_FFT ) # FFT of real signal
FFT_mV = np.abs(2/(totalSamples)*FFT) # FFT magnitude
FFT_dBV = 20*np.log10(FFT_mV/1000)
freqAxis = np.fft.rfftfreq(freqBins_FFT ) # freqBins/2+1
freqAxis_Hz = freqAxis * SAMPLING_FREQUENCY

# Plot FFT
plt.plot(freqAxis_Hz, FFT_dBV)
plt.ylabel('Spectrum magnitude (dBV)')
plt.xlabel('Frequency (Hz)')
plt.grid(True)
plt.show()

# Spectrogram computation
f, t, Sxx = signal.spectrogram(voltageAxis_mV, fs = SAMPLING_FREQUENCY, nperseg = 2**10, nfft = 2**11, scaling = 'spectrum')
plt.pcolormesh(t, f, Sxx, shading='gouraud')
plt.ylabel('Frequency [Hz]')
plt.xlabel('Time [sec]')
##plt.axis([0, 2, 0, 1000])
plt.show()