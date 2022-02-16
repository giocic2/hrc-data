import easygui
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.fft import fftshift

# ANALYSIS SETTINGS
SAMPLING_FREQUENCY = 100e3 # According to "hrc-ps.py" script
FFT_RESOL = 1 # Hz
SMOOTHING_WINDOW = 10 # Hz
FREQUENCY_MIN = -50_000 # Hz (lower limit for doppler centroid estimation)
BANDWIDTH_THRESHOLD = 6 # dB
ITERATIONS = 3

# FFT bins and resolution
freqBins_FFT = int(2**np.ceil(np.log2(abs(SAMPLING_FREQUENCY/2/FFT_RESOL))))
print('FFT resolution: ' + str(SAMPLING_FREQUENCY / freqBins_FFT) + ' Hz')
print('FFT bins: ' + str(freqBins_FFT))
smoothingBins = int(round(SMOOTHING_WINDOW / (SAMPLING_FREQUENCY / freqBins_FFT)))
print('Size of smoothing window (moving average): ' + str(smoothingBins) + ' bins')

done = False
index = 1
while not done:
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
    timeAxis_s = rawSamples_IFI[:,1]
    totalSamples = timeAxis_s.size

    # FFT computation
    complexSignal_mV = np.add(np.asarray(voltageAxis_IFI_mV), 1j*np.asarray(voltageAxis_IFQ_mV))
    complexSignal_mV_win = complexSignal_mV * np.hamming(totalSamples)
    FFT = np.fft.fftshift(np.fft.fft(complexSignal_mV_win, n = freqBins_FFT)) # FFT of complex signal
    FFT_mV = np.abs(1/(totalSamples)*FFT) # FFT magnitude
    FFT_max = np.amax(FFT_mV)
    FFT_dBV = 20*np.log10(FFT_mV/1000)
    freqAxis = np.fft.fftshift(np.fft.fftfreq(freqBins_FFT)) # freqBins+1
    freqAxis_Hz = freqAxis * SAMPLING_FREQUENCY

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
            if startBand <= peakFreq and stopBand >= peakFreq:
                centroidDetected = True
                break

    print('Detected Doppler frequency: {:.1f}'.format(peakFreq) + ' Hz')
    print('Amplitude of this FFT peak (norm.smooth.): {:.1f}'.format(FFT_norm_dB_smooth_max) + ' dB')
    print('Bandwidth threshold: {:.1f}'.format(BANDWIDTH_THRESHOLD) + ' dB')
    print('Bandwidth: {:.1f}'.format(stopBand - startBand) + ' Hz')
    print('Bandwidth starts at {:.1f}'.format(startBand) + ' Hz')
    print('Bandwidth stops at {:.1f}'.format(stopBand) + ' Hz')

    # Plot FFT: normalized and smoothed
    # plt.plot(freqAxis_Hz, FFT_norm_dB)
    plt.plot(freqAxis_Hz, FFT_norm_dB_smooth)
    plt.ylabel('Spectrum magnitude (dB)')
    plt.xlabel('Frequency (Hz)')
    plt.grid(True)

    if index == ITERATIONS:
        done = True
    index += 1
plt.legend(['1','2','3'])
plt.show()