import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', family='arial')
plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-small')
plt.rcParams["figure.figsize"] = (7,5)


freq_fan = 1/2 * np.asarray([14.1, 16.3, 18.6, 20.8, 23.0, 25.0, 27.3, 30.3, 33.3])
freq_IF = np.asarray([196.8, 230.4, 262.5, 289.9, 320.4, 352.5, 379.9, 418.1, 462.3])
# IF peak frequency vs fan rotating frequency
plt.plot(freq_fan, freq_IF, marker='^', markeredgecolor='g', markerfacecolor='none', markersize=12, linestyle='dashed', color='g')
# Plot settings
plt.yscale("linear")
plt.xscale("linear")
plt.ylabel("IF peak frequency (Hz)", fontsize=18)
plt.xlabel("rotation frequency of the fan (Hz)", fontsize=18)
plt.yticks(fontsize=16)
plt.xticks(ticks=np.linspace(start=6, stop=18, num=int((18-6)/2+1), endpoint=True), fontsize=16)
plt.grid(b=True, which='major', linewidth=1)
plt.grid(b=True, which='minor', linewidth=0.3)
# plt.legend([r'- 9 dB', '- 6 dB', '- 4 dB', '- 0 dB'], loc='upper center', fontsize = 18, bbox_to_anchor=(0.5, 1.25), ncol=2, fancybox=True, shadow=True)
plt.tight_layout()
plt.show()