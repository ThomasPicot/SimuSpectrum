import matplotlib.pyplot as plt

from Classes.classspectrum import SpectrumRubidiumD2Line
from Classes.n2_fit_2 import Spectrum_D2, Spectrum_D1

D2 = SpectrumRubidiumD2Line()

"""Temperature = 273
for i in range(273+120, 273+150):"""
# for i in range(0, 30):
plt.plot(D2.detuning()*1e-9, D2.transmission(temp=393, long=0.1, frac=28))
plt.title('Transmission for the D2 line with Rubidium mixture ')
plt.ylabel('Transmission $\mathcal{T}(\Delta)$')
plt.xlabel('Detunning $\Delta = \omega - \omega_i $ [GHz]')
plt.grid()
plt.show()

