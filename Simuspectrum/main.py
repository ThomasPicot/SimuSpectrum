import matplotlib.pyplot as plt

from Classes.classspectrum import SpectrumRubidiumD2Line


sp = SpectrumRubidiumD2Line()
spectrum = sp.transmission(frac=27.8, temp=273+20, long=0.06)
detuning = sp.detuning()
absorption_coefficient = sp.alpha(C_f = 8/31, frac=28, temp=273+20, transition_frequency=196.029e12, deg=8)
print(absorption_coefficient)

plt.plot(detuning, spectrum)
#plt.show()