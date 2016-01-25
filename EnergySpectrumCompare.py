import scipy as np
import matplotlib.pyplot as plt
from scipy.stats import maxwell
from matplotlib.widgets import Slider


fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)
E = np.arange(0.0001,5, 0.0001)
#Initial Values on Slider
a = 0.5 
m2 = 0.5
L=1

k=1.27 #Depends on units used


#Adjust Scale and Location of Energy Spectrum
E_sansOsc=maxwell.pdf(E,0,1) #Energy Spectrum without Oscillations (Maxwell-Boltzmann Distribution)
N=E_sansOsc*(1-(a*(np.sin(m2*k*(L/E))**2))) #Energy Spectrum with Oscillations 

#Plot both curves
l, =plt.plot(E, E_sansOsc, lw=2, color='red')
l1, =plt.plot(E, N, lw=0.5, color='blue')
plt.xlabel("E (MeV)")
plt.ylabel("Number of Entries")
plt.title("Energy Spectrum Comparison")


#Postion of Slider
axcolor = 'lightgoldenrodyellow'

axm2 = plt.axes([0.28, 0.1, 0.65, 0.03], axisbg=axcolor)
axa = plt.axes([0.28, 0.05, 0.65, 0.03], axisbg=axcolor)
axlen = plt.axes([0.28, 0.002, 0.65, 0.03], axisbg=axcolor)

#Attributes of Slider

sm2 = Slider(axm2, 'Mass Square Difference (eV^2)', 0, 10, valinit=m2)
sa = Slider(axa, 'sin^2(2*Mixing Angle)', 0, 1, valinit=a)
slen = Slider(axlen, 'Length (m)', 0, 20, valinit=L)

#Update Graph
def update(val):
    m2 = sm2.val
    a = sa.val
    L= slen.val
    l1.set_ydata(E_sansOsc*(1-(a*(np.sin(m2*k*(L/E))**2))))
    fig.canvas.draw_idle()
sa.on_changed(update)
sm2.on_changed(update)
slen.on_changed(update)

plt.show()
