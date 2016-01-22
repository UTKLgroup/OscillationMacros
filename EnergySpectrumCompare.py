#import math
import scipy as np
import matplotlib.pyplot as plt
from scipy.special import factorial
from scipy.stats import maxwell
from matplotlib.widgets import Slider
#from scipy.special import factorial

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)
E = np.arange(0.0001,5, 0.0001)
med=2
a = 0.5
m2 = 0.5
k=1.27
mu=0.5
L=1

def fact(x):
    return factorial(x)
    




#p_pois=(mu**(np.absolute(E-med))*np.exp(-mu))/fact((np.absolute(E-med)))

#Adjust Scale and Location of Energy Flux
maxw_dist=maxwell.pdf(E,0,1)#Energy Flux without Oscillations (Boltzmann Distribution)
N=maxwell.pdf(E,0,1)*(1-(a*(np.sin(m2*k*(L/E))**2)))#Energy Flux with Oscillations 

#Plot both curves
l, =plt.plot(E, maxw_dist, lw=2, color='red')
l1, =plt.plot(E, N, lw=0.5, color='blue')
plt.xlabel("E (MeV)")
plt.ylabel("Number of Entries")
plt.title("Energy Spectrum Comparison")


#Postion of Slider
axcolor = 'lightgoldenrodyellow'

axamp = plt.axes([0.28, 0.1, 0.65, 0.03], axisbg=axcolor)
axfreq = plt.axes([0.28, 0.05, 0.65, 0.03], axisbg=axcolor)
axlen = plt.axes([0.28, 0.002, 0.65, 0.03], axisbg=axcolor)

#Attributes of Slider

samp = Slider(axamp, 'Mass Square Difference (eV^2)', 0, 10, valinit=m2)
sfreq = Slider(axfreq, 'sin^2(2theta)', 0, 1, valinit=a)
slen = Slider(axlen, 'Length (m)', 0, 20, valinit=L)

#Update Graph
def update(val):
    m2 = samp.val
    a = sfreq.val
    L= slen.val
    l1.set_ydata(maxwell.pdf(E,0,1)*(1-(a*(np.sin(m2*k*(L/E))**2))) )
    fig.canvas.draw_idle()
sfreq.on_changed(update)
samp.on_changed(update)
slen.on_changed(update)



plt.show()
