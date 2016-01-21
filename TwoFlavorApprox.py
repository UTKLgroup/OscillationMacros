import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)
t = np.arange(0,10, 0.0001)
a = 0.5
m2 = 10
k=1.27
#s = a0*np.sin(2*np.pi*f0*t)
s=a*(np.sin(m2*k*(t))**2)
l, =plt.semilogx(t, s, lw=2, color='red')
plt.axis([0, 1, 0, 1])
plt.xlabel("L/E (m/MeV)")
plt.ylabel("Probability of Disappearance")

#Postion of Slider
axcolor = 'lightgoldenrodyellow'
axfreq = plt.axes([0.28, 0.05, 0.65, 0.03], axisbg=axcolor)
axamp = plt.axes([0.28, 0.1, 0.65, 0.03], axisbg=axcolor)

sfreq = Slider(axfreq, 'sin^2(2theta)', 0, 1, valinit=a)
samp = Slider(axamp, 'Mass Square Difference (eV^2)', 0, 100, valinit=m2)

def update(val):
    m2 = samp.val
    a = sfreq.val
    l.set_ydata(a*(np.sin(m2*k*(t))**2))
    fig.canvas.draw_idle()
sfreq.on_changed(update)
samp.on_changed(update)

plt.show()