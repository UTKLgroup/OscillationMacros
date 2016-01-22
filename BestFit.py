
import scipy as np
import matplotlib.pyplot as plt
#from scipy.special import factorial
from scipy.stats import maxwell
#from matplotlib.widgets import Slider
from scipy.stats import chisquare
#from scipy.special import factorial

E = np.arange(0.0001,5, 0.0001)
med=2
k=1.27
mu=0.5
L=1



Data=[]

a_res=50
dm_res=50

a_lim=1
m2_lim=1

da=a_lim/a_res
dm2=m2_lim/dm_res


#p_pois=(mu**(np.absolute(E-med))*np.exp(-mu))/fact((np.absolute(E-med)))

#Adjust Scale and Location of Energy Flux
'''
def chi(mix,massdiff):
    P=mix*(np.sin(massdiff*k*(L/E))**2)
    N=maxwell.pdf(E,0,1)*(1-P)#Energy Flux with Oscillations 
    v=chisquare(N,f_exp=maxwell.pdf(E,0,1))
    return v[0]
'''
ax=np.arange(0,a_lim,da)
m2x=np.arange(0,m2_lim,dm2)

#Best Fit Values
a_bf=0.2
m2_bf=0.2

E_dist=maxwell.pdf(E,0,1)*(1-(a_bf*(np.sin(m2_bf*k*(L/E))**2)))

plt.plot([a_bf],[m2_bf], '*', color='w' )

a=0
while (a<a_lim):
    m2=0
    temp=[] 
    while (m2<m2_lim):
        P=a*(np.sin(m2*k*(L/E))**2) 
        N=maxwell.pdf(E,0,1)*(1-P)#Energy Flux with Oscillations
        v=chisquare(N,f_exp=E_dist)          
        temp.append(v[0])
        m2+=dm2
        m2=round(m2,5)
    Data.append(temp)
    a+=da
    a=round(a,5)
        
arr=np.transpose(np.array(Data))

A,M2=np.meshgrid(ax, m2x)

print(len(ax))
print(len(m2x))
print(len(M2))
print(len(A))
print(len(arr[0]))

#CS=plt.contourf(A,M2,arr, lvl)
'''
CS=plt.contour(A,M2,arr, levels=[2.3], color='b')
CS=plt.contour(A,M2,arr, levels=[4.61], color='r')

lvl=np.arange(0,10,0.05)

CS=plt.contourf(A,M2,arr, levels=lvl)
plt.colorbar(CS, format="%.2f")
'''

cs=plt.contour(A,M2,arr, levels=[2.30,4.61,9.21], colors=['r','b','k'])
lines=[cs.collections[0],cs.collections[1],cs.collections[2]]
labels = ['68% C.L.','90% C.L.','99% C.L.']
plt.legend(lines, labels)




#arr=np.array(Data)
#arrT=np.transpose(arr)
#print (arrT)

#plt.figure()
#plt.contour(arrT[0],arrT[1], arrT[2])



'''
for a in range(0,1+da, da):
    for m2 in range(0, 3+dm2, dm2):
        P=a*(np.sin(m2*k*(L/E))**2) 
        N=maxwell.pdf(E,0,1)*(1-P)#Energy Flux with Oscillations
        v=chisquare(N,f_exp=maxwell.pdf(E,0,1))
        temp=[]        
        temp.append(a)
        temp.append(m2)
        temp.append(v)
        Data.append(temp)
'''     

 
  
  




'''
CS = plt.contour(A, M2, v)



a, m2 = np.mgrid[slice(0, 1 + da, da), slice(0, 3 + dm2, dm2)]

maxw_dist=maxwell.pdf(E,0,1)#Energy Flux without Oscillations (Boltzmann Distribution)
N=maxwell.pdf(E,0,1)*(1-(a*(np.sin(m2*k*(L/E))**2)))#Energy Flux with Oscillations 

chi=chisquare(N,f_exp=maxw_dist)


plt.pcolor(a,m2,chi)
'''

'''
#Plot both curves
l, =plt.plot(E, maxw_dist, lw=2, color='red')
l1, =plt.plot(E, N, lw=0.5, color='blue')
plt.xlabel("E (MeV)")
plt.ylabel("Number of Entries")



#Postion of Slider
axcolor = 'lightgoldenrodyellow'

axamp = plt.axes([0.28, 0.1, 0.65, 0.03], axisbg=axcolor)
axfreq = plt.axes([0.28, 0.05, 0.65, 0.03], axisbg=axcolor)
axlen = plt.axes([0.28, 0.002, 0.65, 0.03], axisbg=axcolor)

#Attributes of Slider

samp = Slider(axamp, 'Mass Square Difference (eV^2)', 0, 10, valinit=m2_bf)
sfreq = Slider(axfreq, 'sin^2(2theta)', 0, 1, valinit=a_bf)


#Update Graph
def update(val):
    m2 = samp.val
    a = sfreq.val
    l1.set_ydata(maxwell.pdf(E,0,1)*(1-(a*(np.sin(m2*k*(L/E))**2))) )
    fig.canvas.draw_idle()
sfreq.on_changed(update)
samp.on_changed(update)
slen.on_changed(update)

'''
plt.xlabel("sin^2(2*mixing angle)")
plt.ylabel("Mass Square Difference (eV^2)")

#plt.xscale("log")
#plt.yscale("log")

plt.title(r'$\chi^2$ | L=1m | Best Fit')
plt.show() 
