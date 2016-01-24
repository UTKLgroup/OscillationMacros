
import scipy as np
import matplotlib.pyplot as plt
from scipy.stats import maxwell
from scipy.stats import chisquare

E = np.arange(0.0001,5, 0.0001) #Energy Spectrum
k=1.27
L=1



Data=[]

a_res=50
dm_res=50

a_lim=1
m2_lim=1

da=a_lim/a_res
dm2=m2_lim/dm_res

ax=np.arange(0,a_lim,da)
m2x=np.arange(0,m2_lim,dm2)

#Best Fit Values
#Sensitivity Plot generated when Best Fit values are (0,0)
a_bf=0.2
m2_bf=0.2

E_sansOsc=maxwell.pdf(E,0,1) #Energy spectrum without Oscillations
E_dist=E_sansOsc*(1-(a_bf*(np.sin(m2_bf*k*(L/E))**2)))

plt.plot([a_bf],[m2_bf], '*', color='w' )

a=0
while (a<a_lim):
    m2=0
    temp=[] 
    while (m2<m2_lim):
        P=a*(np.sin(m2*k*(L/E))**2) 
        N=E_sansOsc*(1-P)#Energy Flux with Oscillations
        v=chisquare(N,f_exp=E_dist)          
        temp.append(v[0])
        m2+=dm2
        m2=round(m2,5)
    Data.append(temp)
    a+=da
    a=round(a,5)
        
arr=np.transpose(np.array(Data))

A,M2=np.meshgrid(ax, m2x)

#Just to check
print(len(ax))
print(len(m2x))
print(len(M2))
print(len(A))
print(len(arr[0]))

cs=plt.contour(A,M2,arr, levels=[2.30,4.61,9.21], colors=['r','b','k'])
lines=[cs.collections[0],cs.collections[1],cs.collections[2]]
labels = ['68% C.L.','90% C.L.','99% C.L.']
plt.legend(lines, labels)

plt.xlabel("sin^2(2*mixing angle)")
plt.ylabel("Mass Square Difference (eV^2)")

#plt.xscale("log")
#plt.yscale("log")

plt.title(r'$\chi^2$ | L=1m | Best Fit')
plt.show() 
