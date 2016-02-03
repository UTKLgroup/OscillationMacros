import scipy as np
import matplotlib.pyplot as plt
from scipy.stats import maxwell
from scipy.stats import chisquare

E = np.arange(0.0001,5, 0.0001) #Energy Spectrum | 50,000 entries
k=1.27 #Constant associated with Units
L=1 #Distance of Observation in meters

Data=[]

#Set resoultion. 
#WARNING: 50x50 as resolution takes about 10 seconds to run. 
#Higher resolution would geometrically increase time of execution.
a_res=50.0
m2_res=50.0

#Set limits in parameter space
a_lim=1
m2_lim=1

da=a_lim/a_res
dm2=m2_lim/m2_res

ax=np.arange(da,a_lim+da,da) 
m2x=np.arange(dm2,m2_lim+dm2,dm2) 

#Best Fit Values
#Sensitivity Plot generated when Best Fit values are (0,0)
a_bf=0.4
m2_bf=0.3

E_sansOsc=maxwell.pdf(E,0,1) #Energy spectrum without Oscillations
E_dist=E_sansOsc*(1-(a_bf*(np.sin(m2_bf*k*(L/E))**2)))#Best fit Energy spectrum generation

plt.plot([a_bf],[m2_bf], '*', color='w' )#Best Fit Point in Parametric Space


a=0+da
while (a<=a_lim):
    m2=0+dm2
    temp=[] 
    while (m2<=m2_lim):
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

#Make Legend
cs=plt.contour(A,M2,arr, levels=[2.30,4.61,9.21], colors=['r','b','k'])#Specific confidence levels correspond to specific chi-sq values
lines=[cs.collections[0],cs.collections[1],cs.collections[2]]
labels = ['68% C.L.','90% C.L.','99% C.L.']
plt.legend(lines, labels)

plt.xlabel("sin^2(2*mixing angle)")
plt.ylabel("Mass Square Difference (eV^2)")

'''
To change the scale to log, simply uncomment the below statements.
To change the scale to linear, simply comment out the below statements.
In logscale, the upper limit of an axis corresponds to the a_lim and m2_lim. 
The lower limit of an axis corresponds to the ratios between a_lim & a_res and m2_lim & m2_res.
'''
#plt.xscale("log")
#plt.yscale("log")

plt.title(r'$\chi^2$ | L=1m | Best Fit')#if needed, change Title Appropriately
plt.show() 
