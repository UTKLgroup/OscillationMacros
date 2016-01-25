import scipy as np
import matplotlib.pyplot as plt
from scipy.stats import maxwell
from scipy.stats import chisquare

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)
E = np.arange(0.0001,5, 0.0001) #Energy Spectrum | 50,000 entries
k=1.27 #Constant associated with Units
L=1 #Distance of Observation in meters

Data=[]

#Set resoultion. 
#WARNING: 100x10 as resolution takes about 10 seconds to run. 
#Higher resolution would geometrically increase time of execution.
a_res=100
dm_res=10

#Set limits in parameter space
a_lim=1
m2_lim=1

da=a_lim/a_res
dm2=m2_lim/dm_res


#Compute chi-square for points in parameter space
ax=np.arange(da,a_lim+da,da)
m2x=np.arange(dm2,m2_lim+dm2,dm2)
a=da
while (a<=a_lim):
    m2=dm2
    temp=[] 
    while (m2<=m2_lim): 
        P=a*(np.sin(m2*k*(L/E))**2) 
        N=maxwell.pdf(E,0,1)*(1-P)#Energy Flux with Oscillations
        v=chisquare(N,f_exp=maxwell.pdf(E,0,1))          
        temp.append(v[0])
        m2+=dm2
        m2=round(m2,5)
    Data.append(temp)
    a+=da
    a=round(a,5)
        
arr=np.transpose(np.array(Data))

A,M2=np.meshgrid(ax, m2x)

#Just to Check
print(len(ax))
print(len(m2x))
print(len(M2))
print(len(A))
print(len(arr[0]))


#Plot Contour Fill map
lvl=np.arange(0,10+0.05,0.05)
CS=plt.contourf(A,M2,arr, levels=lvl)
plt.colorbar(CS, format="%.2f")

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

plt.title(r'$\chi^2$ |L=1m')#if needed, change Title Appropriately
plt.show() 
