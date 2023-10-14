import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt



#Asteroid Observation 1 - Frame 14
#Observation Time: 2023-06-29T06:44:45
#RA: 21:01:13.7616
#DEC: +43:58:38.430
#Asteroid Signal: 7801.3692
#RA/DEC Errors: 0.06142687703345614 0.021513360320091164
#Magnitude: 16.709519447026985

#Asteroid Observation 2 - Frame 4 120s
#Observation Time: 2023-07-04T03:57:41
#RA: 20:58:12.4388
#DEC: +43:45:59.771
#RA/DEC Errors: 0.13348680590019915 0.075093662918268
#Asteroid Signal: 790.5583
#Magnitude: 17.057656772778852

#Asteroid Observation 3 - Frame 4 120s (again)
#Observation Time: 2023-07-10T03:27:32
#RA: 20:53:17.3711
#DEC: +43:06:38.987
#RA/DEC Errors: 0.19016218458407455 0.19174134794752723
#Asteroid Signal: 1183.431
#Magnitude: 17.12977738477312

#Asteroid Observation 4 - Frame 4 180s 
#Observation Time: 2023-07-13T05:36:36
#RA: 20:50:18.9100
#DEC: +42:35:09.456
#RA/DEC Errors: 0.34965953233776287 0.21264682207218444
#Asteroid Signal: 2000.1410
#Magnitude: 16.063445075277222

data = np.loadtxt("D:\\CentroidMagnitudeDetermination\\10_REF_PHOTOMETRY4.txt", dtype = float)
print(data)
y = data[:,0] #Rmag
x = np.log10(data[:,1]) #Log Signal
print(y,x)
m, b = np.polyfit(x, y, 1)
plt.scatter(x,y)
plt.plot(x, m*x + b)
rmag_asteroid = m*np.log10(2000.1410) + b
plt.show()
print(rmag_asteroid)

