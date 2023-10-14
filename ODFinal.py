import odlib as od
import math as m
import numpy as np
import vpython as vp
import scipy as sp
file = open("D:\\OD_code\\PLSSSWORK.txt", "r")
fileList = file.readlines()
file.close()
k = 0.0172020989484
cAU = 173.144643267

#Julian Date Conversions
year = int(fileList[0][0:4])
month = int(fileList[0][5:7])
day = int(fileList[0][8:10]) + int(fileList[0][11:13]) / 24 + int(fileList[0][14:16]) / (24*60) + int(fileList[0][17:19]) / (24*60*60)
t1 = 367*year - (7*(year + (month + 9)//12)) // 4 + (275*month) // 9 + day + 1721013.5

year = int(fileList[6][0:4])
month = int(fileList[6][5:7])
day = int(fileList[6][8:10]) + int(fileList[6][11:13]) / 24 + int(fileList[6][14:16]) / (24*60) + int(fileList[6][17:19]) / (24*60*60)
t2 = 367*year - (7*(year + (month + 9)//12)) // 4 + (275*month) // 9 + day + 1721013.5

year = int(fileList[12][0:4])
month = int(fileList[12][5:7])
day = int(fileList[12][8:10]) + int(fileList[12][11:13]) / 24 + int(fileList[12][14:16]) / (24*60) + int(fileList[12][17:19]) / (24*60*60)
t3 = 367*year - (7*(year + (month + 9)//12)) // 4 + (275*month) // 9 + day + 1721013.5

#Original Julian Dates for Time Adjustments
t10 = t1
t20 = t2
t30 = t3

#RA/DEC extraction
raList1 = fileList[1].split(":")
decList1 = fileList[2].split(":")

raList2 = fileList[7].split(":")
decList2 = fileList[8].split(":")

raList3 = fileList[13].split(":")
decList3 = fileList[14].split(":")

ra1 = (float(raList1[0]) + float(raList1[1])/60 + float(raList1[2])/3600)*15
dec1 = m.copysign((abs(float(decList1[0]) )+ float(decList1[1])/60 + float(decList1[2])/3600), float(decList1[0]))

ra2 = (float(raList2[0]) + float(raList2[1])/60 + float(raList2[2])/3600)*15
dec2 = m.copysign((abs(float(decList2[0]) )+ float(decList2[1])/60 + float(decList2[2])/3600), float(decList2[0]))

ra3 = (float(raList3[0]) + float(raList3[1])/60 + float(raList3[2])/3600)*15
dec3 = m.copysign((abs(float(decList3[0]) )+ float(decList3[1])/60 + float(decList3[2])/3600), float(decList3[0]))

#Rhohat (phat) creation 
phatVec1 = vp.vector( (m.cos(m.radians(ra1)))*m.cos(m.radians(dec1)), m.sin(m.radians(ra1))*m.cos(m.radians(dec1)), m.sin(m.radians(dec1)) )
phatVec2 = vp.vector( (m.cos(m.radians(ra2)))*m.cos(m.radians(dec2)), m.sin(m.radians(ra2))*m.cos(m.radians(dec2)), m.sin(m.radians(dec2)) )
phatVec3 = vp.vector( (m.cos(m.radians(ra3)))*m.cos(m.radians(dec3)), m.sin(m.radians(ra3))*m.cos(m.radians(dec3)), m.sin(m.radians(dec3)) )

#Tau creation
tau = k*(t3 - t1)
tau1 = k*(t1 - t2)
tau3 = k*(t3 - t2)

#R Creation - Sun Vectors
R1Array = list(map(float, fileList[3:6]))
R2Array = list(map(float, fileList[9:12]))
R3Array = list(map(float, fileList[15:18]))

#D Constants
D11 = np.dot(np.cross(R1Array, np.array([phatVec2.x, phatVec2.y, phatVec2.z])), np.array([phatVec3.x, phatVec3.y, phatVec3.z]))
D12 = np.dot(np.cross(R2Array, np.array([phatVec2.x, phatVec2.y, phatVec2.z])), np.array([phatVec3.x, phatVec3.y, phatVec3.z]))
D13 = np.dot(np.cross(R3Array, np.array([phatVec2.x, phatVec2.y, phatVec2.z])), np.array([phatVec3.x, phatVec3.y, phatVec3.z]))

D21 = np.dot(np.cross(np.array([phatVec1.x, phatVec1.y, phatVec1.z]), R1Array), np.array([phatVec3.x, phatVec3.y, phatVec3.z]))
D22 = np.dot(np.cross(np.array([phatVec1.x, phatVec1.y, phatVec1.z]), R2Array), np.array([phatVec3.x, phatVec3.y, phatVec3.z]))
D23 = np.dot(np.cross(np.array([phatVec1.x, phatVec1.y, phatVec1.z]), R3Array), np.array([phatVec3.x, phatVec3.y, phatVec3.z]))

D31 = np.dot(np.array([phatVec1.x, phatVec1.y, phatVec1.z]), np.cross(np.array([phatVec2.x, phatVec2.y, phatVec2.z]), R1Array))
D32 = np.dot(np.array([phatVec1.x, phatVec1.y, phatVec1.z]), np.cross(np.array([phatVec2.x, phatVec2.y, phatVec2.z]), R2Array))
D33 = np.dot(np.array([phatVec1.x, phatVec1.y, phatVec1.z]), np.cross(np.array([phatVec2.x, phatVec2.y, phatVec2.z]), R3Array))



D0 = np.dot(np.array([phatVec1.x, phatVec1.y, phatVec1.z]), np.cross(np.array([phatVec2.x, phatVec2.y, phatVec2.z]), np.array([phatVec3.x, phatVec3.y, phatVec3.z])))
print("asdfads",D0)
#SEL Root testing to find r2mag
r2mag, RHO2VAL = od.lagrange(tau, tau1, tau3, R2Array, np.array([phatVec2.x, phatVec2.y, phatVec2.z]), D0, D21, D22, D23)

#Initial f and g values
F1 = 1 - (1 / (2*(r2mag**3)))*tau1**2
F3 = 1 - (1 / (2*(r2mag**3)))*tau3**2
G1 = tau1 - (1 / (6*(r2mag**3)))*tau1**3
G3 = tau3 - (1 / (6*(r2mag**3)))*tau3**3
#Initial c calculations
c1 = G3 / (F1*G3 - G1*F3)
c3 = -G1 / (F1*G3 - G1*F3)
#Initial rho calculations
p1 = (c1*D11 - D12 + c3*D13) / (c1*D0)
p2 = (c1*D21 - D22 + c3*D23) / (-1*D0)
p3 = (c1*D31 - D32 + c3*D33) / (c3*D0)
#Initial position vector calculations
r1 = p1*np.array([phatVec1.x, phatVec1.y, phatVec1.z]) - np.array([R1Array[0], R1Array[1], R1Array[2]])
r2 = p2*np.array([phatVec2.x, phatVec2.y, phatVec2.z]) - np.array([R2Array[0], R2Array[1], R2Array[2]])
r3 = p3*np.array([phatVec3.x, phatVec3.y, phatVec3.z]) - np.array([R3Array[0], R3Array[1], R3Array[2]])

#Initial velocity vector calculations
d1 = -F3 / (F1*G3 - F3*G1)
d3 = F1 / (F1*G3 - F3*G1)
r2dot = d1*r1 + d3*r3

#Newton DeltaE function
def newton(r2, r2dot, tau, n, a):
    r2mag = (r2[0]**2 +r2[1]**2 +r2[2]**2)**(1/2)
    function = lambda x : x - (1 - (r2mag/a))*m.sin(x) + ((np.dot(r2, r2dot)) / (n*a**2))*(1 - m.cos(x)) - tau*n
    return sp.optimize.newton(function, tau*n, tol = 10e-12)

p2minus1 =  0
i  = 0
while(abs(p2minus1 - p2) >= 4e-12):
    #Iteration counter
    i += 1    

    #Updating previous p2 value
    p2minus1 = p2

    #Light time travel correction and updating
    t1 = t10 - p1 / cAU
    t2 = t20 - p2 / cAU
    t3 = t30 - p3 / cAU

    tau1 = k*(t1 - t2)
    tau3 = k*(t3 - t2)
    
    #Magnitude definitions for local use in while loop
    pos = (r2[0]**2 +r2[1]**2 +r2[2]**2)**(1/2)
    vel = (r2dot[0]**2 +r2dot[1]**2 +r2dot[2]**2)**(1/2)
    
    a = 1 / (2 / pos - vel**2)
    n = (1 / a**3)**(1/2)

    #Delta E calculations
    E1 = newton(r2, r2dot, tau1, n, a)
    E3 = newton(r2, r2dot, tau3, n, a)

    #f and g calculations
    f1 = 1 - (a / pos)*(1 - m.cos(E1))
    f3 = 1 - (a / pos)*(1 - m.cos(E3))
    g1 = tau1 + (1 / n)*(m.sin(E1) - E1)
    g3 = tau3 + (1 / n)*(m.sin(E3) - E3)

    #c calculations
    c1 = g3 / (f1*g3 - g1*f3)
    c3 = -g1 / (f1*g3 - g1*f3)

    #rho magnitudes calculations
    p1 = (c1*D11 - D12 + c3*D13) / (c1*D0)
    p2 = (c1*D21 - D22 + c3*D23) / (-1*D0)
    p3 = (c1*D31 - D32 + c3*D33) / (c3*D0)
    
    #position vector calculation
    r1 = p1*np.array([phatVec1.x, phatVec1.y, phatVec1.z]) - np.array([R1Array[0], R1Array[1], R1Array[2]])
    r2 = p2*np.array([phatVec2.x, phatVec2.y, phatVec2.z]) - np.array([R2Array[0], R2Array[1], R2Array[2]])
    r3 = p3*np.array([phatVec3.x, phatVec3.y, phatVec3.z]) - np.array([R3Array[0], R3Array[1], R3Array[2]]) 

    #velocity vector calculation
    d1 = -f3 / (f1*g3 - f3*g1)
    d3 = f1 / (f1*g3 - f3*g1)
    r2dot = d1*r1 + d3*r3   
    
    print(str(i) + ": change in rho2 = " + str(p2minus1 - p2) + " AU; light-travel time = " + str(round((p2 / cAU)*24*3600, 1)) + " sec")


print("In", i , "iterations, r2 and r2dot converged to")
print("r2 = " ,r2)
print("r2dot = ", r2dot)

print("in cartesian equatorial or ")

#Ecliptic Conversion
eps = m.radians(23.4374)
epsMat = np.array([[1,0,0], 
              [0, m.cos(eps), m.sin(eps)], 
              [0, -m.sin(eps), m.cos(eps)]])
r2 = epsMat@r2
r2dot = epsMat@r2dot

print("r2 = " ,r2)
print("r2dot = ", r2dot)
print("in cartesian ecliptic coordinates")

#-------------------------Orbital Elements-------------------------

print("ORBITAL ELEMENTS")

#Copied functions from odlib, changed couple things to fit format
def angularMomentum():
    h = np.cross(r2, r2dot)
    return h

def a():
    r = (r2[0]**2 +r2[1]**2 +r2[2]**2)**(1/2)
    v = (r2dot[0]**2 +r2dot[1]**2 +r2dot[2]**2)**(1/2)
    a = 1 / (2 / r - v**2)
    return a

def e():
    h = np.linalg.norm(angularMomentum())
    A = a()
    e = (1 - h**2 / A)**(1/2)
    return e

def incline():
    h = angularMomentum()
    i = m.atan(((h[0]**2 + h[1]**2)**(1/2)) / h[2])
    i = i *180 / m.pi
    return i 

def longitude_ascending_node():
    h = angularMomentum()
    hmag = np.linalg.norm(h)
    i = incline() / 180 * m.pi
    sin_omega = (h[0] / (hmag*m.sin(i)))
    cos_omega = (-h[1] / (hmag*m.sin(i)))
    omega = od.return_angle(cos_omega, sin_omega) * m.pi / 180
    omega = omega * 180 / m.pi
    return omega

def true_anomaly():
    A = a()
    E = e()
    h = np.linalg.norm(angularMomentum())
    r = (r2[0]**2 +r2[1]**2 +r2[2]**2)**(1/2)
    cosv = ((A*(1-E**2) / r - 1) / E)
    sinv = (A*(1 - E**2) / (E*h))*((np.dot(r2, r2dot)) / np.linalg.norm(r2))
    v = od.return_angle(cosv, sinv)
    return v

def argument_of_periapsis():
   i = incline() / 180 * m.pi
   r = (r2[0]**2 +r2[1]**2 +r2[2]**2)**(1/2)
   om = longitude_ascending_node() * m.pi / 180
   x = r2[0]
   y = r2[1]
   z = r2[2]
   sinU = (z / (r*m.sin(i)))
   cosU = ((x*m.cos(om)) + y*m.sin(om)) / r
   U = od.return_angle(cosU, sinU) * m.pi / 180
   v = true_anomaly() / 180 * m.pi
   w = U - v
   w %= 2*m.pi
   w = w*180 / m.pi
   return w

def mean_anomaly():
    E = e()
    A = a()
    r = (r2[0]**2 +r2[1]**2 +r2[2]**2)**(1/2)
    E0 = m.acos((1 / E)*(1 - r / A))
    M = E0 - E*m.sin(E0)
    M = M * 180 / m.pi
    if (true_anomaly() < 180):
        return 360 - M
    else:
        return M
    
def julian_date():
    M = mean_anomaly() / 180 * m.pi 
    A = a()
    k = 0.01720209895
    JD = t2 - (1 / ((k**2 / A**3)**(1/2)))*M
    return JD


def n():
    A = a()
    return (k**2 / A**3)**(1/2)

def mean_anomaly_2():
    t0 = t2
    t = 2459419.7916667
    M = m.radians(mean_anomaly())
    N = n()
    return M + N*(t - t0)


#print("Angular Momentum:", angularMomentum())
print("Semi-major Axis: a = ", a(), "au")
print("Eccentricity: e = ", e())
print("Incline: i =", incline(), "degrees")
print("Longitude of Ascending Node: omega = ", longitude_ascending_node(), "degrees")
#print("True Anomaly = ", true_anomaly())
print("Argument of Periapsis: Omega = ", argument_of_periapsis(), "degrees")
print("Mean Anomaly: M = ", mean_anomaly(), "degrees")
print("Mean Anomaly 2: M =",  m.degrees(mean_anomaly_2()), "degrees")

print("E:", m.degrees(od.newtonian_raphson(e(), m.radians(mean_anomaly()))), "degrees at central obs.")
print("n:", m.degrees(n()), "deg/day")
print("JD:", julian_date(), " - JD of last perihelion passage")
print("P:", a()**(3/2), "yrs")

def transform():
    A = a()
    E = e()
    M = mean_anomaly_2()
    om = m.radians(longitude_ascending_node())
    i = m.radians(incline())
    w = m.radians(argument_of_periapsis())
    ta = m.radians(true_anomaly())
    exp = (1-E)*(1 - ((np.linalg.norm(r2)) / A))
    if (ta >= 0 and ta <= m.pi):
        E0 = m.acos(exp)
    else:
        E0 = 2*m.pi - m.acos(exp)
    print(A, E, M, om, i, ta)
    #E0 = od.newtonian_raphson(E, M)
    
    omRot = np.array([[m.cos(om), -1*m.sin(om), 0],
                [m.sin(om), m.cos(om),0],
                [0,0,1]])

    iRot = np.array([[1,0,0],
                [0, m.cos(i),-1*m.sin(i)],
                [0, m.sin(i), m.cos(i)]])
    wRot = np.array([[m.cos(w),-1*m.sin(w),0],
                [m.sin(w), m.cos(w),0],
                [0, 0, 1]])
    xyzVec = np.array([[A*m.cos(E0) - A*E],[A*(1-E**2)**(1/2)*m.sin(E0)],[0]])

    ecliptic_matrix = np.matmul(np.matmul(np.matmul(omRot,iRot), wRot),xyzVec)
    epsilon = m.radians(23.5)
    e_matrix= np.array([[1,0,0],
                        [0, m.cos(epsilon), -m.sin(epsilon)],
                        [0,m.sin(epsilon), m.cos(epsilon)]])
    equatorial_matrix = np.matmul(e_matrix, ecliptic_matrix)

    equatorial_vector = vp.vector(equatorial_matrix[0][0], equatorial_matrix[1][0], equatorial_matrix[2][0])
    p = vp.vector(-2.051294646157535E-01, 9.136123910745069E-01, 3.960409296560626E-01) + equatorial_vector
    p = p / p.mag
    declination = m.asin(p.z)
    racos = p.x / m.cos(declination)
    rasin = p.y / m.cos(declination)
    right_ascension = od.return_angle(racos, rasin)
    return ("Right Ascension:", right_ascension, "Declination:", m.degrees(declination))

print(transform())