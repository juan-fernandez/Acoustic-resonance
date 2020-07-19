import numpy as np
import math
import scipy.special as sci
import matplotlib.pyplot as plt


#Temperature
T_celsius = 20 #°C
T_kelvin = T_celsius + 273.15

#Air dynamic viscosity
visc = 18.27 * 0.000001 * ((291.15 + 120) / (T_kelvin + 120)) * ((T_kelvin / 291.15))**(1.5)
#Adiabatic constant
k = 1.4

#Ideal Gas constant
Ro = 8.314

#Molecular weight (kg/mol)
m = 0.02896 # to be defined as input

#Constant of the gas
r = Ro / m

#Speed of sound
speed = (k * r * T_kelvin)**0.5


#Static pressure (Pa)
Ps = 100000

#density
dens = Ps / (287.05 * T_kelvin)

#Prandtl number

cp = 1006.1 #specific heat
kt = 0.026 #Thermal conductivity [W/m K]
#Pr = cp * visc / kt
Pr = 0.71486 #Value given in matlab



freq_start = 1
freq_end = 2000 # in Hz. to be defined as input
step = 5 # in Hz. to be defined as input
Count = 0 # iteration progress display

num = 2 #from input number of sections is taken

freq = range(freq_start, freq_end, step)


ratio = []
prod_ratio = []
alpha = []
n = []
phi = []


rad = [0.004, 0.002, 0.003]
L = [0.06, 0.01, 0.05]
vol = [0, 0.000002, 0.000009]

for j in freq: #loop for frequencies

    omega = 2 * np.pi * j #angular frequency

    for i in range(len(rad)-1, -1, -1): #loop for tube sections

        alpha.insert(0, ((complex(0, 1)) ** 1.5) * rad[i] * np.sqrt(dens * omega / visc))
        n.insert(0, (1 + ((k - 1) / k) * ((sci.jv(2, (alpha[0] * np.sqrt(Pr))))
                                          / (sci.jv(0, (alpha[0] * np.sqrt(Pr)))))) ** (-1))
        phi.insert(0, ((omega / speed) * np.sqrt((sci.jv(0, alpha[0]))/(sci.jv(2, alpha[0]))) * np.sqrt(k / n[0])))

        if i == len(rad)-1: #for the first section of tube

            ratio.insert(0,(((np.cosh(phi[0]*L[i])) +
                            (n[0]*phi[0]*vol[i]/(k*np.pi*(rad[i]**2))) *
                             np.sinh(phi[0]*L[i]))**(-1)))


        # Up to here we´ve calculated the ratio for an isolated section. We apply until here for the first section

        else:

            ratio.insert(0,((np.cosh(phi[0]*L[i]) + (n[0]*phi[0]*vol[i]/(k*np.pi*(rad[i]**2)))*np.sinh(phi[0]*L[i]) +
            (((rad[i+1])/(rad[i]))**2)*((phi[1])/(phi[0])) * #(L[i]/L[i+1]) * ##check rad and L!!!
            ((sci.jv(0,alpha[0]))/(sci.jv(0,alpha[1]))) *
            ((sci.jv(2, alpha[1])) / (sci.jv(2, alpha[0]))) *
            ((np.sinh(phi[0]*L[i]))/(np.sinh(phi[1]*L[i+1]))) *
            ((np.cosh(phi[1]*L[i+1])) - ratio[0])) **(-1)))

    prod_ratio.append(math.prod(ratio)) #multiply ratios for the different tubes and obtains total ratio for a frequency j
    alpha = []
    n = []
    phi = []
    ratio = [] # empty ratio array after all sections have been computes for a frequency

amp_ratio = [np.absolute(x) for x in prod_ratio] # Take absolute value of the complex ratio
phase_shift = [np.angle(x, deg=False) for x in prod_ratio] # Take phase
phase_shift = np.unwrap(phase_shift) # unwrap phases
phase_shift = np.rad2deg(phase_shift) # from rad to degrees

#plot_1 = plt.figure(1, figsize=(10, 5))
#plt.plot(freq, phase_shift)
#plot_2 = plt.figure(2, figsize=(10,5))
#plt.plot(freq, amp_ratio)
#plt.show()