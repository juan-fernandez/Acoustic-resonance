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

amp_ratio = [2*x for x in freq]
phase_shift = freq #define arrays with same length as freq to start with


# take dimensions from GUI
dimensions_matrix = [[0.004, 0.06, 0], [0.002, 0.01, 0.000002], [0.003, 0.05, 0.000009]]
def calculations():

    ratio = []
    prod_ratio = []
    alpha = []
    n = []
    phi = []

   #print("a")
    #print(len(dimensions_array))

    ts_dimensions_matrix = [[dimensions_matrix[i][j] for i in range(2)] for j in range(3)] # transpose matrix
    rad = [float(ts_dimensions_matrix[0][i].get()) for i in range(2)]
    L = [float(ts_dimensions_matrix[1][i].get()) for i in range(2)]
    vol = [float(ts_dimensions_matrix[2][i].get()) for i in range(2)]

    #print(rad, L, vol)

    #rad = [0.004, 0.002, 0.003]
    #L = [0.06, 0.01, 0.05]
    #vol = [0, 0.000002, 0.000009]

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

        prod_ratio.append(math.prod(ratio)) # multiply ratios for the different tubes and obtains total ratio for a frequency j
        alpha = []
        n = []
        phi = []
        ratio = [] # empty ratio array after all sections have been computes for a frequency

    global amp_ratio
    global phase_shift
    amp_ratio = [np.absolute(x) for x in prod_ratio] # Take absolute value of the complex ratio
    phase_shift = [np.angle(x, deg=False) for x in prod_ratio] # Take phase
    phase_shift = np.unwrap(phase_shift) # unwrap phases
    phase_shift = np.rad2deg(phase_shift) # from rad to degrees

#########################################################################################################
##########################     GUI
#########################################################################################################

import tkinter as tk
from tkinter import ttk
from PIL import Image, ImageTk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.figure import Figure

root = tk.Tk()
root.title("-Acoustic tube Bergh-Tidjeman model -")
root.geometry("1300x800")

wrapper1 = tk.LabelFrame(root, text="Tube geometry")
wrapper2 = tk.LabelFrame(root, text="Parameters")
wrapper3 = ttk.Notebook(root)  # Plots
wrapper4 = tk.LabelFrame(root, text="Tube scheme")


dim_input = tk.Label(wrapper1, anchor="s")
dimensions_matrix = [["a" for i in range(3)] for j in range(2)]
#input box Table
dimensions_array=[]
for row in range(3):
    for column in range(3):
        if row == 0:
            if column == 0:
                heading = tk.Label(dim_input, text="Radius [mm]", fg="black", padx=3, pady=3)
            elif column == 1:
                heading = tk.Label(dim_input, text="Length [mm]", fg="black", padx=3, pady=3)
            elif column == 2:
                heading = tk.Label(dim_input, text="Volume [mm^3]", fg="black", padx=3, pady=3)
            #label.config(font=('Arial', 14))
            heading.grid(row=row, column=column)#, sticky="nsew", padx=1, pady=1)
            dim_input.grid_columnconfigure(column, weight=1, minsize=100)
        else:
            dimensions_matrix[row-1][column] = tk.StringVar()
            cell = ttk.Entry(dim_input, textvariable=dimensions_matrix[row-1][column])
            #dimensions_matrix[row-1][column] = dimensions_matrix[row-1][column].get()
            cell.grid(row=row, column=column, sticky="nsew", padx=1, pady=1)
            #dim_input.grid_columnconfigure(column, weight=1)

dim_input.pack(fill="both")

#Turn input dimension array into matrix
#input_matrix = [dimensions_array[i:i+3] for i in range(0,len(dimensions_array), 3)]
#[[print(x[i]) for i in range(3)] for x in input_matrix] # prints input matrix



#Amplitude ratio tab
tab1 = tk.Frame(wrapper3)
tab1.pack(fill="both")
label = tk.Label(tab1, text="jojojojo")
label.pack()
wrapper3.add(tab1, text="Amplitude ratio")

################## amplitude ratio plot

fig1 = Figure(figsize=(10, 2), dpi=100)
#amp_ratio = 0# amp_ratio
ax = fig1.add_subplot(1,1,1)#.plot(freq, amp_ratio)
ax.set_xlim(0, freq_end)
ax.set_ylim(0, max(amp_ratio))
ax.plot(freq, amp_ratio)
canvas = FigureCanvasTkAgg(fig1, master=tab1)
canvas.draw()
#canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

toolbar = NavigationToolbar2Tk(canvas, tab1)  # Plot´s legend
toolbar.update()
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
def amp_ratio_plot():

    #fig1.clear()
    ax = canvas.figure.axes[0]
    ax.set_xlim(0, freq_end)
    ax.set_ylim(0, 1.1 * max(amp_ratio))
    print(amp_ratio)
    ax.plot(freq, amp_ratio)
    #fig1.add_subplot(1,1,1)#.plot(freq, amp_ratio)
    #canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    canvas.draw()

########################################################################################
#Phase shift tab
tab2 = tk.Frame(wrapper3)
tab2.pack(fill="both")
label = tk.Label(tab2, text="jijiji")
label.pack()
wrapper3.add(tab2, text="Phase shift")

fig2 = Figure(figsize=(10, 2), dpi=100)
ax1 = fig2.add_subplot(1,1,1)#.plot(freq, amp_ratio)
ax1.set_xlim(0, freq_end)
ax1.set_ylim(min(phase_shift), max(phase_shift))
ax1.plot(freq, phase_shift)
canvas1 = FigureCanvasTkAgg(fig2, master=tab2)
canvas1.draw()



######################### phase shift plot
#fig2 = Figure(figsize=(10, 2), dpi=100)
#fig2.add_subplot(111).plot(freq, phase_shift)
#canvas = FigureCanvasTkAgg(fig2, master=tab2)
#canvas.draw()  ###################### issue with plotting amp in phase_shift
#canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

toolbar = NavigationToolbar2Tk(canvas1, tab2)  # Plot´s legend
toolbar.update()
canvas1.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

def phase_shift_plot():

    #fig2.clear()
    ax = canvas1.figure.axes[0]
    ax.set_xlim(0, freq_end)
    ax.set_ylim(min(phase_shift), max(phase_shift))
    print(phase_shift)
    ax.plot(freq, phase_shift)
    #fig1.add_subplot(1,1,1)#.plot(freq, amp_ratio)
    #canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    canvas1.draw()


# Entry parameters
input = tk.Label(wrapper2, anchor="s")
for row in range(5):
    for column in range(2):

        cell = ttk.Entry(input)
        #input.pack()
        cell.grid(row=row, column=column, sticky="nsew", padx=1, pady=1)
        cell.grid_columnconfigure(column, weight=1)
input.pack(fill="both")


#input_cell = tk.StringVar()
#cell = ttk.Entry(dim_input, textvariable=input_cell)
#Data_matrix.append(cell)  # append value into data array
#cell.grid(row=row, column=column, sticky="nsew", padx=1, pady=1)
#dim_input.grid_columnconfigure(column, weight=1)

def print_input():
    [[print(x[i].get()) for i in range(3)] for x in dimensions_matrix]


refresh = tk.Button(root, text="Calculate!", command=lambda:\
    [(print_input(), calculations(), amp_ratio_plot(), phase_shift_plot())]).pack(side=tk.TOP)

wrapper3.pack(side=tk.BOTTOM, fill="both", expand="yes", padx=20, pady=10)
wrapper1.pack(side=tk.LEFT, fill="x", expand="yes", padx=20, pady=10)
wrapper2.pack(side=tk.RIGHT, fill="both", expand="yes", padx=20, pady=10)
wrapper4.pack(side=tk.LEFT, fill="both", expand="yes", padx=20,pady=10)


# insert image in GUI
image_name = "C:\\Users\\Daniel\\Desktop\\tube_image.png"
PIL_image = Image.open(image_name)
PIL_image_small = PIL_image.resize((400,150), Image.ANTIALIAS)

img = ImageTk.PhotoImage(PIL_image_small)
in_frame = tk.Label(wrapper4, image = img)
in_frame.pack(side=tk.TOP, fill="both")


root.mainloop()