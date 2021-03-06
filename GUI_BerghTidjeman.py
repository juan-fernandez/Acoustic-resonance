import numpy as np
import math
import cmath
import scipy.special as sci
import os

def calculations():

    ratio = []
    prod_ratio = []
    alpha = []
    n = []
    phi = []

    # Dimensions input

    ts_dimensions_matrix = [[dimensions_matrix[i][j] for i in range(7)] for j in range(3)] # transpose matrix
    rad = [float(ts_dimensions_matrix[0][i].get()) * 0.001 for i in range(7) if ts_dimensions_matrix[0][i].get() != '']
    for r in rad:
        if r == 0:
            tk.messagebox.showerror(title=None, message="Radius cannot be 0")
            return
    L = [float(ts_dimensions_matrix[1][i].get()) * 0.001 for i in range(7) if ts_dimensions_matrix[1][i].get() != '']
    for l in L:
        if l == 0:
            tk.messagebox.showerror(title=None, message="Length cannot be 0")
            return
    vol = [float(ts_dimensions_matrix[2][i].get()) * 1E-9 for i in range(7) if ts_dimensions_matrix[2][i].get() != '']
    if len(rad) == len(L) == len(vol):
        pass
    else:
        tk.messagebox.showwarning(title=None, message="segment dimensions incomplete")

    # Frequency

    freq_start = int(f_start.get())
    freq_end = int(f_end.get())
    step = int(f_step.get())
    freq = range(freq_start, freq_end, step)

    # Air parameters input

    T_celsius = float(temp.get())
    T_kelvin = 273.15 + T_celsius

    Combo_value = Combo.get()
    if Combo_value == 'Air':
        visc = 18.27 * 0.000001 * ((291.15 + 120) / (T_kelvin + 120)) * ((T_kelvin / 291.15)) ** (1.5)
    elif Combo_value == 'Other gas':
        visc = float(vis.get()) * 1E-6


    Ps = float(Pstatic.get())
    dens = Ps / (287.05 * T_kelvin)
    k = float(adiab.get())
    m = float(mol.get())
    Ro = 8.314
    r = Ro / m
    speed = (k * r * T_kelvin) ** 0.5
    kt = float(therm.get())
    cp = float(Sh.get())
    Pr = cp * visc / kt

    var_message = False
    for j in freq: #loop for frequencies

        omega = 2 * np.pi * j #angular frequency


        for i in range(len(rad)-1, -1, -1): #loop for tube sections


            alpha.insert(0, ((complex(0, 1)) ** 1.5) * rad[i] * np.sqrt(dens * omega / visc))

            n.insert(0, (1 + ((k - 1) / k) * ((sci.jv(2, (alpha[0] * np.sqrt(Pr))))
                                              / (sci.jv(0, (alpha[0] * np.sqrt(Pr)))))) ** (-1))

            phi.insert(0, ((omega / speed) * np.sqrt((sci.jv(0, alpha[0]))/(sci.jv(2, alpha[0]))) * np.sqrt(k / n[0])))

            if cmath.isnan(phi[0]) == True and var_message == False:
                tk.messagebox.showerror(title=None, message="For the radius " + str(rad[i] * 1000)
                                + " this model can only simulate up to " + str(j - step) + " Hz of frequency")
                var_message = True

            if i == len(rad)-1: #for the first section of tube

                ratio.insert(0,(((np.cosh(phi[0]*L[i])) +
                                (n[0]*phi[0]*vol[i]/(k*np.pi*(rad[i]**2))) *
                                 np.sinh(phi[0]*L[i]))**(-1)))



            # Up to here we??ve calculated the ratio for an isolated section. We apply until here for the first section

            else:

                ratio.insert(0,((np.cosh(phi[0]*L[i]) + (n[0]*phi[0]*vol[i]/(k*np.pi*(rad[i]**2)))*np.sinh(phi[0]*L[i]) +
                (((rad[i+1])/(rad[i]))**2)*((phi[1])/(phi[0])) *
                ((sci.jv(0,alpha[0]))/(sci.jv(0,alpha[1]))) *
                ((sci.jv(2, alpha[1])) / (sci.jv(2, alpha[0]))) *
                ((np.sinh(phi[0]*L[i]))/(np.sinh(phi[1]*L[i+1]))) *
                ((np.cosh(phi[1]*L[i+1])) - ratio[0])) **(-1)))


        prod_ratio.append(np.prod(ratio)) # multiply ratios for the different tubes and obtains total ratio for a frequency j
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
import tkinter.font as font
from PIL import Image, ImageTk
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.figure import Figure

root = tk.Tk()
root.title("-Acoustic tube Bergh-Tidjeman model -")
root.geometry("1500x900")

wrapper1 = tk.LabelFrame(root, text="Tube geometry")
wrapper2 = tk.LabelFrame(root, text="Parameters")
wrapper3 = ttk.Notebook(root)  # Plots
wrapper4 = tk.LabelFrame(root, text="Tube scheme")


dim_input = tk.Label(wrapper1, anchor="s")
dimensions_matrix = [["a" for i in range(3)] for j in range(7)]
#input box Table
dimensions_array = []
for row in range(8):
    for column in range(4):
        if row == 0:
            if column == 0:
                heading = tk.Label(dim_input, text="Tube segment", fg="black", padx=3, pady=3)
            elif column == 1:
                heading = tk.Label(dim_input, text="Radius [mm]", fg="black", padx=3, pady=3)
            elif column == 2:
                heading = tk.Label(dim_input, text="Length [mm]", fg="black", padx=3, pady=3)
            elif column == 3:
                heading = tk.Label(dim_input, text="Volume [mm^3]", fg="black", padx=3, pady=3)
            #label.config(font=('Arial', 14))
            heading.grid(row=row, column=column)#, sticky="nsew", padx=1, pady=1)
            dim_input.grid_columnconfigure(column, weight=1, minsize=100)
        elif column == 0 and row!=0:
            segment = tk.Label(dim_input, text="# " + str(row), fg="black", padx=3, pady=3, font=("Calibri", 10))
            segment.grid(row=row, column=column)
        else:
            dimensions_matrix[row-1][column-1] = tk.StringVar()
            cell = ttk.Entry(dim_input, textvariable=dimensions_matrix[row-1][column-1], justify='center', width=12)
            if row == 1 and column == 1:
                cell.insert(0, 0.4)
            elif row == 1 and column == 2:
                cell.insert(0, 200)
            elif row == 1 and column == 3:
                cell.insert(0, 5)
            #dimensions_matrix[row-1][column] = dimensions_matrix[row-1][column].get()
            cell.grid(row=row, column=column, sticky="nsew", padx=1, pady=1)
            #dim_input.grid_columnconfigure(column, weight=1)

dim_input.pack(fill="both")

# Input Air parameters

param_input = tk.Label(wrapper2, anchor="s")


# function activating and deactivating viscosity enry
def callbackFunc(event):
    Combo_value = Combo.get()
    if Combo_value == 'Air':
        vis_heading.config(state='disabled')
        vis.config(state='disabled')
    elif Combo_value == 'Other gas':
        vis_heading.config(state='normal')
        vis.config(state='normal')

# selector between Air and other gas
vlist = ['Air', 'Other gas']
Combo = ttk.Combobox(param_input, values=vlist, state='readonly')
Combo.current(0) # shows air as default
Combo.grid(columnspan=2, row=0, column=0, padx=5, pady=5)
Combo_value = Combo.get()
Combo.bind("<<ComboboxSelected>>", callbackFunc)

# viscosity entry
vis_heading = tk.Label(param_input, text="Gas viscosity at temp [m^2/s (10^-6)]", fg="black", padx=8, pady=3)
vis_heading.config(state='disabled')
vis_heading.grid(columnspan=2, row=0, column=2)
visc = tk.StringVar()
vis = ttk.Entry(param_input, justify='center', textvariable=visc)
vis.insert(0, 10)
vis.config(state='disabled')
vis.grid(columnspan=2, row=0, column=4, sticky="nsew", padx=8, pady=5)


sep = ttk.Separator(param_input, orient='horizontal').grid(row=1, column=0, columnspan=6, sticky="ew", pady=8)

temp_heading = tk.Label(param_input, text="Temperature [??C]", fg="black", padx=8, pady=3)
temp_heading.grid(columnspan=3, row=2, column=0)
T_celsius = tk.StringVar()
temp = ttk.Entry(param_input, justify='center', textvariable=T_celsius)
temp.insert(0, 20.0) # default
temp.grid(columnspan=3, row=3, column=0, sticky="nsew", padx=30, pady=1)

Pstatic_heading = tk.Label(param_input, text="Static pressure [Pa]", fg="black", padx=8, pady=3)
Pstatic_heading.grid(columnspan=3, row=4, column=0)
Ps = tk.StringVar()
Pstatic = ttk.Entry(param_input, justify='center', textvariable=Ps)
Pstatic.insert(0, 100000.0)
Pstatic.grid(columnspan=3, row=5, column=0, sticky="nsew", padx=30, pady=1)

adiab_heading = tk.Label(param_input, text="Adiabatic constant", fg="black", padx=8, pady=3)
adiab_heading.grid(columnspan=3, row=6, column=0)
k = tk.StringVar()
adiab = ttk.Entry(param_input, justify='center', textvariable=k)
adiab.insert(0, 1.4)
adiab.grid(columnspan=3, row=7, column=0, sticky="nsew", padx=30, pady=1)

mol_heading = tk.Label(param_input, text="Molecular weight [kg/mol]", fg="black", padx=8, pady=3)
mol_heading.grid(columnspan=3, row=2, column=3)
m = tk.StringVar()
mol = ttk.Entry(param_input, justify='center', textvariable=m)
mol.insert(0, 0.02896)
mol.grid(columnspan=3, row=3, column=3, sticky="nsew", padx=30, pady=1)

therm_heading = tk.Label(param_input, text="Thermal conductivity [W/(m K)]", fg="black", padx=8, pady=3)
therm_heading.grid(columnspan=3, row=4, column=3)
kt = tk.StringVar()
therm = ttk.Entry(param_input, justify='center', textvariable=kt)
therm.insert(0, 0.026)
therm.grid(columnspan=3, row=5, column=3, sticky="nsew", padx=30, pady=1)

Sh_heading = tk.Label(param_input, text="Specific heat (Cp) [J/(K kg)]", fg="black", padx=8, pady=3)
Sh_heading.grid(columnspan=3, row=6, column=3)
cp = tk.StringVar()
Sh = ttk.Entry(param_input, justify='center', textvariable=cp)
Sh.insert(0, 1006.1)
Sh.grid(columnspan=3, row=7, column=3, sticky="nsew", padx=30, pady=1)


# Input frequency range
sep = ttk.Separator(param_input, orient='horizontal').grid(row=8, column=0, columnspan=6, sticky="ew", pady=8)

f_start_heading = tk.Label(param_input, text="Frequency Start [Hz]", fg="black", padx=2, pady=1)
f_start_heading.grid(columnspan=2, row=9, column=0)
freq_start = tk.StringVar
f_start = ttk.Entry(param_input, justify='center', textvariable=freq_start)
f_start.insert(0, 1)
f_start.grid(columnspan=2, row=10, column=0, sticky="nsew", padx=2, pady=1)

f_end_heading = tk.Label(param_input, text="Frequency End [Hz]", fg="black", padx=2, pady=1)
f_end_heading.grid(columnspan=2, row=9, column=2)
freq_end = tk.StringVar
f_end = ttk.Entry(param_input, justify='center', textvariable=freq_end)
f_end.insert(0, 20000)
f_end.grid(columnspan=2, row=10, column=2, sticky="nsew", padx=2, pady=1)


f_step_heading = tk.Label(param_input, text="Frequency Step [Hz]", fg="black", padx=2, pady=1)
f_step_heading.grid(columnspan=2, row=9, column=4)
step = tk.StringVar
f_step = ttk.Entry(param_input, justify='center', textvariable=step)
f_step.insert(0, 5)
f_step.grid(columnspan=2, row=10, column=4, sticky="nsew", padx=2, pady=0)

#freq_input.pack(fill="y")
param_input.pack(fill="both")

freq_start = int(f_start.get())
freq_end = int(f_end.get())
step = int(f_step.get())
freq = range(freq_start, freq_end, step)


#Amplitude ratio tab
tab1 = tk.Frame(wrapper3)
tab1.pack(fill="both")
label = tk.Label(tab1)
label.pack()
wrapper3.add(tab1, text="Amplitude ratio")

################## amplitude ratio plot

calculations()

fig1 = Figure(figsize=(10, 2), dpi=100)

ax = fig1.add_subplot(1,1,1)
ax.set_xlim(0, freq_end)
ax.set_ylim(0, 1.1 * max(amp_ratio) + 0.1)
ax.set_xlabel("Frequency [Hz]")
ax.set_ylabel("Amplitude ratio")
ax.plot(freq, amp_ratio)
canvas1 = FigureCanvasTkAgg(fig1, master=tab1)
ax.grid(True)
canvas1.draw()

toolbar = NavigationToolbar2Tk(canvas1, tab1)  # Plot??s legend
toolbar.update()
canvas1.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

def amp_ratio_plot():

    freq_start = int(f_start.get())
    freq_end = int(f_end.get())
    step = int(f_step.get())
    freq = range(freq_start, freq_end, step)

    ax = canvas1.figure.axes[0]
    ax.set_xlim(0, freq_end)
    ax.set_ylim(0, 1.1 * max(amp_ratio))
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("Amplitude ratio")
    ax.plot(freq, amp_ratio)
    ax.grid(True)
    canvas1.draw()


########################################################################################
#Phase shift tab

tab2 = tk.Frame(wrapper3)
tab2.pack(fill="both")
label = tk.Label(tab2)
label.pack()
wrapper3.add(tab2, text="Phase shift")

fig2 = Figure(figsize=(10, 2), dpi=100)
ax1 = fig2.add_subplot(1,1,1)
ax1.set_xlim(0, freq_end)
ax1.set_ylim(min(phase_shift)-1, max(phase_shift)+1)
ax1.set_xlabel("Frequency [Hz]")
ax1.set_ylabel("Phase shift [degree]")
ax1.plot(freq, phase_shift)
canvas2 = FigureCanvasTkAgg(fig2, master=tab2)
ax1.grid(True)
canvas2.draw()


toolbar = NavigationToolbar2Tk(canvas2, tab2)  # Plot??s legend
toolbar.update()
canvas2.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

def phase_shift_plot():

    freq_start = int(f_start.get())
    freq_end = int(f_end.get())
    step = int(f_step.get())
    freq = range(freq_start, freq_end, step)

    #fig2.clear()
    ax1 = canvas2.figure.axes[0]
    ax1.set_xlim(0, freq_end)
    ax1.set_ylim(min(phase_shift), max(phase_shift))

    ax1.plot(freq, phase_shift)
    ax1.grid(True)
    canvas2.draw()

def clear_plots():

    ax.clear()
    ax1.clear()
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("Amplitude ratio")
    ax1.set_xlabel("Frequency [Hz]")
    ax1.set_ylabel("Phase shift [degree]")
    ax.grid(True)
    ax1.grid(True)
    canvas1.draw()
    canvas2.draw()


# insert image in GUI
central_input = tk.Label(wrapper4, anchor="s")

image_name = str(os.getcwd()) + "\\tube_image.png"#"C:\\Users\\Daniel\\Desktop\\tube_image.png"
PIL_image = Image.open(image_name)
PIL_image_small = PIL_image.resize((400,110), Image.ANTIALIAS)

img = ImageTk.PhotoImage(PIL_image_small)
in_frame = tk.Label(central_input, image = img)
in_frame.grid(columnspan=1, row=0, column=0)
#in_frame.pack(side=tk.BOTTOM, fill="x", pady=3)
central_input.pack(fill="y")

import csv

name_value = tk.StringVar()
csv_name = tk.Entry(central_input, justify='center', textvariable=name_value, width=32, bd=3)
csv_name.insert(0, "Transfer_function.csv")

csv_name.grid(columnspan=1, sticky='W', row=1, column=0, ipady=3, pady=3, padx=2)
#csv_name.pack(expand=True, fill="x", side=tk.LEFT)

def create_csv_file():

    # field names
    fields = ['Frequency [Hz]', 'Amplitude ratio', 'Phase shift [deg]']
    columns = [freq, amp_ratio, phase_shift]
    rows = [[columns[i][j] for i in range(3)] for j in range(len(freq))]# transpose colums to get rows

    filename = name_value.get()

    # writing to csv file
    with open(filename, 'w', newline='') as csvfile:
        # creating a csv writer object
        csvwriter = csv.writer(csvfile, delimiter=';')
        # writing the fields
        csvwriter.writerow(fields)
        # writing the data rows
        csvwriter.writerows(rows)

refresh_plots = tk.Button(central_input, height=1, width=62, text="Clear plots", command=lambda:\
    [clear_plots()]).grid(columnspan=1, row=3, column=0, pady=3)

write_csv = tk.Button(central_input, height=1, width=30, text="Save last data as csv", font=font.Font(family='Calibri', size=11, weight='bold'),command=lambda:[create_csv_file()])\
    .grid(columnspan=1, sticky='E', row=1, column=0, pady=3, padx=2)

calculate_data = tk.Button(central_input, height=3, width=43, text="Calculate!", font=font.Font(family='Calibri', size=14, weight='bold'), command=lambda:\
    [(calculations(), amp_ratio_plot(), phase_shift_plot())]).grid(columnspan=1, row=2, column=0, pady=2)


wrapper3.pack(side=tk.BOTTOM, fill="both", expand="yes", padx=20, pady=10)
wrapper1.pack(side=tk.LEFT, fill="x", expand="yes", padx=20, pady=10)
wrapper2.pack(side=tk.RIGHT, fill="both", expand="yes", padx=20, pady=10)
wrapper4.pack(side=tk.LEFT, fill="both", expand="yes", padx=20,pady=10)


root.mainloop()