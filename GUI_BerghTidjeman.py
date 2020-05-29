import tkinter as tk
from tkinter import ttk

#create instance
parent = tk.Tk()
# Add a title
parent.title("-Acoustic tube Bergh-Tidjeman model -")

###Canvas with drawn polygon

canvas_width = 1000
canvas_height = 500
canvas_main = tk.Canvas(parent, width=canvas_width, height=canvas_height)
canvas_main.pack()

# create input field

entry1 = tk.Entry (parent)
canvas_main.create_window(100, 40, window=entry1) #position of input field
entry1.insert(0, "a default value") # default value



#### frame with two buttons

frame = tk.Frame(parent)
frame.pack()

def write_text():
    print("Tkinter is easy to create GUI!")

text_disp= tk.Button(frame,
                   text="Hello",
                   command=write_text
                   )

text_disp.pack(side=tk.LEFT)

exit_button = tk.Button(frame,
                   text="Exit",
                   fg="green",
                   command=quit)
exit_button.pack(side=tk.RIGHT)

#### combo box

my_str_var = tk.StringVar()

my_combobox = ttk.Combobox(
    parent, textvariable = my_str_var,
    values=["PHP", "Java", "Python"])

my_combobox.pack()

parent.mainloop()