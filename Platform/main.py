#coding: cp1251
import Tkinter as tk
import ttk, tkFileDialog

def askopenfile():
    tkFileDialog.askopenfile(mode='r')
    
root = tk.Tk()

menubar = tk.Menu(root)

filemenu = tk.Menu(menubar, tearoff=0)
filemenu.add_command(label='Open', command=askopenfile)
filemenu.add_command(label='Save')
filemenu.add_command(label='Close')
filemenu.add_separator()
filemenu.add_command(label='Quit', command=root.quit)

editmenu = tk.Menu(menubar, tearoff=0)
editmenu.add_command(label='Cut')
editmenu.add_command(label='Copy')
editmenu.add_command(label='Paste')

menubar.add_cascade(label='File', menu=filemenu)
menubar.add_cascade(label='Edit', menu=editmenu)
root.config(menu=menubar)


root.mainloop()

