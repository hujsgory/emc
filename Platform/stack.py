# -*- coding: cp1251
import Tkinter as tk
import tkFileDialog

class StackEditor():
    def __init__(self, parent):
        self.parent = parent
        self.initialize()

    def initialize(self):
        self.parent.title('Stack Editor')
        try:
            self.parent.iconbitmap(default='stack.ico')
        except tk.TclError:
            pass
        self.menubar = tk.Menu(self.parent, tearoff=0)
        self.menubar.add_command(label='Open', command=self.openfile)
        self.menubar.add_command(label='Save', command=self.savefile)
        self.menubar.add_command(label='Close', command=self.closefile)
        self.parent.config(menu=self.menubar)
        self.txt = tk.Text(self.parent,font='12',height=7)
        self.txt.grid(column=0, row=0, sticky=tk.E+tk.W)
        self.txt.bind('<Return>', self.redraw)
        self.txt.bind('<Up>', self.redraw)
        self.txt.bind('<Down>', self.redraw)
        self.canvas = tk.Canvas(self.parent)
        self.canvas.grid(column=0, row=1, sticky=tk.E+tk.W)
        self.filename = ''

    def openfile(self):
        fn=tkFileDialog.askopenfile(mode='r')
        if fn:
            self.txt.delete('1.0',tk.END)
            self.txt.insert(tk.INSERT, fn.read())
            self.filename = fn.name
            fn.close()
        
    def savefile(self):
        if self.filename == '':
            fn=tkFileDialog.asksaveasfile()
            if fn:
                fn.write(self.txt.get('1.0',tk.END))
                self.filename = fn.name
                fn.close()
        else:
            open(self.filename,'w').write(self.txt.get('1.0',tk.END))

    def closefile(self):
        self.filename = ''
        self.txt.delete('1.0',tk.END)

    def redraw(self, event):
        pass


if __name__ == "__main__":
    root = tk.Tk()
    app = StackEditor(root)
    root.mainloop()
