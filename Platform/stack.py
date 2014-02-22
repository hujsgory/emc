# -*- coding: cp1251
import Tkinter as tk
import tkFileDialog

class StackEditor():
    def __init__(self, parent):
        self.parent = parent
        self.initialize()

    def initialize(self):
        self.parent.title('Stack Editor')
        self.parent.iconbitmap(default='stack.ico')
        self.menubar = tk.Menu(self.parent, tearoff=0)
        self.menubar.add_command(label='Open', command=self.openfile)
        self.menubar.add_command(label='Save', command=self.savefile)
        self.menubar.add_command(label='Close', command=self.closefile)
        self.parent.config(menu=self.menubar)
        self.txt=tk.Text(self.parent)
        self.txt.pack()
        self.filename = ''

    def openfile(self):
        fn=tkFileDialog.askopenfile(mode='r')
        if fn:
            self.txt.delete('1.0',tk.END)
            self.txt.insert(tk.INSERT, fn.read())
            self.filename=fn.name
            fn.close()
        
    def savefile(self):
        if self.filename == '':
            fn=tkFileDialog.asksaveasfile()
            if fn:
                fn.write(self.txt.get('1.0',tk.END))
                self.filename=fn.name
                fn.close()
        else:
            open(self.filename,'w').write(self.txt.get('1.0',tk.END))

    def closefile(self):
        self.filename=''
        self.txt.delete('1.0',tk.END)


if __name__ == "__main__":
    root = tk.Tk()
    app = StackEditor(root)
    root.mainloop()