from tkinter import ttk

class Combobox(ttk.Combobox):
    def __init__(self, master, **kwargs):
        super().__init__(master, **kwargs)
        self.bindtags(self.bindtags()[:-1])