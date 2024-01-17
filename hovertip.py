from idlelib.tooltip import OnHoverTooltipBase
from tkinter import Label, LEFT, SOLID


class Hovertip(OnHoverTooltipBase):
    def __init__(self, anchor_widget, text, hover_delay=1000):
        super().__init__(anchor_widget, hover_delay=hover_delay)
        self.text = text

    def showcontents(self):
        label = Label(
            self.tipwindow,
            text=self.text,
            justify=LEFT,
            relief=SOLID,
            borderwidth=1,
        )
        label.pack()
