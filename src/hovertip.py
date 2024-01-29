from idlelib.tooltip import OnHoverTooltipBase
from tkinter import Label, LEFT, SOLID


class Hovertip(OnHoverTooltipBase):
    def __init__(self, anchor_widget, text, hover_delay=1000, **kwargs):
        super().__init__(anchor_widget, hover_delay=hover_delay)
        self.text = text
        self.kwargs = kwargs

    def showcontents(self):
        if "justify" not in self.kwargs:
            self.kwargs["justify"] = LEFT
        if "relief" not in self.kwargs:
            self.kwargs["relief"] = SOLID
        if "borderwidth" not in self.kwargs:
            self.kwargs["borderwidth"] = 1
        label = Label(
            self.tipwindow,
            text=self.text,
            **self.kwargs,
        )
        label.pack()

    def enable(self):  # Hack
        self._id1 = self.anchor_widget.bind("<Enter>", self._show_event)
        self._id2 = self.anchor_widget.bind("<Leave>", self._hide_event)
        self._id3 = self.anchor_widget.bind("<Button>", self._hide_event)

    def disable(self):  # Hack
        self.anchor_widget.unbind("<Enter>", self._id1)
        self.anchor_widget.unbind("<Leave>", self._id2)
        self.anchor_widget.unbind("<Button>", self._id3)
        self.hidetip()
