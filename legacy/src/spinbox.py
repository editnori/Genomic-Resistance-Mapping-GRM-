import tkinter as tk
from util import Key, force_insertable_value
from decimal import Decimal


class Spinbox(tk.Spinbox):
    def __init__(self, master, default_value: float = None, **kwargs):
        super().__init__(master, **kwargs)

        if "validate" not in kwargs or kwargs["validate"] == tk.NONE:
            super().configure(validate="key")
        if "validatecommand" not in kwargs:
            super().configure(
                validatecommand=(
                    self.register(self.validate_spinbox),
                    "%P",
                    "%W",
                )
            )

        self.set_default_value(default_value)

        self.bind(
            "<FocusOut>",
            lambda e: (
                force_insertable_value(e.widget.get_default_value(), e.widget)
                if not e.widget.get()
                else None
            ),
        )
        self.bind(
            "<KeyRelease>",
            lambda e: (
                master.focus_set() if e.keycode in [Key.ENTER, Key.ESCAPE] else None
            ),
        )
        self.bind(
            "<KeyPress>",
            lambda e: master.focus_set() if e.keycode == Key.SPACE else None,
        )

        def on_scroll(event):
            if event.widget.focus_get() != event.widget:
                return
            if event.delta > 0:
                self.invoke("buttonup")
            else:
                self.invoke("buttondown")

        self.bind("<MouseWheel>", on_scroll)

        self.bindtags(self.bindtags()[:-1])

    def set_default_value(self, default_value: float, set_value: bool = True):
        if default_value:
            default_value = float(default_value)
            self.default_value = (
                int(default_value)
                if default_value.is_integer()
                else float(default_value)
            )
        else:
            self.default_value = (
                int(super().cget("from"))
                if super().cget("from").is_integer()
                else float(super().cget("from"))
            )

        if set_value:
            force_insertable_value(self.default_value, self)

    def get_default_value(self) -> float:
        return self.default_value

    def validate_spinbox(
        self,
        new_value: float,
        widget_name: str,
    ):
        if not new_value:
            return True

        try:
            float(new_value)

            widget = self.nametowidget(widget_name)
            new_value = Decimal(new_value)
            increment = Decimal(widget.cget("increment"))
            from_value = Decimal(widget.cget("from"))
            to_value = Decimal(widget.cget("to"))
        except ValueError:
            return False

        is_in_increment = (new_value - from_value) % increment == 0

        if (from_value <= new_value <= to_value) and is_in_increment:
            return True

        return False
