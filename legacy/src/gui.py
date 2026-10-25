import ctk
import tkinter as tk
from util import select_directory, select_file


class PathSelector(tk.Frame):
    def __init__(
        self,
        parent,
        label="Path",
        is_dir=False,
        label_args=None,
        entry_args=None,
        button_args=None,
        **kwargs,
    ):
        tk.Frame.__init__(self, parent, **kwargs)

        self._is_focused = False
        self.fg_active = "#1f212a"
        self.border_hover = "#63666f"
        self.border_active = "#284cb8"

        if label_args is None:
            label_args = {
                "text_color": "#ffffff",
            }
        if entry_args is None:
            entry_args = {
                "text_color": "#ffffff",
                "fg_color": "#3c404b",
                "border_color": "#3c404b",
            }
        if button_args is None:
            button_args = {
                "fg_color": "#3c404b",
                "hover_color": "#63666f",
                "text_color": "#ffffff",
                "text": "Browse",
            }

        self.label = ctk.CTkLabel(self, text=label, **label_args)
        self.entry_value = tk.StringVar()
        self.entry = ctk.CTkEntry(
            self,
            **entry_args,
            textvariable=self.entry_value,
        )

        def focus_in(event):
            self._is_focused = True
            self.entry.configure(
                fg_color=self.fg_active,
                border_color=self.border_active,
            )

        self.entry.bind(
            "<FocusIn>",
            focus_in,
        )

        def hover_in(event):
            if not self._is_focused:
                self.entry.configure(border_color=self.border_hover)

        self.entry.bind(
            "<Enter>",
            hover_in,
        )

        def hover_out(event):
            if not self._is_focused:
                self.entry.configure(border_color=entry_args["border_color"])

        self.entry.bind(
            "<Leave>",
            hover_out,
        )

        def focus_out(event):
            self._is_focused = False
            self.entry.configure(
                fg_color=entry_args["fg_color"],
                border_color=entry_args["border_color"],
            )

        self.entry.bind(
            "<FocusOut>",
            focus_out,
        )

        self.button = ctk.CTkButton(self, **button_args)

        self.label.grid(row=0, column=0, sticky=tk.W, padx=(0, 5))
        self.entry.grid(row=0, column=1, sticky=tk.EW)
        self.button.grid(row=0, column=2, sticky=tk.E, padx=(5, 0))

        self.grid_columnconfigure(1, weight=3, uniform="uniform")
        self.grid_columnconfigure(2, weight=1, uniform="uniform")

        self.is_dir = is_dir

        self.button.configure(command=self.browse)

    def get(self):
        return self.entry.get()

    def set(self, value):
        self.entry.delete(0, tk.END)
        self.entry.insert(0, value)

    def browse(self):
        self.focus_set()
        if self.is_dir:
            return self.set(select_directory(default=self.get()))
        else:
            return self.set(select_file(default=self.get()))

    def enable(self):
        self.button.configure(state="normal")
        self.entry.configure(state="normal")

    def disable(self):
        self.button.configure(state="disabled")
        self.entry.configure(state="disabled")

    def focus(self):
        self.entry.focus()

    def on_update(self, callback):
        def callback_wrapper(*args):
            callback()

        self.entry_value.trace_add("write", callback_wrapper)
