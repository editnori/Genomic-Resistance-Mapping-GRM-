import customtkinter as ctk
import tkinter as tk
from util import Tag


class ControlFrame(ctk.CTkFrame):
    def __init__(self, master, **kwargs):
        if "fg_color" not in kwargs:
            kwargs["fg_color"] = "transparent"

        if "font" in kwargs:
            self.font = kwargs.pop("font")
        else:
            self.font = lambda size: ("Century Gothic", size)

        master.grid_columnconfigure(tuple(range(10)), weight=1, uniform="col")
        master.grid_rowconfigure(tuple(range(10)), weight=1, uniform="row")

        super().__init__(master, **kwargs)

        self.control_panel = ctk.CTkFrame(
            master,
            corner_radius=15,
            border_width=2,
        )
        self.control_panel.grid(
            row=2, column=2, sticky=tk.NSEW, rowspan=3, columnspan=2, padx=40, pady=40
        )

        self.control_panel_label = ctk.CTkLabel(
            master=self.control_panel,
            text="Control panel",
            font=self.font(20),
        )
        self.control_panel_label.pack(pady=30)

        self.control_panel_frame = ctk.CTkFrame(
            master=self.control_panel,
            corner_radius=15,
            border_width=2,
        )
        self.control_panel_frame.pack(padx=20, pady=(0, 20), fill=tk.BOTH, expand=True)

        self.cmd_output_frame = ctk.CTkFrame(
            master,
            corner_radius=15,
            border_width=2,
        )

        self.cmd_output_frame.grid(
            row=0, column=6, rowspan=10, columnspan=5, sticky=tk.NSEW, padx=40, pady=40
        )

        self.cmd_output_label = ctk.CTkLabel(
            self.cmd_output_frame,
            text="Console output",
            font=self.font(20),
        )

        self.cmd_output_label.pack(pady=30)

        self.cmd_output = ctk.CTkTextbox(
            self.cmd_output_frame,
            font=self.font(14),
            corner_radius=15,
            border_width=2,
            state=tk.DISABLED,
        )

        self.cmd_output.pack(padx=20, pady=(0, 20), fill=tk.BOTH, expand=True)

        self.cmd_output.tag_config(Tag.ERROR, foreground="red")
        self.cmd_output.tag_config(Tag.SUCCESS, foreground="green")
