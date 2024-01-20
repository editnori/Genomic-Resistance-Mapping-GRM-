from ftplib import FTP
import time
import tkinter as tk
from tkinter import ttk, messagebox, filedialog, font

import threading
import os
from subprocess import CalledProcessError, Popen, PIPE, run
import csv
import random
from enum import Enum
import traceback
import re

import customtkinter as ctk
from PIL import Image
import pandas as pd

from kover import KoverDatasetCreator
from hovertip import Hovertip
from table import Table
from label import Label
from ftp_downloader import (
    FTPDownloadApp,
    DownloadWindow,
    get_last_metadata_update_date,
    threaded,
)


class Page(Enum):
    DATA_COLLECTION_PAGE = 0
    PREPROCESSING_PAGE = 1
    KOVER_LEARN_PAGE = 2
    ANALYSIS_PAGE = 3


class Path(str):
    FOREST_DARK = "ui/forest-dark.tcl"
    RAY = "bin/ray/Ray"
    REMOTE_METADATA = "RELEASE_NOTES/genome_metadata"
    RAY_SURVEYOR = "raysurveyor/"
    IMAGES = "ui/test_images/"
    DATA = "data/"
    CONTIGS = DATA + "contigs/"
    FEATURES = DATA + "features/"


class App(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.setup_window()

        self.setup_style()

        self.load_images(Path.IMAGES)

        self.to_remove()  # remove this line when done

        self.create_navigation_frame()

        self.create_data_collection_page()

        self.create_preprocessing_page()

        self.create_kover_learn_page()

        self.create_analysis_page()

        self.set_page(Page.DATA_COLLECTION_PAGE)

    def setup_window(self):
        self.title("Genome analysis tool")

        self.screen_width = self.winfo_screenwidth()
        self.screen_height = self.winfo_screenheight()

        self.geometry(f"{self.screen_width}x{self.screen_height - 100}")
        self.geometry("+0+0")

        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)

    def setup_style(self):
        self.style = ttk.Style(self)
        self.tk.call("source", Path.FOREST_DARK)
        self.style.theme_use("forest-dark")

        self.custom_font = font.nametofont("TkDefaultFont")
        self.custom_font.configure(size=12)

        self.default_font = lambda size: ("Century Gothic", size)

    def load_images(self, path: str) -> dict[ctk.CTkImage]:
        self.images = {}

        self.images["logo"] = ctk.CTkImage(
            Image.open(path + "CustomTkinter_logo_single.png"),
            size=(25, 25),
        )

        self.images["large_test"] = ctk.CTkImage(
            Image.open(path + "large_test_image.png"),
            size=(500, 150),
        )

        self.images["icon"] = ctk.CTkImage(
            Image.open(path + "image_icon_light.png"), size=(25, 25)
        )

        self.images["home"] = ctk.CTkImage(
            light_image=Image.open(path + "database-dark.png"),
            dark_image=Image.open(path + "database-light1.png"),
            size=(25, 25),
        )

        prepocessing_image = Image.open(path + "preprocessing.png")

        self.images["chat"] = ctk.CTkImage(
            light_image=prepocessing_image,
            dark_image=prepocessing_image,
            size=(25, 25),
        )

        self.images["add_user"] = ctk.CTkImage(
            light_image=Image.open(path + "add_user_dark.png"),
            dark_image=Image.open(path + "add_user_light.png"),
            size=(25, 25),
        )

    def create_navigation_frame(self):
        button_height = 80
        button_text_color = ("gray10", "gray90")
        button_hover_color = ("gray70", "gray30")

        self.navigation_frame = ctk.CTkFrame(self, corner_radius=0)

        self.navigation_frame.grid(row=0, column=0, sticky=tk.NSEW)
        self.navigation_frame.grid_rowconfigure(5, weight=1)

        self.navigation_frame_label = ctk.CTkLabel(
            self.navigation_frame,
            text="Patric",
            image=self.images["logo"],
            padx=10,
            compound="left",
            font=ctk.CTkFont(size=15, weight="bold"),
        )

        self.navigation_frame_label.grid(row=0, column=0, pady=20)

        self.data_collection_frame_button = ctk.CTkButton(
            self.navigation_frame,
            corner_radius=0,
            height=button_height,
            text="Data collection",
            border_spacing=10,
            text_color=button_text_color,
            hover_color=button_hover_color,
            image=self.images["home"],
            anchor="w",
            command=lambda: self.set_page(Page.DATA_COLLECTION_PAGE),
        )

        self.data_collection_frame_button.grid(row=1, column=0, sticky="ew")

        self.preprocessing_frame_button = ctk.CTkButton(
            self.navigation_frame,
            corner_radius=0,
            height=button_height,
            text="Data preprocessing",
            border_spacing=10,
            text_color=button_text_color,
            hover_color=button_hover_color,
            image=self.images["chat"],
            anchor="w",
            command=lambda: self.set_page(Page.PREPROCESSING_PAGE),
        )

        self.preprocessing_frame_button.grid(row=2, column=0, sticky="ew")

        self.kover_frame_button = ctk.CTkButton(
            self.navigation_frame,
            corner_radius=0,
            height=button_height,
            text="Kover learn",
            border_spacing=10,
            text_color=button_text_color,
            hover_color=button_hover_color,
            image=self.images["add_user"],
            anchor="w",
            command=lambda: self.set_page(Page.KOVER_LEARN_PAGE),
        )

        self.kover_frame_button.grid(row=3, column=0, sticky="ew")

        self.analysis_frame_button = ctk.CTkButton(
            self.navigation_frame,
            corner_radius=0,
            height=button_height,
            text="Analysis",
            border_spacing=10,
            text_color=button_text_color,
            hover_color=button_hover_color,
            image=self.images["add_user"],
            anchor="w",
            command=lambda: self.set_page(Page.ANALYSIS_PAGE),
        )

        self.analysis_frame_button.grid(row=4, column=0, sticky="ew")

        self.appearance_mode_menu = ctk.CTkOptionMenu(
            self.navigation_frame,
            values=["Dark", "Light", "System"],
            command=ctk.set_appearance_mode,
        )

        self.appearance_mode_menu.grid(row=6, column=0, pady=15, sticky="s")

    def create_data_collection_page(self):
        self.data_collection_frame = ctk.CTkFrame(self, fg_color="transparent")

        self.data_collection_tab_view = ctk.CTkTabview(self.data_collection_frame)

        self.data_collection_tab_view.pack(fill=tk.BOTH, expand=True)
        self.data_collection_tab_view.add("Genomes")
        self.data_collection_tab_view.add("AMR")

        self.create_genomes_tab()

        self.create_amr_tab()

    def create_genomes_tab(self):
        self.genomes_frame = ctk.CTkFrame(self.data_collection_tab_view.tab("Genomes"))

        self.genomes_frame.place(relx=0.5, rely=0.5, anchor=tk.CENTER)

        self.genome_data_frame = ctk.CTkFrame(
            self.genomes_frame,
            corner_radius=15,
            border_width=2,
        )

        self.genome_data_frame.pack(side=tk.LEFT, padx=20)

        self.genome_data_frame_label = ctk.CTkLabel(
            master=self.genome_data_frame,
            text="Genome Data",
            font=self.default_font(20),
        )

        self.genome_data_frame_label.pack(padx=40, pady=(45, 50), fill=tk.X)

        self.genome_data_frame_checkbox_frame = ctk.CTkFrame(
            master=self.genome_data_frame, fg_color="transparent"
        )

        self.genome_data_frame_checkbox_frame.pack(padx=50, pady=5, fill=tk.X)

        self.genome_data_frame_contig_checkbox = ctk.CTkCheckBox(
            master=self.genome_data_frame_checkbox_frame,
            text="Contigs",
            command=self.genome_data_validate_ui,
        )

        self.genome_data_frame_contig_checkbox.pack(side=tk.LEFT)

        self.genome_data_frame_feature_checkbox = ctk.CTkCheckBox(
            master=self.genome_data_frame_checkbox_frame,
            text="Features",
            command=self.genome_data_validate_ui,
        )

        self.genome_data_frame_feature_checkbox.pack(side=tk.LEFT)

        self.genome_data_frame_load_frame = ctk.CTkFrame(
            master=self.genome_data_frame, fg_color="transparent"
        )

        self.genome_data_frame_load_frame.pack(padx=50, pady=5, fill=tk.X)

        self.genome_data_frame_bulk_checkbox = ctk.CTkCheckBox(
            master=self.genome_data_frame_load_frame,
            text="bulk",
            command=self.toggle_bulk_download,
        )

        self.genome_data_frame_bulk_checkbox.grid(row=0, column=0, sticky=tk.W)

        self.genome_data_frame_entry = ctk.CTkEntry(
            master=self.genome_data_frame_load_frame, placeholder_text="Genome ID"
        )

        self.genome_data_frame_entry.grid(row=0, column=1, sticky=tk.W)

        self.genome_data_frame_entry.bind("<KeyRelease>", self.genome_data_validate_ui)

        Hovertip(
            self.genome_data_frame_entry,
            "Genome ID format:\nnumber.number\n\ne.g. 123.456",
        )

        self.genome_data_frame_bulk_button = ctk.CTkButton(
            master=self.genome_data_frame_load_frame,
            text="Load IDs",
            corner_radius=6,
            command=self.select_genome_data_directory,
        )

        self.genome_data_frame_path_saved = ""

        self.genome_data_frame_path = Label(
            master=self.genome_data_frame,
            text=self.genome_data_frame_path_saved,
            fg_color="transparent",
            width=37,
            anchor="w",
        )

        self.genome_data_frame_path.pack(padx=50, pady=5, fill=tk.X)

        self.genome_data_frame_download_button = ctk.CTkButton(
            master=self.genome_data_frame,
            text="Download",
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            command=self.download_genome_data,
        )

        self.genome_data_frame_download_button.pack(padx=50, pady=(5, 30), fill=tk.X)

        self.genome_data_frame_download_button_hover = Hovertip(
            self.genome_data_frame_download_button,
            "",
        )

        self.genome_data_validate_ui()

        self.genome_metadata_frame = ctk.CTkFrame(
            self.genomes_frame,
            corner_radius=15,
            border_width=2,
        )

        self.genome_metadata_frame.pack(side=tk.LEFT, padx=20)

        self.genome_metadata_frame_label = ctk.CTkLabel(
            master=self.genome_metadata_frame,
            text="Latest metadata for Genomes",
            font=self.default_font(20),
        )

        self.genome_metadata_frame_label.pack(padx=40, pady=(45, 50), fill=tk.X)

        self.genome_metadata_frame_download_button = ctk.CTkButton(
            master=self.genome_metadata_frame,
            text="Download",
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            command=self.download_genome_metadata,
        )

        self.genome_metadata_frame_download_button.pack(
            padx=50, pady=(5, 30), fill=tk.X
        )

    def toggle_bulk_download(self, event=None):
        if self.genome_data_frame_bulk_checkbox.get():
            self.genome_data_frame_entry.grid_remove()
            self.genome_data_frame_bulk_button.grid(row=0, column=1, sticky=tk.W)
            self.genome_data_frame_path.configure(
                text=self.genome_data_frame_path_saved
            )
        else:
            self.genome_data_frame_entry.grid(row=0, column=1, sticky=tk.W)
            self.genome_data_frame_bulk_button.grid_remove()
            self.genome_data_frame_path_saved = self.genome_data_frame_path.cget("text")
            self.genome_data_frame_path.configure(text="")

        self.genome_data_validate_ui()

    def select_genome_data_directory(self):
        self.genome_data_frame_path.configure(
            text=filedialog.askopenfilename(filetypes=[("TSV file", "*.tsv")])
        )

        self.genome_data_validate_ui()

    @threaded
    def download_genome_data(self):
        if self.genome_data_frame_bulk_checkbox.get():
            genome_data_ids = pd.read_table(
                self.genome_data_frame_path.cget("text"),
                usecols=["genome_id", "genome_name"],
                converters={"genome_id": str, "genome_name": str},
            )

            genome_data_ids = genome_data_ids[
                genome_data_ids["genome_id"].str.contains(r"^\d+\.\d+$")
            ].values
        else:
            genome_data_ids = [
                (self.genome_data_frame_entry.get(), self.genome_data_frame_entry.get())
            ]

        directory = filedialog.askdirectory()

        if not directory:
            return

        self.genome_data_frame_contig_checkbox.configure(state=tk.DISABLED)
        self.genome_data_frame_feature_checkbox.configure(state=tk.DISABLED)
        self.genome_data_frame_bulk_checkbox.configure(state=tk.DISABLED)
        self.genome_data_frame_entry.configure(state=tk.DISABLED)
        self.genome_data_frame_download_button.configure(
            text="Cancel", command=self.cancel_genome_data_download
        )
        self.cancel_genome_data_download_boolean = False

        number_downloaded = 0

        total_genomes = len(genome_data_ids) * (
            self.genome_data_frame_contig_checkbox.get()
            + self.genome_data_frame_feature_checkbox.get()
        )

        total_progress_bar = ttk.Progressbar(
            master=self.genome_data_frame, mode="determinate", maximum=total_genomes
        )

        total_progress_bar.pack(padx=50, pady=(5, 0), fill=tk.X)

        total_progress_bar["value"] = 0

        total_progress_label = ctk.CTkLabel(
            master=self.genome_data_frame,
            font=self.default_font(10),
            fg_color="transparent",
            text="",
            anchor="w",
        )

        total_progress_label.pack(padx=50, pady=(0, 20), fill=tk.X)

        for genome_data_id, genome_name in genome_data_ids:
            contig_name = f"{genome_data_id}.fna"
            feature_name = f"{genome_data_id}.PATRIC.features.tab"

            local_contig_directory = os.path.join(directory, Path.CONTIGS, genome_name)
            local_feature_directory = os.path.join(
                directory, Path.FEATURES, genome_name
            )

            local_contig_path = os.path.join(local_contig_directory, contig_name)
            local_feature_path = os.path.join(local_feature_directory, feature_name)

            remote_contig = f"genomes/{genome_data_id}/{genome_data_id}.fna"
            remote_feature = (
                f"genomes/{genome_data_id}/{genome_data_id}.PATRIC.features.tab"
            )

            if self.genome_data_frame_contig_checkbox.get():
                os.makedirs(local_contig_directory, exist_ok=True)

                number_downloaded += 1
                total_progress_label.configure(
                    text=f"Downloading: {contig_name} ({number_downloaded}/{total_genomes})"
                )
                total_progress_bar["value"] = number_downloaded

                with open(local_contig_path, "wb") as local_file:
                    ftp = FTP("ftp.bvbrc.org")
                    ftp.login()
                    ftp.voidcmd("TYPE I")

                    contig_size = ftp.size(remote_contig)

                    contig_progress_bar = ttk.Progressbar(
                        master=self.genome_data_frame,
                        mode="determinate",
                        maximum=contig_size,
                    )

                    contig_size_mb = contig_size / 1_048_576

                    contig_progress_bar.pack(padx=50, pady=(5, 0), fill=tk.X)

                    contig_progress_bar["value"] = 0

                    contig_progress_label = ctk.CTkLabel(
                        master=self.genome_data_frame,
                        font=self.default_font(10),
                        fg_color="transparent",
                        text="",
                        anchor="w",
                    )

                    contig_progress_label.pack(padx=50, pady=(0, 20), fill=tk.X)

                    with ftp.transfercmd(f"RETR {remote_contig}") as conn:
                        bytes_received = 0
                        last_update_time = time.time_ns() / 1_000_000
                        while not self.cancel_genome_data_download_boolean and (
                            data := conn.recv(1024)
                        ):
                            local_file.write(data)
                            bytes_received += len(data)
                            current_time_ms = time.time_ns() / 1_000_000
                            if (current_time_ms - last_update_time) > 100:
                                contig_progress_bar["value"] = bytes_received
                                contig_progress_label.configure(
                                    text=f"Downloaded: {bytes_received / 1_048_576:6.2f} MB / {contig_size_mb:6.2f} MB"
                                )
                                last_update_time = current_time_ms

                    contig_progress_bar.destroy()
                    contig_progress_label.destroy()

                if self.cancel_genome_data_download_boolean:
                    os.remove(local_contig_path)
                    break
            if self.genome_data_frame_feature_checkbox.get():
                os.makedirs(local_feature_directory, exist_ok=True)

                number_downloaded += 1
                total_progress_label.configure(
                    text=f"Downloading: {feature_name} ({number_downloaded}/{total_genomes})"
                )
                total_progress_bar["value"] = number_downloaded

                with open(local_feature_path, "wb") as local_file:
                    ftp = FTP("ftp.bvbrc.org")
                    ftp.login()
                    ftp.voidcmd("TYPE I")

                    feature_size = ftp.size(remote_feature)

                    feature_progress_bar = ttk.Progressbar(
                        master=self.genome_data_frame,
                        mode="determinate",
                        maximum=feature_size,
                    )

                    feature_size_mb = feature_size / 1_048_576

                    feature_progress_bar.pack(padx=50, pady=(5, 0), fill=tk.X)

                    feature_progress_bar["value"] = 0

                    feature_progress_label = ctk.CTkLabel(
                        master=self.genome_data_frame,
                        font=self.default_font(10),
                        fg_color="transparent",
                        text="",
                        anchor="w",
                    )

                    feature_progress_label.pack(padx=50, pady=(0, 20), fill=tk.X)

                    with ftp.transfercmd(f"RETR {remote_feature}") as conn:
                        bytes_received = 0
                        last_update_time = time.time_ns() / 1_000_000
                        while not self.cancel_genome_data_download_boolean and (
                            data := conn.recv(1024)
                        ):
                            local_file.write(data)
                            bytes_received += len(data)
                            if (current_time_ms - last_update_time) > 100:
                                feature_progress_bar["value"] = bytes_received
                                feature_progress_label.configure(
                                    text=f"Downloaded: {bytes_received / 1_048_576:6.2f} MB / {feature_size_mb:6.2f} MB"
                                )
                                last_update_time = current_time_ms

                    feature_progress_bar.destroy()
                    feature_progress_label.destroy()

                if self.cancel_genome_data_download_boolean:
                    os.remove(local_feature_path)
                    break
        self.genome_data_frame_contig_checkbox.configure(state=tk.NORMAL)
        self.genome_data_frame_feature_checkbox.configure(state=tk.NORMAL)
        self.genome_data_frame_bulk_checkbox.configure(state=tk.NORMAL)
        self.genome_data_frame_entry.configure(state=tk.NORMAL)
        self.genome_data_frame_download_button.configure(
            text="Download", command=self.download_genome_data
        )

        total_progress_bar.destroy()
        total_progress_label.destroy()

    def cancel_genome_data_download(self):
        if messagebox.askyesno(
            "Confirmation", "Are you sure you want to cancel the download?"
        ):
            self.cancel_genome_data_download_boolean = True

    @threaded
    def download_genome_metadata(self):
        directory = filedialog.askdirectory()

        if not directory:
            return

        self.cancel_genome_metadata_download_boolean = False

        self.genome_metadata_frame_download_button.configure(
            text="Cancel", command=self.cancel_genome_metadata_download
        )

        path = os.path.join(directory, "genome_metadata")

        with open(path, "wb") as local_file:
            ftp = FTP("ftp.bvbrc.org")
            ftp.login()
            ftp.voidcmd("TYPE I")

            metadata_size = ftp.size(Path.REMOTE_METADATA)

            progress_bar = ttk.Progressbar(
                master=self.genome_metadata_frame,
                mode="determinate",
                maximum=metadata_size,
            )

            metadata_size_mb = metadata_size / 1_048_576

            progress_bar.pack(padx=50, pady=(5, 0), fill=tk.X)

            progress_bar["value"] = 0

            progress_label = ctk.CTkLabel(
                master=self.genome_metadata_frame,
                font=self.default_font(10),
                fg_color="transparent",
                text="",
                anchor="w",
            )

            progress_label.pack(padx=50, pady=(0, 20), fill=tk.X)

            with ftp.transfercmd(f"RETR {Path.REMOTE_METADATA}") as conn:
                bytes_received = 0
                last_update_time = time.time_ns() / 1_000_000
                while not self.cancel_genome_metadata_download_boolean and (
                    data := conn.recv(1024)
                ):
                    local_file.write(data)
                    bytes_received += len(data)
                    current_time_ms = time.time_ns() / 1_000_000
                    if (current_time_ms - last_update_time) > 100:
                        progress_bar["value"] = bytes_received
                        progress_label.configure(
                            text=f"Downloaded: {bytes_received / 1_048_576:6.2f} MB / {metadata_size_mb:6.2f} MB"
                        )
                        last_update_time = current_time_ms

            progress_bar.destroy()
            progress_label.destroy()

        if self.cancel_genome_metadata_download_boolean:
            os.remove(path)

        self.genome_metadata_frame_download_button.configure(
            text="Download", command=self.download_genome_metadata
        )

    def cancel_genome_metadata_download(self):
        if messagebox.askyesno(
            "Confirmation", "Are you sure you want to cancel the download?"
        ):
            self.cancel_genome_metadata_download_boolean = True

    def genome_data_validate_ui(self, event=None):
        self.genome_data_frame_download_button_hover.text = ""

        failed = False

        if not (
            self.genome_data_frame_contig_checkbox.get()
            or self.genome_data_frame_feature_checkbox.get()
        ):
            self.genome_data_frame_download_button_hover.text += (
                "• select at least one data type to download.\n"
            )
            failed = True
        if self.genome_data_frame_bulk_checkbox.get():
            if not self.genome_data_frame_path.cget("text"):
                self.genome_data_frame_download_button_hover.text += (
                    "• select genome ids tsv file.\n"
                )
                failed = True
        else:
            if not self.genome_data_frame_entry.get():
                self.genome_data_frame_entry.configure(border_color="#565B5E")
                self.genome_data_frame_download_button_hover.text += (
                    "• input genome id.\n"
                )
                failed = True
            else:
                if (
                    re.match(r"^\d+\.\d+$", self.genome_data_frame_entry.get())
                    is not None
                ):
                    self.genome_data_frame_entry.configure(border_color="green")
                else:
                    self.genome_data_frame_entry.configure(border_color="red")
                    self.genome_data_frame_download_button_hover.text += (
                        "• input a valid genome id.\n"
                    )
                    failed = True

        self.genome_data_frame_download_button_hover.text = (
            self.genome_data_frame_download_button_hover.text.strip("\n")
        )

        if failed:
            self.genome_data_frame_download_button_hover.enable()
            self.genome_data_frame_download_button.configure(state=tk.DISABLED)
        else:
            self.genome_data_frame_download_button_hover.disable()
            self.genome_data_frame_download_button.configure(state=tk.NORMAL)

    def create_amr_tab(self):
        frame4 = ctk.CTkFrame(
            self.data_collection_tab_view.tab("AMR"),
            width=1000,
            height=400,
            corner_radius=15,
            border_width=2,
        )
        frame4.place(x=50, y=470)

        metadata_amr_label = ctk.CTkLabel(
            master=frame4, text="Latest metadata for AMR", font=self.default_font(20)
        )
        metadata_amr_label.place(x=50, y=45)

        # Create custom button
        self.update_date = ctk.CTkLabel(
            master=frame4, text="", font=self.default_font(12), fg_color="transparent"
        )
        self.update_date.place(x=50, y=90)
        dirbtn3 = ctk.CTkButton(
            master=frame4,
            width=150,
            text="Select directory",
            corner_radius=6,
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            command=self,
        )
        dirbtn3.place(x=50, y=120)

        self.viewupdate_date = ctk.CTkButton(
            master=frame4,
            width=150,
            text="View last update date",
            corner_radius=6,
            command=lambda: get_last_metadata_update_date(self.update_date),
        )
        self.viewupdate_date.place(x=250, y=120)

        self.download_button3 = ctk.CTkButton(
            master=frame4, width=150, text="Download", corner_radius=6, command=self
        )
        self.download_button3.place(x=450, y=120)
        self.cancel_button3 = ctk.CTkButton(
            master=frame4, width=150, text="Cancel", corner_radius=6, command=self
        )
        self.cancel_button3.place(x=650, y=120)

        self.progress_bar3 = ttk.Progressbar(
            master=frame4, length=800, mode="determinate"
        )
        self.progress_bar3.place(x=50, y=190)

        self.size_label3 = ctk.CTkLabel(
            master=frame4, text="", font=self.default_font(10), fg_color="transparent"
        )
        self.size_label3.place(x=50, y=200)

        self.download_window2 = DownloadWindow(
            frame4,
            self.download_button3,
            self.cancel_button3,
            dirbtn3,
            self.size_label3,
            self.progress_bar3,
            self.size_label3,
        )

        self.download_app = FTPDownloadApp(
            self.download_window2, "RELEASE_NOTES/PATRIC_genomes_AMR.txt"
        )
        # logic

        # creating AMR metadata frame
        amr_frame = ctk.CTkFrame(
            self.data_collection_tab_view.tab("AMR"),
            width=500,
            height=400,
            corner_radius=15,
            border_width=2,
        )

        amr_frame.place(x=50, y=45)

        list_amr_label = ctk.CTkLabel(
            master=amr_frame,
            text="List available AMR datasets",
            font=self.default_font(20),
        )
        list_amr_label.place(x=50, y=20)
        self.download_button4 = ctk.CTkButton(
            master=amr_frame,
            text="load amr list",
            corner_radius=6,
            command=self.load_amr_data,
        )
        self.download_button4.place(x=50, y=60)
        self.species_filter = ttk.Combobox(master=amr_frame, state=tk.DISABLED)
        self.species_filter.bind("<<ComboboxSelected>>", self.update_table)
        self.species_filter.bind(
            "<Control-BackSpace>", lambda e: self.species_filter.set("")
        )
        self.species_filter.bind("<KeyRelease>", self.update_table)
        self.species_filter.bind(
            "<FocusOut>", lambda _: self.species_filter.selection_clear()
        )

        self.species_filter.place(x=280, y=60)

        self.amr_list_filter_checkbox = ctk.CTkCheckBox(
            master=amr_frame,
            text="Phenotype count ≥ 25",
            command=self.on_amr_list_filter_check,
            state=tk.DISABLED,
        )

        self.amr_list_filter_checkbox.place(x=50, y=100)

        Hovertip(
            self.amr_list_filter_checkbox,
            "Only show species with ≥ 25\nphenotypes (resistant and susceptible)",
        )

        self.total_label = tk.Label(master=amr_frame, text="Total: not loaded")
        self.total_label.place(x=280, y=100)
        columns = ["Species", "Antibiotic"]
        self.amr_list_table = Table(
            master=amr_frame,
            columns=columns,
            show="headings",
            selectmode="browse",
            height=9,
        )

        self.species_filter.bind("<Return>", lambda _: self.amr_list_table.focus_set())

        self.amr_list_table.place(x=50, y=130)
        for col in columns:
            self.amr_list_table.column(col, width=185)
            self.amr_list_table.heading(col, text=col)

        self.amr_list_table.bind("<Double-1>", self.on_table_select)
        self.amr_list_table.bind("<Return>", self.on_table_select)
        self.amr_list_table.bind("<space>", self.on_table_select)

        self.amr_list = pd.DataFrame()

        frame6 = ctk.CTkFrame(
            self.data_collection_tab_view.tab("AMR"),
            width=800,
            height=400,
            corner_radius=15,
            border_width=2,
        )
        frame6.place(x=580, y=45)
        # Create input fields for antibiotic and species
        full_amr_label = ctk.CTkLabel(
            master=frame6,
            text="Get amr data by species and antibiotic",
            font=self.default_font(20),
        )
        full_amr_label.place(x=50, y=20)
        self.antibiotic_selection = ttk.Combobox(
            master=frame6, state=tk.DISABLED, width=18
        )
        self.antibiotic_selection.bind(
            "<<ComboboxSelected>>", self.on_antibiotic_select
        )
        self.antibiotic_selection.place(x=50, y=50)

        self.species_selection = ttk.Combobox(
            master=frame6, state=tk.DISABLED, width=18
        )

        self.species_selection.bind("<<ComboboxSelected>>", self.update_amr_full)
        self.species_selection.place(x=220, y=50)

        self.drop_intermediate_checkbox = ctk.CTkCheckBox(
            master=frame6,
            text="Drop Intermediate",
            command=self.update_amr_full,
            state=tk.DISABLED,
        )
        self.drop_intermediate_checkbox.place(x=50, y=90)

        self.numeric_phenotypes_checkbox = ctk.CTkCheckBox(
            master=frame6,
            text="Numeric Phenotypes",
            command=self.update_amr_full,
            state=tk.DISABLED,
        )

        self.numeric_phenotypes_checkbox.place(x=220, y=90)

        Hovertip(
            self.numeric_phenotypes_checkbox, "0: Resistant\n1: Susceptible\n2: Other"
        )

        # Create a button to select the AMR metadata file
        self.save_table_button = ctk.CTkButton(
            master=frame6,
            text="Export to .tsv",
            command=self.save_to_tsv,
            state=tk.DISABLED,
        )
        self.save_table_button.place(x=400, y=50)
        self.totaldata = tk.Label(master=frame6, text="Total phenotypes:")
        self.totaldata.place(x=400, y=90)
        self.totalresistance_label = tk.Label(
            master=frame6, text="Total resistance phenotypes available:"
        )
        self.totalresistance_label.place(x=50, y=120)
        self.totalsusceptible_label = tk.Label(
            master=frame6, text="Total susceptible phenotypes available:"
        )
        self.totalsusceptible_label.place(x=400, y=120)

        # Create a table to display the results
        columns = [
            "Genome ID",
            "Genome Name",
            "Phenotype",
            "Measurements",
        ]
        self.results_table = Table(
            master=frame6,
            columns=columns,
            show="headings",
            height=8,
        )

        self.species_selection.bind(
            "<Return>", lambda _: self.results_table.focus_set()
        )
        self.antibiotic_selection.bind(
            "<Return>", lambda _: self.results_table.focus_set()
        )

        for col in columns:
            self.results_table.column(col, width=150, anchor=tk.CENTER)
            self.results_table.heading(col, text=col)

        self.results_table.place(x=50, y=150)

    def create_preprocessing_page(self):
        self.preprocessing_frame = ctk.CTkFrame(
            self, corner_radius=0, fg_color="transparent"
        )

        preprocessing_tab_view = ctk.CTkTabview(
            self.preprocessing_frame,
            width=self.screen_width - 230,
            height=self.screen_height - 150,
        )
        preprocessing_tab_view.grid(row=0, column=0, padx=(20, 0), pady=(20, 0))
        preprocessing_tab_view.add("preprocessing")
        self.controlpanel = ctk.CTkFrame(
            preprocessing_tab_view.tab("preprocessing"),
            width=550,
            height=500,
            corner_radius=15,
            border_width=2,
        )
        self.controlpanel.place(x=50, y=20)
        l2 = ctk.CTkLabel(
            master=self.controlpanel, text="Control panel", font=self.default_font(20)
        )
        l2.place(x=50, y=45)

        self.datasetbtn = ctk.CTkButton(
            master=self.controlpanel,
            width=220,
            text="Pick dataset",
            corner_radius=6,
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            command=self.browse_dataset,
        )
        self.datasetbtn.place(x=50, y=150)

        self.entry1 = ctk.CTkEntry(
            master=self.controlpanel,
            width=380,
            placeholder_text="Dataset path",
            textvariable=self.dataset_folder,
        )
        self.entry1.place(x=50, y=200)

        self.outputbtn = ctk.CTkButton(
            master=self.controlpanel,
            width=220,
            text="Pick output directory",
            corner_radius=6,
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            command=self.browse_output_dir,
        )
        self.outputbtn.place(x=50, y=250)

        self.entry2 = ctk.CTkEntry(
            master=self.controlpanel,
            width=380,
            placeholder_text="Output directory path",
            textvariable=self.output_dir,
        )
        self.entry2.place(x=50, y=300)

        l3 = ctk.CTkLabel(
            master=self.controlpanel,
            text="Enter kmer length",
            font=self.default_font(15),
        )
        l3.place(x=50, y=350)
        self.kmer_length = ctk.CTkEntry(
            master=self.controlpanel, width=100, placeholder_text="kmer length"
        )
        self.kmer_length.place(x=50, y=400)
        l4 = ctk.CTkLabel(
            master=self.controlpanel, text="Pick kmer tool", font=self.default_font(15)
        )
        l4.place(x=300, y=350)

        self.kmer_tool = tk.StringVar(value="Ray Surveyor")

        self.kmertool_filter = ttk.Combobox(
            master=self.controlpanel,
            textvariable=self.kmer_tool,
            values=["Ray Surveyor", "DSK"],
            state="readonly",
        )

        self.kmertool_filter.place(x=300, y=400)

        self.runbtn = ctk.CTkButton(
            master=self.controlpanel,
            width=150,
            text="Run " + self.kmer_tool.get(),
            corner_radius=6,
            command=self.run_preprocessing,
        )
        self.runbtn.place(x=150, y=450)

        self.kmer_tool.trace_add("write", self.update_button_text)

        outputscreen = ctk.CTkFrame(
            preprocessing_tab_view.tab("preprocessing"),
            width=650,
            height=650,
            corner_radius=15,
            border_width=2,
        )

        outputscreen.place(x=800, y=20)
        outputscreen.grid_propagate(False)
        outputscreen.grid_rowconfigure(0, weight=1)
        outputscreen.grid_columnconfigure(0, weight=1)

        self.cmd_output = tk.Text(
            master=outputscreen,
            height=650,
            width=650,
            state="disabled",
            font=self.custom_font,
        )

        self.scrollbar = tk.Scrollbar(outputscreen, command=self.cmd_output.yview)
        self.scrollbar.grid(row=0, column=1, sticky=tk.NSEW)
        self.cmd_output["yscrollcommand"] = self.scrollbar.set
        self.cmd_output.grid(row=0, column=0, sticky=tk.NSEW, padx=2, pady=2)

    def create_kover_learn_page(self):
        self.kover_frame = ctk.CTkFrame(self, corner_radius=0, fg_color="transparent")
        kover_tab_view = ctk.CTkTabview(
            self.kover_frame,
            width=self.screen_width - 230,
            height=self.screen_height - 150,
        )
        kover_tab_view.grid(row=0, column=0, padx=(20, 0), pady=(20, 0))
        kover_tab_view.add("Create dataset")
        kover_tab_view.add("split dataset")
        kover_tab_view.add("kover learn")

        create_dataset_frame = ctk.CTkFrame(
            kover_tab_view.tab("Create dataset"),
            width=700,
            height=self.screen_height - 230,
            corner_radius=15,
            border_width=2,
        )
        create_dataset_frame.place(x=50, y=20)

        l1 = ctk.CTkLabel(
            master=create_dataset_frame,
            text="Create Dataset",
            font=self.default_font(20),
        )
        l1.place(x=50, y=45)

        # Selection option for dataset type
        dataset_type_label = ctk.CTkLabel(
            master=create_dataset_frame,
            text="Select Dataset Type",
            font=self.default_font(15),
        )
        dataset_type_label.place(x=50, y=100)
        self.dataset_type_var1 = tk.StringVar(value="contigs")  # Set a default value
        dataset_type_options = ["reads", "contigs", "kmer matrix"]
        dataset_type_menu = ttk.Combobox(
            master=create_dataset_frame,
            textvariable=self.dataset_type_var1,
            values=dataset_type_options,
            state="readonly",
        )
        dataset_type_menu.place(x=50, y=150)
        dataset_type_menu.bind("<<ComboboxSelected>>", self.on_dataset_type_selected)
        # Buttons to pick dataset, output directory, phenotype description, phenotype metadata
        self.pickdataset1_btn = ctk.CTkButton(
            master=create_dataset_frame,
            width=220,
            text="Pick Dataset",
            corner_radius=6,
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            command=self.browse_dataset,
        )
        self.pickdataset1_btn.place(x=50, y=200)

        self.pickdataset1_entry = ctk.CTkEntry(
            master=create_dataset_frame,
            width=380,
            placeholder_text="Dataset path",
            textvariable=self.dataset_folder,
        )
        self.pickdataset1_entry.place(x=50, y=250)

        self.pickoutput_btn = ctk.CTkButton(
            master=create_dataset_frame,
            width=220,
            text="Pick Output Directory",
            corner_radius=6,
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            command=self.browse_output_dir,
        )
        self.pickoutput_btn.place(x=50, y=300)

        self.pickoutput_entry = ctk.CTkEntry(
            master=create_dataset_frame,
            width=380,
            placeholder_text="Output directory path",
            textvariable=self.output_dir,
        )
        self.pickoutput_entry.place(x=50, y=350)

        self.phenotype_desc_btn = ctk.CTkButton(
            master=create_dataset_frame,
            width=220,
            text="Pick Phenotype Description",
            corner_radius=6,
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            command=self.pick_desctsv_file,
        )
        self.phenotype_desc_btn.place(x=50, y=400)

        self.phenotype_desc_entry = ctk.CTkEntry(
            master=create_dataset_frame,
            width=380,
            placeholder_text="Phenotype description path",
            textvariable=self.desctsv_file_path,
        )
        self.phenotype_desc_entry.place(x=50, y=450)

        # Kmer entry with max kmer length 128
        self.kmer_label = ctk.CTkLabel(
            master=create_dataset_frame,
            text="Enter Kmer Length (max 128)",
            font=self.default_font(15),
        )
        self.kmer_label.place(x=50, y=600)
        self.kmer_length_var = tk.StringVar(value="31")
        self.kmer_length_spinbox = tk.Spinbox(
            master=create_dataset_frame,
            from_=0,
            to=128,
            width=5,
            textvariable=self.kmer_length_var,
            wrap=True,
            fg="black",
            buttonbackground="white",
        )
        self.kmer_length_spinbox.place(x=50, y=630)

        # Compression level input (0 to 9)
        compression_label = ctk.CTkLabel(
            master=create_dataset_frame,
            text="Enter Compression Level (0-9)",
            font=self.default_font(15),
        )
        compression_label.place(x=300, y=600)
        self.compression_var = tk.StringVar(value="4")

        self.compression_spinbox = tk.Spinbox(
            master=create_dataset_frame,
            from_=0,
            to=9,
            width=5,
            textvariable=self.compression_var,
            wrap=True,
            fg="black",
            buttonbackground="white",
        )
        self.compression_spinbox.place(x=300, y=630)

        # Create button to initiate the dataset creation
        self.create_dataset_btn = ctk.CTkButton(
            master=create_dataset_frame,
            width=150,
            text="Create Dataset",
            corner_radius=6,
            command=self.create_dataset_thread,
        )
        self.create_dataset_btn.place(x=500, y=630)

        phenotype_metadata_btn = ctk.CTkButton(
            master=create_dataset_frame,
            width=220,
            text="Pick Phenotype Metadata",
            corner_radius=6,
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            command=self.pick_metatsv_file,
        )
        phenotype_metadata_btn.place(x=50, y=500)

        self.phenotype_metadata_entry = ctk.CTkEntry(
            master=create_dataset_frame,
            width=380,
            placeholder_text="Phenotype metadata path",
            textvariable=self.metatsv_file_path,
        )
        self.phenotype_metadata_entry.place(x=50, y=550)

        outputscreen = ctk.CTkFrame(
            kover_tab_view.tab("Create dataset"),
            width=650,
            height=650,
            corner_radius=15,
            border_width=2,
        )

        outputscreen.place(x=800, y=20)
        outputscreen.grid_propagate(False)
        outputscreen.grid_rowconfigure(0, weight=1)
        outputscreen.grid_columnconfigure(0, weight=1)

        self.cmd_output1 = tk.Text(
            master=outputscreen,
            height=650,
            width=650,
            state="disabled",
            font=self.custom_font,
        )

        self.scrollbar1 = tk.Scrollbar(outputscreen, command=self.cmd_output1.yview)
        self.scrollbar1.grid(row=0, column=1, sticky=tk.NSEW)
        self.cmd_output1["yscrollcommand"] = self.scrollbar1.set
        self.cmd_output1.grid(row=0, column=0, sticky=tk.NSEW, padx=2, pady=2)

        create_dataset_frame = ctk.CTkFrame(
            kover_tab_view.tab("split dataset"),
            width=700,
            height=self.screen_height - 230,
            corner_radius=15,
            border_width=2,
        )
        create_dataset_frame.place(x=50, y=20)
        l1 = ctk.CTkLabel(
            master=create_dataset_frame,
            text="Split Dataset",
            font=self.default_font(20),
        )
        l1.place(x=50, y=45)

        dataset_btn = ctk.CTkButton(
            master=create_dataset_frame,
            width=220,
            text="Pick kover Dataset",
            corner_radius=6,
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            command=self.pickkover,
        )
        dataset_btn.place(x=50, y=150)

        dataset_entry = ctk.CTkEntry(
            master=create_dataset_frame,
            width=380,
            placeholder_text="Dataset path",
            textvariable=self.selected_kover,
        )
        dataset_entry.place(x=50, y=200)

        output_btn = ctk.CTkButton(
            master=create_dataset_frame,
            width=220,
            text="Pick Output Directory",
            corner_radius=6,
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            command=self.browse_output_dir,
        )
        output_btn.place(x=50, y=250)

        output_entry = ctk.CTkEntry(
            master=create_dataset_frame,
            width=380,
            placeholder_text="Output directory path",
            textvariable=self.output_dir,
        )
        output_entry.place(x=50, y=300)

        l1 = ctk.CTkLabel(
            master=create_dataset_frame,
            text="Train size(%)",
            font=self.default_font(16),
        )
        l1.place(x=50, y=350)

        self.trainsize_spinbox = tk.Spinbox(
            master=create_dataset_frame,
            from_=1,
            to=100,
            width=5,
            textvariable=self.compression_var,
            wrap=True,
        )

        self.trainsize_spinbox.place(x=170, y=355)

        trainlabel = ctk.CTkLabel(
            master=create_dataset_frame,
            text="Use train ids and tests",
            font=self.default_font(16),
        )
        trainlabel.place(x=250, y=350)
        self.trainvar = tk.StringVar(value="no")  # Set a default value
        train_options = ["yes", "no"]
        train_options_menu = ttk.Combobox(
            master=create_dataset_frame,
            textvariable=self.trainvar,
            values=train_options,
            width=10,
            state="readonly",
        )
        train_options_menu.place(x=420, y=350)

        train_options_menu.bind("<<ComboboxSelected>>", self.toggle_train_test_widgets)

        self.foldvar = tk.StringVar(value="2")
        foldlabel = ctk.CTkLabel(
            master=create_dataset_frame, text="Folds", font=self.default_font(16)
        )
        foldlabel.place(x=50, y=400)
        foldentry = ctk.CTkEntry(
            master=create_dataset_frame,
            width=80,
            placeholder_text="Enter number of folds",
            textvariable=self.foldvar,
        )

        foldentry.place(x=100, y=400)

        randomseedlabel = ctk.CTkLabel(
            master=create_dataset_frame, text="Random seed", font=self.default_font(16)
        )
        randomseedlabel.place(x=200, y=400)
        self.randomseedentry = ctk.CTkEntry(
            master=create_dataset_frame,
            width=100,
            placeholder_text="random seed",
        )

        self.randomseedentry.place(x=330, y=400)
        randomseed_btn = ctk.CTkButton(
            master=create_dataset_frame,
            width=150,
            text="Generate random seed",
            command=lambda: self.generate_random_seed(self.randomseedentry),
        )
        randomseed_btn.place(x=450, y=400)

        self.trainids_btn = ctk.CTkButton(
            master=create_dataset_frame,
            width=220,
            text="Pick Train ids file",
            corner_radius=6,
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            command=self.browse_dataset,
        )

        self.trainids_entry = ctk.CTkEntry(
            master=create_dataset_frame,
            width=200,
            placeholder_text="train ids file path",
            textvariable=self.dataset_folder,
        )

        self.testids_btn = ctk.CTkButton(
            master=create_dataset_frame,
            width=220,
            text="Pick Test ids file",
            corner_radius=6,
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            command=self.browse_dataset,
        )

        self.testids_entry = ctk.CTkEntry(
            master=create_dataset_frame,
            width=200,
            placeholder_text="train ids file path",
            textvariable=self.dataset_folder,
        )

        self.uniqueidlabel = ctk.CTkLabel(
            master=create_dataset_frame,
            text="Unique split id",
            font=self.default_font(16),
        )
        self.uniqueidlabel.place(x=50, y=450)
        self.uniqueidentry = ctk.CTkEntry(
            master=create_dataset_frame,
            width=100,
            placeholder_text="Eg. split1",
        )

        self.uniqueidentry.place(x=180, y=450)
        self.split_btn = ctk.CTkButton(
            master=create_dataset_frame,
            width=150,
            text="Split dataset",
            command=self.split_dataset_thread,
        )
        self.split_btn.place(x=200, y=630)

        outputscreen = ctk.CTkFrame(
            kover_tab_view.tab("split dataset"),
            width=650,
            height=650,
            corner_radius=15,
            border_width=2,
        )

        outputscreen.place(x=800, y=20)
        outputscreen.grid_propagate(False)
        outputscreen.grid_rowconfigure(0, weight=1)
        outputscreen.grid_columnconfigure(0, weight=1)

        self.cmd_output2 = tk.Text(
            master=outputscreen,
            height=650,
            width=650,
            state="disabled",
            font=self.custom_font,
        )

        self.scrollbar2 = tk.Scrollbar(outputscreen, command=self.cmd_output2.yview)
        self.scrollbar2.grid(row=0, column=1, sticky=tk.NSEW)
        self.cmd_output2["yscrollcommand"] = self.scrollbar2.set
        self.cmd_output2.grid(row=0, column=0, sticky=tk.NSEW, padx=2, pady=2)

        create_dataset_frame = ctk.CTkFrame(
            kover_tab_view.tab("kover learn"),
            width=700,
            height=self.screen_height - 230,
            corner_radius=15,
            border_width=2,
        )
        create_dataset_frame.place(x=50, y=20)

        l1 = ctk.CTkLabel(
            master=create_dataset_frame, text="Kover learn", font=self.default_font(20)
        )
        l1.place(x=50, y=45)

        # Selection option for dataset type
        dataset_type_label = ctk.CTkLabel(
            master=create_dataset_frame,
            text="Select kover learn Type",
            font=self.default_font(15),
        )
        dataset_type_label.place(x=50, y=100)
        self.learn_type_var = tk.StringVar(value="SCM")  # Set a default value
        learn_type_options = ["SCM", "CART"]
        learn_type_menu = ttk.Combobox(
            master=create_dataset_frame,
            textvariable=self.learn_type_var,
            values=learn_type_options,
            state="readonly",
        )
        learn_type_menu.place(x=50, y=150)
        learn_type_menu.bind("<<ComboboxSelected>>", self.on_learn_type_selected)
        # Buttons to pick dataset, output directory, phenotype description, phenotype metadata
        self.pickdataset_btn = ctk.CTkButton(
            master=create_dataset_frame,
            width=220,
            text="Pick Dataset",
            corner_radius=6,
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            command=self.pickkover,
        )
        self.pickdataset_btn.place(x=50, y=200)

        self.pickdataset_entry = ctk.CTkEntry(
            master=create_dataset_frame,
            width=380,
            placeholder_text="Dataset path",
            textvariable=self.selected_kover,
        )
        self.pickdataset_entry.place(x=280, y=200)

        self.pickoutput_btn = ctk.CTkButton(
            master=create_dataset_frame,
            width=220,
            text="Pick Output Directory",
            corner_radius=6,
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            command=self.browse_output_dir,
        )
        self.pickoutput_btn.place(x=50, y=250)

        self.pickoutput_entry = ctk.CTkEntry(
            master=create_dataset_frame,
            width=380,
            placeholder_text="Output directory path",
            textvariable=self.output_dir,
        )
        self.pickoutput_entry.place(x=280, y=250)

        self.hyperparameterlabel = ctk.CTkLabel(
            master=create_dataset_frame,
            text="Hyper-parameter:",
            font=self.default_font(12),
        )
        self.hyperparameterlabel.place(x=50, y=300)
        self.hyperparameterentry = ctk.CTkEntry(
            master=create_dataset_frame,
            width=300,
            placeholder_text="Enter hyperparamters eg: 0.1 0.178 0.316",
        )
        self.hyperparameterentry.place(x=180, y=300)

        self.modellabel = ctk.CTkLabel(
            master=create_dataset_frame, text="Model type:", font=self.default_font(12)
        )
        self.modellabel.place(x=50, y=350)
        self.model_type_var = tk.StringVar(value="conjunction")  # Set a default value
        self.model_type_options = ["conjunction", "disjunction"]
        self.model_type_menu = ttk.Combobox(
            master=create_dataset_frame,
            textvariable=self.model_type_var,
            values=self.model_type_options,
            width=10,
            state="readonly",
        )
        self.model_type_menu.place(x=130, y=350)

        splitidentifierlabel = ctk.CTkLabel(
            master=create_dataset_frame,
            text="Split identifier",
            font=self.default_font(12),
        )
        splitidentifierlabel.place(x=280, y=350)
        self.splitidentifierentry = ctk.CTkEntry(
            master=create_dataset_frame,
            width=200,
            placeholder_text="Enter split identifier in the dataset",
        )

        self.splitidentifierentry.place(x=370, y=350)

        self.maxrulelabel = ctk.CTkLabel(
            master=create_dataset_frame,
            text="Maximum rules :",
            font=self.default_font(12),
        )
        self.maxrulelabel.place(x=50, y=400)
        self.maxrules_var = tk.StringVar(value=1000)
        self.maxrules_spinbox = tk.Spinbox(
            master=create_dataset_frame,
            from_=1,
            to=10000,
            width=5,
            textvariable=self.maxrules_var,
            wrap=True,
        )

        self.maxrules_spinbox.place(x=150, y=400)

        self.maxequivrulelabel = ctk.CTkLabel(
            master=create_dataset_frame,
            text="Maximum equivalent rules :",
            font=self.default_font(12),
        )
        self.maxequivrulelabel.place(x=280, y=400)
        self.maxequivrules_var = tk.StringVar(value=1000)
        self.maxequivrules_spinbox = tk.Spinbox(
            master=create_dataset_frame,
            from_=1,
            to=10000,
            width=5,
            textvariable=self.maxequivrules_var,
            wrap=True,
        )

        self.maxequivrules_spinbox.place(x=460, y=400)

        hptypelabel = ctk.CTkLabel(
            master=create_dataset_frame, text="Hp choice :", font=self.default_font(12)
        )
        hptypelabel.place(x=50, y=450)
        self.hp_type_var = tk.StringVar(value="none")  # Set a default value
        self.hp_type_options = ["bound", "cv", "none"]
        self.hp_type_menu = ttk.Combobox(
            master=create_dataset_frame,
            textvariable=self.hp_type_var,
            values=self.hp_type_options,
            width=10,
            state="readonly",
        )
        self.hp_type_menu.place(x=130, y=450)

        self.maxboundsize = ctk.CTkLabel(
            master=create_dataset_frame,
            text="Max bound genome size :",
            font=self.default_font(12),
        )

        self.maxboundsize_var = tk.StringVar()
        self.maxboundsize_spinbox = tk.Spinbox(
            master=create_dataset_frame,
            from_=1,
            to=10000,
            width=5,
            textvariable=self.maxboundsize_var,
            wrap=True,
        )

        self.hp_type_menu.bind("<<ComboboxSelected>>", self.toggle_max_bound_widgets)

        randomseedlabel = ctk.CTkLabel(
            master=create_dataset_frame, text="Random seed", font=self.default_font(12)
        )
        randomseedlabel.place(x=50, y=500)
        self.randomseedentry2 = ctk.CTkEntry(
            master=create_dataset_frame,
            width=100,
            placeholder_text="random seed",
        )

        self.randomseedentry2.place(x=180, y=500)
        randomseed_btn2 = ctk.CTkButton(
            master=create_dataset_frame,
            width=150,
            text="Generate random seed",
            command=lambda: self.generate_random_seed(self.randomseedentry2),
        )
        randomseed_btn2.place(x=300, y=500)

        self.koverlearn_btn = ctk.CTkButton(
            master=create_dataset_frame,
            width=150,
            text="kover learn",
            command=self.kover_learn_thread,
        )
        self.koverlearn_btn.place(x=200, y=630)

        outputscreen = ctk.CTkFrame(
            kover_tab_view.tab("kover learn"),
            width=650,
            height=650,
            corner_radius=15,
            border_width=2,
        )

        outputscreen.place(x=800, y=20)
        outputscreen.grid_propagate(False)
        outputscreen.grid_rowconfigure(0, weight=1)
        outputscreen.grid_columnconfigure(0, weight=1)

        self.cmd_output3 = tk.Text(
            master=outputscreen,
            height=650,
            width=650,
            state="disabled",
            font=self.custom_font,
        )

        self.scrollbar3 = tk.Scrollbar(outputscreen, command=self.cmd_output3.yview)
        self.scrollbar3.grid(row=0, column=1, sticky=tk.NSEW)
        self.cmd_output3["yscrollcommand"] = self.scrollbar3.set
        self.cmd_output3.grid(row=0, column=0, sticky=tk.NSEW, padx=2, pady=2)

    def create_analysis_page(self):
        self.analysis_frame = ctk.CTkFrame(
            self, corner_radius=0, fg_color="transparent"
        )

    def set_page(self, page: Page):
        page_frame = {
            Page.DATA_COLLECTION_PAGE: (
                self.data_collection_frame,
                self.data_collection_frame_button,
            ),
            Page.PREPROCESSING_PAGE: (
                self.preprocessing_frame,
                self.preprocessing_frame_button,
            ),
            Page.KOVER_LEARN_PAGE: (self.kover_frame, self.kover_frame_button),
            Page.ANALYSIS_PAGE: (self.analysis_frame, self.analysis_frame_button),
        }

        for frame, button in page_frame.values():
            frame.grid_forget()
            button.configure(fg_color="transparent")

        page_frame[page][0].grid(row=0, column=1, sticky=tk.NSEW, padx=20, pady=(0, 20))
        page_frame[page][1].configure(fg_color=("gray75", "gray25"))

    def run_preprocessing(self):
        self.cmd_output["state"] = tk.NORMAL
        self.cmd_output.delete("1.0", tk.END)
        if self.kmer_tool.get() == "Ray Surveyor":
            self.run_ray_surveyor()
        elif self.kmer_tool.get() == "DSK":
            self.run_dsk()

    def run_ray_surveyor(self):
        dataset_folder = self.dataset_folder.get()
        output_dir = self.output_dir.get()
        output_dir = output_dir.replace("\\", "/")
        output_dir = output_dir.replace("C:", "/mnt/c")
        kmer_length = self.kmer_length.get()

        if (
            not dataset_folder
            or not self.ray_surveyor_dir
            or not output_dir
            or not kmer_length
        ):
            self.update_cmd_output(
                "Please fill in all required fields.\n", self.cmd_output
            )
            return

        input_files = [
            f"{dataset_folder}/{file}" for file in os.listdir(dataset_folder)
        ]

        if not input_files:
            self.update_cmd_output(
                "No .fna files found in the selected dataset folder.", self.cmd_output
            )
            return

        os.makedirs(output_dir, exist_ok=True)
        self.generate_survey_conf(
            self.ray_surveyor_dir, input_files, kmer_length, output_dir
        )

        def run_ray_surveyor_command():
            ray_surveyor_command1 = f"cd {self.ray_surveyor_dir}"
            ray_surveyor_command2 = (
                f'wsl mpiexec -n 2 "{self.ray_executable}" survey1.conf'
            )

            try:
                self.update_cmd_output("Running Ray Surveyor...", self.cmd_output)

                process = Popen(
                    ray_surveyor_command1,
                    shell=True,
                    stdout=PIPE,
                    stderr=PIPE,
                    universal_newlines=True,
                )
                self.display_process_output(process, self.cmd_output)

                process = Popen(
                    ray_surveyor_command2,
                    shell=True,
                    cwd=self.ray_surveyor_dir,
                    stdout=PIPE,
                    stderr=PIPE,
                    universal_newlines=True,
                )
                self.display_process_output(process, self.cmd_output)

                self.update_cmd_output(
                    "Ray Surveyor completed successfully.\n check the output at:"
                    + output_dir,
                    self.cmd_output,
                )
            except CalledProcessError as e:
                self.update_cmd_output(
                    "Ray Surveyor encountered an error.", self.cmd_output, "error"
                )
                self.update_cmd_output(e.stdout, self.cmd_output, "error")
                self.update_cmd_output(e.stderr, self.cmd_output, "error")

        cmd_thread = threading.Thread(target=run_ray_surveyor_command)
        cmd_thread.start()

    def run_dsk(self):
        kmer_length = self.kmer_length.get()
        # dsk_dir = self.dsk_dir.get() # ERROR
        dataset_folder = self.dataset_folder.get()
        dataset_folder = dataset_folder.replace("\\", "/")
        dataset_folder = dataset_folder.replace("C:", "/mnt/c")
        output_dir = self.output_dir.get()

        # if not dsk_dir or not dataset_folder or not output_dir:
        # self.update_cmd_output("Please fill in all required fields for DSK.")
        # return

        list_reads_file = os.path.join(output_dir, "list_reads")
        list_reads_file = list_reads_file.replace("\\", "/")
        dsk_output_dir = os.path.join(output_dir, "dsk_output")
        dsk_output_dir = dsk_output_dir.replace("\\", "/")
        dsk_output_dir = dsk_output_dir.replace("C:", "/mnt/c")

        print(dsk_output_dir)

        try:
            # # Run the 'ls' command to list files in the dataset folder and save to list_reads
            # ls_command1 = f"cd { output_dir}"

            # run(ls_command1, shell=True, check=True,)
            ls_command = f" wsl ls -1 {dataset_folder}/* > list_reads_file"
            run(
                ls_command,
                shell=True,
                check=True,
            )

            # Run the 'dsk' command
            dsk_command = f"wsl bin/dsk/dsk -file list_reads_file -out-dir {dsk_output_dir} -kmer-size {kmer_length}"
            self.update_cmd_output("Running DSK...", self.cmd_output)

            process = Popen(
                dsk_command,
                shell=True,
                stdout=PIPE,
                stderr=PIPE,
                universal_newlines=True,
            )
            self.display_process_output(process, self.cmd_output)

            self.update_cmd_output(
                "DSK completed successfully. Output stored in:", dsk_output_dir
            )
        except CalledProcessError as e:
            self.update_cmd_output(
                "DSK encountered an error.", self.cmd_output, "error"
            )
            self.update_cmd_output(e.stdout, self.cmd_output, "error")
            self.update_cmd_output(e.stderr, self.cmd_output, "error")

    def koverlearn_error(self):
        if not self.pickdataset_entry.get():
            self.display_error_message("Please select kover dataset.")
            return False
        if not self.pickoutput_entry.get():
            self.display_error_message("Please select output directory path.")
            return False
        if not self.hyperparameterentry.get():
            self.display_error_message("Please enter hyperparameters")
            return False
        if not self.splitidentifierentry.get():
            self.display_error_message("Please enter split identifier")
            return False
        if not self.randomseedentry2.get():
            self.display_error_message("Please enter random seed")
            return False

        return True

    def koverlearn(self):
        if self.koverlearn_error():
            self.koverlearn_btn.configure(text="learning", state=tk.DISABLED)
            koverdataset_path1 = self.selected_kover.get()
            koverdataset_path1 = koverdataset_path1.replace("\\", "/")
            koverdataset_path1 = koverdataset_path1.replace("C:", "/mnt/c")
            output_path1 = self.output_dir.get()
            output_path1 = output_path1.replace("\\", "/")
            output_path1 = output_path1.replace("C:", "/mnt/c")
            hyperparameter = self.hyperparameterentry.get().split()
            model_type = self.model_type_var.get()
            splitid_identifier = self.splitidentifierentry.get()
            max_rules = self.maxrules_var.get()
            max_equiv_rules = self.maxequivrules_var.get()
            random_seed = self.randomseedentry2.get()
            hp_choice = self.hp_type_var.get()
            max_bound = self.maxboundsize_var.get()
            try:
                if self.learn_type_var.get() == "SCM":
                    command = KoverDatasetCreator.kover_learn_scm(
                        self,
                        koverdataset_path1,
                        splitid_identifier,
                        model_type,
                        hyperparameter,
                        max_rules,
                        max_equiv_rules,
                        hp_choice,
                        max_bound,
                        random_seed,
                        output_path1,
                    )

                else:
                    # Handle other dataset types if needed
                    command = KoverDatasetCreator.kover_learn_cart(
                        self,
                        koverdataset_path1,
                        splitid_identifier,
                        model_type,
                        hyperparameter,
                        max_rules,
                        max_equiv_rules,
                        hp_choice,
                        max_bound,
                        output_path1,
                    )

                process = Popen(
                    command,
                    shell=True,
                    stdout=PIPE,
                    stderr=PIPE,
                    universal_newlines=True,
                )
                self.display_process_output(process, self.cmd_output3)

            except Exception as e:
                # Handle any exceptions that occur during the dataset creation process
                print(f"An error occurred: {e}")

            finally:
                # Enable the button after the dataset creation process is complete
                self.koverlearn_btn.configure(text="learn", state=tk.NORMAL)

    def split_error_check(self):
        # Check if dataset entry is filled
        if not self.selected_kover.get():
            self.display_error_message("Please select kover dataset.")
            return False

        # Check if output entry is filled
        if not self.output_dir.get():
            self.display_error_message("Please select an output directory path.")
            return False

        # Check if trainids and testids are required and filled
        if self.trainvar.get() == "yes":
            if not self.trainids_entry.get():
                self.display_error_message("Please select a train IDs file.")
                return False

        if self.trainvar.get() == "yes":
            if not self.testids_entry.get():
                self.display_error_message("Please select a test IDs file.")
                return False

        # Check if random seed entry is filled
        if not self.randomseedentry.get():
            self.display_error_message("Please enter a random seed.")
            return False
        if not self.uniqueidentry.get():
            self.display_error_message("Please enter a unique splitid.")
            return False

        # Additional checks as needed...

        # If all checks passed, return True
        return True

    def display_error_message(self, message):
        # Use messagebox.showerror to display an error message
        messagebox.showerror("Error", message)

    def split_dataset(self):
        if self.split_error_check():
            # Get values from GUI variables and save them in respective variables
            self.split_btn.configure(text="splitting", state=tk.DISABLED)
            koverdataset_path = self.selected_kover.get()
            koverdataset_path = koverdataset_path.replace("\\", "/")
            koverdataset_path = koverdataset_path.replace("C:", "/mnt/c")
            output_path = self.output_dir.get()
            output_path = output_path.replace("\\", "/")
            output_path = output_path.replace("C:", "/mnt/c")
            train_size = float(self.trainsize_spinbox.get()) / 100.0
            use_train_ids = self.trainvar.get()
            folds = self.foldvar.get()
            random_seed = (
                int(self.randomseedentry.get()) if self.randomseedentry.get() else None
            )
            train_ids_path = (
                self.trainids_entry.get() if use_train_ids == "yes" else None
            )
            test_ids_path = self.testids_entry.get() if use_train_ids == "yes" else None
            unique_split_id = self.uniqueidentry.get()
            splitoutput = os.path.join(output_path, "split_DATASET.kover")
            if use_train_ids == "no":
                command = KoverDatasetCreator.split_dataset(
                    self,
                    koverdataset_path,
                    unique_split_id,
                    train_size,
                    folds,
                    random_seed,
                )
                process = Popen(
                    command,
                    shell=True,
                    stdout=PIPE,
                    stderr=PIPE,
                    universal_newlines=True,
                )
                self.display_process_output(process, self.cmd_output2)
                self.split_btn.configure(text="split dataset", state=tk.NORMAL)

    def errorcheck_create_dataset(self):
        # Check if necessary fields are filled correctly
        self.dataset_path = self.pickdataset1_entry.get()
        self.output_path = self.pickoutput_entry.get()
        self.kmer_length_str = self.kmer_length_var.get()
        self.compression_str = self.compression_var.get()
        self.phenotype_desc_path = self.phenotype_desc_entry.get()
        self.phenotype_metadata_path = self.phenotype_metadata_entry.get()

        try:
            # Check if dataset path is provided
            if not self.dataset_path:
                messagebox.showerror("Error", "Please provide a dataset path.")

            # Check if output path is provided
            if not self.output_path:
                messagebox.showerror(
                    "Error", "Please provide an output directory path."
                )

            # Check if kmer length is an integer and within the allowed range
            kmer_length = int(self.kmer_length_str)
            if not 0 <= kmer_length <= 128:
                messagebox.showerror(
                    "Error", "Kmer length should be between 0 and 128."
                )

            # Check if compression is an integer and within the allowed range
            self.compression = int(self.compression_str)
            if not 0 <= self.compression <= 9:
                messagebox.showerror(
                    "Error", "Compression level should be between 0 and 9."
                )

            # Check if phenotype description path is provided
            if not self.phenotype_desc_path:
                messagebox.showerror(
                    "Error", "Please provide a phenotype description file."
                )

            # Check if phenotype metadata path is provided
            if not self.phenotype_metadata_path:
                messagebox.showerror(
                    "Error", "Please provide a phenotype metadata file."
                )

            # If all checks pass, return True
            return True

        except ValueError as e:
            # Display an error message
            # print(f"Error: {e}")
            # Return False if any check fails
            return False

    def create_dataset_thread(self):
        thread = threading.Thread(target=self.create_dataset)
        thread.start()

    def split_dataset_thread(self):
        thread = threading.Thread(target=self.split_dataset)
        thread.start()

    def kover_learn_thread(self):
        self.cmd_output3["state"] = tk.NORMAL
        self.cmd_output3.delete("1.0", tk.END)
        thread = threading.Thread(target=self.koverlearn)
        thread.start()

    def create_dataset(self):
        if self.errorcheck_create_dataset():
            self.create_dataset_btn.configure(text="creating", state=tk.DISABLED)
            output = os.path.join(self.output_path, "DATASET.kover")
            output = output.replace("\\", "/")
            output = output.replace("C:", "/mnt/c")
            self.phenotype_desc_path = self.phenotype_desc_path.replace("\\", "/")
            self.phenotype_desc_path = self.phenotype_desc_path.replace("C:", "/mnt/c")
            self.phenotype_metadata_path = self.phenotype_metadata_path.replace(
                "\\", "/"
            )
            self.phenotype_metadata_path = self.phenotype_metadata_path.replace(
                "C:", "/mnt/c"
            )
            try:
                if self.dataset_type_var1.get() == "contigs":
                    dataset = KoverDatasetCreator.contigs_parser(
                        self, self.dataset_path, "out.tsv"
                    )
                    command = KoverDatasetCreator.create_from_contigs(
                        self,
                        "out.tsv",
                        self.phenotype_desc_path,
                        self.phenotype_metadata_path,
                        output,
                        self.kmer_length_var.get(),
                        self.compression,
                    )
                elif self.dataset_type_var1.get() == "kmer matrix":
                    dataset = self.selected_kmer_matrix.get()
                    dataset = dataset.replace("\\", "/")
                    dataset = dataset.replace("C:", "/mnt/c")
                    command = KoverDatasetCreator.create_from_tsv(
                        self,
                        dataset,
                        self.phenotype_desc_path,
                        self.phenotype_metadata_path,
                        output,
                    )
                else:
                    # Handle other dataset types if needed
                    pass

                process = Popen(
                    command,
                    shell=True,
                    stdout=PIPE,
                    stderr=PIPE,
                    universal_newlines=True,
                )
                self.display_process_output(process, self.cmd_output1)
                # self.display_process_output("Done!", self.cmd_output1)
                self.create_dataset_btn.configure(
                    text="create Dataset", state=tk.NORMAL
                )
            except Exception:
                # Handle any exceptions that occur during the dataset creation process
                traceback.print_exc()

            finally:
                # Enable the button after the dataset creation process is complete
                self.create_dataset_btn.configure(
                    text="create Dataset", state=tk.NORMAL
                )

    def generate_random_seed(self, randomseedentry):
        # Generate a random seed and place it in randomseedentry
        random_seed = random.randint(1, 10000)
        randomseedentry.delete(0, tk.END)  # Clear any existing value
        randomseedentry.insert(0, str(random_seed))  # Place the generated random seed

    def pick_metatsv_file(self):
        # Open a file dialog for selecting TSV files
        tsv_file_path = filedialog.askopenfilename(
            filetypes=[("TSV Files", "*.tsv")], title="Select a TSV File"
        )
        self.metatsv_file_path.set(tsv_file_path)

    def pick_desctsv_file(self):
        # Open a file dialog for selecting TSV files
        tsv_file_path = filedialog.askopenfilename(
            filetypes=[("TSV Files", "*.tsv")], title="Select a TSV File"
        )
        self.desctsv_file_path.set(tsv_file_path)

    def toggle_train_test_widgets(self, event):
        selected_option = self.trainvar.get()
        if selected_option == "yes":
            self.trainids_btn.place(x=50, y=450)
            self.trainids_entry.place(x=300, y=450)
            self.testids_btn.place(x=50, y=500)
            self.testids_entry.place(x=300, y=500)
            self.uniqueidlabel.place(x=50, y=550)
            self.uniqueidentry.place(x=180, y=550)
        if selected_option == "no":
            self.trainids_btn.place_forget()
            self.trainids_entry.place_forget()
            self.testids_btn.place_forget()
            self.testids_entry.place_forget()
            self.uniqueidlabel.place(x=50, y=450)
            self.uniqueidentry.place(x=180, y=450)

    def update_button_text(self, *args):
        self.runbtn["text"] = "Run " + self.kmer_tool.get()

    def browse_dataset(self):
        dataset_folder = filedialog.askdirectory()
        self.dataset_folder.set(dataset_folder)

    def browse_output_dir(self):
        output_dir = filedialog.askdirectory()
        self.output_dir.set(output_dir)

    @threaded
    def load_amr_data(self):
        self.download_button4.configure(text="Loading..", state=tk.DISABLED)
        self.amr_list_filter_checkbox.configure(state=tk.DISABLED)
        self.species_filter.configure(state=tk.DISABLED)
        self.species_selection.configure(state=tk.DISABLED)
        self.antibiotic_selection.configure(state=tk.DISABLED)
        self.drop_intermediate_checkbox.configure(state=tk.DISABLED)
        self.numeric_phenotypes_checkbox.configure(state=tk.DISABLED)
        self.save_table_button.configure(state=tk.DISABLED)
        amr_metadata_file = filedialog.askopenfilename(
            filetypes=[("AMR Text Files", "*.txt")]
        )
        if amr_metadata_file:
            self.total_label.config(text="Total: ...")

            self.amr_full = pd.read_table(
                amr_metadata_file,
                usecols=[
                    "genome_id",
                    "genome_name",
                    "antibiotic",
                    "resistant_phenotype",
                    "measurement",
                    "measurement_unit",
                ],
                converters={
                    "genome_id": float,
                    "genome_name": lambda x: " ".join(x.lower().split()[:2])
                    .replace("[", "")
                    .replace("]", ""),
                    "antibiotic": str,
                    "resistant_phenotype": str,
                    "measurement": str,
                    "measurement_unit": str,
                },
            )

            self.amr_full.drop_duplicates(inplace=True)

            self.amr_full = self.amr_full[
                (self.amr_full["genome_id"] != "")
                & (self.amr_full["genome_name"] != "")
                & (self.amr_full["antibiotic"] != "")
                & (self.amr_full["resistant_phenotype"] != "")
                & (self.amr_full["measurement"] != "")
                & (self.amr_full["measurement_unit"] != "")
            ]

            self.amr_full = self.amr_full[self.amr_full["measurement_unit"] != "mm"]

            self.amr_full["measurement"] += self.amr_full["measurement_unit"]

            self.amr_list_unfiltered = self.amr_full[
                ["genome_name", "antibiotic"]
            ].drop_duplicates()

            self.amr_list_filtered = (
                self.amr_full.groupby(["genome_name", "antibiotic"])
                .filter(
                    lambda x: (x["resistant_phenotype"] == "Resistant").sum() >= 25
                    and (x["resistant_phenotype"] == "Susceptible").sum() >= 25
                )[["genome_name", "antibiotic"]]
                .drop_duplicates()
            )

            species_list = self.amr_list_unfiltered["genome_name"].unique()
            species_list.sort()
            species_list = ["All"] + species_list.tolist()

            antibiotic_list = self.amr_list_unfiltered["antibiotic"].unique()
            antibiotic_list.sort()
            antibiotic_list = ["All"] + antibiotic_list.tolist()

            self.species_filter["values"] = species_list
            self.species_filter.current(0)

            self.species_selection["values"] = species_list
            self.species_selection.current(0)

            self.antibiotic_selection["values"] = antibiotic_list
            self.antibiotic_selection.current(0)

            if self.amr_list_filter_checkbox.get():
                self.amr_list = self.amr_list_filtered
            else:
                self.amr_list = self.amr_list_unfiltered

            self.update_table()
            self.update_amr_full()

        if hasattr(self, "amr_full"):
            self.amr_list_filter_checkbox.configure(state=tk.NORMAL)
            self.species_filter.configure(state=tk.NORMAL)
            self.species_selection.configure(state="readonly")
            self.antibiotic_selection.configure(state="readonly")
            self.drop_intermediate_checkbox.configure(state=tk.NORMAL)
            self.numeric_phenotypes_checkbox.configure(state=tk.NORMAL)
            self.save_table_button.configure(state=tk.NORMAL)

        self.download_button4.configure(text="load amr list", state=tk.NORMAL)

    def update_table(self, event=None):
        selected_species = self.species_filter.get()
        self.species_filter.set(selected_species)
        self.amr_list_table.load_start = 0

        if len(self.amr_list) > 0:
            if not selected_species or selected_species == "All":
                self.amr_list_table.filtered_data = self.amr_list
            else:
                self.amr_list_table.filtered_data = self.amr_list[
                    self.amr_list["genome_name"].str.contains(
                        pat=selected_species, case=False, regex=False
                    )
                    | self.amr_list["antibiotic"].str.contains(
                        pat=selected_species, case=False, regex=False
                    )
                ]

            self.amr_list_table.reset_view()

            self.total_label.config(
                text=f"Total: {len(self.amr_list_table.filtered_data)}"
            )

    def on_amr_list_filter_check(self, event=None):
        self.amr_list_filter_checkbox.configure(state=tk.DISABLED)
        self.species_filter.configure(state=tk.DISABLED)

        if self.amr_list_filter_checkbox.get():
            self.amr_list = self.amr_list_filtered
        else:
            self.amr_list = self.amr_list_unfiltered

        species_list = self.amr_list["genome_name"].unique()
        species_list.sort()
        species_list = ["All"] + species_list.tolist()

        self.species_filter["values"] = species_list
        self.species_filter.current(0)

        self.update_table()

        self.amr_list_filter_checkbox.configure(state=tk.NORMAL)
        self.species_filter.configure(state=tk.NORMAL)

    def on_antibiotic_select(self, event=None):
        antibiotic = self.antibiotic_selection.get()

        if antibiotic == "All":
            filtered = self.amr_list["genome_name"].unique()
            filtered.sort()
            self.species_selection["values"] = ["All"] + filtered.tolist()
        else:
            filtered = self.amr_list[self.amr_list["antibiotic"] == antibiotic][
                "genome_name"
            ].unique()
            filtered.sort()
            self.species_selection["values"] = ["All"] + filtered.tolist()

        self.species_selection.current(0)

        self.update_amr_full()

    def on_table_select(self, event=None):
        selected_item = self.amr_list_table.focus()

        if not selected_item:
            return

        selected_species = self.amr_list_table.item(selected_item, "values")[0]
        selected_antibiotics = self.amr_list_table.item(selected_item, "values")[1]

        self.species_selection.set(selected_species)
        self.antibiotic_selection.set(selected_antibiotics)

        self.update_amr_full()

    def phenotype_mask(self, data_frame: pd.DataFrame):
        numeric_phenotypes = data_frame.copy()

        numeric_phenotypes.loc[
            ~numeric_phenotypes["resistant_phenotype"].isin(
                ["Resistant", "Susceptible"]
            ),
            "resistant_phenotype",
        ] = 2

        numeric_phenotypes.loc[
            numeric_phenotypes["resistant_phenotype"] == "Susceptible",
            "resistant_phenotype",
        ] = 0

        numeric_phenotypes.loc[
            numeric_phenotypes["resistant_phenotype"] == "Resistant",
            "resistant_phenotype",
        ] = 1

        return numeric_phenotypes

    @threaded
    def update_amr_full(self, event=None):
        antibiotic = self.antibiotic_selection.get()
        species = self.species_selection.get()

        self.antibiotic_selection.selection_clear()
        self.species_selection.selection_clear()

        if len(self.amr_full) > 0:
            self.totalresistance_label.configure(
                text="Total resistance phenotypes available: ..."
            )
            self.totalsusceptible_label.configure(
                text="Total susceptible phenotypes available: ..."
            )
            self.totaldata.configure(text="Total Phenotypes: ...")
            self.save_table_button.configure(state=tk.DISABLED)
            self.species_selection.configure(state=tk.DISABLED)
            self.antibiotic_selection.configure(state=tk.DISABLED)

            try:
                self.results_table.filtered_data = self.amr_full

                if antibiotic == "All" and species == "All":
                    pass
                elif antibiotic == "All":
                    self.results_table.filtered_data = self.amr_full[
                        self.amr_full["genome_name"] == species
                    ]
                elif species == "All":
                    self.results_table.filtered_data = self.amr_full[
                        self.amr_full["antibiotic"] == antibiotic
                    ]
                else:
                    self.results_table.filtered_data = self.amr_full[
                        (self.amr_full["antibiotic"] == antibiotic)
                        & (self.amr_full["genome_name"] == species)
                    ]

                self.results_table.filtered_data = self.results_table.filtered_data[
                    ["genome_id", "genome_name", "resistant_phenotype", "measurement"]
                ]

                if self.drop_intermediate_checkbox.get():
                    self.results_table.filtered_data = self.results_table.filtered_data[
                        (
                            self.results_table.filtered_data["resistant_phenotype"]
                            == "Resistant"
                        )
                        | (
                            self.results_table.filtered_data["resistant_phenotype"]
                            == "Susceptible"
                        )
                    ]

                total_resistant = len(
                    self.results_table.filtered_data[
                        self.results_table.filtered_data["resistant_phenotype"]
                        == "Resistant"
                    ]
                )

                total_susceptible = len(
                    self.results_table.filtered_data[
                        self.results_table.filtered_data["resistant_phenotype"]
                        == "Susceptible"
                    ]
                )

                if self.numeric_phenotypes_checkbox.get():
                    self.results_table.filtered_data = self.phenotype_mask(
                        self.results_table.filtered_data
                    )

                total = len(self.results_table.filtered_data)

                self.totalresistance_label.configure(
                    text=f"Total resistance phenotypes available: {total_resistant}"
                )
                self.totalsusceptible_label.configure(
                    text=f"Total susceptible phenotypes available: {total_susceptible}"
                )

                self.results_table.reset_view()

                self.totaldata.configure(text=f"Total Phenotypes: {total}")
            except Exception:
                traceback.print_exc()
                messagebox.showerror("Error", "Error while reading the metadata file")
            finally:
                self.save_table_button.configure(state=tk.NORMAL)
                self.species_selection.configure(state="readonly")
                self.antibiotic_selection.configure(state="readonly")

    def save_description_to_tsv(self, species: str, antibiotics: str, file_path: str):
        with open(file_path, "w", newline="") as tsv_file:
            tsv_writer = csv.writer(tsv_file)
            tsv_writer.writerow([f"Species: {species}"])
            tsv_writer.writerow([f"Antibiotics: {antibiotics}"])

    def save_all_table_to_tsv(self, table: Table, file_path: str):
        with open(file_path, "w", newline="") as tsv_file:
            columns = table.filtered_data.columns.tolist()

            tsv_writer = csv.writer(tsv_file, delimiter="\t")

            tsv_writer.writerow(columns)

            tsv_writer.writerows(table.filtered_data.values)

    def save_phenotype_metadata_to_tsv(self, table: Table, file_path: str):
        with open(file_path, "w", newline="") as tsv_file:
            tsv_writer = csv.writer(tsv_file, delimiter="\t")
            tsv_writer.writerows(
                table.filtered_data.iloc[:, [0, 2]].drop_duplicates().values
            )

    def save_id_name_to_tsv(self, table: Table, file_path: str):
        with open(file_path, "w", newline="") as tsv_file:
            tsv_writer = csv.writer(tsv_file, delimiter="\t")
            columns = table.filtered_data.columns[:2].tolist()
            tsv_writer.writerow(columns)
            tsv_writer.writerows(table.filtered_data.iloc[:, [0, 1]].values)

    @threaded
    def save_to_tsv(self):
        species_text = self.species_selection.get().replace(" ", "-").replace("/", "-")
        antibiotic_text = (
            self.antibiotic_selection.get().replace(" ", "-").replace("/", "-")
        )

        if species_text == "" or antibiotic_text == "":
            messagebox.showerror("Error", "Please load amr data first.")
            return

        self.save_table_button.configure(state=tk.DISABLED, text="Saving...")

        selected_directory = filedialog.askdirectory()

        if not selected_directory:
            self.save_table_button.configure(state=tk.NORMAL, text="Export to .tsv")
            return

        file_name = f"{species_text}_{antibiotic_text}.tsv"
        species = species_text
        antibiotics = antibiotic_text

        tsv_folder_path = os.path.join(selected_directory, Path.DATA, species_text)
        tsv_folder_path = os.path.join(tsv_folder_path, antibiotic_text)

        os.makedirs(tsv_folder_path, exist_ok=True)

        base_name, extension = os.path.splitext(file_name)

        phenometafile = f"{base_name}_phenotype_metadata{extension}"
        file_path = os.path.join(tsv_folder_path, phenometafile)
        self.save_phenotype_metadata_to_tsv(self.results_table, file_path)

        descfilename = f"{base_name}_description{extension}"
        file_path = os.path.join(tsv_folder_path, descfilename)
        self.save_description_to_tsv(species, antibiotics, file_path)

        allfilename = f"{base_name}_full{extension}"
        file_path = os.path.join(tsv_folder_path, allfilename)
        self.save_all_table_to_tsv(self.results_table, file_path)

        id_name_filename = f"{base_name}_id_name{extension}"
        file_path = os.path.join(tsv_folder_path, id_name_filename)
        self.save_id_name_to_tsv(self.results_table, file_path)

        self.save_table_button.configure(state=tk.NORMAL, text="Export to .tsv")

    def select_directory(self):
        # Use the filedialog to select a directory and set it in the FTPDownloadApp
        self.selected_directory = filedialog.askdirectory()
        return self.selected_directory
        # self.download_app.selected_directory = self.selected_directory

    def generate_survey_conf(
        self, ray_surveyor_dir, input_files, kmer_length, output_dir
    ):
        survey_conf_path = os.path.join(ray_surveyor_dir, "survey1.conf")
        with open(survey_conf_path, "w") as survey_conf_file:
            survey_conf_file.write(f"-k {kmer_length}\n")
            survey_conf_file.write("-run-surveyor\n")
            survey_conf_file.write(f"-output {output_dir}\n")
            survey_conf_file.write("-write-kmer-matrix\n")
            survey_conf_file.write("\n")

            for input_file in input_files:
                file_name = os.path.splitext(os.path.basename(input_file))[0]
                input_file = input_file.replace("\\", "/")
                input_file = input_file.replace("C:", "/mnt/c")
                survey_conf_file.write(
                    f"-read-sample-assembly {file_name} {input_file}\n"
                )

    def display_process_output(self, process: Popen, output_target=None):
        for line in process.stdout:
            self.update_cmd_output(line, output_target)
        for line in process.stderr:
            self.update_cmd_output(line, output_target, "error")

    def update_cmd_output(self, message: str, output_target: tk.Text, tag: str = None):
        output_target["state"] = tk.NORMAL

        if tag:
            output_target.insert(tk.END, f"{tag}: {message}")
        else:
            output_target.insert(tk.END, message)

        output_target.see(tk.END)
        output_target["state"] = tk.DISABLED
        output_target.update_idletasks()

    def pickkover(self):
        file_path = filedialog.askopenfilename(
            title="Select a .kover file", filetypes=[("Kover Files", "*.kover")]
        )
        if file_path:
            # Do something with the selected file path (e.g., store it in a variable)
            self.selected_kover.set(file_path)

    def pickkmer_matrix(self):
        file_path = filedialog.askopenfilename(
            title="Select a kmer matrix file", filetypes=[("Kover Files", "*.tsv")]
        )
        if file_path:
            # Do something with the selected file path (e.g., store it in a variable)
            self.selected_kmer_matrix.set(file_path)

    def on_learn_type_selected(self, event):
        self.update_kover_learn_gui()

    def update_kover_learn_gui(self):
        selected_type = self.learn_type_var.get()
        if selected_type == "CART":
            self.hyperparameterlabel.configure(text="Class importance:")
            self.modellabel.configure(text="Criterion:")
            self.model_type_var.set("gini")  # Set the default value to "gini"
            self.model_type_options = ["gini", "crossentropy"]
            self.model_type_menu["values"] = self.model_type_options
            self.maxrulelabel.configure(text="max depth")
            self.maxrules_var.set(10)
            self.hp_type_var.set("cv")  # Set a default value
            self.hp_type_options = ["bound", "cv"]
            self.hp_type_menu["values"] = self.hp_type_options
            self.maxequivrulelabel.configure(text="min samples split")
            self.maxequivrules_var.set(2)
        else:
            self.hyperparameterlabel.configure(text="Hyper-parameter:")
            self.modellabel.configure(text="Model type:")
            self.model_type_var.set("conjunction")  # Set the default value to "gini"
            self.model_type_options = ["conjunction", "disjunction"]
            self.model_type_menu["values"] = self.model_type_options
            self.maxrulelabel.configure(text="Maximum rules")
            self.maxrules_var.set(10)
            self.hp_type_var.set("none")  # Set a default value
            self.hp_type_options = ["bound", "cv", "none"]
            self.hp_type_menu["values"] = self.hp_type_options
            self.maxequivrulelabel.configure(text="Maximum equivalent rules :")
            self.maxequivrules_var.set(10)

    def on_dataset_type_selected(self, event):
        self.update_pickdataset_button()

    def update_pickdataset_button(self):
        selected_type = self.dataset_type_var1.get()
        if selected_type == "kmer matrix":
            self.pickdataset1_btn.configure(
                text="Pick Kmer Matrix file", command=self.pickkmer_matrix
            )
            self.pickdataset1_entry.configure(textvariable=self.selected_kmer_matrix)
        elif selected_type == "contigs":
            self.pickdataset1_btn.configure(
                text="Pick Contigs Folder", command=self.browse_dataset
            )
            self.pickdataset1_entry.configure(textvariable=self.dataset_folder)
        elif selected_type == "reads":
            self.pickdataset1_btn.configure(
                text="Pick Reads Folder", command=self.browse_dataset
            )
            self.pickdataset1_entry.configure(textvariable=self.dataset_folder)

    def toggle_max_bound_widgets(self, event):
        selected_type = self.hp_type_var.get()

        if selected_type == "bound":
            self.maxboundsize_spinbox.place(x=450, y=450)
            self.maxboundsize.place(x=280, y=450)
        elif selected_type == "cv":
            pass
        else:
            self.maxboundsize_spinbox.place_forget()
            self.maxboundsize.place_forget()

    def to_remove(self):
        self.bind_all("<ButtonPress>", lambda event: print("debug:", event.widget))
        self.dataset_folder = tk.StringVar()
        self.output_dir = tk.StringVar()
        self.desctsv_file_path = tk.StringVar()
        self.metatsv_file_path = tk.StringVar()
        self.selected_kover = tk.StringVar()
