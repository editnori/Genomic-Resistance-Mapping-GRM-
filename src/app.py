from ftplib import FTP
import time
import tkinter as tk
from tkinter import ttk, messagebox, font

import os
import pathlib as pl
from subprocess import CalledProcessError, Popen
import csv
import random
from enum import Enum
import traceback
import re

import ctk
from PIL import Image
import pandas as pd

import util
from util import threaded, Tag
from spinbox import Spinbox
from combobox import Combobox
from controlframe import ControlFrame
import kover
from hovertip import Hovertip
from table import Table
from label import Label
from ftp_downloader import FTPDownloadApp, DownloadWindow, get_last_metadata_update_date

from tkwebview2.tkwebview2 import WebView2, have_runtime, install_runtime

from pythonnet import load

load("coreclr")

from clr import AddReference

AddReference("System.Threading")

from System.Threading import Thread, ApartmentState, ThreadStart

class Page(Enum):
    DATA_COLLECTION_PAGE = 0
    PREPROCESSING_PAGE = 1
    KOVER_LEARN_PAGE = 2
    ANALYSIS_PAGE = 3


class Path(str):
    ROOT = pl.Path(__file__).parent.parent
    REMOTE_METADATA = "RELEASE_NOTES/genome_metadata"
    FOREST_DARK = os.path.join(ROOT, "ui/forest-dark.tcl")
    RAY = os.path.join(ROOT, "bin/ray/Ray")
    DSK = os.path.join(ROOT, "bin/dsk/dsk")
    KOVER = os.path.join(ROOT, "bin/kover/kover")
    IMAGES = os.path.join(ROOT, "ui/test_images/")
    TEMP = os.path.join(ROOT, ".temp/")
    DATA = os.path.join(ROOT, "data/")
    CONTIGS = "contigs/"
    FEATURES = "features/"


class App(ctk.CTk):
    def __init__(self):
        super().__init__()

        os.makedirs(Path.DATA, exist_ok=True)

        self.setup_window()

        self.setup_style()

        self.load_images(Path.IMAGES)

        self.create_navigation_frame()

        self.create_data_collection_page()

        self.create_preprocessing_page()

        self.create_kover_learn_page()

        self.create_analysis_page()

        self.set_page(Page.DATA_COLLECTION_PAGE)

    def setup_window(self):
        if not have_runtime():
            install_runtime()

        self.call("encoding", "system", "utf-8")

        self.title("Genome analysis tool")

        self.screen_width = self.winfo_screenwidth()
        self.screen_height = self.winfo_screenheight()

        self.geometry("1200x800")

        self.protocol("WM_DELETE_WINDOW", self.end_app)

        self.minsize(1280, 720)

        self.after(0, lambda: self.state("zoomed"))

        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)

        def on_click(event=None):
            try:
                event.widget.focus_set()
            except Exception:
                pass

        self.bind_all("<Button-1>", on_click)

    def end_app(self):
        if messagebox.askokcancel("Quit", "Do you want to quit?"):
            self.destroy()

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

        kover_image = Image.open(path + "kover.png")

        self.images["kover"] = ctk.CTkImage(
            light_image=kover_image,
            dark_image=kover_image,
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

        self.data_collection_frame_button.grid(row=1, column=0, sticky=tk.EW)

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

        self.preprocessing_frame_button.grid(row=2, column=0, sticky=tk.EW)

        self.kover_frame_button = ctk.CTkButton(
            self.navigation_frame,
            corner_radius=0,
            height=button_height,
            text="Kover learn",
            border_spacing=10,
            text_color=button_text_color,
            hover_color=button_hover_color,
            image=self.images["kover"],
            anchor="w",
            command=lambda: self.set_page(Page.KOVER_LEARN_PAGE),
        )

        self.kover_frame_button.grid(row=3, column=0, sticky=tk.EW)

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

        self.analysis_frame_button.grid(row=4, column=0, sticky=tk.EW)

        self.appearance_mode_menu = ctk.CTkOptionMenu(
            self.navigation_frame,
            values=["Dark", "Light", "System"],
            command=ctk.set_appearance_mode,
        )

        self.appearance_mode_menu.grid(row=6, column=0, pady=15, sticky="s")

        ctk.set_appearance_mode("Dark")

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

        self.total_progress_bar = ttk.Progressbar(
            master=self.genome_data_frame, mode="determinate"
        )

        self.total_progress_label = ctk.CTkLabel(
            master=self.genome_data_frame,
            font=self.default_font(10),
            fg_color="transparent",
            text="",
            anchor="w",
        )

        self.file_progress_bar = ttk.Progressbar(
            master=self.genome_data_frame, mode="determinate"
        )

        self.file_progress_label = ctk.CTkLabel(
            master=self.genome_data_frame,
            font=self.default_font(10),
            fg_color="transparent",
            text="",
            anchor="w",
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
            self.genome_data_frame_path.config(text="")

        self.genome_data_validate_ui()

    def select_genome_data_directory(self):
        self.genome_data_frame_path.configure(
            text=util.select_file(filetypes=[("TSV file", "*.tsv")])
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

        genome_names = [
            util.sanitize_filename(genome_name) for _, genome_name in genome_data_ids
        ]

        directory = util.select_directory()

        if not directory:
            return

        contigs_directory = os.path.join(directory, Path.CONTIGS)
        features_directory = os.path.join(directory, Path.FEATURES)

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
        self.total_progress_bar.pack(padx=50, fill=tk.X)
        self.total_progress_label.pack(padx=50, pady=(0, 20), fill=tk.X)
        self.file_progress_bar.pack(padx=50, fill=tk.X)
        self.file_progress_label.pack(padx=50, pady=(0, 20), fill=tk.X)

        self.total_progress_bar.configure(maximum=total_genomes)
        self.total_progress_bar["value"] = 0

        for genome_data_id, genome_name in genome_data_ids:
            genome_name = util.sanitize_filename(genome_name)
            contig_name = f"{genome_data_id}.fna"
            feature_name = f"{genome_data_id}.PATRIC.features.tab"

            local_contig_directory = os.path.join(contigs_directory, genome_name)
            local_feature_directory = os.path.join(features_directory, genome_name)

            local_contig_path = os.path.join(local_contig_directory, contig_name)
            local_feature_path = os.path.join(local_feature_directory, feature_name)

            remote_contig = f"genomes/{genome_data_id}/{genome_data_id}.fna"
            remote_feature = (
                f"genomes/{genome_data_id}/{genome_data_id}.PATRIC.features.tab"
            )

            if self.genome_data_frame_contig_checkbox.get():
                os.makedirs(local_contig_directory, exist_ok=True)

                number_downloaded += 1
                self.total_progress_label.configure(
                    text=f"Downloading: {contig_name} ({number_downloaded}/{total_genomes})"
                )
                self.total_progress_bar["value"] = number_downloaded
                try:
                    with open(local_contig_path, "wb") as local_file:
                        ftp = FTP("ftp.bvbrc.org")
                        ftp.login()
                        ftp.voidcmd("TYPE I")

                        contig_size = ftp.size(remote_contig)

                        self.file_progress_bar.configure(maximum=contig_size)

                        contig_size_mb = contig_size / 1_048_576

                        self.file_progress_bar["value"] = 0

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
                                    self.file_progress_bar["value"] = bytes_received
                                    self.file_progress_label.configure(
                                        text=f"Downloaded: {bytes_received / 1_048_576:6.2f} MB / {contig_size_mb:6.2f} MB"
                                    )
                                    last_update_time = current_time_ms

                        self.file_progress_bar["value"] = 0
                        self.file_progress_label.configure(text="")
                except Exception as e:
                    messagebox.showerror(
                        "Error", f"An error occurred while downloading contigs\n\n{e}"
                    )
                    try:
                        os.remove(local_contig_path)
                        if len(os.listdir(local_contig_directory)) == 0:
                            os.rmdir(local_contig_directory)
                    except Exception:
                        pass

                if self.cancel_genome_data_download_boolean:
                    try:
                        os.remove(local_contig_path)
                        if len(os.listdir(local_contig_directory)) == 0:
                            os.rmdir(local_contig_directory)
                    except Exception:
                        pass
                    break

            if self.genome_data_frame_feature_checkbox.get():
                os.makedirs(local_feature_directory, exist_ok=True)

                number_downloaded += 1
                self.total_progress_label.configure(
                    text=f"Downloading: {feature_name} ({number_downloaded}/{total_genomes})"
                )
                self.total_progress_bar["value"] = number_downloaded
                try:
                    with open(local_feature_path, "wb") as local_file:
                        ftp = FTP("ftp.bvbrc.org")
                        ftp.login()
                        ftp.voidcmd("TYPE I")

                        feature_size = ftp.size(remote_feature)

                        self.file_progress_bar.configure(maximum=feature_size)

                        feature_size_mb = feature_size / 1_048_576

                        self.file_progress_bar["value"] = 0

                        with ftp.transfercmd(f"RETR {remote_feature}") as conn:
                            bytes_received = 0
                            last_update_time = time.time_ns() / 1_000_000
                            while not self.cancel_genome_data_download_boolean and (
                                data := conn.recv(1024)
                            ):
                                local_file.write(data)
                                bytes_received += len(data)
                                if (current_time_ms - last_update_time) > 100:
                                    self.file_progress_bar["value"] = bytes_received
                                    self.file_progress_label.configure(
                                        text=f"Downloaded: {bytes_received / 1_048_576:6.2f} MB / {feature_size_mb:6.2f} MB"
                                    )
                                    last_update_time = current_time_ms

                        self.file_progress_bar["value"] = 0
                        self.file_progress_label.configure(text="")
                except Exception as e:
                    messagebox.showerror(
                        "Error", f"An error occurred while downloading features\n\n{e}"
                    )
                    try:
                        os.remove(local_feature_path)
                        if len(os.listdir(local_feature_directory)) == 0:
                            os.rmdir(local_feature_directory)
                    except Exception:
                        pass

                if self.cancel_genome_data_download_boolean:
                    try:
                        os.remove(local_feature_path)
                        if len(os.listdir(local_feature_directory)) == 0:
                            os.rmdir(local_feature_directory)
                    except Exception:
                        pass
                    break

        [
            kover.create_contigs_path_tsv(contigs_directory, genome_name)
            for genome_name in genome_names
        ]

        self.genome_data_frame_contig_checkbox.configure(state=tk.NORMAL)
        self.genome_data_frame_feature_checkbox.configure(state=tk.NORMAL)
        self.genome_data_frame_bulk_checkbox.configure(state=tk.NORMAL)
        self.genome_data_frame_entry.configure(state=tk.NORMAL)
        self.genome_data_frame_download_button.configure(
            text="Download", command=self.download_genome_data
        )

        self.total_progress_bar["value"] = 0
        self.total_progress_label.configure(text="")

        self.total_progress_bar.pack_forget()
        self.total_progress_label.pack_forget()
        self.file_progress_bar.pack_forget()
        self.file_progress_label.pack_forget()

    def cancel_genome_data_download(self):
        if messagebox.askyesno(
            "Confirmation", "Are you sure you want to cancel the download?"
        ):
            self.cancel_genome_data_download_boolean = True

    @threaded
    def download_genome_metadata(self):
        directory = util.select_directory()

        if not directory:
            return

        self.cancel_genome_metadata_download_boolean = False

        self.genome_metadata_frame_download_button.configure(
            text="Cancel", command=self.cancel_genome_metadata_download
        )

        path = os.path.join(directory, "genome_metadata")
        try:
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
        except Exception as e:
            messagebox.showerror(
                "Error", f"An error occurred while downloading the metadata\n\n{e}"
            )
            try:
                os.remove(path)
            except Exception:
                pass

        if self.cancel_genome_metadata_download_boolean:
            try:
                os.remove(path)
            except Exception:
                pass

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
        self.data_collection_tab_view.tab("AMR").grid_rowconfigure(
            tuple(range(10)), weight=1, uniform="row"
        )
        self.data_collection_tab_view.tab("AMR").grid_columnconfigure(
            tuple(range(10)), weight=1, uniform="col"
        )
        self.data_collection_main_frame = ctk.CTkScrollableFrame(
            self.data_collection_tab_view.tab("AMR"),
            width=1000,
            height=600,
        )

        self.data_collection_main_frame.pack(fill=tk.BOTH, expand=True)

        frame4 = ctk.CTkFrame(
            self.data_collection_main_frame,
            width=1000,
            height=400,
            corner_radius=15,
            border_width=2,
        )
        frame4.grid(row=1, column=0, columnspan=2, sticky=tk.W, padx=50)

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
            self.data_collection_main_frame,
            width=500,
            height=400,
            corner_radius=15,
            border_width=2,
        )

        amr_frame.grid(row=0, column=0, padx=50, pady=50)

        list_amr_label = ctk.CTkLabel(
            master=amr_frame,
            text="List available AMR datasets",
            font=self.default_font(20),
        )
        list_amr_label.place(x=50, y=20)
        self.download_button4 = ctk.CTkButton(
            master=amr_frame,
            text="load amr list",
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            font=self.default_font(12),
            command=self.load_amr_data,
        )
        self.download_button4.place(x=50, y=60)
        self.species_filter = Combobox(master=amr_frame, state=tk.DISABLED)
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
            text="Phenotype count ≥ 50",
            command=self.on_amr_list_filter_check,
            state=tk.DISABLED,
        )

        self.amr_list_filter_checkbox.place(x=50, y=100)

        Hovertip(
            self.amr_list_filter_checkbox,
            "Only show species with ≥ 50\nphenotypes (resistant and susceptible)",
        )

        self.total_label = ctk.CTkLabel(master=amr_frame, text="Total: not loaded")
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
            self.data_collection_main_frame,
            width=800,
            height=400,
            corner_radius=15,
            border_width=2,
        )
        frame6.grid(row=0, column=1, pady=50)
        # Create input fields for antibiotic and species
        full_amr_label = ctk.CTkLabel(
            master=frame6,
            text="Get amr data by species and antibiotic",
            font=self.default_font(20),
        )
        full_amr_label.place(x=50, y=20)
        self.antibiotic_selection = Combobox(master=frame6, state=tk.DISABLED, width=18)
        self.antibiotic_selection.bind(
            "<<ComboboxSelected>>", self.on_antibiotic_select
        )
        self.antibiotic_selection.place(x=50, y=50)

        self.species_selection = Combobox(master=frame6, state=tk.DISABLED, width=18)

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

        self.filter_contradictions_checkbox = ctk.CTkCheckBox(
            master=frame6,
            text="Filter Contradictions",
            command=self.update_amr_full,
            state=tk.DISABLED,
        )

        self.filter_contradictions_checkbox.place(x=390, y=90)

        Hovertip(
            self.numeric_phenotypes_checkbox, "0: Resistant\n1: Susceptible\n2: Other"
        )

        # Create a button to select the AMR metadata file
        self.save_table_button = ctk.CTkButton(
            master=frame6,
            text="Export to .tsv",
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            font=self.default_font(12),
            command=self.save_to_tsv,
            state=tk.DISABLED,
        )
        self.save_table_button.place(x=400, y=50)
        self.totaldata = ctk.CTkLabel(master=frame6, text="Total:")
        self.totaldata.place(x=50, y=120)
        self.totalresistance_label = ctk.CTkLabel(
            master=frame6, text="Total Resistant:"
        )
        self.totalresistance_label.place(x=220, y=120)
        self.totalsusceptible_label = ctk.CTkLabel(
            master=frame6, text="Total Susceptible:"
        )
        self.totalsusceptible_label.place(x=390, y=120)

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
        self.preprocessing_frame = ctk.CTkFrame(self, fg_color="transparent")

        self.preprocessing_frame_tab_view = ctk.CTkTabview(self.preprocessing_frame)

        self.preprocessing_frame_tab_view.pack(fill=tk.BOTH, expand=True)

        self.preprocessing_frame_tab_view.add("Preprocessing")
        self.preprocessing_frame_control = ControlFrame(
            self.preprocessing_frame_tab_view.tab("Preprocessing")
        )

        self.preprocessing_frame_control.control_panel.configure(width=200, height=200)

        self.preprocessing_frame_dataset_path = Label(
            master=self.preprocessing_frame_control.control_panel,
            fg_color="transparent",
            width=27,
            anchor="w",
        )

        self.preprocessing_frame_control_panel_dataset_button = ctk.CTkButton(
            master=self.preprocessing_frame_control.control_panel,
            text="Pick dataset",
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            command=lambda: (
                self.preprocessing_frame_dataset_path.configure(
                    text=util.select_directory()
                ),
                self.preprocessing_validate_ui(),
            ),
        )

        self.preprocessing_frame_control_panel_dataset_button.pack(
            anchor=tk.W, padx=20, pady=(30, 5)
        )

        self.preprocessing_frame_dataset_path.pack(anchor=tk.W, padx=20, pady=5)

        self.preprocessing_frame_control_panel_kmer_size_entry = ctk.CTkEntry(
            master=self.preprocessing_frame_control.control_panel,
            placeholder_text="K-mer size (Default: 31)",
            width=150,
        )
        self.preprocessing_frame_control_panel_kmer_size_entry.pack(
            anchor=tk.W, padx=20, pady=5
        )

        self.preprocessing_frame_control_panel_kmer_size_entry.bind(
            "<KeyRelease>", self.preprocessing_validate_ui
        )

        Hovertip(
            self.preprocessing_frame_control_panel_kmer_size_entry,
            "K-mer size format:\nodd_number\n\ne.g. 21",
        )

        self.kmer_tools = ["Ray Surveyor", "DSK"]

        self.preprocessing_frame_control_panel_kmer_tool_selector = Combobox(
            master=self.preprocessing_frame_control.control_panel,
            values=self.kmer_tools,
            state="readonly",
        )

        self.preprocessing_frame_control_panel_kmer_tool_selector.current(0)

        self.preprocessing_frame_control_panel_kmer_tool_selector.pack(
            anchor=tk.W, padx=20, pady=5
        )

        self.preprocessing_frame_control_panel_run_button = ctk.CTkButton(
            master=self.preprocessing_frame_control.control_panel,
            text=f"Run {self.preprocessing_frame_control_panel_kmer_tool_selector.get()}",
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            font=self.default_font(12),
            command=self.run_preprocessing,
        )
        self.preprocessing_frame_control_panel_run_button.pack(
            anchor=tk.W, padx=20, pady=5
        )

        self.preprocessing_frame_control_panel_kmer_tool_selector.bind(
            "<<ComboboxSelected>>",
            lambda e: self.preprocessing_frame_control_panel_run_button.configure(
                text=(e.widget.selection_clear(), f"Run {e.widget.get()}")[1]
            ),
        )

        self.preprocessing_frame_run_button_hover = Hovertip(
            self.preprocessing_frame_control_panel_run_button,
            "",
        )

        self.preprocessing_validate_ui()

    @threaded
    def run_preprocessing(self):
        output_directory = util.select_directory(title="Select output directory")

        if not output_directory:
            return

        self.preprocessing_frame_control_panel_dataset_button.configure(
            state=tk.DISABLED
        )
        self.preprocessing_frame_control_panel_kmer_size_entry.configure(
            state=tk.DISABLED
        )
        self.preprocessing_frame_control_panel_kmer_tool_selector.configure(
            state=tk.DISABLED
        )

        self.preprocessing_frame_control_panel_run_button.configure(
            text="Loading...", state=tk.DISABLED
        )

        self.clear_cmd_output(self.preprocessing_frame_control.cmd_output)

        if (
            self.preprocessing_frame_control_panel_kmer_tool_selector.get()
            == self.kmer_tools[0]
        ):
            self.run_ray_surveyor(output_directory)
        elif (
            self.preprocessing_frame_control_panel_kmer_tool_selector.get()
            == self.kmer_tools[1]
        ):
            self.run_dsk(output_directory)

        self.preprocessing_frame_control_panel_dataset_button.configure(state=tk.NORMAL)
        self.preprocessing_frame_control_panel_kmer_size_entry.configure(
            state=tk.NORMAL
        )
        self.preprocessing_frame_control_panel_kmer_tool_selector.configure(
            state="readonly"
        )
        self.preprocessing_frame_control_panel_run_button.configure(
            text=f"Run {self.preprocessing_frame_control_panel_kmer_tool_selector.get()}",
            command=self.run_preprocessing,
        )

    def cancel_process(self, process: Popen, output: ctk.CTkTextbox):
        if messagebox.askyesno(
            "Confirmation", "Are you sure you want to cancel this process?"
        ):
            try:
                process.kill()
            except Exception:
                pass

    def run_ray_surveyor(self, output_directory: str):
        dataset_folder = self.preprocessing_frame_dataset_path.cget("text")

        input_files = [
            f"{dataset_folder}/{file}"
            for file in os.listdir(dataset_folder)
            if file.endswith(".fna")
        ]

        if not input_files:
            util.update_cmd_output(
                "No .fna files found in the selected dataset folder.",
                self.preprocessing_frame_control.cmd_output,
            )
            return

        kmer_size = self.preprocessing_frame_control_panel_kmer_size_entry.get()

        if not kmer_size:
            kmer_size = 31

        config_path = self.generate_survey_conf(
            input_files, kmer_size, output_directory
        )

        try:
            util.update_cmd_output(
                "Running Ray Surveyor...\n", self.preprocessing_frame_control.cmd_output
            )

            ray_surveyor_command = f'mpiexec -n 4 "{util.to_linux_path(Path.RAY)}" "{util.to_linux_path(config_path)}"'

            process = util.run_bash_command(ray_surveyor_command, Path.TEMP)

            self.preprocessing_frame_control_panel_run_button.configure(
                text="Cancel",
                state=tk.NORMAL,
                command=lambda: self.cancel_process(
                    process, self.preprocessing_frame_control.cmd_output
                ),
            )

            util.display_process_output(
                process, self.preprocessing_frame_control.cmd_output
            )

            if process.poll() == 0:
                util.update_cmd_output(
                    f"\nRay Surveyor completed successfully.\n\nOutput stored in: {output_directory}",
                    self.preprocessing_frame_control.cmd_output,
                    Tag.SUCCESS,
                )
            else:
                util.update_cmd_output(
                    "\n\nProcess cancelled.",
                    self.preprocessing_frame_control.cmd_output,
                    Tag.ERROR,
                )
        except CalledProcessError as e:
            util.update_cmd_output(
                "Ray Surveyor encountered an error.",
                self.preprocessing_frame_control.cmd_output,
                Tag.ERROR,
            )
            util.update_cmd_output(
                e.stdout, self.preprocessing_frame_control.cmd_output, Tag.ERROR
            )
            util.update_cmd_output(
                e.stderr, self.preprocessing_frame_control.cmd_output, Tag.ERROR
            )
        finally:
            try:
                os.remove(config_path)
            except Exception:
                pass

    def run_dsk(self, output_directory: str):
        dataset_folder = self.preprocessing_frame_dataset_path.cget("text")

        kmer_size = self.preprocessing_frame_control_panel_kmer_size_entry.get()

        if not kmer_size:
            kmer_size = 31

        config_path = f"{output_directory}/dsk_output"

        try:
            util.update_cmd_output(
                "Running DSK...\n", self.preprocessing_frame_control.cmd_output
            )

            ls_command = f'ls -1 {util.to_linux_path(dataset_folder)}/*.fna > "{util.to_linux_path(config_path)}"'
            dsk_command = f'"{util.to_linux_path(Path.DSK)}" -file "{util.to_linux_path(config_path)}" -out-dir "{util.to_linux_path(output_directory)}" -kmer-size {kmer_size}'

            process = util.run_bash_command(f"{ls_command}\n{dsk_command}", Path.TEMP)

            self.preprocessing_frame_control_panel_run_button.configure(
                text="Cancel",
                state=tk.NORMAL,
                command=lambda: self.cancel_process(
                    process, self.preprocessing_frame_control.cmd_output
                ),
            )

            util.display_process_output(
                process, self.preprocessing_frame_control.cmd_output
            )

            if process.poll() == 0:
                util.update_cmd_output(
                    f"\nDSK completed successfully.\n\nOutput stored in: {output_directory}",
                    self.preprocessing_frame_control.cmd_output,
                    Tag.SUCCESS,
                )
            else:
                util.update_cmd_output(
                    "\n\nProcess cancelled.",
                    self.preprocessing_frame_control.cmd_output,
                    Tag.ERROR,
                )
        except CalledProcessError as e:
            util.update_cmd_output(
                "DSK encountered an error.",
                self.preprocessing_frame_control.cmd_output,
                Tag.ERROR,
            )
            util.update_cmd_output(
                e.stdout, self.preprocessing_frame_control.cmd_output, Tag.ERROR
            )
            util.update_cmd_output(
                e.stderr, self.preprocessing_frame_control.cmd_output, Tag.ERROR
            )
        finally:
            try:
                os.remove(config_path)
            except Exception:
                pass

    def preprocessing_validate_ui(self, event=None):
        self.preprocessing_frame_run_button_hover.text = ""

        failed = False

        if not self.preprocessing_frame_dataset_path.cget("text"):
            self.preprocessing_frame_run_button_hover.text += (
                "• select dataset directory.\n"
            )
            failed = True
        if not self.preprocessing_frame_control_panel_kmer_size_entry.get():
            self.preprocessing_frame_control_panel_kmer_size_entry.configure(
                border_color="#565B5E"
            )
        else:
            if (
                re.match(
                    r"^\d*[13579]$",
                    self.preprocessing_frame_control_panel_kmer_size_entry.get(),
                )
                is not None
            ):
                self.preprocessing_frame_control_panel_kmer_size_entry.configure(
                    border_color="green"
                )
            else:
                self.preprocessing_frame_control_panel_kmer_size_entry.configure(
                    border_color="red"
                )
                self.preprocessing_frame_run_button_hover.text += (
                    "• input a valid K-mer size.\n"
                )
                failed = True

        self.preprocessing_frame_run_button_hover.text = (
            self.preprocessing_frame_run_button_hover.text.strip("\n")
        )

        if failed:
            self.preprocessing_frame_run_button_hover.enable()
            self.preprocessing_frame_control_panel_run_button.configure(
                state=tk.DISABLED
            )
        else:
            self.preprocessing_frame_run_button_hover.disable()
            self.preprocessing_frame_control_panel_run_button.configure(state=tk.NORMAL)

    def create_kover_learn_page(self):
        self.kover_frame = ctk.CTkFrame(self, fg_color="transparent")

        self.kover_frame_tab_view = ctk.CTkTabview(self.kover_frame)

        self.kover_frame_tab_view.pack(fill=tk.BOTH, expand=True)

        self.kover_frame_tab_view.add("Create dataset")
        self.dataset_creation_frame = ControlFrame(
            self.kover_frame_tab_view.tab("Create dataset")
        )

        self.dataset_creation_frame.control_panel.grid(
            row=2,
            column=1,
            sticky=tk.NSEW,
            rowspan=3,
            columnspan=2,
            padx=40,
            pady=40,
        )

        self.dataset_creation_frame.cmd_output_frame.grid(
            row=0, column=3, rowspan=10, columnspan=7, sticky=tk.NSEW, padx=40, pady=40
        )

        self.dataset_creation_frame.control_panel.grid_rowconfigure(
            tuple(range(10)), pad=20
        )
        self.dataset_creation_frame.control_panel.grid_columnconfigure(
            tuple(range(2)), weight=1, uniform="column"
        )

        self.dataset_type = ["contigs", "kmer matrix"]

        self.dataset_creation_control_panel_dataset_type_selector = Combobox(
            master=self.dataset_creation_frame.control_panel,
            values=self.dataset_type,
            state="readonly",
        )

        self.dataset_creation_control_panel_dataset_type_selector.current(0)

        self.dataset_creation_control_panel_dataset_type_selector.grid(
            row=0, column=0, sticky=tk.W, padx=20, pady=20
        )

        self.dataset_creation_control_panel_dataset_type_selector.bind(
            "<<ComboboxSelected>>",
            lambda e: (e.widget.selection_clear(), self.dataset_creation_validate_ui()),
        )

        self.dataset_creation_frame_dataset_path = ctk.CTkEntry(
            master=self.dataset_creation_frame.control_panel,
            fg_color="transparent",
            state=tk.DISABLED,
        )

        self.dataset_creation_control_panel_dataset_button = ctk.CTkButton(
            master=self.dataset_creation_frame.control_panel,
            text="Pick dataset",
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            font=self.default_font(12),
            command=lambda: self.update_entry(
                self.dataset_creation_frame_dataset_path,
                util.select_file(
                    filetypes=[("TSV Files", "*.tsv")],
                    title="Select Dataset File",
                ),
                self.dataset_creation_validate_ui,
            ),
        )

        self.dataset_creation_control_panel_dataset_button.grid(
            row=1, column=0, sticky=tk.W, padx=(20, 50)
        )

        self.dataset_creation_frame_dataset_path.grid(
            row=1, column=1, sticky=tk.EW, padx=20
        )

        self.dataset_creation_frame_description_path = ctk.CTkEntry(
            master=self.dataset_creation_frame.control_panel,
            fg_color="transparent",
            state=tk.DISABLED,
        )

        self.dataset_creation_control_panel_description_button = ctk.CTkButton(
            master=self.dataset_creation_frame.control_panel,
            text="Pick Phenotype Description",
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            font=self.default_font(12),
            command=lambda: self.update_entry(
                self.dataset_creation_frame_description_path,
                util.select_file(
                    filetypes=[("TSV Files", "*.tsv")],
                    title="Select Phenotype Description File",
                ),
                self.dataset_creation_validate_ui,
            ),
        )

        self.dataset_creation_control_panel_description_button.grid(
            row=2, column=0, sticky=tk.W, padx=(20, 50), pady=0
        )

        self.dataset_creation_frame_description_path.grid(
            row=2, column=1, sticky=tk.EW, padx=20, pady=0
        )

        self.dataset_creation_frame_metadata_path = ctk.CTkEntry(
            master=self.dataset_creation_frame.control_panel,
            fg_color="transparent",
            state=tk.DISABLED,
        )

        self.dataset_creation_control_panel_metadata_button = ctk.CTkButton(
            master=self.dataset_creation_frame.control_panel,
            text="Pick Phenotype Metadata",
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            font=self.default_font(12),
            command=lambda: self.update_entry(
                self.dataset_creation_frame_metadata_path,
                util.select_file(
                    filetypes=[("TSV Files", "*.tsv")],
                    title="Select Phenotype Metadata File",
                ),
                self.dataset_creation_validate_ui,
            ),
        )

        self.dataset_creation_control_panel_metadata_button.grid(
            row=3, column=0, sticky=tk.W, padx=(20, 50), pady=0
        )

        self.dataset_creation_frame_metadata_path.grid(
            row=3, column=1, sticky=tk.EW, padx=20, pady=0
        )

        self.dataset_creation_frame_temp_path = ctk.CTkEntry(
            master=self.dataset_creation_frame.control_panel,
            fg_color="transparent",
            state=tk.DISABLED,
        )

        self.dataset_creation_control_panel_temp_button = ctk.CTkButton(
            master=self.dataset_creation_frame.control_panel,
            text="Pick Temp Directory",
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            font=self.default_font(12),
            command=lambda: self.update_entry(
                self.dataset_creation_frame_temp_path,
                util.select_directory(
                    title="Select Temp Directory",
                ),
            ),
        )

        self.dataset_creation_control_panel_temp_button.grid(
            row=4, column=0, sticky=tk.W, padx=(20, 50), pady=0
        )

        self.dataset_creation_frame_temp_path.grid(
            row=4, column=1, sticky=tk.EW, padx=20, pady=0
        )

        self.dataset_creation_frame_kmer_size_label = ctk.CTkLabel(
            master=self.dataset_creation_frame.control_panel,
            text="K-mer size",
            font=self.default_font(15),
        )
        self.dataset_creation_frame_kmer_size_label.grid(
            row=5, column=0, sticky=tk.W, padx=20
        )

        self.dataset_creation_frame_kmer_size_spinbox = Spinbox(
            master=self.dataset_creation_frame.control_panel,
            from_=3,
            to=128,
            increment=2,
            wrap=True,
            buttonbackground="#2b2b2b",
            disabledbackground="#595959",
            font=self.default_font(10),
        )
        self.dataset_creation_frame_kmer_size_spinbox.grid(
            row=6, column=0, sticky=tk.W, padx=20, pady=0
        )

        self.dataset_creation_frame_kmer_size_spinbox.set_default_value(31)

        self.dataset_creation_frame_compression_label = ctk.CTkLabel(
            master=self.dataset_creation_frame.control_panel,
            text="Compression Level",
            font=self.default_font(15),
        )
        self.dataset_creation_frame_compression_label.grid(
            row=5, column=1, sticky=tk.W, padx=20
        )

        self.dataset_creation_frame_compression_spinbox = Spinbox(
            master=self.dataset_creation_frame.control_panel,
            from_=0,
            to=9,
            wrap=True,
            buttonbackground="#2b2b2b",
            font=self.default_font(10),
        )
        self.dataset_creation_frame_compression_spinbox.grid(
            row=6, column=1, sticky=tk.W, padx=20
        )
        self.dataset_creation_frame_compression_spinbox.set_default_value(4)

        self.dataset_creation_control_panel_kmer_min_abundance_label = ctk.CTkLabel(
            master=self.dataset_creation_frame.control_panel,
            text="Minimum K-mer Abundance",
            font=self.default_font(15),
        )

        # self.dataset_creation_control_panel_kmer_min_abundance_label.grid(
        #     row=7, column=0, sticky=tk.W, padx=20, pady=(20, 0)
        # )

        self.dataset_creation_control_panel_kmer_min_abundance_spinbox = Spinbox(
            master=self.dataset_creation_frame.control_panel,
            from_=1,
            to=100,
            wrap=True,
            buttonbackground="#2b2b2b",
            font=self.default_font(10),
        )

        # self.dataset_creation_control_panel_kmer_min_abundance_spinbox.grid(
        #     row=8, column=0, sticky=tk.W, padx=20, pady=(20, 0)
        # )

        self.dataset_creation_control_panel_kmer_min_abundance_spinbox.set_default_value(
            1
        )

        self.dataset_creation_control_panel_singleton_kmer_checkbox = ctk.CTkCheckBox(
            master=self.dataset_creation_frame.control_panel,
            text="Singleton K-mers",
        )

        self.dataset_creation_control_panel_singleton_kmer_checkbox.grid(
            row=9, column=0, sticky=tk.W, padx=20, pady=(20, 0)
        )

        self.dataset_creation_control_panel_cpu_label = ctk.CTkLabel(
            master=self.dataset_creation_frame.control_panel,
            text="CPUs",
            font=self.default_font(15),
        )

        self.dataset_creation_control_panel_cpu_label.grid(
            row=7, column=1, sticky=tk.W, padx=20, pady=(20, 0)
        )

        self.dataset_creation_control_panel_cpu_spinbox = Spinbox(
            master=self.dataset_creation_frame.control_panel,
            from_=1,
            to=64,
            wrap=True,
            buttonbackground="#2b2b2b",
            disabledbackground="#595959",
            font=self.default_font(10),
        )

        self.dataset_creation_control_panel_cpu_spinbox.grid(
            row=8, column=1, sticky=tk.W, padx=20, pady=(20, 0)
        )

        self.dataset_creation_control_panel_cpu_spinbox.set_default_value(4)

        self.dataset_creation_frame_create_dataset_button = ctk.CTkButton(
            master=self.dataset_creation_frame.control_panel,
            text="Create Dataset",
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            font=self.default_font(12),
            command=self.create_dataset,
            state=tk.DISABLED,
        )

        self.dataset_creation_frame_create_dataset_button.grid(
            row=10, column=0, sticky=tk.W, padx=20, pady=(20, 0)
        )

        self.dataset_creation_frame_create_dataset_button_hover = Hovertip(
            self.dataset_creation_frame_create_dataset_button,
            "",
        )

        self.dataset_creation_validate_ui()

        self.kover_frame_tab_view.add("Split dataset")
        self.dataset_split_frame = ControlFrame(
            self.kover_frame_tab_view.tab("Split dataset")
        )

        self.dataset_split_frame.control_panel.grid(
            row=2,
            column=1,
            sticky=tk.NSEW,
            rowspan=3,
            columnspan=2,
            padx=40,
            pady=40,
        )

        self.dataset_split_frame.cmd_output_frame.grid(
            row=0, column=3, rowspan=10, columnspan=7, sticky=tk.NSEW, padx=40, pady=40
        )

        self.dataset_split_frame.control_panel.grid_rowconfigure(
            tuple(range(6)), pad=20
        )
        self.dataset_split_frame.control_panel.grid_columnconfigure(
            tuple(range(2)), weight=1, uniform="column"
        )

        self.dataset_split_frame_dataset_path = ctk.CTkEntry(
            master=self.dataset_split_frame.control_panel,
            fg_color="transparent",
            state=tk.DISABLED,
        )

        self.dataset_split_control_panel_dataset_button = ctk.CTkButton(
            master=self.dataset_split_frame.control_panel,
            text="Pick dataset",
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            font=self.default_font(12),
            command=lambda: self.update_entry(
                self.dataset_split_frame_dataset_path,
                util.select_file(
                    filetypes=[("Kover Files", "*.kover")],
                    title="Select Dataset File",
                ),
                self.dataset_split_validate_ui,
            ),
        )

        self.dataset_split_control_panel_dataset_button.grid(
            row=1, column=0, sticky=tk.W, padx=(20, 50), pady=(20, 0)
        )

        self.dataset_split_frame_dataset_path.grid(
            row=1, column=1, sticky=tk.EW, padx=20, pady=(20, 0)
        )

        self.dataset_split_frame_train_ids_path = ctk.CTkEntry(
            master=self.dataset_split_frame.control_panel,
            fg_color="transparent",
            state=tk.DISABLED,
        )

        self.dataset_split_control_panel_train_ids_button = ctk.CTkButton(
            master=self.dataset_split_frame.control_panel,
            text="Pick Train IDs",
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            font=self.default_font(12),
            command=lambda: self.update_entry(
                self.dataset_split_frame_train_ids_path,
                util.select_file(
                    filetypes=[("All Files", "*.*")],
                    title="Select Train IDs File",
                ),
                self.dataset_split_validate_ui,
            ),
        )

        self.dataset_split_control_panel_train_ids_button.grid(
            row=2, column=0, sticky=tk.W, padx=(20, 50), pady=0
        )

        self.dataset_split_frame_train_ids_path.grid(
            row=2, column=1, sticky=tk.EW, padx=20, pady=0
        )

        self.dataset_split_frame_test_ids_path = ctk.CTkEntry(
            master=self.dataset_split_frame.control_panel,
            fg_color="transparent",
            state=tk.DISABLED,
        )

        self.dataset_split_control_panel_test_ids_button = ctk.CTkButton(
            master=self.dataset_split_frame.control_panel,
            text="Pick Test IDs",
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            font=self.default_font(12),
            command=lambda: self.update_entry(
                self.dataset_split_frame_test_ids_path,
                util.select_file(
                    filetypes=[("All Files", "*.*")],
                    title="Select Test IDs File",
                ),
                self.dataset_split_validate_ui,
            ),
        )

        self.dataset_split_control_panel_test_ids_button.grid(
            row=3, column=0, sticky=tk.W, padx=(20, 50), pady=0
        )

        self.dataset_split_frame_test_ids_path.grid(
            row=3, column=1, sticky=tk.EW, padx=20, pady=0
        )

        self.dataset_split_frame_kmer_size_label = ctk.CTkLabel(
            master=self.dataset_split_frame.control_panel,
            text="Train size (%)",
            font=self.default_font(15),
        )
        self.dataset_split_frame_kmer_size_label.grid(
            row=4, column=0, sticky=tk.W, padx=20
        )

        self.dataset_split_frame_train_size_spinbox = Spinbox(
            master=self.dataset_split_frame.control_panel,
            from_=0,
            to=100,
            increment=1,
            wrap=True,
            buttonbackground="#2b2b2b",
            disabledbackground="#595959",
            font=self.default_font(10),
        )
        self.dataset_split_frame_train_size_spinbox.grid(
            row=5, column=0, sticky=tk.W, padx=20, pady=0
        )

        self.dataset_split_frame_train_size_spinbox.set_default_value(50)

        self.dataset_split_control_panel_cpu_label = ctk.CTkLabel(
            master=self.dataset_split_frame.control_panel,
            text="CPUs",
            font=self.default_font(15),
        )

        self.dataset_split_control_panel_cpu_label.grid(
            row=4, column=1, sticky=tk.W, padx=20
        )

        self.dataset_split_control_panel_cpu_spinbox = Spinbox(
            master=self.dataset_split_frame.control_panel,
            from_=1,
            to=64,
            wrap=True,
            buttonbackground="#2b2b2b",
            disabledbackground="#595959",
            font=self.default_font(10),
        )

        self.dataset_split_control_panel_cpu_spinbox.grid(
            row=5, column=1, sticky=tk.W, padx=20
        )

        self.dataset_split_control_panel_cpu_spinbox.set_default_value(4)

        self.dataset_split_frame_fold_label = ctk.CTkLabel(
            master=self.dataset_split_frame.control_panel,
            text="Folds",
            font=self.default_font(15),
        )
        self.dataset_split_frame_fold_label.grid(
            row=6, column=1, sticky=tk.W, padx=20, pady=(20, 0)
        )

        self.dataset_split_frame_fold_spinbox = Spinbox(
            master=self.dataset_split_frame.control_panel,
            from_=2,
            to=100,
            wrap=True,
            buttonbackground="#2b2b2b",
            disabledbackground="#595959",
            font=self.default_font(10),
        )
        self.dataset_split_frame_fold_spinbox.grid(
            row=7, column=1, sticky=tk.W, padx=20, pady=(20, 0)
        )
        self.dataset_split_frame_fold_spinbox.set_default_value(2)

        self.dataset_split_frame_fold_spinbox.configure(state=tk.DISABLED)

        self.dataset_split_control_panel_checkboxes = ctk.CTkFrame(
            master=self.dataset_split_frame.control_panel,
            fg_color="transparent",
        )

        self.dataset_split_control_panel_checkboxes.grid(
            row=6, column=0, sticky=tk.W, padx=20, pady=(20, 0)
        )

        self.dataset_split_control_panel_fold_checkbox = ctk.CTkCheckBox(
            master=self.dataset_split_control_panel_checkboxes,
            text="Folds",
            command=self.dataset_split_validate_ui,
        )

        self.dataset_split_control_panel_fold_checkbox.grid(row=0, column=1)

        self.dataset_split_control_panel_IDs_checkbox = ctk.CTkCheckBox(
            master=self.dataset_split_control_panel_checkboxes,
            text="Select IDs",
            command=self.dataset_split_validate_ui,
        )

        self.dataset_split_control_panel_IDs_checkbox.grid(row=0, column=0)

        self.dataset_split_control_panel_id_entry = ctk.CTkEntry(
            master=self.dataset_split_frame.control_panel,
            placeholder_text="Enter ID",
        )

        self.dataset_split_control_panel_id_entry.grid(
            row=7, column=0, sticky=tk.W, padx=20, pady=(20, 0)
        )

        self.dataset_split_control_panel_id_entry.bind(
            "<KeyRelease>", self.dataset_split_validate_ui
        )

        self.dataset_split_control_panel_seed_entry = ctk.CTkEntry(
            master=self.dataset_split_frame.control_panel,
            validate="key",
            validatecommand=(
                self.register(
                    lambda new_value: (
                        True if not new_value or new_value.isdigit() else False
                    )
                ),
                "%P",
            ),
        )

        self.dataset_split_control_panel_seed_entry.grid(
            row=8, column=1, sticky=tk.W, padx=20, pady=(20, 0)
        )

        self.dataset_split_control_panel_seed_entry.bind(
            "<FocusOut>",
            lambda _: (
                util.force_insertable_value(
                    0, self.dataset_split_control_panel_seed_entry
                )
                if not self.dataset_split_control_panel_seed_entry.get()
                else None
            ),
        )

        util.force_insertable_value(0, self.dataset_split_control_panel_seed_entry)

        self.dataset_split_control_panel_seed_button = ctk.CTkButton(
            master=self.dataset_split_frame.control_panel,
            text="Random seed",
            font=self.default_font(12),
            command=lambda: self.generate_random_seed(
                self.dataset_split_control_panel_seed_entry
            ),
        )

        self.dataset_split_control_panel_seed_button.grid(
            row=8, column=0, sticky=tk.W, padx=20, pady=(20, 0)
        )

        self.dataset_split_frame_split_dataset_button = ctk.CTkButton(
            master=self.dataset_split_frame.control_panel,
            text="Split Dataset",
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            font=self.default_font(12),
            command=self.split_dataset,
        )

        self.dataset_split_frame_split_dataset_button.grid(
            row=9, column=0, sticky=tk.W, padx=20, pady=(20, 0)
        )

        self.dataset_split_frame_split_dataset_button_hover = Hovertip(
            self.dataset_split_frame_split_dataset_button,
            "",
        )

        self.dataset_split_validate_ui()

        self.kover_frame_tab_view.add("View dataset")
        self.dataset_view_frame = ControlFrame(
            self.kover_frame_tab_view.tab("View dataset")
        )

        self.dataset_view_frame.control_panel.configure(height=220)

        self.dataset_view_frame.control_panel.grid(
            row=2,
            column=1,
            sticky=tk.NSEW,
            rowspan=3,
            columnspan=2,
            padx=40,
            pady=40,
        )
        self.dataset_view_frame.cmd_output_frame.grid(
            row=0, column=3, rowspan=10, columnspan=8, sticky=tk.NSEW, padx=40, pady=40
        )

        self.dataset_view_frame_a_checkbox = ctk.CTkCheckBox(
            master=self.dataset_view_frame.control_panel,
            text="All",
            text_color="red",
        )
        self.dataset_view_frame_a_checkbox.grid(
            row=0, column=0, sticky=tk.W, padx=20, pady=(20, 0)
        )
        Hovertip(
            self.dataset_view_frame_a_checkbox,
            "Warning:\n\nincludes kmers",
            fg="red",
        )
        self.dataset_view_frame_genome_type_checkbox = ctk.CTkCheckBox(
            master=self.dataset_view_frame.control_panel, text="Genome Type"
        )
        self.dataset_view_frame_genome_type_checkbox.select()
        self.dataset_view_frame_genome_type_checkbox.grid(
            row=0, column=1, sticky=tk.W, padx=20, pady=(20, 0)
        )
        self.dataset_view_frame_genome_source_checkbox = ctk.CTkCheckBox(
            master=self.dataset_view_frame.control_panel, text="Genome Source"
        )
        self.dataset_view_frame_genome_source_checkbox.select()
        self.dataset_view_frame_genome_source_checkbox.grid(
            row=0, column=2, sticky=tk.W, padx=20, pady=(20, 0)
        )
        self.dataset_view_frame_genome_ids_checkbox = ctk.CTkCheckBox(
            master=self.dataset_view_frame.control_panel, text="Genome IDs"
        )
        self.dataset_view_frame_genome_ids_checkbox.select()
        self.dataset_view_frame_genome_ids_checkbox.grid(
            row=1, column=0, sticky=tk.W, padx=20, pady=(20, 0)
        )
        self.dataset_view_frame_genome_count_checkbox = ctk.CTkCheckBox(
            master=self.dataset_view_frame.control_panel, text="Genome Count"
        )
        self.dataset_view_frame_genome_count_checkbox.select()
        self.dataset_view_frame_genome_count_checkbox.grid(
            row=1, column=1, sticky=tk.W, padx=20, pady=(20, 0)
        )
        self.dataset_view_frame_kmers_checkbox = ctk.CTkCheckBox(
            master=self.dataset_view_frame.control_panel,
            text="Kmers",
            text_color="red",
        )
        self.dataset_view_frame_kmers_checkbox.grid(
            row=1, column=2, sticky=tk.W, padx=20, pady=(20, 0)
        )
        Hovertip(
            self.dataset_view_frame_kmers_checkbox,
            "Warning:\n\nprints all kmers which could be in millions",
            fg="red",
        )
        self.dataset_view_frame_kmer_len_checkbox = ctk.CTkCheckBox(
            master=self.dataset_view_frame.control_panel, text="Kmer Length"
        )
        self.dataset_view_frame_kmer_len_checkbox.select()
        self.dataset_view_frame_kmer_len_checkbox.grid(
            row=2, column=0, sticky=tk.W, padx=20, pady=(20, 0)
        )
        self.dataset_view_frame_kmer_count_checkbox = ctk.CTkCheckBox(
            master=self.dataset_view_frame.control_panel, text="Kmer Count"
        )
        self.dataset_view_frame_kmer_count_checkbox.select()
        self.dataset_view_frame_kmer_count_checkbox.grid(
            row=2, column=1, sticky=tk.W, padx=20, pady=(20, 0)
        )
        self.dataset_view_frame_phenotype_description_checkbox = ctk.CTkCheckBox(
            master=self.dataset_view_frame.control_panel,
            text="Phenotype Description",
        )
        self.dataset_view_frame_phenotype_description_checkbox.select()
        self.dataset_view_frame_phenotype_description_checkbox.grid(
            row=2, column=2, sticky=tk.W, padx=20, pady=(20, 0)
        )
        self.dataset_view_frame_phenotype_metadata_checkbox = ctk.CTkCheckBox(
            master=self.dataset_view_frame.control_panel,
            text="Phenotype Metadata",
        )
        self.dataset_view_frame_phenotype_metadata_checkbox.select()
        self.dataset_view_frame_phenotype_metadata_checkbox.grid(
            row=3, column=0, sticky=tk.W, padx=20, pady=(20, 0)
        )
        self.dataset_view_frame_phenotype_tags_checkbox = ctk.CTkCheckBox(
            master=self.dataset_view_frame.control_panel, text="Phenotype Tags"
        )
        self.dataset_view_frame_phenotype_tags_checkbox.select()
        self.dataset_view_frame_phenotype_tags_checkbox.grid(
            row=3, column=1, sticky=tk.W, padx=20, pady=(20, 0)
        )
        self.dataset_view_frame_splits_checkbox = ctk.CTkCheckBox(
            master=self.dataset_view_frame.control_panel, text="Splits"
        )
        self.dataset_view_frame_splits_checkbox.select()
        self.dataset_view_frame_splits_checkbox.grid(
            row=3, column=2, sticky=tk.W, padx=20, pady=(20, 0)
        )
        self.dataset_view_frame_uuid_checkbox = ctk.CTkCheckBox(
            master=self.dataset_view_frame.control_panel, text="UUID"
        )
        self.dataset_view_frame_uuid_checkbox.select()
        self.dataset_view_frame_uuid_checkbox.grid(
            row=4, column=0, sticky=tk.W, padx=20, pady=(20, 0)
        )
        self.dataset_view_frame_compression_checkbox = ctk.CTkCheckBox(
            master=self.dataset_view_frame.control_panel, text="Compression"
        )
        self.dataset_view_frame_compression_checkbox.select()
        self.dataset_view_frame_compression_checkbox.grid(
            row=4, column=1, sticky=tk.W, padx=20, pady=(20, 0)
        )
        self.dataset_view_frame_classification_type_checkbox = ctk.CTkCheckBox(
            master=self.dataset_view_frame.control_panel,
            text="Classification Type",
        )
        self.dataset_view_frame_classification_type_checkbox.select()
        self.dataset_view_frame_classification_type_checkbox.grid(
            row=4, column=2, sticky=tk.W, padx=20, pady=(20, 0)
        )

        self.dataset_view_frame_view_dataset_button = ctk.CTkButton(
            master=self.dataset_view_frame.control_panel,
            text="View Dataset",
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            font=self.default_font(12),
            command=self.view_dataset,
        )

        self.dataset_view_frame_view_dataset_button.grid(
            row=5, column=0, sticky=tk.W, padx=20, pady=(100, 0)
        )

        self.kover_frame_tab_view.add("Kover learn")
        self.kover_learn_frame = ControlFrame(
            self.kover_frame_tab_view.tab("Kover learn")
        )

        self.kover_learn_frame.control_panel.grid_rowconfigure(tuple(range(13)), pad=20)
        self.kover_learn_frame.control_panel.grid_columnconfigure(
            tuple(range(2)), weight=1, uniform="column", minsize=300
        )

        self.kover_models = ["SCM", "CART"]

        self.kover_learn_frame_control_panel_kover_models_selector = Combobox(
            master=self.kover_learn_frame.control_panel,
            state=tk.DISABLED,
            values=self.kover_models,
        )

        self.kover_learn_frame_control_panel_kover_models_selector.current(0)

        self.kover_learn_frame_control_panel_kover_models_selector.bind(
            "<<ComboboxSelected>>", self.on_kover_model_selected
        )

        self.kover_learn_frame_dataset_path = ctk.CTkEntry(
            master=self.kover_learn_frame.control_panel,
            fg_color="transparent",
            state="readonly",
        )

        self.kover_learn_control_panel_dataset_button = ctk.CTkButton(
            master=self.kover_learn_frame.control_panel,
            text="Pick dataset",
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            font=self.default_font(12),
            command=self.kover_select_dataset,
        )

        self.kover_learn_frame_control_panel_split_selector = Combobox(
            master=self.kover_learn_frame.control_panel,
            state=tk.DISABLED,
        )
        self.kover_learn_frame_control_panel_split_selector.set("Split ID")

        self.kover_learn_frame_control_panel_split_selector.bind(
            "<<ComboboxSelected>>", lambda e: e.widget.selection_clear()
        )

        self.kover_model_types = (
            kover.ModelType.CONJUNCTION,
            kover.ModelType.DISJUNCTION,
            kover.ModelType.BOTH,
        )
        self.kover_criteria = (
            kover.Criterion.GINI,
            kover.Criterion.CROSS_ENTROPY,
        )
        self.kover_learn_frame_control_panel_model_criterion_selector = Combobox(
            master=self.kover_learn_frame.control_panel, state="readonly"
        )

        self.kover_learn_frame_control_panel_model_criterion_selector.bind(
            "<<ComboboxSelected>>", lambda e: e.widget.selection_clear()
        )

        self.kover_learn_control_panel_max_rules_max_depth_label = ctk.CTkLabel(
            master=self.kover_learn_frame.control_panel,
            font=self.default_font(15),
            text="",
        )

        self.kover_learn_control_panel_max_rules_max_depth_spinbox = Spinbox(
            master=self.kover_learn_frame.control_panel,
            to=100,
            wrap=True,
            buttonbackground="#2b2b2b",
            disabledbackground="#595959",
            font=self.default_font(10),
        )
        self.kover_learn_control_panel_max_equiv_rules_min_samples_split_label = (
            ctk.CTkLabel(
                master=self.kover_learn_frame.control_panel,
                font=self.default_font(15),
                text="",
            )
        )

        self.kover_learn_control_panel_max_equiv_rules_min_samples_split_spinbox = (
            Spinbox(
                master=self.kover_learn_frame.control_panel,
                wrap=True,
                buttonbackground="#2b2b2b",
                disabledbackground="#595959",
                font=self.default_font(10),
            )
        )

        self.kover_learn_control_panel_p_class_importance_entry = ctk.CTkEntry(
            master=self.kover_learn_frame.control_panel,
        )

        self.kover_learn_control_panel_p_class_importance_entry.bind(
            "<KeyRelease>", self.kover_learn_validate_ui
        )

        self.kover_learn_frame_kmer_blacklist_path = ctk.CTkEntry(
            master=self.kover_learn_frame.control_panel,
            fg_color="transparent",
            state="readonly",
        )

        self.kover_learn_control_panel_kmer_blacklist_button = ctk.CTkButton(
            master=self.kover_learn_frame.control_panel,
            text="Pick Kmer Blacklist",
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            font=self.default_font(12),
            command=lambda: self.update_entry(
                self.kover_learn_frame_kmer_blacklist_path,
                util.select_file(
                    filetypes=[("Fasta Files", "*.fa")],
                    title="Select Kmer Blacklist File",
                ),
            ),
        )

        self.kover_learn_frame_control_panel_hp_selector = Combobox(
            master=self.kover_learn_frame.control_panel,
            state="readonly",
            values=kover.HpChoice.BOUND,
        )

        self.kover_learn_frame_control_panel_hp_selector.current(0)

        def on_hp_selected(event=None):
            event.widget.selection_clear()
            if event.widget.get() == kover.HpChoice.BOUND:
                self.kover_learn_control_panel_bound_max_genome_size_spinbox.configure(
                    state=tk.NORMAL
                )
            else:
                self.kover_learn_control_panel_bound_max_genome_size_spinbox.configure(
                    state=tk.DISABLED
                )

        self.kover_learn_frame_control_panel_hp_selector.bind(
            "<<ComboboxSelected>>", on_hp_selected
        )

        self.kover_learn_control_panel_bound_max_genome_size_label = ctk.CTkLabel(
            master=self.kover_learn_frame.control_panel,
            text="Max genome size (0 for Default)",
            font=self.default_font(15),
            anchor=tk.W,
        )

        self.kover_learn_control_panel_bound_max_genome_size_spinbox = Spinbox(
            master=self.kover_learn_frame.control_panel,
            from_=0,
            to=100,
            wrap=True,
            buttonbackground="#2b2b2b",
            disabledbackground="#595959",
            font=self.default_font(10),
        )

        self.kover_learn_control_panel_bound_max_genome_size_spinbox.set_default_value(
            0
        )

        self.kover_learn_control_panel_seed_entry = ctk.CTkEntry(
            master=self.kover_learn_frame.control_panel,
            validate="key",
            validatecommand=(
                self.register(
                    lambda new_value: (
                        True if not new_value or new_value.isdigit() else False
                    )
                ),
                "%P",
            ),
        )

        self.kover_learn_control_panel_seed_button = ctk.CTkButton(
            master=self.kover_learn_frame.control_panel,
            text="Random seed",
            font=self.default_font(12),
            command=lambda: self.generate_random_seed(
                self.kover_learn_control_panel_seed_entry
            ),
        )

        self.kover_learn_control_panel_cpu_label = ctk.CTkLabel(
            master=self.kover_learn_frame.control_panel,
            text="CPUs",
            font=self.default_font(15),
        )

        self.kover_learn_control_panel_cpu_spinbox = Spinbox(
            master=self.kover_learn_frame.control_panel,
            from_=1,
            to=64,
            wrap=True,
            buttonbackground="#2b2b2b",
            disabledbackground="#595959",
            font=self.default_font(10),
        )

        self.kover_learn_control_panel_cpu_spinbox.set_default_value(4)

        self.kover_learn_frame_Initiate_button = ctk.CTkButton(
            master=self.kover_learn_frame.control_panel,
            text="Learn",
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            font=self.default_font(12),
            command=self.Initiate_kover_learn,
        )

        self.kover_learn_frame_Initiate_button_hover = Hovertip(
            self.kover_learn_frame_Initiate_button,
            "",
        )

        self.kover_learn_frame.control_panel.grid(
            row=2,
            column=1,
            sticky=tk.NSEW,
            rowspan=3,
            columnspan=2,
            padx=40,
            pady=40,
        )
        self.kover_learn_frame.cmd_output_frame.grid(
            row=0, column=3, rowspan=10, columnspan=8, sticky=tk.NSEW, padx=40, pady=40
        )

        self.kover_learn_control_panel_dataset_button.grid(
            row=0, column=0, sticky=tk.W, padx=(20, 50), pady=(20, 0)
        )
        self.kover_learn_frame_dataset_path.grid(
            row=0, column=1, sticky=tk.EW, padx=20, pady=(20, 0)
        )
        self.kover_learn_control_panel_kmer_blacklist_button.grid(
            row=1, column=0, sticky=tk.W, padx=(20, 50), pady=(20, 0)
        )
        self.kover_learn_frame_kmer_blacklist_path.grid(
            row=1, column=1, sticky=tk.EW, padx=20, pady=(20, 0)
        )
        self.kover_learn_frame_control_panel_kover_models_selector.grid(
            row=2, column=0, sticky=tk.W, padx=(20, 50), pady=(20, 0)
        )
        self.kover_learn_control_panel_p_class_importance_entry.grid(
            row=2, column=1, sticky=tk.EW, padx=20, pady=(20, 0)
        )
        self.kover_learn_frame_control_panel_hp_selector.grid(
            row=3, column=0, sticky=tk.W, padx=(20, 50), pady=(20, 0)
        )
        self.kover_learn_frame_control_panel_model_criterion_selector.grid(
            row=3, column=1, sticky=tk.W, padx=20, pady=(20, 0)
        )
        self.kover_learn_control_panel_bound_max_genome_size_label.grid(
            row=4, column=0, sticky=tk.W, padx=20, pady=(20, 0)
        )
        self.kover_learn_control_panel_max_rules_max_depth_label.grid(
            row=4, column=1, sticky=tk.W, padx=20, pady=(20, 0)
        )
        self.kover_learn_control_panel_bound_max_genome_size_spinbox.grid(
            row=5, column=0, sticky=tk.W, padx=20
        )
        self.kover_learn_control_panel_max_rules_max_depth_spinbox.grid(
            row=5, column=1, sticky=tk.W, padx=20
        )
        self.kover_learn_control_panel_cpu_label.grid(
            row=6, column=0, sticky=tk.W, padx=20
        )
        self.kover_learn_control_panel_max_equiv_rules_min_samples_split_label.grid(
            row=6, column=1, sticky=tk.W, padx=20, pady=(20, 0)
        )
        self.kover_learn_control_panel_cpu_spinbox.grid(
            row=7, column=0, sticky=tk.W, padx=20
        )
        self.kover_learn_control_panel_max_equiv_rules_min_samples_split_spinbox.grid(
            row=7, column=1, sticky=tk.W, padx=20
        )
        self.kover_learn_control_panel_seed_button.grid(
            row=8, column=0, sticky=tk.W, padx=20, pady=(20, 0)
        )
        self.kover_learn_control_panel_seed_entry.grid(
            row=8, column=1, sticky=tk.W, padx=20, pady=(20, 0)
        )
        self.kover_learn_frame_control_panel_split_selector.grid(
            row=9, column=0, sticky=tk.W, padx=(20, 50), pady=(20, 0)
        )
        self.kover_learn_frame_Initiate_button.grid(
            row=9, column=1, sticky=tk.W, padx=20, pady=(20, 0)
        )

        self.kover_learn_validate_ui()

    def kover_learn_validate_ui(self, event=None):
        self.kover_learn_frame_Initiate_button_hover.text = ""
        failed = False

        if not self.kover_learn_control_panel_p_class_importance_entry.get():
            self.kover_learn_control_panel_p_class_importance_entry.configure(
                border_color="#565B5E"
            )
        elif (
            re.match(
                r"^(\d+.\d+ )*(\d+.\d+){1}$",
                self.kover_learn_control_panel_p_class_importance_entry.get(),
            )
            is not None
        ):
            self.kover_learn_control_panel_p_class_importance_entry.configure(
                border_color="green"
            )
        else:
            failed = True
            self.kover_learn_control_panel_p_class_importance_entry.configure(
                border_color="red"
            )
            self.kover_learn_frame_Initiate_button_hover.text += (
                "• Invalid Hyper Parameter.\n"
            )

        if not self.kover_learn_frame_dataset_path.get():
            failed = True
            self.kover_learn_frame_Initiate_button_hover.text += "• Invalid Split ID.\n"

        self.kover_learn_frame_Initiate_button_hover.text = (
            self.kover_learn_frame_Initiate_button_hover.text.strip("\n")
        )

        if not failed:
            self.kover_learn_frame_Initiate_button_hover.disable()
            self.kover_learn_frame_Initiate_button.configure(state=tk.NORMAL)
        else:
            self.kover_learn_frame_Initiate_button_hover.enable()
            self.kover_learn_frame_Initiate_button.configure(state=tk.DISABLED)

    @threaded
    def Initiate_kover_learn(self):
        try:
            output_directory = util.select_directory(title="Select Output Directory")

            if not output_directory:
                return

            selected_dataset = self.kover_learn_frame_dataset_path.get()
            blacklist = self.kover_learn_frame_kmer_blacklist_path.get()
            split = self.kover_learn_frame_control_panel_split_selector.get()
            max_genome_size = (
                self.kover_learn_control_panel_bound_max_genome_size_spinbox.get()
            )
            cpu = self.kover_learn_control_panel_cpu_spinbox.get()
            seed = self.kover_learn_control_panel_seed_entry.get()
            hp_choice = self.kover_learn_frame_control_panel_hp_selector.get()
            p_class_importance = (
                self.kover_learn_control_panel_p_class_importance_entry.get()
            )

            learn_type = (
                self.kover_learn_frame_control_panel_kover_models_selector.get()
            )
            if learn_type == "SCM":
                model_type = (
                    self.kover_learn_frame_control_panel_model_criterion_selector.get()
                )
                max_rules = (
                    self.kover_learn_control_panel_max_rules_max_depth_spinbox.get()
                )
                max_equiv_rules = (
                    self.kover_learn_control_panel_max_equiv_rules_min_samples_split_spinbox.get()
                )
                seed = self.kover_learn_control_panel_seed_entry.get()

                command = kover.scm_command(
                    Path.KOVER,
                    selected_dataset,
                    split,
                    model_type,
                    p_class_importance,
                    max_rules,
                    max_equiv_rules,
                    blacklist,
                    hp_choice,
                    max_genome_size,
                    seed,
                    cpu,
                    output_directory,
                    True,
                    False,
                )

            elif learn_type == "CART":
                criterion = (
                    self.kover_learn_frame_control_panel_model_criterion_selector.get()
                )
                max_depth = (
                    self.kover_learn_control_panel_max_rules_max_depth_spinbox.get()
                )
                min_samples_split = (
                    self.kover_learn_control_panel_max_equiv_rules_min_samples_split_spinbox.get()
                )
                cpu = self.kover_learn_control_panel_cpu_spinbox.get()

                command = kover.tree_command(
                    Path.KOVER,
                    selected_dataset,
                    split,
                    criterion,
                    max_depth,
                    min_samples_split,
                    p_class_importance,
                    blacklist,
                    hp_choice,
                    max_genome_size,
                    cpu,
                    output_directory,
                    True,
                    False,
                )

            process = util.run_bash_command(command, Path.TEMP)

            self.kover_learn_frame_Initiate_button.configure(
                text="Cancel",
                command=lambda: self.cancel_process(
                    process, self.kover_learn_frame.cmd_output
                ),
            )

            self.clear_cmd_output(self.kover_learn_frame.cmd_output)

            util.update_cmd_output(
                "Initializing Kover learn...\n\n",
                self.kover_learn_frame.cmd_output,
                Tag.SYSTEM,
            )

            util.display_process_output(process, self.kover_learn_frame.cmd_output)

            if process.returncode == 0:
                util.update_cmd_output(
                    "\nKover learn completed successfully.",
                    self.kover_learn_frame.cmd_output,
                    Tag.SUCCESS,
                )
            else:
                util.update_cmd_output(
                    "\nKover learn failed.",
                    self.kover_learn_frame.cmd_output,
                    Tag.ERROR,
                )

        except Exception as e:
            messagebox.showerror("Error", e)
            traceback.print_exc()
        finally:
            self.kover_learn_frame_Initiate_button.configure(
                text="Learn", command=self.Initiate_kover_learn
            )

    def on_kover_model_selected(self, event=None):
        self.kover_learn_frame_control_panel_kover_models_selector.selection_clear()

        selected = self.kover_learn_frame_control_panel_kover_models_selector.get()

        self.kover_learn_control_panel_p_class_importance_entry.delete(0, tk.END)

        values = [kover.HpChoice.BOUND, kover.HpChoice.CV, kover.HpChoice.NONE]
        if (
            self.fold_count[self.kover_learn_frame_control_panel_split_selector.get()]
            == 0
        ):
            values.remove(kover.HpChoice.CV)
        if selected == "SCM":
            self.kover_learn_control_panel_seed_button.configure(state=tk.NORMAL)
            self.kover_learn_control_panel_seed_entry.configure(state=tk.NORMAL)
            self.kover_learn_frame_control_panel_hp_selector.configure(values=values)
            self.kover_learn_frame_control_panel_model_criterion_selector.configure(
                values=self.kover_model_types
            )
            self.kover_learn_frame_control_panel_model_criterion_selector.current(2)
            self.kover_learn_control_panel_max_rules_max_depth_label.configure(
                text="Max rules"
            )
            self.kover_learn_control_panel_max_rules_max_depth_spinbox.configure(
                from_=0
            )
            self.kover_learn_control_panel_max_rules_max_depth_spinbox.set_default_value(
                10
            )
            self.kover_learn_control_panel_max_equiv_rules_min_samples_split_label.configure(
                text="Max equivalent rules (0 for Default)"
            )
            self.kover_learn_control_panel_max_equiv_rules_min_samples_split_spinbox.configure(
                from_=0, to=10000
            )
            self.kover_learn_control_panel_max_equiv_rules_min_samples_split_spinbox.set_default_value(
                5
            )
            self.kover_learn_control_panel_p_class_importance_entry.configure(
                placeholder_text="Trade-off parameter (Ex: 0.1 0.178 0.316 0.562 1.0 ...)"
            )

        elif selected == "CART":
            values.remove(kover.HpChoice.NONE)
            self.kover_learn_control_panel_seed_button.configure(state=tk.DISABLED)
            self.kover_learn_control_panel_seed_entry.configure(state=tk.DISABLED)
            self.kover_learn_frame_control_panel_hp_selector.configure(values=values)
            self.kover_learn_frame_control_panel_model_criterion_selector.configure(
                values=self.kover_criteria
            )
            self.kover_learn_frame_control_panel_model_criterion_selector.current(0)
            self.kover_learn_control_panel_max_rules_max_depth_label.configure(
                text="Max depth (Default 10)"
            )
            self.kover_learn_control_panel_max_rules_max_depth_spinbox.configure(
                from_=1
            )
            self.kover_learn_control_panel_max_rules_max_depth_spinbox.set_default_value(
                10
            )
            self.kover_learn_control_panel_max_equiv_rules_min_samples_split_label.configure(
                text="Min samples split (Default 2)"
            )
            self.kover_learn_control_panel_max_equiv_rules_min_samples_split_spinbox.configure(
                from_=2, to=100
            )
            self.kover_learn_control_panel_max_equiv_rules_min_samples_split_spinbox.set_default_value(
                2
            )
            self.kover_learn_control_panel_p_class_importance_entry.configure(
                placeholder_text="Class importance (Ex: 0.25 0.5 0.75 1.0)"
            )
            if (
                self.kover_learn_frame_control_panel_hp_selector.get()
                == kover.HpChoice.NONE
            ):
                self.kover_learn_frame_control_panel_hp_selector.current(1)
                self.kover_learn_control_panel_bound_max_genome_size_spinbox.configure(
                    state=tk.DISABLED
                )

    @threaded
    def kover_select_dataset(self):
        dataset_path = util.select_file(
            filetypes=[("Kover Files", "*.kover")], title="Select Dataset File"
        )

        if not dataset_path:
            return

        self.clear_cmd_output(self.kover_learn_frame.cmd_output)

        self.kover_learn_frame_control_panel_split_selector.configure(state=tk.DISABLED)
        self.kover_learn_frame_control_panel_split_selector.set("Split ID")
        self.kover_learn_frame_control_panel_split_selector.unbind(
            "<<ComboboxSelected>>"
        )
        self.kover_learn_control_panel_dataset_button.configure(state=tk.DISABLED)

        util.update_cmd_output(
            "Fetching splits from dataset...\n\n",
            self.kover_learn_frame.cmd_output,
            Tag.SYSTEM,
        )

        command = kover.info_command(Path.KOVER, dataset_path, splits=True)
        process = util.run_bash_command(command, Path.TEMP)
        self.fold_count = {
            s.split()[0]: int(re.search(r"(Folds: )(\d+)", s).group(2))
            for s in process.stdout.read().splitlines()[1:]
        }
        split_count = len(self.fold_count)
        if split_count == 0:
            util.update_cmd_output(
                "No splits in datatset.\n\n",
                self.kover_learn_frame.cmd_output,
                Tag.ERROR,
            )
        else:
            self.kover_learn_frame_control_panel_split_selector.configure(
                values=tuple(self.fold_count.keys()), state="readonly"
            )
            self.kover_learn_frame_control_panel_split_selector.current(0)
            util.update_cmd_output(
                f"Fetched {split_count} split{'s' if split_count > 1 else ''}!\n\n{'ID':<10}Folds\n\n",
                self.kover_learn_frame.cmd_output,
                Tag.SUCCESS,
            )
            util.update_cmd_output(
                "\n".join(
                    [f"{key:<10}{value}" for key, value in self.fold_count.items()]
                ),
                self.kover_learn_frame.cmd_output,
                Tag.SUCCESS,
            )
            self.update_entry(self.kover_learn_frame_dataset_path, dataset_path)

            def on_split_selected(event=None):
                event.widget.selection_clear()
                values = [
                    kover.HpChoice.BOUND,
                    kover.HpChoice.CV,
                    kover.HpChoice.NONE,
                ]
                if (
                    self.kover_learn_frame_control_panel_kover_models_selector.get()
                    == "CART"
                ):
                    values.remove(kover.HpChoice.NONE)
                if self.fold_count[event.widget.get()] == 0:
                    values.remove(kover.HpChoice.CV)
                    self.kover_learn_frame_control_panel_hp_selector.configure(
                        values=values
                    )
                    self.kover_learn_frame_control_panel_hp_selector.current(0)
                else:
                    self.kover_learn_frame_control_panel_hp_selector.configure(
                        values=values
                    )
                    self.kover_learn_frame_control_panel_hp_selector.current(0)

            self.kover_learn_frame_control_panel_split_selector.bind(
                "<<ComboboxSelected>>", on_split_selected
            )

        self.kover_learn_control_panel_dataset_button.configure(state=tk.NORMAL)
        self.kover_learn_frame_control_panel_kover_models_selector.configure(
            state="readonly"
        )
        self.kover_learn_validate_ui()
        self.on_kover_model_selected()

    def dataset_split_validate_ui(self, event=None):
        self.dataset_split_frame_split_dataset_button_hover.text = ""

        failed = False

        if self.dataset_split_control_panel_IDs_checkbox.get():
            self.dataset_split_control_panel_train_ids_button.configure(state=tk.NORMAL)
            self.dataset_split_control_panel_test_ids_button.configure(state=tk.NORMAL)
            self.dataset_split_frame_train_size_spinbox.configure(state=tk.DISABLED)
        else:
            self.dataset_split_control_panel_train_ids_button.configure(
                state=tk.DISABLED
            )
            self.dataset_split_control_panel_test_ids_button.configure(
                state=tk.DISABLED
            )
            self.dataset_split_frame_train_size_spinbox.configure(state=tk.NORMAL)

        if self.dataset_split_control_panel_fold_checkbox.get():
            self.dataset_split_frame_fold_spinbox.configure(state=tk.NORMAL)
        else:
            self.dataset_split_frame_fold_spinbox.configure(state=tk.DISABLED)

        if not self.dataset_split_frame_dataset_path.get():
            self.dataset_split_frame_split_dataset_button_hover.text += (
                "• select dataset directory.\n"
            )
            failed = True
        if self.dataset_split_control_panel_IDs_checkbox.get() and not (
            self.dataset_split_frame_train_ids_path.get()
            and self.dataset_split_frame_test_ids_path.get()
        ):
            self.dataset_split_frame_split_dataset_button_hover.text += "• select either both train IDs and test IDs files or disable IDs path.\n"
            failed = True
        if not self.dataset_split_control_panel_id_entry.get():
            self.dataset_split_frame_split_dataset_button_hover.text += "• enter ID.\n"
            failed = True

        self.dataset_split_frame_split_dataset_button_hover.text = (
            self.dataset_split_frame_split_dataset_button_hover.text.strip("\n")
        )

        if failed:
            self.dataset_split_frame_split_dataset_button_hover.enable()
            self.dataset_split_frame_split_dataset_button.configure(state=tk.DISABLED)
        else:
            self.dataset_split_frame_split_dataset_button_hover.disable()
            self.dataset_split_frame_split_dataset_button.configure(state=tk.NORMAL)

    def update_entry(self, entry: ctk.CTkEntry, value: str, after: callable = None):
        entry.configure(state=tk.NORMAL)
        entry.delete(0, tk.END)
        if value:
            entry.insert(0, value)
        entry.configure(state=tk.DISABLED)

        if after:
            after()

    def dataset_creation_validate_ui(self):
        failed = False
        self.dataset_creation_frame_create_dataset_button_hover.text = ""

        if (
            self.dataset_creation_control_panel_dataset_type_selector.get()
            == self.dataset_type[1]
        ):
            self.dataset_creation_control_panel_singleton_kmer_checkbox.configure(
                state=tk.DISABLED
            )
            self.dataset_creation_frame_kmer_size_spinbox.configure(state=tk.DISABLED)
            self.dataset_creation_control_panel_cpu_spinbox.configure(state=tk.DISABLED)
        else:
            self.dataset_creation_control_panel_singleton_kmer_checkbox.configure(
                state=tk.NORMAL
            )
            self.dataset_creation_frame_kmer_size_spinbox.configure(state=tk.NORMAL)
            self.dataset_creation_control_panel_cpu_spinbox.configure(state=tk.NORMAL)

        if not self.dataset_creation_frame_dataset_path.get():
            self.dataset_creation_frame_create_dataset_button_hover.text += (
                "• select dataset directory.\n"
            )
            failed = True
        if not self.dataset_creation_frame_description_path.get():
            self.dataset_creation_frame_create_dataset_button_hover.text += (
                "• select phenotype description file.\n"
            )
            failed = True
        if not self.dataset_creation_frame_metadata_path.get():
            self.dataset_creation_frame_create_dataset_button_hover.text += (
                "• select phenotype metadata file.\n"
            )
            failed = True

        self.dataset_creation_frame_create_dataset_button_hover.text = (
            self.dataset_creation_frame_create_dataset_button_hover.text.strip("\n")
        )

        if failed:
            self.dataset_creation_frame_create_dataset_button_hover.enable()
            self.dataset_creation_frame_create_dataset_button.configure(
                state=tk.DISABLED
            )
        else:
            self.dataset_creation_frame_create_dataset_button_hover.disable()
            self.dataset_creation_frame_create_dataset_button.configure(state=tk.NORMAL)

    def create_analysis_page(self):
        self.analysis_frame = ctk.CTkFrame(
            self, corner_radius=0, fg_color="transparent"
        )

        self.webview = WebView2(self.analysis_frame, 500, 500)
        self.webview.pack(fill="both", expand=True, padx=100, pady=100)
        self.webview.load_url("http://127.0.0.1:5500/data/kover-amr-platform-gh-pages/index.html")

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

    @threaded
    def split_dataset(self):
        try:
            train_path = None
            test_path = None
            train_size = None
            split_id = self.dataset_split_control_panel_id_entry.get()
            folds = 0
            if self.dataset_split_control_panel_fold_checkbox.get():
                folds = self.dataset_split_frame_fold_spinbox.get()

            if self.dataset_split_control_panel_IDs_checkbox.get():
                train_path = self.dataset_split_frame_train_ids_path.get()
                test_path = self.dataset_split_frame_test_ids_path.get()
            else:
                train_size = (
                    float(self.dataset_split_frame_train_size_spinbox.get()) / 100
                )

            command = kover.split_command(
                Path.KOVER,
                self.dataset_split_frame_dataset_path.get(),
                split_id,
                train_size,
                train_path,
                test_path,
                folds,
                self.dataset_split_control_panel_seed_entry.get(),
                True,
                False,
            )

            process = util.run_bash_command(command, Path.TEMP)

            self.dataset_split_frame_split_dataset_button.configure(
                text="Cancel",
                command=lambda: self.cancel_process(
                    process, self.dataset_split_frame.cmd_output
                ),
            )

            self.clear_cmd_output(self.dataset_split_frame.cmd_output)

            util.update_cmd_output(
                "Processing dataset split request...\n\n",
                self.dataset_split_frame.cmd_output,
                Tag.SYSTEM,
            )

            util.display_process_output(process, self.dataset_split_frame.cmd_output)

            if process.returncode == 0:
                util.update_cmd_output(
                    "\nDataset split completed successfully.",
                    self.dataset_split_frame.cmd_output,
                    Tag.SUCCESS,
                )
            else:
                util.update_cmd_output(
                    "\nDataset split failed.",
                    self.dataset_split_frame.cmd_output,
                    Tag.ERROR,
                )

        except Exception as e:
            messagebox.showerror("Error", e)
            traceback.print_exc()
        finally:
            self.dataset_split_frame_split_dataset_button.configure(
                text="Split Dataset", command=self.split_dataset
            )

    @threaded
    def view_dataset(self):
        try:
            dataset = util.select_file(
                filetypes=[("Kover Files", "*.kover")], title="Select Dataset File"
            )

            if not dataset:
                return

            command = kover.info_command(
                Path.KOVER,
                dataset,
                self.dataset_view_frame_a_checkbox.get(),
                self.dataset_view_frame_genome_type_checkbox.get(),
                self.dataset_view_frame_genome_source_checkbox.get(),
                self.dataset_view_frame_genome_ids_checkbox.get(),
                self.dataset_view_frame_genome_count_checkbox.get(),
                self.dataset_view_frame_kmers_checkbox.get(),
                self.dataset_view_frame_kmer_len_checkbox.get(),
                self.dataset_view_frame_kmer_count_checkbox.get(),
                self.dataset_view_frame_phenotype_description_checkbox.get(),
                self.dataset_view_frame_phenotype_metadata_checkbox.get(),
                self.dataset_view_frame_phenotype_tags_checkbox.get(),
                self.dataset_view_frame_splits_checkbox.get(),
                self.dataset_view_frame_uuid_checkbox.get(),
                self.dataset_view_frame_compression_checkbox.get(),
                self.dataset_view_frame_classification_type_checkbox.get(),
            )

            process = util.run_bash_command(command, Path.TEMP)

            self.dataset_view_frame_view_dataset_button.configure(
                text="Cancel",
                command=lambda: self.cancel_process(
                    process, self.dataset_view_frame.cmd_output
                ),
            )

            self.clear_cmd_output(self.dataset_view_frame.cmd_output)

            util.update_cmd_output(
                "Processing dataset view request...\n\n",
                self.dataset_view_frame.cmd_output,
                Tag.SYSTEM,
            )

            util.display_process_output(
                process, self.dataset_view_frame.cmd_output, 1, 100
            )

            if process.returncode == 0:
                util.update_cmd_output(
                    "\n\nDataset view completed successfully.",
                    self.dataset_view_frame.cmd_output,
                    Tag.SUCCESS,
                )
            else:
                util.update_cmd_output(
                    "\n\nDataset view failed.",
                    self.dataset_view_frame.cmd_output,
                    Tag.ERROR,
                )

        except Exception as e:
            messagebox.showerror("Error", e)
            traceback.print_exc()
        finally:
            self.dataset_view_frame_view_dataset_button.configure(
                text="View Dataset", command=self.view_dataset
            )

    @threaded
    def create_dataset(self):
        try:
            output_path = util.select_directory(title="Select Output Directory")

            if not output_path:
                return

            output_path = os.path.join(output_path, "DATASET.kover")

            selected_source = (
                self.dataset_creation_control_panel_dataset_type_selector.get()
            )

            if selected_source == self.dataset_type[0]:
                source = kover.Source.CONTIGS
            elif selected_source == self.dataset_type[1]:
                source = kover.Source.K_MER_MATREX
            else:
                return

            command = kover.create_command(
                Path.KOVER,
                source,
                self.dataset_creation_frame_dataset_path.get(),
                output_path,
                self.dataset_creation_frame_description_path.get(),
                self.dataset_creation_frame_metadata_path.get(),
                self.dataset_creation_frame_kmer_size_spinbox.get(),
                self.dataset_creation_control_panel_kmer_min_abundance_spinbox.get(),
                self.dataset_creation_control_panel_singleton_kmer_checkbox.get(),
                self.dataset_creation_control_panel_cpu_spinbox.get(),
                self.dataset_creation_frame_compression_spinbox.get(),
                self.dataset_creation_frame_temp_path.get(),
                True,
                False,
            )

            process = util.run_bash_command(command, Path.TEMP)

            self.dataset_creation_frame_create_dataset_button.configure(
                text="Cancel",
                command=lambda: self.cancel_process(
                    process, self.dataset_creation_frame.cmd_output
                ),
            )

            self.clear_cmd_output(self.dataset_creation_frame.cmd_output)

            util.update_cmd_output(
                "Processing dataset creation request...\n\n",
                self.dataset_creation_frame.cmd_output,
                Tag.SYSTEM,
            )

            util.display_process_output(process, self.dataset_creation_frame.cmd_output)

            if process.returncode == 0:
                util.update_cmd_output(
                    "\nDataset creation completed successfully.",
                    self.dataset_creation_frame.cmd_output,
                    Tag.SUCCESS,
                )
            else:
                util.update_cmd_output(
                    "\nDataset creation failed.",
                    self.dataset_creation_frame.cmd_output,
                    Tag.ERROR,
                )

        except Exception as e:
            messagebox.showerror("Error", e)
            traceback.print_exc()
        finally:
            self.dataset_creation_frame_create_dataset_button.configure(
                text="Create Dataset", command=self.create_dataset
            )

    def generate_random_seed(self, seed_entry: ctk.CTkEntry):
        util.force_insertable_value(random.randint(1, 10000), seed_entry)

    @threaded
    def load_amr_data(self):
        self.download_button4.configure(text="Loading..", state=tk.DISABLED)
        self.amr_list_filter_checkbox.configure(state=tk.DISABLED)
        self.species_filter.configure(state=tk.DISABLED)
        self.species_selection.configure(state=tk.DISABLED)
        self.antibiotic_selection.configure(state=tk.DISABLED)
        self.drop_intermediate_checkbox.configure(state=tk.DISABLED)
        self.numeric_phenotypes_checkbox.configure(state=tk.DISABLED)
        self.filter_contradictions_checkbox.configure(state=tk.DISABLED)
        self.save_table_button.configure(state=tk.DISABLED)
        amr_metadata_file = util.select_file(filetypes=[("AMR Text Files", "*.txt")])
        if amr_metadata_file:
            self.total_label.configure(text="Total: ...")

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
                    "genome_id": str,
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
                    lambda x: (x["resistant_phenotype"] == "Resistant").sum() >= 50
                    and (x["resistant_phenotype"] == "Susceptible").sum() >= 50
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
            self.filter_contradictions_checkbox.configure(state=tk.NORMAL)
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

            self.total_label.configure(
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
            self.totalresistance_label.configure(text="Total Resistant: ...")
            self.totalsusceptible_label.configure(text="Total Susceptible: ...")
            self.totaldata.configure(text="Total: ...")
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

                if self.filter_contradictions_checkbox.get():
                    self.results_table.filtered_data = (
                        self.results_table.filtered_data.groupby("genome_id")
                        .filter(
                            lambda x: not (
                                len(x) > 1
                                and len(x["resistant_phenotype"].unique()) > 1
                            )
                        )
                        .reset_index(drop=True)
                    )

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
                    text=f"Total Resistant: {total_resistant}"
                )
                self.totalsusceptible_label.configure(
                    text=f"Total Susceptible: {total_susceptible}"
                )

                self.results_table.reset_view()

                self.totaldata.configure(text=f"Total: {total}")
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
        species = util.sanitize_filename(self.species_selection.get())
        antibiotics = util.sanitize_filename(self.antibiotic_selection.get())

        if species == "" or antibiotics == "":
            messagebox.showerror("Error", "Please load amr data first.")
            return

        self.save_table_button.configure(state=tk.DISABLED, text="Saving...")

        selected_directory = util.select_directory()

        if not selected_directory:
            self.save_table_button.configure(state=tk.NORMAL, text="Export to .tsv")
            return

        file_name = f"{species}_{antibiotics}.tsv"

        tsv_folder_path = os.path.join(selected_directory, species, antibiotics)

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

    def generate_survey_conf(
        self,
        input_files: list[str],
        kmer_size: int | str,
        output_dir: str,
    ) -> str:
        survey_conf_path = os.path.join(output_dir, "survey.conf")

        with open(survey_conf_path, "w") as survey_conf_file:
            survey_conf_file.write(f"-k {kmer_size}\n")
            survey_conf_file.write("-run-surveyor\n")
            survey_conf_file.write(
                f"-output {util.to_linux_path(output_dir)}/survey.res\n"
            )
            survey_conf_file.write("-write-kmer-matrix\n")

            for input_file in input_files:
                file_name = pl.Path(input_file).stem
                input_file = util.to_linux_path(input_file)
                survey_conf_file.write(
                    f"-read-sample-assembly {file_name} {input_file}\n"
                )

        return survey_conf_path

    def clear_cmd_output(self, output_target: ctk.CTkTextbox):
        output_target.configure(state=tk.NORMAL)
        output_target.delete("1.0", tk.END)
        output_target.configure(state=tk.DISABLED)
        output_target.update_idletasks()
