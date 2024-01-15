# from asyncio import sleep
import time
import tkinter as tk
from tkinter import ttk, messagebox, filedialog, font
import threading
import os
from subprocess import CalledProcessError, Popen, PIPE, run
from datetime import datetime
import csv
import random
from enum import Enum
import traceback

import customtkinter
from PIL import Image
import pandas as pd
import numpy as np

from kover import KoverDatasetCreator
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


class App(customtkinter.CTk):
    def __init__(self):
        super().__init__()

        self.ray_executable = "../bin/ray/Ray"
        self.ray_surveyor_dir = "raysurveyor"

        self.selected_directory = None
        self.genome_id = ""
        self.local_path = None
        self.amr_full = pd.DataFrame()
        customtkinter.set_appearance_mode("Dark")
        self.kmer_tool = tk.StringVar()
        self.kmer_tool.set("Ray Surveyor")
        self.dataset_folder = tk.StringVar()
        self.output_dir = tk.StringVar()
        self.metatsv_file_path = tk.StringVar()
        self.desctsv_file_path = tk.StringVar()
        self.selected_kover = tk.StringVar()
        self.selected_kmer_matrix = tk.StringVar()
        self.title("Genome analysis tool")
        # Get the screen resolution
        # Get the screen resolution
        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()

        # Set the window size to the screen resolution
        self.geometry(f"{screen_width}x{screen_height-100}")

        # Position the window at the top left corner of the screen
        self.geometry("+0+0")

        # set grid layout 1x2
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)

        # load images with light and dark mode image
        image_path = os.path.join(
            os.path.dirname(os.path.realpath(__file__)), "ui/test_images"
        )
        self.logo_image = customtkinter.CTkImage(
            Image.open(os.path.join(image_path, "CustomTkinter_logo_single.png")),
            size=(26, 26),
        )
        self.large_test_image = customtkinter.CTkImage(
            Image.open(os.path.join(image_path, "large_test_image.png")),
            size=(500, 150),
        )
        self.image_icon_image = customtkinter.CTkImage(
            Image.open(os.path.join(image_path, "image_icon_light.png")), size=(20, 20)
        )
        self.home_image = customtkinter.CTkImage(
            light_image=Image.open(os.path.join(image_path, "database-dark.png")),
            dark_image=Image.open(os.path.join(image_path, "database-light1.png")),
            size=(30, 30),
        )
        self.chat_image = customtkinter.CTkImage(
            light_image=Image.open(os.path.join(image_path, "preprocessing.png")),
            dark_image=Image.open(os.path.join(image_path, "preprocessing.png")),
            size=(30, 30),
        )
        self.add_user_image = customtkinter.CTkImage(
            light_image=Image.open(os.path.join(image_path, "add_user_dark.png")),
            dark_image=Image.open(os.path.join(image_path, "add_user_light.png")),
            size=(20, 20),
        )

        # create navigation frame
        self.navigation_frame = customtkinter.CTkFrame(self, corner_radius=0)
        self.navigation_frame.grid(row=0, column=0, sticky=tk.NSEW)
        self.navigation_frame.grid_rowconfigure(5, weight=1)

        self.navigation_frame_label = customtkinter.CTkLabel(
            self.navigation_frame,
            text="  Patric ",
            image=self.logo_image,
            compound="left",
            font=customtkinter.CTkFont(size=15, weight="bold"),
        )
        self.navigation_frame_label.grid(row=0, column=0, padx=20, pady=20)

        self.data_collection_frame_button = customtkinter.CTkButton(
            self.navigation_frame,
            corner_radius=0,
            height=80,
            border_spacing=10,
            text="Data collection",
            fg_color="transparent",
            text_color=("gray10", "gray90"),
            hover_color=("gray70", "gray30"),
            image=self.home_image,
            anchor="w",
            command=lambda: self.select_frame_by_name(Page.DATA_COLLECTION_PAGE),
        )
        self.data_collection_frame_button.grid(row=1, column=0, sticky="ew")

        self.preprocessing_frame_button = customtkinter.CTkButton(
            self.navigation_frame,
            corner_radius=0,
            height=80,
            border_spacing=10,
            text="Data preprocessing",
            fg_color="transparent",
            text_color=("gray10", "gray90"),
            hover_color=("gray70", "gray30"),
            image=self.chat_image,
            anchor="w",
            command=lambda: self.select_frame_by_name(Page.PREPROCESSING_PAGE),
        )
        self.preprocessing_frame_button.grid(row=2, column=0, sticky="ew")

        self.kover_frame_button = customtkinter.CTkButton(
            self.navigation_frame,
            corner_radius=0,
            height=80,
            border_spacing=10,
            text="Kover learn",
            fg_color="transparent",
            text_color=("gray10", "gray90"),
            hover_color=("gray70", "gray30"),
            image=self.add_user_image,
            anchor="w",
            command=lambda: self.select_frame_by_name(Page.KOVER_LEARN_PAGE),
        )
        self.kover_frame_button.grid(row=3, column=0, sticky="ew")

        self.analysis_frame_button = customtkinter.CTkButton(
            self.navigation_frame,
            corner_radius=0,
            height=80,
            border_spacing=10,
            text="Analysis",
            fg_color="transparent",
            text_color=("gray10", "gray90"),
            hover_color=("gray70", "gray30"),
            image=self.add_user_image,
            anchor="w",
            command=lambda: self.select_frame_by_name(Page.ANALYSIS_PAGE),
        )
        self.analysis_frame_button.grid(row=4, column=0, sticky="ew")

        self.appearance_mode_menu = customtkinter.CTkOptionMenu(
            self.navigation_frame,
            values=["Dark", "Light", "System"],
            command=customtkinter.set_appearance_mode,
        )
        self.appearance_mode_menu.grid(row=6, column=0, padx=20, pady=20, sticky="s")

        # create home frame
        self.data_collection_frame = customtkinter.CTkFrame(
            self, corner_radius=0, fg_color="transparent"
        )
        tab_view = customtkinter.CTkTabview(
            self.data_collection_frame,
            width=screen_width - 230,
            height=screen_height - 150,
        )
        tab_view.grid(row=0, column=0, padx=(20, 0), pady=(20, 0))
        tab_view.add("Genomes")
        tab_view.add("AMR")

        self.contigframe = customtkinter.CTkFrame(
            tab_view.tab("Genomes"),
            width=320,
            height=400,
            corner_radius=15,
            border_width=2,
        )
        self.contigframe.place(x=50, y=45)

        l2 = customtkinter.CTkLabel(
            master=self.contigframe,
            text="Contigs for specific Genome",
            font=("Century Gothic", 20),
        )
        l2.place(x=50, y=45)

        self.entry1 = customtkinter.CTkEntry(
            master=self.contigframe, width=220, placeholder_text="Genome Id"
        )
        self.entry1.place(x=50, y=110)

        # Create custom button
        self.dirbtn = customtkinter.CTkButton(
            master=self.contigframe,
            width=220,
            text="Select directory",
            corner_radius=6,
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            command=self.select_directory,
        )
        self.dirbtn.place(x=50, y=165)
        self.loadtsv_button = customtkinter.CTkButton(
            master=self.contigframe,
            width=220,
            text="Load bulk Genomes ids",
            corner_radius=6,
            command=self.select_tsv_file,
        )
        self.loadtsv_button.place(x=50, y=200)

        self.download_button = customtkinter.CTkButton(
            master=self.contigframe,
            width=220,
            text="Download",
            corner_radius=6,
            command=self.download_contigs,
        )
        self.download_button.place(x=50, y=245)
        self.cancel_button = customtkinter.CTkButton(
            master=self.contigframe,
            width=220,
            text="Cancel",
            corner_radius=6,
            command=self.cancelcontigs,
            state=tk.DISABLED,
        )
        self.cancel_button.place(x=50, y=290)

        self.progress_bar = ttk.Progressbar(
            master=self.contigframe, length=200, mode="determinate"
        )
        self.progress_bar.place(x=50, y=340)

        self.size_label = customtkinter.CTkLabel(
            master=self.contigframe,
            text="",
            font=("Century Gothic", 10),
            fg_color="transparent",
        )
        self.size_label.place(x=50, y=370)

        # creating specific genome frame
        self.genomeframe = customtkinter.CTkFrame(
            tab_view.tab("Genomes"),
            width=320,
            height=400,
            corner_radius=15,
            border_width=2,
        )
        self.genomeframe.place(x=500, y=45)

        l2 = customtkinter.CTkLabel(
            master=self.genomeframe,
            text="Features for specific Genome",
            font=("Century Gothic", 20),
        )
        l2.place(x=50, y=45)

        self.entry2 = customtkinter.CTkEntry(
            master=self.genomeframe, width=220, placeholder_text="Genome Id"
        )
        self.entry2.place(x=50, y=110)

        # Create custom button
        dirbtn1 = customtkinter.CTkButton(
            master=self.genomeframe,
            width=220,
            text="Select directory",
            corner_radius=6,
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            command=self.select_directory,
        )
        dirbtn1.place(x=50, y=165)
        self.loadtsv_button1 = customtkinter.CTkButton(
            master=self.genomeframe,
            width=220,
            text="Load bulk Genomes ids",
            corner_radius=6,
            command=self.select_tsv_file,
        )
        self.loadtsv_button1.place(x=50, y=200)

        self.download_button1 = customtkinter.CTkButton(
            master=self.genomeframe,
            width=220,
            text="Download",
            corner_radius=6,
            command=self.downloadfeatures,
        )
        self.download_button1.place(x=50, y=245)
        self.cancel_button1 = customtkinter.CTkButton(
            master=self.genomeframe,
            width=220,
            text="Cancel",
            corner_radius=6,
            command=self.cancelfeatures,
            state=tk.DISABLED,
        )
        self.cancel_button1.place(x=50, y=290)

        self.progress_bar1 = ttk.Progressbar(
            master=self.genomeframe, length=220, mode="determinate"
        )
        self.progress_bar1.place(x=50, y=340)

        self.size_label1 = customtkinter.CTkLabel(
            master=self.genomeframe,
            text="",
            font=("Century Gothic", 10),
            fg_color="transparent",
        )
        self.size_label1.place(x=50, y=370)

        # creating genome metadata frame
        metadataframe = customtkinter.CTkFrame(
            tab_view.tab("Genomes"),
            width=320,
            height=400,
            corner_radius=15,
            border_width=2,
        )
        metadataframe.place(x=900, y=45)

        l2 = customtkinter.CTkLabel(
            master=metadataframe,
            text="Latest metadata for Genomes",
            font=("Century Gothic", 20),
        )
        l2.place(x=50, y=45)

        # Create custom button
        dirbtn2 = customtkinter.CTkButton(
            master=metadataframe,
            width=220,
            text="Select directory",
            corner_radius=6,
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            command=self,
        )
        dirbtn2.place(x=50, y=165)

        self.download_button2 = customtkinter.CTkButton(
            master=metadataframe,
            width=220,
            text="Download",
            corner_radius=6,
            command=self,
        )
        self.download_button2.place(x=50, y=240)
        self.cancel_button2 = customtkinter.CTkButton(
            master=metadataframe,
            width=220,
            text="Cancel",
            corner_radius=6,
            command=self,
        )
        self.cancel_button2.place(x=50, y=290)

        self.progress_bar2 = ttk.Progressbar(
            master=metadataframe, length=220, mode="determinate"
        )
        self.progress_bar2.place(x=50, y=340)

        self.size_label2 = customtkinter.CTkLabel(
            master=metadataframe,
            text="",
            font=("Century Gothic", 10),
            fg_color="transparent",
        )
        self.size_label2.place(x=50, y=370)

        self.download_window = DownloadWindow(
            metadataframe,
            download_button=self.download_button2,
            cancel_button=self.cancel_button2,
            select_path_button=dirbtn2,
            progress_bar=self.progress_bar2,
            size_label=self.size_label2,
            path_label=self.size_label2,
        )

        self.download_app = FTPDownloadApp(
            self.download_window, "RELEASE_NOTES/genome_metadata"
        )

        # logic

        # creating AMR metadata frame
        frame4 = customtkinter.CTkFrame(
            tab_view.tab("AMR"),
            width=1000,
            height=400,
            corner_radius=15,
            border_width=2,
        )
        frame4.place(x=50, y=470)

        l2 = customtkinter.CTkLabel(
            master=frame4, text="Latest metadata for AMR", font=("Century Gothic", 20)
        )
        l2.place(x=50, y=45)

        # Create custom button
        self.update_date = customtkinter.CTkLabel(
            master=frame4, text="", font=("Century Gothic", 12), fg_color="transparent"
        )
        self.update_date.place(x=50, y=90)
        dirbtn3 = customtkinter.CTkButton(
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

        self.viewupdate_date = customtkinter.CTkButton(
            master=frame4,
            width=150,
            text="View last update date",
            corner_radius=6,
            command=lambda: get_last_metadata_update_date(self.update_date),
        )
        self.viewupdate_date.place(x=250, y=120)

        self.download_button3 = customtkinter.CTkButton(
            master=frame4, width=150, text="Download", corner_radius=6, command=self
        )
        self.download_button3.place(x=450, y=120)
        self.cancel_button3 = customtkinter.CTkButton(
            master=frame4, width=150, text="Cancel", corner_radius=6, command=self
        )
        self.cancel_button3.place(x=650, y=120)

        self.progress_bar3 = ttk.Progressbar(
            master=frame4, length=800, mode="determinate"
        )
        self.progress_bar3.place(x=50, y=190)

        self.size_label3 = customtkinter.CTkLabel(
            master=frame4, text="", font=("Century Gothic", 10), fg_color="transparent"
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
        frame5 = customtkinter.CTkFrame(
            tab_view.tab("AMR"),
            width=500,
            height=400,
            corner_radius=15,
            border_width=2,
        )

        frame5.place(x=50, y=45)
        # Create a style
        style = ttk.Style(frame5)

        # Import the tcl file
        frame5.tk.call("source", "ui/forest-dark.tcl")

        # Set the theme with the theme_use method
        style.theme_use("forest-dark")

        l2 = customtkinter.CTkLabel(
            master=frame5,
            text="List available AMR datasets",
            font=("Century Gothic", 20),
        )
        l2.place(x=50, y=20)
        self.download_button4 = customtkinter.CTkButton(
            master=frame5,
            text="load amr list",
            corner_radius=6,
            command=self.load_amr_data,
        )
        self.download_button4.place(x=50, y=60)
        self.species_filter = ttk.Combobox(master=frame5, state="readonly")
        self.species_filter.bind("<<ComboboxSelected>>", self.update_table)
        self.species_filter.bind("<KeyRelease>", self.update_table)
        self.species_filter.place(x=280, y=60)
        self.total_label = tk.Label(master=frame5, text="Total: not loaded")
        self.total_label.place(x=280, y=100)
        columns = ["Species", "Antibiotic"]
        self.amr_list_table = ttk.Treeview(
            master=frame5,
            columns=columns,
            show="headings",
            style="Treeview.Treeview",
            selectmode="browse",
            height=9,
        )
        self.amr_list_table.place(x=50, y=130)
        for index, col in enumerate(columns):
            self.amr_list_table.heading(
                col,
                text=col,
                command=lambda i=index: self.sort_amr_list_table(i),
            )
            self.amr_list_table.column(col, width=185)
        self.amr_list_table.bind("<Double-1>", self.on_table_double_click)
        self.amr_list_table.bind("<MouseWheel>", self.on_table_scroll)

        self.amr_list_columns_sort_order = list(range(len(columns)))
        self.sort_amr_list_ascending = False

        self.filtered_data = pd.DataFrame()
        self.amr_list = pd.DataFrame()
        self.start = 0

        frame6 = customtkinter.CTkFrame(
            tab_view.tab("AMR"),
            width=800,
            height=400,
            corner_radius=15,
            border_width=2,
        )
        frame6.place(x=580, y=45)
        # Create input fields for antibiotic and species
        l2 = customtkinter.CTkLabel(
            master=frame6,
            text="Get amr data by species and antibiotic",
            font=("Century Gothic", 20),
        )
        l2.place(x=50, y=20)
        self.antibiotic_entry = customtkinter.CTkEntry(
            master=frame6, width=150, placeholder_text="Antibiotic"
        )
        self.antibiotic_entry.bind("<KeyRelease>", self.update_amr_full)
        self.antibiotic_entry.place(x=50, y=50)

        self.species_entry = customtkinter.CTkEntry(
            master=frame6, width=150, placeholder_text="Species(optional)"
        )
        self.species_entry.bind("<KeyRelease>", self.update_amr_full)
        self.species_entry.place(x=220, y=50)

        self.drop_intermediate = tk.BooleanVar(value=False)
        self.drop_intermediate_checkbox = customtkinter.CTkCheckBox(
            master=frame6, text="Drop Intermediate", variable=self.drop_intermediate
        )
        self.drop_intermediate_checkbox.place(x=50, y=90)

        # Create a button to select the AMR metadata file
        self.save_table_button = customtkinter.CTkButton(
            master=frame6, text="Download data", command=self.save_to_tsv
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
            "Genome Name",
            "Genome ID",
            "Phenotype",
            "Measurements",
        ]
        self.results_table = ttk.Treeview(
            master=frame6,
            columns=columns,
            show="headings",
            height=8,
        )
        for index, col in enumerate(columns):
            self.results_table.heading(
                col,
                text=col,
                command=lambda i=index: self.sort_results_table(i),
            )
            self.results_table.column(col, width=165)

        self.results_table.place(x=50, y=150)
        self.amr_full_columns_sort_order = list(range(len(columns)))
        self.sort_amr_full_ascending = False

        # create second frame
        self.preprocessing_frame = customtkinter.CTkFrame(
            self, corner_radius=0, fg_color="transparent"
        )

        tab_view = customtkinter.CTkTabview(
            self.preprocessing_frame,
            width=screen_width - 230,
            height=screen_height - 150,
        )
        tab_view.grid(row=0, column=0, padx=(20, 0), pady=(20, 0))
        tab_view.add("preprocessing")
        self.controlpanel = customtkinter.CTkFrame(
            tab_view.tab("preprocessing"),
            width=550,
            height=500,
            corner_radius=15,
            border_width=2,
        )
        self.controlpanel.place(x=50, y=20)
        l2 = customtkinter.CTkLabel(
            master=self.controlpanel, text="Control panel", font=("Century Gothic", 20)
        )
        l2.place(x=50, y=45)

        self.datasetbtn = customtkinter.CTkButton(
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

        self.entry1 = customtkinter.CTkEntry(
            master=self.controlpanel,
            width=380,
            placeholder_text="Dataset path",
            textvariable=self.dataset_folder,
        )
        self.entry1.place(x=50, y=200)

        self.outputbtn = customtkinter.CTkButton(
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

        self.entry2 = customtkinter.CTkEntry(
            master=self.controlpanel,
            width=380,
            placeholder_text="Output directory path",
            textvariable=self.output_dir,
        )
        self.entry2.place(x=50, y=300)

        l3 = customtkinter.CTkLabel(
            master=self.controlpanel,
            text="Enter kmer length",
            font=("Century Gothic", 15),
        )
        l3.place(x=50, y=350)
        self.kmer_length = customtkinter.CTkEntry(
            master=self.controlpanel, width=100, placeholder_text="kmer length"
        )
        self.kmer_length.place(x=50, y=400)
        l4 = customtkinter.CTkLabel(
            master=self.controlpanel, text="Pick kmer tool", font=("Century Gothic", 15)
        )
        l4.place(x=300, y=350)
        self.kmertool_filter = ttk.Combobox(
            master=self.controlpanel,
            textvariable=self.kmer_tool,
            values=["Ray Surveyor", "DSK"],
            state="readonly",
        )
        self.kmertool_filter.place(x=300, y=400)

        self.runbtn = customtkinter.CTkButton(
            master=self.controlpanel,
            width=150,
            text="Run " + self.kmer_tool.get(),
            corner_radius=6,
            command=self.run_preprocessing,
        )
        self.runbtn.place(x=150, y=450)

        self.kmer_tool.trace_add("write", self.update_button_text)

        custom_font = font.nametofont("TkDefaultFont")
        custom_font.configure(size=12)

        outputscreen = customtkinter.CTkFrame(
            tab_view.tab("preprocessing"),
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
            font=custom_font,
        )

        self.scrollbar = tk.Scrollbar(outputscreen, command=self.cmd_output.yview)
        self.scrollbar.grid(row=0, column=1, sticky=tk.NSEW)
        self.cmd_output["yscrollcommand"] = self.scrollbar.set
        self.cmd_output.grid(row=0, column=0, sticky=tk.NSEW, padx=2, pady=2)

        # create third frame
        self.kover_frame = customtkinter.CTkFrame(
            self, corner_radius=0, fg_color="transparent"
        )
        tab_view = customtkinter.CTkTabview(
            self.kover_frame, width=screen_width - 230, height=screen_height - 150
        )
        tab_view.grid(row=0, column=0, padx=(20, 0), pady=(20, 0))
        tab_view.add("Create dataset")
        tab_view.add("split dataset")
        tab_view.add("kover learn")

        create_dataset_frame = customtkinter.CTkFrame(
            tab_view.tab("Create dataset"),
            width=700,
            height=screen_height - 230,
            corner_radius=15,
            border_width=2,
        )
        create_dataset_frame.place(x=50, y=20)

        l1 = customtkinter.CTkLabel(
            master=create_dataset_frame,
            text="Create Dataset",
            font=("Century Gothic", 20),
        )
        l1.place(x=50, y=45)

        # Selection option for dataset type
        dataset_type_label = customtkinter.CTkLabel(
            master=create_dataset_frame,
            text="Select Dataset Type",
            font=("Century Gothic", 15),
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
        self.pickdataset1_btn = customtkinter.CTkButton(
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

        self.pickdataset1_entry = customtkinter.CTkEntry(
            master=create_dataset_frame,
            width=380,
            placeholder_text="Dataset path",
            textvariable=self.dataset_folder,
        )
        self.pickdataset1_entry.place(x=50, y=250)

        self.pickoutput_btn = customtkinter.CTkButton(
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

        self.pickoutput_entry = customtkinter.CTkEntry(
            master=create_dataset_frame,
            width=380,
            placeholder_text="Output directory path",
            textvariable=self.output_dir,
        )
        self.pickoutput_entry.place(x=50, y=350)

        self.phenotype_desc_btn = customtkinter.CTkButton(
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

        self.phenotype_desc_entry = customtkinter.CTkEntry(
            master=create_dataset_frame,
            width=380,
            placeholder_text="Phenotype description path",
            textvariable=self.desctsv_file_path,
        )
        self.phenotype_desc_entry.place(x=50, y=450)

        # Kmer entry with max kmer length 128
        self.kmer_label = customtkinter.CTkLabel(
            master=create_dataset_frame,
            text="Enter Kmer Length (max 128)",
            font=("Century Gothic", 15),
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
        compression_label = customtkinter.CTkLabel(
            master=create_dataset_frame,
            text="Enter Compression Level (0-9)",
            font=("Century Gothic", 15),
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
        self.create_dataset_btn = customtkinter.CTkButton(
            master=create_dataset_frame,
            width=150,
            text="Create Dataset",
            corner_radius=6,
            command=self.create_dataset_thread,
        )
        self.create_dataset_btn.place(x=500, y=630)

        phenotype_metadata_btn = customtkinter.CTkButton(
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

        self.phenotype_metadata_entry = customtkinter.CTkEntry(
            master=create_dataset_frame,
            width=380,
            placeholder_text="Phenotype metadata path",
            textvariable=self.metatsv_file_path,
        )
        self.phenotype_metadata_entry.place(x=50, y=550)

        outputscreen = customtkinter.CTkFrame(
            tab_view.tab("Create dataset"),
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
            font=custom_font,
        )

        self.scrollbar1 = tk.Scrollbar(outputscreen, command=self.cmd_output1.yview)
        self.scrollbar1.grid(row=0, column=1, sticky=tk.NSEW)
        self.cmd_output1["yscrollcommand"] = self.scrollbar1.set
        self.cmd_output1.grid(row=0, column=0, sticky=tk.NSEW, padx=2, pady=2)

        create_dataset_frame = customtkinter.CTkFrame(
            tab_view.tab("split dataset"),
            width=700,
            height=screen_height - 230,
            corner_radius=15,
            border_width=2,
        )
        create_dataset_frame.place(x=50, y=20)
        l1 = customtkinter.CTkLabel(
            master=create_dataset_frame,
            text="Split Dataset",
            font=("Century Gothic", 20),
        )
        l1.place(x=50, y=45)

        dataset_btn = customtkinter.CTkButton(
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

        dataset_entry = customtkinter.CTkEntry(
            master=create_dataset_frame,
            width=380,
            placeholder_text="Dataset path",
            textvariable=self.selected_kover,
        )
        dataset_entry.place(x=50, y=200)

        output_btn = customtkinter.CTkButton(
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

        output_entry = customtkinter.CTkEntry(
            master=create_dataset_frame,
            width=380,
            placeholder_text="Output directory path",
            textvariable=self.output_dir,
        )
        output_entry.place(x=50, y=300)

        l1 = customtkinter.CTkLabel(
            master=create_dataset_frame,
            text="Train size(%)",
            font=("Century Gothic", 16),
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

        trainlabel = customtkinter.CTkLabel(
            master=create_dataset_frame,
            text="Use train ids and tests",
            font=("Century Gothic", 16),
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
        foldlabel = customtkinter.CTkLabel(
            master=create_dataset_frame, text="Folds", font=("Century Gothic", 16)
        )
        foldlabel.place(x=50, y=400)
        foldentry = customtkinter.CTkEntry(
            master=create_dataset_frame,
            width=80,
            placeholder_text="Enter number of folds",
            textvariable=self.foldvar,
        )

        foldentry.place(x=100, y=400)

        randomseedlabel = customtkinter.CTkLabel(
            master=create_dataset_frame, text="Random seed", font=("Century Gothic", 16)
        )
        randomseedlabel.place(x=200, y=400)
        self.randomseedentry = customtkinter.CTkEntry(
            master=create_dataset_frame,
            width=100,
            placeholder_text="random seed",
        )

        self.randomseedentry.place(x=330, y=400)
        randomseed_btn = customtkinter.CTkButton(
            master=create_dataset_frame,
            width=150,
            text="Generate random seed",
            command=lambda: self.generate_random_seed(self.randomseedentry),
        )
        randomseed_btn.place(x=450, y=400)

        self.trainids_btn = customtkinter.CTkButton(
            master=create_dataset_frame,
            width=220,
            text="Pick Train ids file",
            corner_radius=6,
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            command=self.browse_dataset,
        )

        self.trainids_entry = customtkinter.CTkEntry(
            master=create_dataset_frame,
            width=200,
            placeholder_text="train ids file path",
            textvariable=self.dataset_folder,
        )

        self.testids_btn = customtkinter.CTkButton(
            master=create_dataset_frame,
            width=220,
            text="Pick Test ids file",
            corner_radius=6,
            fg_color="transparent",
            border_width=1,
            border_color="#FFCC70",
            command=self.browse_dataset,
        )

        self.testids_entry = customtkinter.CTkEntry(
            master=create_dataset_frame,
            width=200,
            placeholder_text="train ids file path",
            textvariable=self.dataset_folder,
        )

        self.uniqueidlabel = customtkinter.CTkLabel(
            master=create_dataset_frame,
            text="Unique split id",
            font=("Century Gothic", 16),
        )
        self.uniqueidlabel.place(x=50, y=450)
        self.uniqueidentry = customtkinter.CTkEntry(
            master=create_dataset_frame,
            width=100,
            placeholder_text="Eg. split1",
        )

        self.uniqueidentry.place(x=180, y=450)
        self.split_btn = customtkinter.CTkButton(
            master=create_dataset_frame,
            width=150,
            text="Split dataset",
            command=self.split_dataset_thread,
        )
        self.split_btn.place(x=200, y=630)

        outputscreen = customtkinter.CTkFrame(
            tab_view.tab("split dataset"),
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
            font=custom_font,
        )

        self.scrollbar2 = tk.Scrollbar(outputscreen, command=self.cmd_output2.yview)
        self.scrollbar2.grid(row=0, column=1, sticky=tk.NSEW)
        self.cmd_output2["yscrollcommand"] = self.scrollbar2.set
        self.cmd_output2.grid(row=0, column=0, sticky=tk.NSEW, padx=2, pady=2)

        ##KOVER LEARN

        create_dataset_frame = customtkinter.CTkFrame(
            tab_view.tab("kover learn"),
            width=700,
            height=screen_height - 230,
            corner_radius=15,
            border_width=2,
        )
        create_dataset_frame.place(x=50, y=20)

        l1 = customtkinter.CTkLabel(
            master=create_dataset_frame, text="Kover learn", font=("Century Gothic", 20)
        )
        l1.place(x=50, y=45)

        # Selection option for dataset type
        dataset_type_label = customtkinter.CTkLabel(
            master=create_dataset_frame,
            text="Select kover learn Type",
            font=("Century Gothic", 15),
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
        self.pickdataset_btn = customtkinter.CTkButton(
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

        self.pickdataset_entry = customtkinter.CTkEntry(
            master=create_dataset_frame,
            width=380,
            placeholder_text="Dataset path",
            textvariable=self.selected_kover,
        )
        self.pickdataset_entry.place(x=280, y=200)

        self.pickoutput_btn = customtkinter.CTkButton(
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

        self.pickoutput_entry = customtkinter.CTkEntry(
            master=create_dataset_frame,
            width=380,
            placeholder_text="Output directory path",
            textvariable=self.output_dir,
        )
        self.pickoutput_entry.place(x=280, y=250)

        self.hyperparameterlabel = customtkinter.CTkLabel(
            master=create_dataset_frame,
            text="Hyper-parameter:",
            font=("Century Gothic", 12),
        )
        self.hyperparameterlabel.place(x=50, y=300)
        self.hyperparameterentry = customtkinter.CTkEntry(
            master=create_dataset_frame,
            width=300,
            placeholder_text="Enter hyperparamters eg: 0.1 0.178 0.316",
        )
        self.hyperparameterentry.place(x=180, y=300)

        self.modellabel = customtkinter.CTkLabel(
            master=create_dataset_frame, text="Model type:", font=("Century Gothic", 12)
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

        splitidentifierlabel = customtkinter.CTkLabel(
            master=create_dataset_frame,
            text="Split identifier",
            font=("Century Gothic", 12),
        )
        splitidentifierlabel.place(x=280, y=350)
        self.splitidentifierentry = customtkinter.CTkEntry(
            master=create_dataset_frame,
            width=200,
            placeholder_text="Enter split identifier in the dataset",
        )

        self.splitidentifierentry.place(x=370, y=350)

        self.maxrulelabel = customtkinter.CTkLabel(
            master=create_dataset_frame,
            text="Maximum rules :",
            font=("Century Gothic", 12),
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

        self.maxequivrulelabel = customtkinter.CTkLabel(
            master=create_dataset_frame,
            text="Maximum equivalent rules :",
            font=("Century Gothic", 12),
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

        hptypelabel = customtkinter.CTkLabel(
            master=create_dataset_frame, text="Hp choice :", font=("Century Gothic", 12)
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

        self.maxboundsize = customtkinter.CTkLabel(
            master=create_dataset_frame,
            text="Max bound genome size :",
            font=("Century Gothic", 12),
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

        randomseedlabel = customtkinter.CTkLabel(
            master=create_dataset_frame, text="Random seed", font=("Century Gothic", 12)
        )
        randomseedlabel.place(x=50, y=500)
        self.randomseedentry2 = customtkinter.CTkEntry(
            master=create_dataset_frame,
            width=100,
            placeholder_text="random seed",
        )

        self.randomseedentry2.place(x=180, y=500)
        randomseed_btn2 = customtkinter.CTkButton(
            master=create_dataset_frame,
            width=150,
            text="Generate random seed",
            command=lambda: self.generate_random_seed(self.randomseedentry2),
        )
        randomseed_btn2.place(x=300, y=500)

        self.koverlearn_btn = customtkinter.CTkButton(
            master=create_dataset_frame,
            width=150,
            text="kover learn",
            command=self.kover_learn_thread,
        )
        self.koverlearn_btn.place(x=200, y=630)

        outputscreen = customtkinter.CTkFrame(
            tab_view.tab("kover learn"),
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
            font=custom_font,
        )

        self.scrollbar3 = tk.Scrollbar(outputscreen, command=self.cmd_output3.yview)
        self.scrollbar3.grid(row=0, column=1, sticky=tk.NSEW)
        self.cmd_output3["yscrollcommand"] = self.scrollbar3.set
        self.cmd_output3.grid(row=0, column=0, sticky=tk.NSEW, padx=2, pady=2)

        # create fourth frame
        self.analysis_frame = customtkinter.CTkFrame(
            self, corner_radius=0, fg_color="transparent"
        )

        # select default frame
        self.select_frame_by_name(Page.DATA_COLLECTION_PAGE)

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
            except Exception as e:
                # Handle any exceptions that occur during the dataset creation process
                print(e.with_traceback())

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

    def select_tsv_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("TSV files", "*.tsv")])

        if file_path:
            self.file_name = os.path.basename(file_path)
            self.selected_file_path = file_path
            self.download_button["state"] = tk.NORMAL

    def create_contigs_folder(self):
        contigs_folder = "contigs"
        file_name_parts = self.file_name.split("_")
        datafolder = "data"
        if len(file_name_parts) > 1:
            tsv_folder_path = os.path.join(
                self.selected_directory, datafolder, file_name_parts[0]
            )
            file_name_parts[1] = file_name_parts[1].split(".")[0]
            tsv_folder_path = os.path.join(
                tsv_folder_path, file_name_parts[1], contigs_folder
            )
        else:
            tsv_folder_path = os.path.join(
                self.selected_directory,
                self.datafolder,
                *tsv_file_name_parts,
                contigs_folder,
            )

        os.makedirs(tsv_folder_path, exist_ok=True)

        contigs_folder = tsv_folder_path
        return contigs_folder

    def create_features_folder(self):
        features_folder = "features"
        file_name_parts = self.file_name.split("_")
        datafolder = "data"
        if len(file_name_parts) > 1:
            tsv_folder_path = os.path.join(
                self.selected_directory, datafolder, file_name_parts[0]
            )
            file_name_parts[1] = file_name_parts[1].split(".")[0]
            tsv_folder_path = os.path.join(
                tsv_folder_path, file_name_parts[1], features_folder
            )
        else:
            tsv_folder_path = os.path.join(
                self.selected_directory,
                datafolder,
                *tsv_file_name_parts,
                features_folder,
            )

        os.makedirs(tsv_folder_path, exist_ok=True)

        contigs_folder = tsv_folder_path
        return contigs_folder

    def get_genome_ids_from_file(self):
        genome_ids = []
        with open(self.selected_file_path, "r") as tsv_file:
            for line in tsv_file:
                columns = line.strip().split("\t")
                if len(columns) > 1:
                    genome_id = columns[1].replace(" ", "")  # Remove spaces
                    genome_ids.append(genome_id)
        return genome_ids

    def download_contigs(self):
        if hasattr(self, "selected_file_path") and self.selected_file_path:
            genome_ids = self.get_genome_ids_from_file()
        else:
            genome_ids = []

        if self.entry1.get():
            genome_ids.append("Genome id")
            genome_ids.append(self.entry1.get())

        if not genome_ids:
            messagebox.showerror("Error", "No genome IDs found.")
            return

        if not self.selected_directory:
            messagebox.showerror("Error", "Please select a directory for the download.")
            return

        self.download_button.configure(text="Downloading..", state=tk.DISABLED)
        self.cancel_button.configure(state=tk.NORMAL)
        self.progress_bar["value"] = 0
        self.size_label.configure(text="")
        self.downloading = True
        contigs_folder = self.create_contigs_folder()

        def download_genomes():
            total_genomes = len(genome_ids[1:])  # Start from the second index
            for index, genome_id in enumerate(genome_ids[1:]):
                if self.cancel_download:
                    break

                self.remote_path = f"genomes/{genome_id}/{genome_id}.fna"
                self.download_ftp(contigs_folder)
                progress_percentage = (index + 1) / total_genomes * 100
                self.progress_bar["value"] = progress_percentage
                self.size_label.configure(
                    text=f"Downloaded {index + 1}/{total_genomes} genomes"
                )

            self.download_button.configure(text="Download", state=tk.NORMAL)
            self.cancel_button.configure(state=tk.DISABLED)
            self.progress_bar[
                "value"
            ] = 100  # Set the progress bar to 100 when downloads are completed
            self.cancel_download = False
            self.downloading = False

        self.cancel_download = False

        # Create a separate thread for downloading
        self.download_thread = threading.Thread(target=download_genomes)
        self.download_thread.start()

    def cancelcontigs(self):
        if self.download_thread and self.download_thread.is_alive():
            if messagebox.askyesno(
                "Confirmation", "Are you sure you want to cancel the download?"
            ):
                self.cancel_download = True
                self.cancel_button.configure(state=tk.DISABLED)
                self.download_button.configure(text="Download", state=tk.NORMAL)
                self.progress_bar["value"] = 0
                self.size_label.configure(text="")

    def downloadfeatures(self):
        if hasattr(self, "selected_file_path") and self.selected_file_path:
            genome_ids = self.get_genome_ids_from_file()

        else:
            genome_ids = []

        if self.entry2.get():
            genome_ids.append("Genome id")
            genome_ids.append(self.entry2.get())

        if not genome_ids:
            messagebox.showerror("Error", "No genome IDs found.")
            return

        if not self.selected_directory:
            messagebox.showerror("Error", "Please select a directory for the download.")
            return

        self.download_button1.configure(text="Downloading..", state=tk.DISABLED)
        self.cancel_button1.configure(state=tk.NORMAL)
        self.progress_bar1["value"] = 0
        self.size_label1.configure(text="")
        self.downloading = True
        features_folder = self.create_features_folder()

        def download_features():
            total_genomes = len(genome_ids[1:])  # Start from the second index
            for index, genome_id in enumerate(genome_ids[1:]):
                if self.cancel_download:
                    break

                self.remote_path = (
                    "genomes/" + genome_id + "/" + genome_id + ".PATRIC.features.tab"
                )
                self.download_ftp(features_folder)
                progress_percentage = (index + 1) / total_genomes * 100
                self.progress_bar1["value"] = progress_percentage
                self.size_label1.configure(
                    text=f"Downloaded {index + 1}/{total_genomes} genomes"
                )

            self.download_button1.configure(text="Download", state=tk.NORMAL)
            self.cancel_button1.configure(state=tk.DISABLED)
            self.progress_bar1[
                "value"
            ] = 100  # Set the progress bar to 100 when downloads are completed
            self.cancel_download = False
            self.downloading = False

        self.cancel_download = False

        # Create a separate thread for downloading
        self.download_thread = threading.Thread(target=download_features)
        self.download_thread.start()

    def cancelfeatures(self):
        if self.download_thread and self.download_thread.is_alive():
            if messagebox.askyesno(
                "Confirmation", "Are you sure you want to cancel the download?"
            ):
                self.cancel_download = True
                self.cancel_button1.configure(state=tk.DISABLED)
                self.download_button1.configure(text="Download", state=tk.NORMAL)
                self.progress_bar1["value"] = 0
                self.size_label1.configure(text="")

    @threaded
    def load_amr_data(self):
        self.download_button4.configure(text="Loading..", state=tk.DISABLED)
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
                    "measurement_sign",
                    "measurement_value",
                    "measurement_unit",
                ],
                converters={
                    "genome_id": str,
                    "genome_name": lambda x: " ".join(x.lower().split()[:2]),
                },
            )

            self.amr_full = self.amr_full.drop_duplicates().dropna()
            self.amr_full = self.amr_full[self.amr_full["genome_name"] != ""]

            self.amr_list = pd.read_table(
                amr_metadata_file,
                usecols=["genome_name", "antibiotic"],
                converters={
                    "genome_name": lambda x: " ".join(x.lower().split()[:2]),
                },
            )

            self.amr_list = self.amr_list.drop_duplicates().dropna()
            self.amr_list = self.amr_list[self.amr_list["genome_name"] != ""]

            species_list = self.amr_list["genome_name"].unique()
            species_list.sort()
            self.species_filter["values"] = species_list.tolist()
            self.species_filter.configure(state=tk.DISABLED)
            self.update_table()
            self.species_filter.configure(state=tk.NORMAL)
        self.download_button4.configure(text="load amr list", state=tk.NORMAL)

    def update_table(self, event=None):
        selected_species = self.species_filter.get()
        self.start = 0
        if len(self.amr_list) > 0:
            if not selected_species:
                self.filtered_data = self.amr_list
            else:
                self.filtered_data = self.amr_list[
                    self.amr_list["genome_name"].str.contains(
                        pat=selected_species, case=False, regex=False
                    )
                    | self.amr_list["antibiotic"].str.contains(
                        pat=selected_species, case=False, regex=True
                    )
                ]

            data = self.filtered_data.values

            self.amr_list_table.delete(*self.amr_list_table.get_children())

            for name, antibiotic in data[:10]:
                self.amr_list_table.insert(
                    "",
                    tk.END,
                    values=(name, antibiotic),
                )

            self.total_label.config(text=f"Total: {len(self.filtered_data)}")

    def on_table_scroll(self, event):
        size = 10

        data = self.filtered_data.values

        if event.delta < 0 and self.start < len(data) - size:
            self.start += 1
            self.amr_list_table.delete(self.amr_list_table.get_children()[0])
            self.amr_list_table.insert(
                "",
                tk.END,
                values=(data[self.start + size - 1][0], data[self.start + size - 1][1]),
            )
        elif event.delta > 0 and self.start > 0:
            self.start -= 1
            self.amr_list_table.delete(self.amr_list_table.get_children()[-1])
            self.amr_list_table.insert(
                "",
                0,
                values=(data[self.start][0], data[self.start][1]),
            )

    def on_table_double_click(self, event):
        selected_item = self.amr_list_table.selection()

        if not selected_item:
            return

        item = selected_item  # Get the selected item
        selected_species: str = self.amr_list_table.item(item, "values")[0]
        selected_antibiotics: str = self.amr_list_table.item(item, "values")[1]
        selected_species = selected_species.replace("[", "").replace("]", "")
        # Update the entry box with the selected species

        self.species_entry.delete(0, tk.END)  # Clear the current entry
        self.antibiotic_entry.delete(0, tk.END)  # Clear the current entry

        self.species_entry.insert(0, selected_species)
        self.antibiotic_entry.insert(0, selected_antibiotics)

        self.update_amr_full()

        # You can now display all the antibiotics for the selected species
        self.display_antibiotics_for_species(selected_species)

    def display_antibiotics_for_species(self, selected_species):
        pass
        # Here, you can filter the dataset for the selected_species and display the antibiotics
        # You can use a new table or a label to display the antibiotics
        # Example:
        # antibiotics_for_species = self.get_antibiotics_for_species(selected_species)
        # print(antibiotics_for_species)
        # self.display_antibiotics(antibiotics_for_species)

    def get_antibiotics_for_species(self, species):
        pass
        # Filter the dataset to get antibiotics for the selected species
        # print(self.amr_list)
        # filtered_data = self.amr_list[
        #     self.amr_list["genome_name"].str.strip("").str.lower() == species
        # ]
        # antibiotics = filtered_data["antibiotic"].unique()
        # print(filtered_data)
        # print(antibiotics)

        # return filtered_data["antibiotic"].unique()

    def _remove_duplicates(self, data):
        # Keep only one measurement for the same genome, antibiotic and phenotype
        data = data.drop_duplicates(
            subset=["genome_id", "antibiotic", "resistant_phenotype"], keep="first"
        )
        # Drop all genomes/antibiotic combinations for which we have contradictory measurements
        data = data.drop_duplicates(subset=["genome_id", "antibiotic"], keep=False)
        return data

    def sort_amr_list_table(self, column_index: int):
        if len(self.filtered_data) == 0:
            return
        self.update_amr_list_column_sort(column_index)
        self.filtered_data.sort_values(
            by=[self.amr_list.columns[i] for i in self.amr_list_columns_sort_order],
            inplace=True,
            ascending=[self.sort_amr_list_ascending] * len(self.amr_list.columns),
        )
        self.update_table()

    def sort_results_table(self, column_index: int):
        if len(self.amr_full) == 0:
            return
        self.update_amr_full_column_sort(column_index)
        self.amr_full.sort_values(
            by=[self.amr_full.columns[i] for i in self.amr_full_columns_sort_order],
            inplace=True,
            ascending=[self.sort_amr_full_ascending] * len(self.amr_full.columns),
        )
        self.update_amr_full()

    def update_amr_full_column_sort(self, column_index: int):
        if self.amr_full_columns_sort_order[0] == column_index:
            self.sort_amr_full_ascending = not self.sort_amr_full_ascending
        else:
            self.sort_amr_full_ascending = True
        self.amr_full_columns_sort_order.remove(column_index)
        self.amr_full_columns_sort_order.insert(0, column_index)

    def update_amr_list_column_sort(self, column_index: int):
        if self.amr_list_columns_sort_order[0] == column_index:
            self.sort_amr_list_ascending = not self.sort_amr_list_ascending
        else:
            self.sort_amr_list_ascending = True
        self.amr_list_columns_sort_order.remove(column_index)
        self.amr_list_columns_sort_order.insert(0, column_index)

    def update_amr_full(self, event=None):
        antibiotic = self.antibiotic_entry.get()
        species = self.species_entry.get()

        if (antibiotic or species) and len(self.amr_full) > 0:
            try:
                if antibiotic:
                    amr = self.amr_full.loc[self.amr_full.antibiotic == antibiotic]

                if species:
                    species = species.replace("[", "").replace("]", "").lower()
                    amr = self.amr_full.loc[self.amr_full.amr["genome_name"] == species]

                amr = self._remove_duplicates(amr)
                amr = amr.loc[amr.resistant_phenotype != "Not defined"]

                if self.drop_intermediate.get():
                    amr = amr.loc[amr.resistant_phenotype != "Intermediate"]

                numeric_phenotypes = np.zeros(amr.shape[0], dtype=np.uint8)
                numeric_phenotypes[amr.resistant_phenotype.values == "Resistant"] = 1
                numeric_phenotypes[
                    amr.resistant_phenotype.values == "Non-susceptible"
                ] = 1
                numeric_phenotypes[
                    amr.resistant_phenotype.values == "Nonsusceptible"
                ] = 1
                numeric_phenotypes[amr.resistant_phenotype.values == "Intermediate"] = 2
                numeric_phenotypes[
                    amr.resistant_phenotype.values == "Susceptible-dose dependent"
                ] = 2

                total_resistant = np.sum(numeric_phenotypes == 1)
                total_susceptible = np.sum(numeric_phenotypes == 0)
                total = np.sum(numeric_phenotypes)

                self.totalresistance_label.configure(
                    text=f"Total resistance phenotypes available:{total_resistant}"
                )
                self.totalsusceptible_label.configure(
                    text=f"Total susceptible phenotypes available:{total_susceptible}"
                )

                amr["Measurement"] = (
                    amr["measurement_sign"].astype(str)
                    + amr["measurement_value"].astype(str)
                    + amr["measurement_unit"].astype(str)
                )

                self.results_table.delete(*self.results_table.get_children())
                for name, id, phenotype, measurement in zip(
                    amr["genome_name"].values,
                    amr["genome_id"].values,
                    numeric_phenotypes,
                    amr["Measurement"].values,
                ):
                    self.results_table.insert(
                        "",
                        tk.END,
                        values=(
                            "               " + name,
                            "                       " + id + "                       ",
                            phenotype,
                            "              " + measurement,
                        ),
                    )

                selected_items = self.results_table.get_children()
                number_of_rows = len(selected_items)
                self.totaldata.configure(text=f"Total Phenotypes:{number_of_rows}")
            except Exception:
                traceback.print_exc()
                messagebox.showerror("Error", "Error while reading the metadata file")

    def save_description_to_tsv(self, species, antibiotics, file_path):
        with open(file_path, "w", newline="") as tsv_file:
            tsv_writer = csv.writer(tsv_file)
            tsv_writer.writerow("Species: " + species)
            tsv_writer.writerow("Antibiotics: " + antibiotics)

    def save_all_table_to_tsv(self, treeview, file_path):
        # Get the column headings from the treeview
        columns = treeview["columns"]

        # Open the TSV file for writing
        with open(file_path, "w", newline="") as tsv_file:
            tsv_writer = csv.writer(tsv_file, delimiter="\t")

            # Write the column headings as the first row
            tsv_writer.writerow(columns)

            # Iterate through the treeview items and write their values to the TSV file
            for item in treeview.get_children():
                row_data = [treeview.item(item, "values")]
                tsv_writer.writerows(row_data)

    def save_table_to_tsv(self, treeview, file_path):
        # Open the TSV file for writing
        with open(file_path, "w", newline="") as tsv_file:
            tsv_writer = csv.writer(tsv_file, delimiter="\t")

            # Iterate through the treeview items and write the values of the second and third columns to the TSV file
            for item in treeview.get_children():
                values = treeview.item(item, "values")
                if values:
                    # Only include the second and third values in the row_data
                    row_data = [values[1], values[2]]
                    tsv_writer.writerow(row_data)

    def save_to_tsv(self):
        if not self.selected_directory:
            self.selected_directory = self.select_directory()
            print(self.selected_directory)
        species_text = self.species_entry.get().replace(" ", "-")
        antibiotic_text = self.antibiotic_entry.get()

        if species_text and antibiotic_text:
            file_name = f"{species_text}_{antibiotic_text}.tsv"
            species = species_text
            antibiotics = antibiotic_text
        elif species_text:
            file_name = f"{species_text}.tsv"
            species = species_text
            antibiotics = ""
        elif antibiotic_text:
            file_name = f"{antibiotic_text}.tsv"
            species = ""
            antibiotics = antibiotic_text
        else:
            # Prompt the user to enter either antibiotic or species
            # You can display an error message or take appropriate action here
            return

        file_name_parts = file_name.split("_")
        self.datafolder = "data"
        if len(file_name_parts) > 1:
            tsv_folder_path = os.path.join(
                self.selected_directory, self.datafolder, file_name_parts[0]
            )
            file_name_parts[1] = file_name_parts[1].split(".")[0]
            tsv_folder_path = os.path.join(tsv_folder_path, file_name_parts[1])
        else:
            file_name = file_name.split(".")[0]
            tsv_folder_path = os.path.join(
                self.selected_directory, self.datafolder, file_name
            )

        os.makedirs(tsv_folder_path, exist_ok=True)

        base_name, extension = os.path.splitext(file_name)
        phenometafile = f"{base_name}_phenotype_metadata{extension}"
        file_path = os.path.join(tsv_folder_path, phenometafile)
        self.save_table_to_tsv(self.results_table, file_path)
        descfilename = f"{base_name}_description{extension}"
        file_path = os.path.join(tsv_folder_path, descfilename)
        self.save_description_to_tsv(species, antibiotics, file_path)
        allfilename = f"{base_name}_full{extension}"
        file_path = os.path.join(tsv_folder_path, allfilename)
        self.save_all_table_to_tsv(self.results_table, file_path)
        # Schedule a function to change it back to the original text after 2000 milliseconds (2 seconds)

    def reset_button_text(self):
        # Change the button text back to the original text
        self.save_table_button.configure(text="Download data")

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

    def select_frame_by_name(self, page: Page):
        self.data_collection_frame.grid_forget()
        self.preprocessing_frame.grid_forget()
        self.kover_frame.grid_forget()
        self.analysis_frame.grid_forget()

        self.data_collection_frame_button.configure(fg_color="transparent")
        self.preprocessing_frame_button.configure(fg_color="transparent")
        self.kover_frame_button.configure(fg_color="transparent")
        self.analysis_frame_button.configure(fg_color="transparent")

        if page == Page.DATA_COLLECTION_PAGE:
            self.data_collection_frame.grid(row=0, column=1, sticky=tk.NSEW)
            self.data_collection_frame_button.configure(fg_color=("gray75", "gray25"))
        elif page == Page.PREPROCESSING_PAGE:
            self.preprocessing_frame.grid(row=0, column=1, sticky=tk.NSEW)
            self.preprocessing_frame_button.configure(fg_color=("gray75", "gray25"))
        elif page == Page.KOVER_LEARN_PAGE:
            self.kover_frame.grid(row=0, column=1, sticky=tk.NSEW)
            self.kover_frame_button.configure(fg_color=("gray75", "gray25"))
        elif page == Page.ANALYSIS_PAGE:
            self.analysis_frame.grid(row=0, column=1, sticky=tk.NSEW)
            self.analysis_frame_button.configure(fg_color=("gray75", "gray25"))
