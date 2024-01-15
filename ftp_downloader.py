from datetime import datetime
import time
import tkinter as tk
from tkinter import messagebox, filedialog
from ftplib import FTP
from threading import Thread
from concurrent.futures import Future
import os
from tkinter.ttk import Progressbar
import traceback
from typing import Optional

from customtkinter import CTkButton, CTkLabel, CTkFrame


def threaded(fn):
    """@threaded decorator from https://stackoverflow.com/a/19846691"""

    def call_with_future(fn, future, args, kwargs):
        try:
            result = fn(*args, **kwargs)
            future.set_result(result)
        except Exception as exc:
            future.set_exception(exc)

    def wrapper(*args, **kwargs):
        future = Future()
        Thread(target=call_with_future, args=(fn, future, args, kwargs)).start()
        return future

    return wrapper


@threaded
def get_last_metadata_update_date(label: CTkLabel):
    label.configure(text="Last Update Time: Loading...")
    try:
        ftps = FTP("ftp.bvbrc.org")
        ftps.login()
        mod_time = ftps.sendcmd(
            "MDTM " + "/RELEASE_NOTES/PATRIC_genomes_AMR.txt"
        ).split()[1]
        last_update_time = datetime.strptime(mod_time, "%Y%m%d%H%M%S")
        label.configure(text=f"Last Update Time: {last_update_time}")
    except Exception:
        print(traceback.format_exc())  # LOGGER


class DownloadWindow:
    def __init__(
        self,
        frame: CTkFrame,
        download_button: Optional[CTkButton],
        cancel_button: Optional[CTkButton],
        select_path_button: Optional[CTkButton],
        path_label: Optional[CTkLabel],
        progress_bar: Optional[Progressbar],
        size_label: Optional[CTkLabel],
    ):
        self._capacity: int = 1
        self._is_downloading: bool = False
        self._last_update_time: float = 0

        self._frame = frame
        self._download_button = download_button
        self._cancel_button = cancel_button
        self._select_path_button = select_path_button
        self._path_label = path_label
        self._progress_bar = progress_bar
        self._size_label = size_label

    def set_downloading(self, downloading: bool):
        if downloading:
            if self._download_button:
                self._download_button.configure(state=tk.DISABLED)
            if self._cancel_button:
                self._cancel_button.configure(state=tk.NORMAL)
            if self._select_path_button:
                self._select_path_button.configure(state=tk.DISABLED)
            self._is_downloading = True
        else:
            if self._download_button:
                self._download_button.configure(state=tk.NORMAL)
            if self._cancel_button:
                self._cancel_button.configure(state=tk.DISABLED)
            if self._select_path_button:
                self._select_path_button.configure(state=tk.NORMAL)
            self._is_downloading = False

    def is_downloading(self) -> bool:
        return self._is_downloading

    def set_path(self, path: str):
        if self._path_label:
            self._path_label.configure(text=path)

    def get_path(self) -> Optional[str]:
        if self._path_label and self._path_label.cget("text"):
            return self._path_label.cget("text")
        return None

    def set_capacity(self, capacity: int):
        self._capacity = capacity

    def set_progress_label(self, size: int):
        current_time_ms = time.time_ns() / 1_000_000
        if (current_time_ms - self._last_update_time) > 100:
            if self._size_label:
                self._size_label.configure(
                    text=f"Downloading: {size / 1_048_576:.2f} MB / {self._capacity / 1_048_576:.2f} MB"
                )
            if self._progress_bar:
                self._progress_bar["value"] = size / self._capacity * 100
            self._last_update_time = current_time_ms
        elif size == self._capacity:
            if self._size_label:
                self._size_label.configure(
                    text=f"Downloading: complete ({size / 1_048_576:.2f} MB)"
                )
            if self._progress_bar:
                self._progress_bar["value"] = 0
            self._last_update_time = current_time_ms

    def set_buttons_commands(
        self,
        download_command=None,
        cancel_command=None,
        select_path_command=None,
    ):
        if self._download_button:
            self._download_button.configure(command=download_command)
        if self._cancel_button:
            self._cancel_button.configure(command=cancel_command)
        if self._select_path_button:
            self._select_path_button.configure(command=select_path_command)


class FTPDownloadApp:
    def __init__(self, download_window: DownloadWindow, remote_path: str):
        self.local_path: Optional[str] = None
        self.selected_directory: Optional[str] = None
        self.download_thread: Optional[Thread] = None

        self.download_window = download_window
        self.remote_path = remote_path

        self.download_window.set_downloading(False)

        self.download_window.set_buttons_commands(
            download_command=self.download,
            cancel_command=self.cancel,
            select_path_command=self.select_directory,
        )

    @threaded
    def download(self):
        try:
            self.download_window.set_downloading(True)
            filename = self.remote_path.rsplit("/", 2)[-1]

            if not self.download_window.get_path():
                raise IOError("No download path selected.")
            else:
                selected_path = self.download_window.get_path()

            local_path = os.path.join(selected_path, filename)

            with open(local_path, "wb") as local_file:
                ftp = FTP("ftp.bvbrc.org")
                ftp.login()
                ftp.voidcmd("TYPE I")
                self.download_window.set_capacity(ftp.size(self.remote_path))
                with ftp.transfercmd(f"RETR {self.remote_path}") as conn:
                    bytes_received = 0
                    while self.download_window.is_downloading() and (
                        data := conn.recv(1024)
                    ):
                        local_file.write(data)
                        bytes_received += len(data)
                        self.download_window.set_progress_label(bytes_received)
        except IOError:
            messagebox.showerror("Error", "Please change the download path.")
            print(traceback.format_exc())  # LOGGER
            self.download_window.set_downloading(False)
        except Exception:
            print(traceback.format_exc())  # LOGGER
            self.download_window.set_downloading(False)
        finally:
            if ftp:
                ftp.quit()
            self.download_window.set_downloading(False)

    def cancel(self):
        if messagebox.askyesno(
            "Confirmation", "Are you sure you want to cancel the download?"
        ):
            self.download_window.set_downloading(False)

    def select_directory(self):
        self.download_window.set_path(filedialog.askdirectory())
