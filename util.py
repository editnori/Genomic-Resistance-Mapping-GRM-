from subprocess import PIPE, Popen
from tkinter import messagebox
from threading import Thread
from traceback import print_exc
import os
import re
from concurrent.futures import Future
from typing import IO, Optional

from customtkinter import filedialog


CRLF = b"\r\n"
LF = b"\n"


def select_directory(title: Optional[str] = None) -> Optional[str]:
    selected_directory: str = filedialog.askdirectory(title=title)

    if not selected_directory:
        return None

    if not os.path.isdir(selected_directory):
        messagebox.showerror(
            "Error", "Directory is invalid.\n\nPlease select a valid directory."
        )
        return None

    if (
        selected_directory.find(" ") != -1
        or selected_directory.find("(") != -1
        or selected_directory.find(")") != -1
    ):
        messagebox.showerror(
            "Error",
            "Directory path is invalid.\n\nPlease select a directory path without spaces or parentheses.",
        )
        return None

    return selected_directory


def threaded(fn) -> Future:
    """@threaded decorator from https://stackoverflow.com/a/19846691"""

    def call_with_future(fn, future, args, kwargs):
        try:
            result = fn(*args, **kwargs)
            future.set_result(result)
        except Exception as exc:
            print_exc()  # LOGGER
            future.set_exception(exc)

    def wrapper(*args, **kwargs):
        future = Future()
        Thread(target=call_with_future, args=(fn, future, args, kwargs)).start()
        return future

    return wrapper


@threaded
def read_output(io: IO, *tags: str) -> list[str]:
    lines = []

    for line in io:
        lines.append(line)

    return lines


def to_linux_path(path: str) -> str:
    path = os.path.abspath(path)
    path = path.replace(path[0], path[0].lower(), 1)
    path = re.sub(r"(\w)(:\\)", r"/mnt/\1/", path)
    path = path.replace("\\", "/")
    return path


def run_bash_command(command: str, temp_directory: str) -> Optional[Popen]:
    try:
        temp_file = os.path.join(temp_directory, "tmp.sh")
        with open(temp_file, "wb") as bash_file:
            bash_file.write(
                f'#!/bin/bash\n{command}\nrm "$0"'.encode("UTF-8").replace(CRLF, LF)
            )  # scary rm command
    except Exception as e:
        messagebox.showerror(
            "Error", f"An error occurred while creating the bash file.\n\n{e}"
        )

        return None

    return Popen(
        f'wsl -e "{to_linux_path(temp_file)}"',
        shell=True,
        universal_newlines=True,
        stdout=PIPE,
        stderr=PIPE,
    )
