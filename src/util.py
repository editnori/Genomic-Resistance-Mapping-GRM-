from subprocess import PIPE, Popen
from tkinter import messagebox, NORMAL, DISABLED, END, Spinbox
from customtkinter import CTkEntry, CTkTextbox
from threading import Thread
from traceback import print_exc
import os
import re
from concurrent.futures import Future
from typing import IO, Iterable, Optional
from queue import Queue, Empty

from customtkinter import filedialog


CRLF = b"\r\n"
LF = b"\n"


class Tag(str):
    ERROR = "error"
    SUCCESS = "success"
    NORMAL = None


class Key(int):
    ENTER = 13
    SPACE = 32
    ESCAPE = 27


def select_directory(title: Optional[str] = "Select Folder") -> Optional[str]:
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


def select_file(
    filetypes: Optional[Iterable[str]] = (("All Files", "*.*"),),
    title: Optional[str] = "Open",
) -> Optional[str]:
    selected_file: str = filedialog.askopenfilename(filetypes=filetypes, title=title)

    if not selected_file:
        return None

    if not os.path.isfile(selected_file):
        messagebox.showerror(
            "Error", "File path is invalid.\n\nPlease select a valid file path."
        )
        return None

    if (
        selected_file.find(" ") != -1
        or selected_file.find("(") != -1
        or selected_file.find(")") != -1
    ):
        messagebox.showerror(
            "Error",
            "File path is invalid.\n\nPlease select a file path without spaces or parentheses.",
        )
        return None

    return selected_file


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


def to_linux_path(path: str) -> str:
    path = os.path.abspath(path)
    path = path.replace(path[0], path[0].lower(), 1)
    path = re.sub(r"(\w)(:\\)", r"/mnt/\1/", path)
    path = path.replace("\\", "/")
    return path


def try_pass_except(func, *args, **kwargs):
    try:
        func(*args, **kwargs)
    except:
        pass


def run_bash_command(command: str) -> Optional[Popen]:
    temp_file = "tmp.sh"

    try:
        with open(temp_file, "wb") as bash_file:
            bash_file.write(
                f'#!/bin/bash\n{command}\nrm "$0"'.encode("UTF-8").replace(CRLF, LF)
            )  # scary rm command
    except Exception as e:
        messagebox.showerror(
            "Error", f"An error occurred while creating the bash file.\n\n{e}"
        )

        return None

    process = Popen(
        f'wsl -e "{to_linux_path(temp_file)}"',
        universal_newlines=True,
        stdout=PIPE,
        stderr=PIPE,
    )

    Thread(
        target=lambda: (process.wait(), try_pass_except(os.remove, "tmp.sh"))
    ).start()

    return process


@threaded
def enqueue_output(out: IO, queue: Queue, tag: Optional[str] = None):
    """Credit: https://stackoverflow.com/a/4896288"""
    for line in out:
        queue.put((tag, line))
    out.close()


def sanitize_filename(filename: str) -> str:
    return filename.replace("/", "-").replace(" ", "-").replace(".", "-")


def force_insertable_value(
    new_value: float,
    spinbox: Spinbox | CTkEntry,
):
    validation = spinbox.cget("validate")
    spinbox.configure(validate="none")
    spinbox.delete(0, END)
    spinbox.insert(0, new_value)
    spinbox.configure(validate=validation)


def update_cmd_output(message: str, output_target: CTkTextbox, *tags: str):
    output_target.configure(state=NORMAL)
    is_at_end = output_target.yview()[1] > 0.95
    output_target.insert(END, message, tags)
    if is_at_end:
        output_target.see(END)
    output_target.configure(state=DISABLED)
    output_target.update_idletasks()


def display_process_output(process: Popen, output_target: CTkTextbox = None):
    messages = Queue()

    enqueue_output(process.stdout, messages, Tag.NORMAL)
    enqueue_output(process.stderr, messages, Tag.ERROR)

    has_messages = True

    while (process.poll() is None) or has_messages:
        try:
            has_messages = True

            tag, message = messages.get(timeout=0.1)

            if output_target:
                update_cmd_output(message, output_target, tag)
            else:
                print(f"{tag}: {message}", end="")
        except Empty:
            has_messages = False
