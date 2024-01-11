
import tkinter as tk
from tkinter import messagebox, filedialog
from ftplib import FTP
import threading
import os

class FTPDownloadApp:
    def __init__(self, root, remote_path, download_button, cancel_button, select_dir_button, progress_bar, size_label, local_path_label):
        self.root = root
        #self.root.title("FTP File Downloader")
        
        self.ftp = None
        self.remote_path = remote_path
        self.local_path = None
        self.selected_directory = None
        
        
        
        self.download_button = download_button
        self.cancel_button = cancel_button
        self.select_dir_button = select_dir_button
        self.progress_bar = progress_bar
        self.size_label = size_label
        self.local_path_label = local_path_label
        
        self.download_button.configure(command=self.download)
        self.cancel_button.configure(command=self.cancel, state=tk.DISABLED)
        self.select_dir_button.configure(command=self.select_directory)
        
        self.download_thread = None
        
    def download(self):
       
        if not self.selected_directory:
            messagebox.showerror("Error", "Please select a directory for the download.")
            return
        
        self.download_button.configure(state=tk.DISABLED)
        self.cancel_button.configure(state=tk.NORMAL)
        
        def download_ftp():
            try:
                self.ftp = FTP("ftp.bvbrc.org")
                self.ftp.login()
                self.ftp.voidcmd('TYPE I')
                total_size = self.ftp.size(self.remote_path)
                total_size_mb = total_size / (1024 * 1024)
                self.size_label.configure(text=f"Downloading: 0 MB / {total_size_mb:.2f} MB")
                filename = self.remote_path.split("/")[-1]
                self.local_path = os.path.join(self.selected_directory, filename)
                
                with self.ftp.transfercmd("RETR " + self.remote_path) as conn:
                    with open(self.local_path, "wb") as local_file:
                        bytes_received = 0
                        while True:
                            if self.cancel_download:
                                break
                            data = conn.recv(1024)
                            if not data:
                                break
                            local_file.write(data)
                            bytes_received += len(data)
                            self.root.update_idletasks()
                            self.progress_bar["value"] = (bytes_received / total_size) * 100
                            self.size_label.configure(text=f"Downloading: {bytes_received / (1024 * 1024):.2f} MB / {total_size_mb:.2f} MB")
                
                if not self.cancel_download:
                    self.local_path_label.configure(text=self.local_path)
            except Exception as e:
                print("Error:", e)
            finally:
                if self.ftp:
                    self.ftp.quit()
                self.download_button.configure(state=tk.NORMAL)
                self.cancel_button.configure(state=tk.DISABLED)
                self.cancel_download = False

        self.cancel_download = False
        self.download_thread = threading.Thread(target=download_ftp)
        self.download_thread.start()
    
    def cancel(self):
        if self.download_thread and self.download_thread.is_alive():
            if messagebox.askyesno("Confirmation", "Are you sure you want to cancel the download?"):
                self.cancel_download = True
                self.cancel_button.configure(state=tk.DISABLED)
                self.download_button.configure(state=tk.NORMAL)
    
    def select_directory(self):
        self.selected_directory = filedialog.askdirectory()


