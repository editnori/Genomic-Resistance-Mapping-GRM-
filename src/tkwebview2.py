from tkinter import Frame
from tkinter import filedialog
import ctypes
import os
from threading import Semaphore
import clr
from webview.window import Window
from webview.platforms.edgechromium import EdgeChrome

clr.AddReference("System.Windows.Forms")
clr.AddReference("System.Threading")
from System.Windows.Forms import Control
from System.Threading import SynchronizationContext, SendOrPostCallback

user32 = ctypes.windll.user32


class WebView2(Frame):
    def __init__(self, parent, width: int, height: int, url: str = "", **kw):
        Frame.__init__(self, parent, width=width, height=height, **kw)
        control = Control()
        window = Window(
            "master",
            str(id(self)),
            url=None,
            html=None,
            js_api=None,
            width=width,
            height=height,
            x=None,
            y=None,
            resizable=True,
            fullscreen=False,
            min_size=(200, 100),
            hidden=False,
            frameless=False,
            easy_drag=True,
            minimized=False,
            on_top=False,
            confirm_close=False,
            background_color="#FFFFFF",
            transparent=False,
            text_select=True,
            localization=None,
            zoomable=True,
            draggable=True,
            vibrancy=False,
        )
        self.window = window
        self.web_view = EdgeChrome(control, window, None)
        self.control = control
        self.web = self.web_view.web_view
        self.width = width
        self.height = height
        self.parent = parent
        self.chwnd = int(str(self.control.Handle))
        user32.SetParent(self.chwnd, self.winfo_id())
        user32.MoveWindow(self.chwnd, 0, 0, width, height, True)
        self.loaded = window.events.loaded
        self.__go_bind()
        if url != "":
            self.load_url(url)
        self.core = None
        self.web.CoreWebView2InitializationCompleted += self.__load_core

    def __go_bind(self):
        self.bind("<Destroy>", lambda event: self.web.Dispose())
        self.bind("<Configure>", self.__resize_webview)

    def __resize_webview(self, event):
        user32.MoveWindow(
            self.chwnd, 0, 0, self.winfo_width(), self.winfo_height(), True
        )

    def __load_core(self, sender, _):
        self.core = sender.CoreWebView2
        self.core.NewWindowRequested -= self.web_view.on_new_window_request
        settings = sender.CoreWebView2.Settings
        settings.AreDefaultContextMenusEnabled = True
        settings.AreDevToolsEnabled = True

    def __download_file(self, sender, args):
        def UpdateProgress(download):
            def state(sender, e):
                s = download.State.ToString()
                if s == 0:
                    return
                elif s == "Interrupted":
                    print(args.DownloadOperation.InterruptReason)
                elif s == "Completed":
                    print("Complete")

            download.StateChanged += state

        def __dd(e):
            path = filedialog.askdirectory(initialdir=args.ResultFilePath, title="下载")
            fname = os.path.basename(args.DownloadOperation.ResultFilePath)
            if path == None:
                args.Cancel = True
            else:
                self.core.OpenDefaultDownloadDialog()
                args.ResultFilePath = path + "/" + fname
                UpdateProgress(args.DownloadOperation)
                print("ok")

        SynchronizationContext.Current.Post(SendOrPostCallback(__dd), None)

    def get_url(self):
        return self.web_view.get_current_url()

    def evaluate_js(self, script, callback=None):
        js_r = []
        semaphore = Semaphore(1)
        if callback != None:
            return self.web_view.evaluate_js(script, semaphore, js_r, callback)
        else:
            return self.web_view.evaluate_js(script, semaphore, js_r)

    def load_css(self, css):
        self.web_view.load_css(css)

    def load_html(self, content, base_uri=None):
        self.web_view.load_html(content, base_uri)

    def load_url(self, url):
        self.web_view.load_url(url)

    def reload(self):
        self.core.Reload()

    def event_core_completed(self, command):
        self.web.CoreWebView2InitializationCompleted += command


def have_runtime():
    from webview.platforms.winforms import _is_chromium

    return _is_chromium()


def install_runtime():
    from urllib import request
    import subprocess
    import os

    url = r"https://go.microsoft.com/fwlink/p/?LinkId=2124703"
    path = os.getcwd() + "\\webview2runtimesetup.exe"
    unit = request.urlopen(url).read()
    with open(path, mode="wb") as uf:
        uf.write(unit)
    cmd = path
    p = subprocess.Popen(cmd, shell=True)
    return_code = p.wait()
    os.remove(path)
    return return_code
