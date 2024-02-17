from app import App, Thread, ApartmentState, ThreadStart
import os


def main():
    app = App()
    app.mainloop()


t = Thread(ThreadStart(main))
t.ApartmentState = ApartmentState.STA
t.Start()
t.Join()
os._exit(0)