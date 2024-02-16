from app import App
from System.Threading import Thread, ApartmentState, ThreadStart


def main():
    app = App()
    app.mainloop()


t = Thread(ThreadStart(main))
t.ApartmentState = ApartmentState.STA
t.Start()
t.Join()
main()
