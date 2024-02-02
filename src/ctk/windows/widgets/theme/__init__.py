from .theme_manager import ThemeManager

# load default blue theme
try:
    ThemeManager.load_theme("blue")
except FileNotFoundError as err:
    raise FileNotFoundError(f"{err}\nThe .json theme file for ctk could not be found.\n" +
                            f"If packaging with pyinstaller was used, have a look at the wiki:\n" +
                            f"https://github.com/TomSchimansky/ctk/wiki/Packaging#windows-pyinstaller-auto-py-to-exe")
