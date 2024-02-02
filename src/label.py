from typing import Any, Optional, Tuple, Union

from ctk import CTkLabel, CTkImage, CTkFont
from hovertip import Hovertip


class Label(CTkLabel):
    """wrapper for CTkLabel that adds a hovertip and truncates text if it is too long
    credit to https://stackoverflow.com/a/51144221"""

    def __init__(
        self,
        master: Any,
        width: int = 0,
        height: int = 28,
        corner_radius: Optional[int] = None,
        bg_color: Union[str, Tuple[str, str]] = "transparent",
        fg_color: Optional[Union[str, Tuple[str, str]]] = None,
        text_color: Optional[Union[str, Tuple[str, str]]] = None,
        text_color_disabled: Optional[Union[str, Tuple[str, str]]] = None,
        text: str = "",
        font: Optional[Union[tuple, CTkFont]] = None,
        image: Union[CTkImage, None] = None,
        compound: str = "center",
        anchor: str = "center",  # label anchor: center, n, e, s, w
        wraplength: int = 0,
        **kwargs
    ):
        super().__init__(
            master,
            width,
            height,
            corner_radius,
            bg_color,
            fg_color,
            text_color,
            text_color_disabled,
            text,
            font,
            image,
            compound,
            anchor,
            wraplength,
            **kwargs
        )

        self.width = width

        self.hovertip = Hovertip(self, text)

        self.hovertip.disable()

        self.configure(text=text)

    def configure_text(self, text: str) -> str:
        self.hovertip.text = text

        if text and hasattr(self, "width") and len(text) > self.width:
            self.hovertip.enable()
            text = text[: self.width - 3] + "..."
        else:
            self.hovertip.disable()
        return text

    def configure(self, require_redraw=False, **kwargs):
        if "text" in kwargs:
            kwargs["text"] = self.configure_text(kwargs["text"])
        super().configure(require_redraw, **kwargs)

    def cget(self, option):
        if option == "text":
            return self.hovertip.text
        return super().cget(option)
