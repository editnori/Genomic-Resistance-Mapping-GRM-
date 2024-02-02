from tkinter import ttk, END
import pandas as pd


class Key(int):
    UP = 38
    DOWN = 40


class Table(ttk.Treeview):
    def __init__(self, master, *args, **kwargs):
        super().__init__(master, *args, **kwargs)

        self.tag_configure("all", background="#313131")

        self.load_start: int = 0
        self.load_size: int = kwargs.get("height", 0) + 1
        self.columns: tuple[str] = tuple(kwargs.get("columns", tuple()))
        self.columns_order: list[str] = list(range(len(self.columns)))
        self.sort_ascending = False
        self.filtered_data: pd.DataFrame = pd.DataFrame()

        self._mouse_pressed_before_focus = False

        post_class_id = f"post-class-{int(id(self))}"

        bindtags = list(self.bindtags()[:-1])
        bindtags.insert(2, post_class_id)
        self.bindtags(bindtags)

        self.event_add("<<maneuver>>", "<MouseWheel>", "<Down>", "<Up>")
        self.event_add("<<focus>>", "<ButtonPress>", "<FocusIn>")

        self.bind_class(post_class_id, "<<maneuver>>", self._on_maneuver)

        self.bind("<ButtonPress>", self._mouse_pressed)
        self.bind("<FocusIn>", self._select_first)
        self.bind("<FocusOut>", self.selection_clear)

    def _mouse_pressed(self, event=None):
        self._mouse_pressed_before_focus = True

    def _select_first(self, event=None):
        if self._mouse_pressed_before_focus:
            self._mouse_pressed_before_focus = False
            return

        children = super().get_children()
        if len(children) > 0:
            super().selection_set(children[0])
            super().focus(children[0])

    def selection_clear(self, event=None):
        for i in super().selection():
            super().selection_remove(i)

    def reset_view(self):
        super().delete(*super().get_children())

        self.load_start = 0

        if len(self.filtered_data) == 0:
            return

        for items in self.filtered_data.values[: self.load_size]:
            super().insert(
                "",
                END,
                values=list(items),
                tags=("all",),
            )

    def sort_by_column(self, column_index: int):
        if len(self.filtered_data) == 0:
            return

        self.sort_ascending = (
            not self.sort_ascending if self.columns_order[0] == column_index else True
        )

        self.columns_order.remove(column_index)
        self.columns_order.insert(0, column_index)

        columns = [self.filtered_data.columns[i] for i in self.columns_order]

        self.filtered_data = self.filtered_data.sort_values(
            by=columns,
            ascending=[self.sort_ascending] * len(columns),
        )

        self.reset_view()

    def heading(self, column, option=None, **kw):
        super().heading(
            column,
            option,
            command=lambda: self.sort_by_column(self.columns.index(column)),
            **kw,
        )

    def _on_maneuver(self, event):
        children = super().get_children()

        if len(children) == 0:
            return

        selected_index = self.index(super().focus())

        if (event.keycode == Key.DOWN or event.keycode == Key.UP) and (
            0 < selected_index < self.load_size - 1
        ):
            return

        if (event.delta < 0 or event.keycode == Key.DOWN) and self.load_start < len(
            self.filtered_data
        ) - self.load_size:
            self.load_start += 1

            super().delete(children[0])
            super().insert(
                "",
                END,
                values=(
                    list(
                        self.filtered_data.values[self.load_start + self.load_size - 1]
                    )
                ),
                tags=("all",),
            )

        elif (event.delta > 0 or event.keycode == Key.UP) and self.load_start > 0:
            self.load_start -= 1

            super().delete(children[-1])
            super().insert(
                "",
                0,
                values=list(self.filtered_data.values[self.load_start]),
                tags=("all",),
            )
