"""Utility functions and objects."""

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor
from PySide6.QtWidgets import QLabel, QColorDialog


class QRightLabel(QLabel):
    """Label, but default right aligned.

    Parameters
    ----------
    text : str
        Text label to display.
    parent : QWidget or None, optional
        Parent widget. Default is ``None``.
    """

    def __init__(self, text: str, parent=None):
        super().__init__(text, parent=parent)
        self.setAlignment(Qt.AlignRight)


class QColorPatch(QLabel):
    """A patch of color.

    Parameters
    ----------
    r : float
        Red channel value in ``[0, 1]``.
    g : float
        Green channel value in ``[0, 1]``.
    b : float
        Blue channel value in ``[0, 1]``.
    parent : QWidget or None, optional
        Parent widget. Default is ``None``.
    """

    colorChanged = Signal()

    def __init__(self, r: float, g: float, b: float, parent=None):
        super().__init__("          ", parent=parent)

        self.parent_widget=parent

        self.rgb = r, g, b
        self.setColor(r, g, b)

    def setColor(self, r: float, g: float, b: float):
        """Set the color of this patch.

        Parameters
        ----------
        r : float
            Red channel value in ``[0, 1]``.
        g : float
            Green channel value in ``[0, 1]``.
        b : float
            Blue channel value in ``[0, 1]``.
        """
        self.rgb = r, g, b

        r = int(r * 255)
        g = int(g * 255)
        b = int(b * 255)

        self.setStyleSheet(f"""
            QLabel {{
                background-color: rgb({r}, {g}, {b});
                color: white;
                border: 2px solid black;
            }}
            
            QLabel:hover {{
                border: 2px solid white;
            }}
        """) #noqa: W293

    def mousePressEvent(self, event):
        """Qt override for mouse click."""
        if event.button() == Qt.LeftButton:
            color = QColorDialog.getColor(
                initial=QColor.fromRgbF(*self.rgb),
                parent=self.parent_widget,
                title="Select Color"
            )
            if color.isValid():
                self.setColor(color.redF(), color.greenF(), color.blueF())
                self.colorChanged.emit()
            else:
                pass
