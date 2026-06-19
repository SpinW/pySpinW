""" Utility functions and objects """

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor
from PySide6.QtWidgets import QLabel, QColorDialog


class QRightLabel(QLabel):
    """ Label, but default right aligned"""

    def __init__(self, test: str, parent=None):
        super().__init__(test, parent=parent)
        self.setAlignment(Qt.AlignRight)


class QColorPatch(QLabel):
    """ A patch of color"""

    colorChanged = Signal()

    def __init__(self, r: float, g: float, b: float, parent=None):
        super().__init__("          ", parent=parent)

        self.parent_widget=parent

        self.rgb = r, g, b
        self.setColor(r, g, b)

    def setColor(self, r: float, g: float, b: float):
        """ Set the color of this patch """
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
        """ Qt override for mouse click """
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
