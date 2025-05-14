from PySide6.QtCore import Qt
from PySide6.QtWidgets import QLabel


class QRightLabel(QLabel):
    """ Right aligned label"""
    def __init__(self, text: str):
        super().__init__(text)
        self.setAlignment(Qt.AlignRight)