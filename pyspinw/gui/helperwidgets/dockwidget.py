""" DockWidget Base Class """

from PySide6.QtGui import QCloseEvent, Qt
from PySide6.QtWidgets import QDockWidget


class SpinWDockWidget(QDockWidget):
    """ Base class for spinW DockWidgets"""

    def __init__(self, parent=None):
        super().__init__(parent=parent)

        self.setAllowedAreas(Qt.DockWidgetArea.LeftDockWidgetArea | Qt.DockWidgetArea.RightDockWidgetArea)


    def closeEvent(self, event: QCloseEvent):
        """Qt Override: hide the window, don't close it"""
        self.setVisible(False)
        event.ignore()
