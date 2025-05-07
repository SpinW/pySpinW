from PySide6.QtGui import QCloseEvent, Qt
from PySide6.QtWidgets import QDockWidget


class SpinWDockWidget(QDockWidget):

    def __init__(self, parent=None):
        super().__init__(parent=parent)

        self.setAllowedAreas(Qt.DockWidgetArea.LeftDockWidgetArea | Qt.DockWidgetArea.RightDockWidgetArea)


    def closeEvent(self, event: QCloseEvent):
        self.setVisible(False)
        event.ignore()