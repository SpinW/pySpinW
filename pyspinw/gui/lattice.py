from vispy.app.qt import QWidget

class LatticeParameters(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.a =