from PySide6.QtWidgets import QWidget


def layout_debug(widget: QWidget):
    """ Set the stylesheet on the given widget to make it easy to debug the layout"""
    widget.setStyleSheet("""
                        * {
                            border: 1px solid red;
                            background-color: rgba(255, 0, 0, 0.05);
                        }
                        QLabel {
                            border: 1px solid blue;
                            background-color: rgba(0, 0, 255, 0.1);
                        }
                        QLineEdit {
                            border: 1px solid orange;
                            background-color: rgba(255, 165, 0, 0.1);
                        }
                        QPushButton {
                            border: 1px solid green;
                            background-color: rgba(0, 255, 0, 0.1);
                        }
                        """)
