"""ActionLabel - a label with a link, for consistent look and feel"""

from PySide6.QtCore import Signal, Qt
from PySide6.QtWidgets import QLabel


class ActionLabel(QLabel):
    """ Label that prompts for an optional action - NOT VISIBLE BY DEFAULT"""

    action = Signal()

    def __init__(self, main_text: str, action_text: str = "click"):
        super().__init__(main_text + f' <a href="do_action"><font color="orange">{action_text}</font></a>')

        self.setStyleSheet("color: red;")

        self.setWordWrap(True)

        self.setOpenExternalLinks(False)  # Don't open in browser
        self.setTextInteractionFlags(self.textInteractionFlags() |
                                           Qt.LinksAccessibleByMouse)

        self.linkActivated.connect(self._on_link_activated)
        self.setVisible(False)

    def _on_link_activated(self, href: str):
        self.action.emit()
