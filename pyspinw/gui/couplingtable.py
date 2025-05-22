class CouplingTable(QTableWidget):
    def __init__(self):
        self._couplings = []

        self.setRowCount(0)
        self.setColumnCount(14)

        self.setHorizontalHeaderLabels(["Name",
                                        "Site 1",
                                        "Site 2",
                                        "Type",
                                        ])