""" Loading of icons with dark mode or light mode """
from importlib import resources

from PySide6.QtCore import QByteArray, QBuffer, QRectF
from PySide6.QtGui import QIcon, Qt, QGuiApplication
from PySide6.QtSvg import QSvgRenderer
from PySide6.QtGui import QPixmap, QPainter

def pixmap_from_package(svg_package: str, svg_name: str, size: int) -> QPixmap:
    # Get the path to the SVG file inside the package
    svg_path = resources.files(svg_package) / svg_name

    # Load SVG directly from file path (QtSvg accepts pathlib.Path)
    renderer = QSvgRenderer(str(svg_path))
    if not renderer.isValid():
        raise Exception(f"SVG is invalid: {svg_path}")


    # Render SVG to pixmap
    pixmap = QPixmap(size, size)
    pixmap.fill(Qt.transparent)

    painter = QPainter(pixmap)
    renderer.render(painter, QRectF(0, 0, size, size))
    painter.end()

    return pixmap


def icon_from_package(svg_package: str, svg_name: str, size: int) -> QIcon:
    """
    Load an SVG icon from a Python package using importlib.resources
    and return a QIcon of the given size.

    :param svg_package: package path as string, e.g., 'pyspinw.gui.icons'
    :param svg_name: filename
    :param size: icon size in pixels
    :return: QIcon
    """

    pixmap = pixmap_from_package(svg_package, svg_name, size)


    # Create QIcon
    icon = QIcon()
    icon.addPixmap(pixmap)
    return icon

def load_icon(name: str, size: int=128) -> QIcon:
    """ Load an icon, choose based on dark/light mode """

    if QGuiApplication.styleHints().colorScheme() == Qt.ColorScheme.Dark:
        name += "-dark"


    return icon_from_package("pyspinw.gui.icons.svg", name + ".svg", size)


def load_pixmap(name: str, size: int=128) -> QPixmap:
    """ Load an icon, choose based on dark/light mode """

    if QGuiApplication.styleHints().colorScheme() == Qt.ColorScheme.Dark:
        name += "-dark"


    return pixmap_from_package("pyspinw.gui.icons.svg", name + ".svg", size)

