""" Loading of icons with dark mode or light mode """
from importlib import resources

from PySide6.QtCore import QByteArray, QBuffer, QRectF
from PySide6.QtGui import QIcon, Qt, QGuiApplication
from PySide6.QtSvg import QSvgRenderer
from PySide6.QtGui import QPixmap, QPainter

def pixmap_from_package(svg_package: str, svg_name: str, size: int) -> QPixmap:
    """ Load an svg file using importlib and convert to a pixmap"""
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
    """Load an svg icon from a package and return a QIcon  """
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

def png_icon(name):
    """ Load a png icon from the icons directory"""
    # TODO: Can this be made robust against weird packaging stuff???
    with resources.as_file(resources.files("pyspinw.gui.icons") / f"{name}.png") as icon_path:
        return QIcon(str(icon_path))
