"""Loading of icons with dark mode or light mode."""
from importlib import resources

from PySide6.QtCore import QByteArray, QBuffer, QRectF
from PySide6.QtGui import QIcon, Qt, QGuiApplication
from PySide6.QtSvg import QSvgRenderer
from PySide6.QtGui import QPixmap, QPainter

def pixmap_from_package(svg_package: str, svg_name: str, size: int) -> QPixmap:
    """Load an svg file using importlib and convert to a pixmap.

    Parameters
    ----------
    svg_package : str
        Name of the package containing the SVG resource.
    svg_name : str
        Filename of the SVG file inside the package.
    size : int
        Width and height of the output pixmap in pixels.

    Returns
    -------
    QPixmap
        Rendered pixmap from the SVG file.

    Raises
    ------
    Exception
        If the SVG file is invalid.
    """
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
    """Load an svg icon from a package and return a QIcon.

    Parameters
    ----------
    svg_package : str
        Name of the package containing the SVG resource.
    svg_name : str
        Filename of the SVG file inside the package.
    size : int
        Width and height of the icon pixmap in pixels.

    Returns
    -------
    QIcon
        Loaded icon containing the rendered pixmap.
    """
    pixmap = pixmap_from_package(svg_package, svg_name, size)


    # Create QIcon
    icon = QIcon()
    icon.addPixmap(pixmap)
    return icon

def load_icon(name: str, size: int=128) -> QIcon:
    """Load an icon, choose based on dark/light mode.

    Parameters
    ----------
    name : str
        Base name of the SVG icon without extension.
    size : int, optional
        Width and height of the icon pixmap in pixels. Default is ``128``.

    Returns
    -------
    QIcon
        Loaded icon appropriate for the current color scheme.
    """
    if QGuiApplication.styleHints().colorScheme() == Qt.ColorScheme.Dark:
        name += "-dark"


    return icon_from_package("pyspinw.gui.icons.svg", name + ".svg", size)


def load_pixmap(name: str, size: int=128) -> QPixmap:
    """Load an icon, choose based on dark/light mode.

    Parameters
    ----------
    name : str
        Base name of the SVG icon without extension.
    size : int, optional
        Width and height of the pixmap in pixels. Default is ``128``.

    Returns
    -------
    QPixmap
        Loaded pixmap appropriate for the current color scheme.
    """
    if QGuiApplication.styleHints().colorScheme() == Qt.ColorScheme.Dark:
        name += "-dark"


    return pixmap_from_package("pyspinw.gui.icons.svg", name + ".svg", size)

def png_icon(name):
    """Load a png icon from the icons directory.

    Parameters
    ----------
    name : str
        Base name of the PNG icon without the ``.png`` extension.

    Returns
    -------
    QIcon
        Loaded PNG icon.
    """
    with resources.as_file(resources.files("pyspinw.gui.icons") / f"{name}.png") as icon_path:
        return QIcon(str(icon_path))
