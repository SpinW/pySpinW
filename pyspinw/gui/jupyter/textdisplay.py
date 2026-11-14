"""Jupyter ipytree-based table and tree display for sites and exchanges."""

from collections import defaultdict
from collections.abc import Callable
import html

import ipytree
import ipywidgets as widgets
import numpy as np

from pyspinw.anisotropy import Anisotropy
from pyspinw.gui.rendermodel import RenderExchange, RenderModel
from pyspinw.gui.textdisplay_format import format_triple, format_spin_data, format_anisotropies, format_g

_TREE_STYLE_HTML = """
<style>
.pyspinw-tree-container .jstree-default .jstree-wholerow-clicked,
.pyspinw-tree-container .jstree-default .jstree-wholerow-clicked.jstree-wholerow-hovered {
    background: #DD9911 !important;
    background-color: #DD9911 !important;
    color: black !important;
}
.pyspinw-tree-container .jstree-default .jstree-clicked {
    color: black !important;
}
.pyspinw-tree-container .jstree-anchor {
    width: calc(100% - 24px);
}
.pyspinw-tree-row {
    display: inline-flex;
    width: 100%;
    font-family: monospace;
    font-size: 12px;
    align-items: center;
}
.pyspinw-tree-cell {
    display: inline-block;
    flex: 1 1 0;
    min-width: 70px;
    padding: 0px 4px;
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
}
.pyspinw-tree-header {
    display: flex;
    width: 100%;
    box-sizing: border-box;
    padding-left: 24px;
    border-bottom: 1px solid #666;
    font-family: sans-serif;
    font-size: 12px;
    font-weight: bold;
}
.pyspinw-tree-header-cell {
    flex: 1 1 0;
    min-width: 70px;
    padding: 3px 4px;
    text-align: left;
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
}
</style>
"""


def _format_row_html(columns: list[str], is_bold: bool = False) -> str:
    """Format a sequence of column strings as an inline-flex HTML row for ``ipytree.Node``."""
    cells = "".join(f'<span class="pyspinw-tree-cell">{html.escape(str(col))}</span>' for col in columns)
    weight = "font-weight:bold;" if is_bold else "font-weight:normal;"
    return f'<span class="pyspinw-tree-row" style="{weight}">{cells}</span>'


def _format_header_html(headers: list[str]) -> str:
    """Format column header labels as an aligned HTML header bar."""
    cells = "".join(f'<span class="pyspinw-tree-header-cell">{html.escape(str(h))}</span>' for h in headers)
    return f'<div class="pyspinw-tree-header">{cells}</div>'


class DisplayItem(ipytree.Node):
    """Tree node item representing a single row in a ``ParameterTable``.

    Mirrors :class:`pyspinw.gui.textdisplay.DisplayItem` using :class:`ipytree.Node`.

    Parameters
    ----------
    columns : list of str or str
        Column strings (or single label string) for this row.
    ids : list of int
        Render IDs associated with this item.
    opened : bool, optional
        Whether the node is initially expanded. Default is ``False``.
    """

    def __init__(self, columns: list[str] | str, ids: list[int], opened: bool = False):
        if isinstance(columns, str):
            self.columns = [columns]
        else:
            self.columns = [str(c) for c in columns]

        self.ids = ids
        self.parent_item: DisplayItem | None = None
        self._is_bold: bool = False

        super().__init__(
            name=_format_row_html(self.columns, is_bold=False),
            nodes=(),
            opened=opened,
            show_icon=False,
        )

    def append_row(self, child: "DisplayItem"):
        """Append a child ``DisplayItem`` row to this node.

        Parameters
        ----------
        child : DisplayItem
            Child row node to attach under this item.
        """
        child.parent_item = self
        self.add_node(child)

    def set_bold(self, is_bold: bool):
        """Set whether this row is rendered with a bold font.

        Parameters
        ----------
        is_bold : bool
            Whether the row text should be bold.
        """
        if self._is_bold != is_bold:
            self._is_bold = is_bold
            self.name = _format_row_html(self.columns, is_bold=is_bold)


class ParameterTable(ipytree.Tree):
    """Tree widget for the text-based site viewer and exchange viewer.

    Mirrors :class:`pyspinw.gui.textdisplay.ParameterTable` using :class:`ipytree.Tree`.

    Parameters
    ----------
    headers : list of str or None, optional
        Column header labels for the table. Default is ``None``.
    **kwargs
        Additional keyword arguments forwarded to :class:`ipytree.Tree`.
    """

    def __init__(self, headers: list[str] | None = None, **kwargs):
        if "layout" not in kwargs:
            kwargs["layout"] = widgets.Layout(width="100%", max_height="220px", overflow="auto", margin="0px")
        super().__init__(multiple_selection=True, animation=100, **kwargs)
        self.add_class("pyspinw-tree-container")

        self.headers: list[str] = list(headers) if headers is not None else []
        self.header_widget = widgets.HTML(
            value=_format_header_html(self.headers) if self.headers else "",
            layout=widgets.Layout(width="100%", margin="0px"),
        )

        self.hover_ids: list[int] = []
        self._hover_callbacks: list[Callable[[], None]] = []
        self._selection_callbacks: list[Callable[[], None]] = []
        self._observed_nodes: set[str] = set()
        self._block_signals: bool = False

        self.observe(self._on_nodes_changed, names="nodes")
        self.observe(self._on_selected_nodes_changed, names="selected_nodes")

    def set_horizontal_header_labels(self, headers: list[str]):
        """Set the column header labels for this parameter table.

        Parameters
        ----------
        headers : list of str
            Column header names.
        """
        self.headers = list(headers)
        self.header_widget.value = _format_header_html(self.headers)

    def iter_all_items(self) -> list[DisplayItem]:
        """Return a flat list of all ``DisplayItem`` nodes in this tree.

        Returns
        -------
        list of DisplayItem
            All top-level and child ``DisplayItem`` instances.
        """
        result: list[DisplayItem] = []

        def _collect(nodes):
            for node in nodes:
                if isinstance(node, DisplayItem):
                    result.append(node)
                if node.nodes:
                    _collect(node.nodes)

        _collect(self.nodes)
        return result

    def _register_node_observers(self, node: ipytree.Node):
        """Attach ``selected`` and ``nodes`` observers recursively to ``node``."""
        if node._id not in self._observed_nodes:
            self._observed_nodes.add(node._id)
            node.observe(self._on_node_selected_changed, names="selected")
            node.observe(self._on_child_nodes_changed, names="nodes")
        for child in node.nodes:
            self._register_node_observers(child)

    def _on_nodes_changed(self, _change):
        for node in self.nodes:
            self._register_node_observers(node)

    def _on_child_nodes_changed(self, change):
        owner = change["owner"]
        for child in owner.nodes:
            self._register_node_observers(child)

    def _on_node_selected_changed(self, _change):
        if self._block_signals:
            return

        selected = tuple(item for item in self.iter_all_items() if item.selected)
        self._block_signals = True
        try:
            self.set_trait("selected_nodes", selected)
        finally:
            self._block_signals = False

        for cb in self._selection_callbacks:
            cb()

    def _on_selected_nodes_changed(self, _change):
        if self._block_signals:
            return

        selected_set = {node._id for node in self.selected_nodes}
        self._block_signals = True
        try:
            for item in self.iter_all_items():
                should_select = item._id in selected_set
                if item.selected != should_select:
                    item.selected = should_select
        finally:
            self._block_signals = False

        for cb in self._selection_callbacks:
            cb()

    def select_items(self, items: list[DisplayItem]):
        """Select the given ``DisplayItem`` rows while blocking selection signals.

        Parameters
        ----------
        items : list of DisplayItem
            Items to mark as selected; all other items are deselected.
        """
        selected_ids = {item._id for item in items}
        self._block_signals = True
        try:
            for item in self.iter_all_items():
                should_select = item._id in selected_ids
                if item.selected != should_select:
                    item.selected = should_select
            self.set_trait("selected_nodes", tuple(items))
        finally:
            self._block_signals = False

    def on_hover_changed(self, callback: Callable[[], None]):
        """Register a callback invoked when hovered items change.

        Parameters
        ----------
        callback : callable
            Zero-argument function.
        """
        self._hover_callbacks.append(callback)

    def on_selection_changed(self, callback: Callable[[], None]):
        """Register a callback invoked when selected items change.

        Parameters
        ----------
        callback : callable
            Zero-argument function.
        """
        self._selection_callbacks.append(callback)


class JupyterTextDisplay(widgets.VBox):
    """Interactive text display of sites and exchanges using ``ipytree.Tree``.

    Mirrors :class:`pyspinw.gui.textdisplay.TextDisplay` with two :class:`ParameterTable`
    tree widgets for sites and exchanges.

    Parameters
    ----------
    render_model : RenderModel
        Render model containing the sites, exchanges, and anisotropies to display.
    """

    SITE_HEADERS = ["Name", "Lattice Pos.", "Cartesian Pos.", "Moment", "Anisotropies", "g-Factor"]
    EXCHANGE_HEADERS = ["Name", "Offset", "Type", "Site 1", "Site 2", "Parameters"]

    def __init__(self, render_model: RenderModel):
        self.render_model = render_model

        self._hover_callbacks: list[Callable[[], None]] = []
        self._selection_callbacks: list[Callable[[], None]] = []

        #
        # Sites and anisotropies
        #
        site_uid_to_anisotropies: dict[int, list[Anisotropy]] = defaultdict(list)

        for anisotropy in render_model.expanded.anisotropies:
            site_uid_to_anisotropies[anisotropy.site.unique_id].append(anisotropy)

        for anisotropy in render_model.hamiltonian.anisotropies:
            site_uid_to_anisotropies[anisotropy.site.unique_id].append(anisotropy)

        self.site_tree = ParameterTable(headers=self.SITE_HEADERS)

        small_cell = render_model.hamiltonian.structure.unit_cell
        big_cell = render_model.expanded.structure.unit_cell

        # Mapping from render ID to [unexpanded_item, expanded_item]
        self.site_render_id_to_row_item: dict[int, list[DisplayItem]] = defaultdict(list)

        # Mapping from expanded site UID to offset
        site_expanded_uid_to_offset: dict[int, tuple[int, int, int] | None] = {}

        #
        # Fill out site data
        #
        site_root_nodes: list[DisplayItem] = []
        for unexpanded_index, site in enumerate(render_model.original_sites):
            parent_ids: list[int] = []

            unexpanded_row = DisplayItem(
                [
                    site.name,
                    format_triple(site.ijk),
                    format_triple(small_cell.lattice_units_to_cartesian(site.ijk)),
                    format_spin_data(site.spin_data),
                    format_anisotropies(site_uid_to_anisotropies[site.unique_id]),
                    format_g(site.g),
                ],
                parent_ids,
            )

            for expanded_index in render_model.original_index_to_expanded[unexpanded_index]:
                expanded_render_site = render_model.sites[expanded_index]
                expanded_site = expanded_render_site.site

                ids = [expanded_render_site.render_id]
                parent_ids += ids

                offset = expanded_render_site.offset
                site_expanded_uid_to_offset[expanded_site.unique_id] = offset

                expanded_row = DisplayItem(
                    [
                        f"{expanded_site.name} {offset}",
                        format_triple(expanded_site.ijk),
                        format_triple(big_cell.lattice_units_to_cartesian(expanded_site.ijk)),
                        format_triple(expanded_site.base_spin),
                        format_anisotropies(site_uid_to_anisotropies[expanded_site.unique_id]),
                        format_g(expanded_site.g),
                    ],
                    ids=ids,
                )

                unexpanded_row.append_row(expanded_row)

                for render_id in ids:
                    self.site_render_id_to_row_item[render_id] = [unexpanded_row, expanded_row]

            site_root_nodes.append(unexpanded_row)

        self.site_tree.nodes = tuple(site_root_nodes)

        #
        # Exchanges
        #
        self.exchange_render_id_to_row_item: dict[int, list[DisplayItem]] = defaultdict(list)

        self.exchanges_original_to_expanded: dict[int, list[int]] = defaultdict(list)
        for expanded_index, unexpanded_index in enumerate(render_model.exchange_mapping):
            self.exchanges_original_to_expanded[unexpanded_index].append(expanded_index)

        self.exchange_tree = ParameterTable(headers=self.EXCHANGE_HEADERS)

        exchange_root_nodes: list[DisplayItem] = []
        for unexpanded_index, unexpanded_exchange in enumerate(render_model.hamiltonian.exchanges):
            parent_ids = []

            unexpanded_row = DisplayItem(
                [
                    unexpanded_exchange.name,
                    str(unexpanded_exchange.cell_offset.as_tuple),
                    unexpanded_exchange.exchange_type,
                    unexpanded_exchange.site_1.name,
                    unexpanded_exchange.site_2.name,
                    unexpanded_exchange.parameter_string,
                ],
                parent_ids,
            )

            for expanded_index in self.exchanges_original_to_expanded[unexpanded_index]:
                expanded_exchange: RenderExchange = render_model.exchanges[expanded_index]

                ids = [expanded_exchange.render_id]
                parent_ids += ids

                site_1_offset = site_expanded_uid_to_offset.get(expanded_exchange.exchange.site_1.unique_id)
                site_1_name = f"{expanded_exchange.exchange.site_1.name} {site_1_offset}"
                site_2_offset = site_expanded_uid_to_offset.get(expanded_exchange.exchange.site_2.unique_id)
                site_2_name = f"{expanded_exchange.exchange.site_2.name} {site_2_offset}"

                expanded_row = DisplayItem(
                    [
                        expanded_exchange.exchange.name,
                        str(expanded_exchange.exchange.cell_offset.as_tuple),
                        expanded_exchange.exchange.exchange_type,
                        site_1_name,
                        site_2_name,
                        expanded_exchange.exchange.parameter_string,
                    ],
                    ids,
                )

                unexpanded_row.append_row(expanded_row)

                for render_id in ids:
                    self.exchange_render_id_to_row_item[render_id] = [unexpanded_row, expanded_row]

            exchange_root_nodes.append(unexpanded_row)

        self.exchange_tree.nodes = tuple(exchange_root_nodes)

        #
        # Layout components
        #
        style_widget = widgets.HTML(value=_TREE_STYLE_HTML, layout=widgets.Layout(display="none"))
        sites_box = widgets.VBox(
            [self.site_tree.header_widget, self.site_tree],
            layout=widgets.Layout(width="100%", margin="0px 0px 6px 0px", border="1px solid #444"),
        )
        exchanges_box = widgets.VBox(
            [self.exchange_tree.header_widget, self.exchange_tree],
            layout=widgets.Layout(width="100%", margin="0px", border="1px solid #444"),
        )

        super().__init__(
            children=[style_widget, sites_box, exchanges_box],
            layout=widgets.Layout(width="520px", padding="0px 6px"),
        )

        #
        # Hovering and selection state
        #
        self.hover_ids: list[int] = []
        self.current_selection: list[int] = []
        self.site_tree.on_hover_changed(self.on_hover_changed)
        self.exchange_tree.on_hover_changed(self.on_hover_changed)
        self.site_tree.on_selection_changed(self.on_selection_changed)
        self.exchange_tree.on_selection_changed(self.on_selection_changed)

    def on_hover_changed(self, callback: Callable[[], None] | None = None):
        """Handle tree hover changes or register a callback invoked when hovered items change.

        Parameters
        ----------
        callback : callable or None, optional
            If provided, registers ``callback`` to be invoked when hovered items change.
            If ``None``, updates ``self.hover_ids`` from the trees and notifies callbacks.
        """
        if callback is not None:
            self._hover_callbacks.append(callback)
            return

        self.hover_ids = self.site_tree.hover_ids + self.exchange_tree.hover_ids
        for cb in self._hover_callbacks:
            cb()

    def on_selection_changed(self, callback: Callable[[], None] | None = None):
        """Handle tree selection changes or register a callback invoked when selected items change.

        Parameters
        ----------
        callback : callable or None, optional
            If provided, registers ``callback`` to be invoked when selected items change.
            If ``None``, collects selected ``render_ids`` from both trees and notifies callbacks.
        """
        if callback is not None:
            self._selection_callbacks.append(callback)
            return

        render_ids: list[int] = []

        # Sites
        for node in self.site_tree.selected_nodes:
            if isinstance(node, DisplayItem):
                render_ids += node.ids

        # Couplings
        for node in self.exchange_tree.selected_nodes:
            if isinstance(node, DisplayItem):
                render_ids += node.ids

        # Remove duplicates while preserving order
        self.current_selection = list(dict.fromkeys(render_ids))

        for cb in self._selection_callbacks:
            cb()

    def set_row_boldness(self, item: DisplayItem, is_bold: bool):
        """Set the boldness of the font for a given row item.

        Parameters
        ----------
        item : DisplayItem
            The row item to update.
        is_bold : bool
            Whether the row font should be bold.
        """
        item.set_bold(is_bold)

    def set_hover(self, render_ids: list[int]):
        """Set the hover state for specified ``render_ids``.

        Parameters
        ----------
        render_ids : list of int
            Render IDs of currently hovered elements.
        """
        self.hover_ids = list(render_ids)

        for items in self.site_render_id_to_row_item.values():
            for item in items:
                self.set_row_boldness(item, False)

        for render_id in render_ids:
            for item in self.site_render_id_to_row_item[render_id]:
                self.set_row_boldness(item, True)

        for items in self.exchange_render_id_to_row_item.values():
            for item in items:
                self.set_row_boldness(item, False)

        for render_id in render_ids:
            for item in self.exchange_render_id_to_row_item[render_id]:
                self.set_row_boldness(item, True)

    def set_selection(self, render_ids: list[int]):
        """Set the selection based on ``render_ids`` without emitting recursive selection signals.

        Parameters
        ----------
        render_ids : list of int
            Render IDs of elements to select in the tree views.
        """
        self.current_selection = list(dict.fromkeys(render_ids))

        #
        # Sites
        #
        selected_site_items: list[DisplayItem] = []
        for render_id in render_ids:
            for item in self.site_render_id_to_row_item[render_id]:
                parent = item.parent_item
                if parent is not None:
                    selected_site_items.append(item)
                    parent.opened = True

        self.site_tree.select_items(selected_site_items)

        #
        # Couplings
        #
        selected_exchange_items: list[DisplayItem] = []
        for render_id in render_ids:
            for item in self.exchange_render_id_to_row_item[render_id]:
                parent = item.parent_item
                if parent is not None:
                    selected_exchange_items.append(item)
                    parent.opened = True

        self.exchange_tree.select_items(selected_exchange_items)
