"""Jupyter ipywidgets text-based table and tree display for sites and exchanges."""

from collections import defaultdict
from collections.abc import Callable
from dataclasses import dataclass, field
import html

import ipywidgets as widgets
import numpy as np

from pyspinw.anisotropy import Anisotropy
from pyspinw.gui.rendermodel import RenderExchange, RenderModel


def format_triple(triple) -> str:
    """Format a 3-element numeric sequence as ``"(a, b, c)"``.

    Parameters
    ----------
    triple : tuple or array_like of float
        Sequence of three numeric values.

    Returns
    -------
    str
        Formatted string representation ``"(a, b, c)"``.
    """
    a, b, c = triple
    return f"({a:.3g}, {b:.3g}, {c:.3g})"


def format_spin_data(spin_data: np.ndarray) -> str:
    """Format a 2D array of spin vectors as comma-separated triples.

    Parameters
    ----------
    spin_data : numpy.ndarray
        2D array of shape ``(n, 3)`` containing spin vectors.

    Returns
    -------
    str
        Comma-separated string of formatted spin triples.
    """
    return ", ".join([format_triple(spin_data[i, :]) for i in range(spin_data.shape[0])])


def format_anisotropies(anisotropies: list[Anisotropy]) -> str:
    """Format a list of site anisotropies for table display.

    Parameters
    ----------
    anisotropies : list of Anisotropy
        List of anisotropy terms associated with a site.

    Returns
    -------
    str
        Semicolon-separated parameter strings, or ``"None"`` if empty.
    """
    if anisotropies:
        return "; ".join([a.parameter_string for a in anisotropies])
    return "None"


def format_g(g_matrix: np.ndarray) -> str:
    """Format a 3x3 g-tensor matrix into a compact string.

    Parameters
    ----------
    g_matrix : (3, 3) numpy.ndarray
        3x3 g-tensor matrix.

    Returns
    -------
    str
        Compact string representation of the g-tensor.
    """
    diag = np.diagonal(g_matrix)
    if np.all(g_matrix == np.diag(diag)):
        a, b, c = diag
        if a == b and b == c:
            return f"{a:.3g}"
        return f"{a:.3g}, {b:.3g}, {c:.3g}"

    m = g_matrix.reshape(-1)
    return (
        f"[[{m[0]:.3g}, {m[1]:.3g}, {m[2]:.3g}],"
        f" [{m[3]:.3g}, {m[4]:.3g}, {m[5]:.3g}],"
        f" [{m[6]:.3g}, {m[7]:.3g}, {m[8]:.3g}]]"
    )


@dataclass
class TableRowNode:
    """Hierarchical row node for a site or exchange entry.

    Parameters
    ----------
    label : str
        Identifier label shown in the selection list.
    columns : list of str
        Formatted column strings for this row.
    ids : list of int
        Render IDs associated with this row.
    children : list of TableRowNode, optional
        Expanded supercell child rows.
    """

    label: str
    columns: list[str]
    ids: list[int]
    children: list["TableRowNode"] = field(default_factory=list)


class JupyterTextDisplay(widgets.VBox):
    """Interactive text display of sites and exchanges using ``ipywidgets``.

    Replicates the hierarchical Sites and Exchanges tables from ``TextDisplay``
    with selection and hover synchronization.

    Parameters
    ----------
    render_model : RenderModel
        Render model containing the sites, exchanges, and anisotropies to display.
    """

    SITE_HEADERS = ["Name", "Lattice Pos.", "Cartesian Pos.", "Moment", "Anisotropies", "g-Factor"]
    EXCHANGE_HEADERS = ["Name", "Offset", "Type", "Site 1", "Site 2", "Parameters"]

    def __init__(self, render_model: RenderModel):
        self.render_model = render_model
        self.hover_ids: list[int] = []
        self.current_selection: list[int] = []

        self._hover_callbacks: list[Callable[[], None]] = []
        self._selection_callbacks: list[Callable[[], None]] = []
        self._block_selection_events = False

        # Build data models for sites and exchanges
        self.site_nodes: list[TableRowNode] = []
        self.exchange_nodes: list[TableRowNode] = []
        self._site_option_to_ids: dict[str, list[int]] = {}
        self._exchange_option_to_ids: dict[str, list[int]] = {}
        self._site_render_id_to_option: dict[int, str] = {}
        self._exchange_render_id_to_option: dict[int, str] = {}

        self._build_data(render_model)

        # Interactive multi-select widgets
        site_options = list(self._site_option_to_ids.keys())
        exchange_options = list(self._exchange_option_to_ids.keys())

        self.site_selector = widgets.SelectMultiple(
            options=site_options,
            value=(),
            rows=min(8, max(3, len(site_options))),
            description="Select Sites:",
            style={"description_width": "95px"},
            layout=widgets.Layout(width="98%"),
        )
        self.exchange_selector = widgets.SelectMultiple(
            options=exchange_options,
            value=(),
            rows=min(8, max(3, len(exchange_options))),
            description="Exchanges:",
            style={"description_width": "95px"},
            layout=widgets.Layout(width="98%"),
        )
        self.clear_selection_button = widgets.Button(
            description="Clear Selection",
            tooltip="Clear all selected sites and exchanges",
            layout=widgets.Layout(width="130px", height="28px", margin="2px 0px"),
        )

        self.site_selector.observe(self._on_selector_changed, names="value")
        self.exchange_selector.observe(self._on_selector_changed, names="value")
        self.clear_selection_button.on_click(self._on_clear_selection_clicked)

        # HTML Table Views
        self.site_table_html = widgets.HTML(
            value=self._render_table_html(self.SITE_HEADERS, self.site_nodes),
            layout=widgets.Layout(width="100%", max_height="240px", overflow="auto"),
        )
        self.exchange_table_html = widgets.HTML(
            value=self._render_table_html(self.EXCHANGE_HEADERS, self.exchange_nodes),
            layout=widgets.Layout(width="100%", max_height="240px", overflow="auto"),
        )

        sites_box = widgets.VBox(
            [
                widgets.HTML("<b>Sites</b>"),
                self.site_table_html,
                self.site_selector,
            ],
            layout=widgets.Layout(width="100%", margin="0px 0px 6px 0px"),
        )
        exchanges_box = widgets.VBox(
            [
                widgets.HTML("<b>Exchanges</b>"),
                self.exchange_table_html,
                self.exchange_selector,
            ],
            layout=widgets.Layout(width="100%", margin="0px 0px 4px 0px"),
        )

        super().__init__(
            children=[sites_box, exchanges_box, self.clear_selection_button],
            layout=widgets.Layout(width="480px", padding="0px 6px"),
        )

    def _build_data(self, render_model: RenderModel):
        """Populate ``self.site_nodes`` and ``self.exchange_nodes`` from ``render_model``."""
        site_uid_to_anisotropies: dict[int, list[Anisotropy]] = defaultdict(list)
        for anisotropy in render_model.expanded.anisotropies:
            site_uid_to_anisotropies[anisotropy.site.unique_id].append(anisotropy)
        for anisotropy in render_model.hamiltonian.anisotropies:
            site_uid_to_anisotropies[anisotropy.site.unique_id].append(anisotropy)

        small_cell = render_model.hamiltonian.structure.unit_cell
        big_cell = render_model.expanded.structure.unit_cell

        site_expanded_uid_to_offset: dict[int, tuple[int, int, int] | None] = {}

        # 1. Sites
        for unexpanded_index, site in enumerate(render_model.original_sites):
            parent_ids: list[int] = []
            children: list[TableRowNode] = []

            for expanded_index in render_model.original_index_to_expanded[unexpanded_index]:
                expanded_render_site = render_model.sites[expanded_index]
                expanded_site = expanded_render_site.site
                ids = [expanded_render_site.render_id]
                parent_ids.extend(ids)

                offset = expanded_render_site.offset
                site_expanded_uid_to_offset[expanded_site.unique_id] = offset
                child_name = f"{expanded_site.name} {offset}"
                child_cols = [
                    child_name,
                    format_triple(expanded_site.ijk),
                    format_triple(big_cell.lattice_units_to_cartesian(expanded_site.ijk)),
                    format_triple(expanded_site.base_spin),
                    format_anisotropies(site_uid_to_anisotropies[expanded_site.unique_id]),
                    format_g(expanded_site.g),
                ]
                child_label = f"  ↳ {child_name} (#{expanded_render_site.render_id})"
                child_node = TableRowNode(label=child_label, columns=child_cols, ids=ids)
                children.append(child_node)

                self._site_option_to_ids[child_label] = ids
                self._site_render_id_to_option[expanded_render_site.render_id] = child_label

            parent_cols = [
                site.name,
                format_triple(site.ijk),
                format_triple(small_cell.lattice_units_to_cartesian(site.ijk)),
                format_spin_data(site.spin_data),
                format_anisotropies(site_uid_to_anisotropies[site.unique_id]),
                format_g(site.g),
            ]
            parent_label = f"{site.name} [all {len(children)}]"
            parent_node = TableRowNode(
                label=parent_label,
                columns=parent_cols,
                ids=list(parent_ids),
                children=children,
            )
            self.site_nodes.append(parent_node)
            self._site_option_to_ids[parent_label] = list(parent_ids)

        # 2. Exchanges
        exchanges_original_to_expanded: dict[int, list[int]] = defaultdict(list)
        for expanded_index, unexpanded_index in enumerate(render_model.exchange_mapping):
            exchanges_original_to_expanded[unexpanded_index].append(expanded_index)

        for unexpanded_index, unexpanded_exchange in enumerate(render_model.hamiltonian.exchanges):
            parent_ids = []
            children = []

            for expanded_index in exchanges_original_to_expanded[unexpanded_index]:
                expanded_exchange: RenderExchange = render_model.exchanges[expanded_index]
                ids = [expanded_exchange.render_id]
                parent_ids.extend(ids)

                site_1_offset = site_expanded_uid_to_offset.get(expanded_exchange.exchange.site_1.unique_id)
                site_1_name = f"{expanded_exchange.exchange.site_1.name} {site_1_offset}"
                site_2_offset = site_expanded_uid_to_offset.get(expanded_exchange.exchange.site_2.unique_id)
                site_2_name = f"{expanded_exchange.exchange.site_2.name} {site_2_offset}"

                child_cols = [
                    expanded_exchange.exchange.name,
                    str(expanded_exchange.exchange.cell_offset.as_tuple),
                    expanded_exchange.exchange.exchange_type,
                    site_1_name,
                    site_2_name,
                    expanded_exchange.exchange.parameter_string,
                ]
                child_label = (
                    f"  ↳ {expanded_exchange.exchange.name}: {site_1_name} → {site_2_name} "
                    f"(#{expanded_exchange.render_id})"
                )
                child_node = TableRowNode(label=child_label, columns=child_cols, ids=ids)
                children.append(child_node)

                self._exchange_option_to_ids[child_label] = ids
                self._exchange_render_id_to_option[expanded_exchange.render_id] = child_label

            parent_cols = [
                unexpanded_exchange.name,
                str(unexpanded_exchange.cell_offset.as_tuple),
                unexpanded_exchange.exchange_type,
                unexpanded_exchange.site_1.name,
                unexpanded_exchange.site_2.name,
                unexpanded_exchange.parameter_string,
            ]
            parent_label = (
                f"{unexpanded_exchange.name} "
                f"({unexpanded_exchange.site_1.name}→{unexpanded_exchange.site_2.name}) "
                f"[all {len(children)}]"
            )
            parent_node = TableRowNode(
                label=parent_label,
                columns=parent_cols,
                ids=list(parent_ids),
                children=children,
            )
            self.exchange_nodes.append(parent_node)
            self._exchange_option_to_ids[parent_label] = list(parent_ids)

    def _render_table_html(self, headers: list[str], nodes: list[TableRowNode]) -> str:
        """Render hierarchical ``TableRowNode`` items as an HTML table with highlight styling.

        Parameters
        ----------
        headers : list of str
            Column headers.
        nodes : list of TableRowNode
            Top-level row nodes.

        Returns
        -------
        str
            HTML table string.
        """
        hover_set = set(self.hover_ids)
        selection_set = set(self.current_selection)

        th_cells = "".join(
            f'<th style="padding:3px 6px; border-bottom:1px solid #666; text-align:left; font-size:12px;">'
            f"{html.escape(str(h))}</th>"
            for h in headers
        )

        rows_html: list[str] = []
        for parent in nodes:
            parent_hovered = bool(hover_set.intersection(parent.ids))
            parent_selected = bool(parent.ids) and all(rid in selection_set for rid in parent.ids)
            any_child_selected = bool(selection_set.intersection(parent.ids))

            parent_style = "font-size:12px;"
            if parent_selected:
                parent_style += "background-color:#DD9911; color:black;"
            if parent_hovered:
                parent_style += "font-weight:bold;"

            tds = "".join(
                f'<td style="padding:2px 6px; border-bottom:1px solid #333;">{html.escape(str(col))}</td>'
                for col in parent.columns
            )
            rows_html.append(f'<tr style="{parent_style}">{tds}</tr>')

            if parent.children:
                child_rows: list[str] = []
                for child in parent.children:
                    child_hovered = bool(hover_set.intersection(child.ids))
                    child_selected = bool(selection_set.intersection(child.ids))
                    child_style = "font-size:11px;"
                    if child_selected:
                        child_style += "background-color:#DD9911; color:black;"
                    if child_hovered:
                        child_style += "font-weight:bold;"

                    child_tds = "".join(
                        f'<td style="padding:2px 6px; border-bottom:1px dotted #333;">'
                        f"{'&nbsp;&nbsp;↳ ' if idx == 0 else ''}{html.escape(str(col))}</td>"
                        for idx, col in enumerate(child.columns)
                    )
                    child_rows.append(f'<tr style="{child_style}">{child_tds}</tr>')

                open_attr = " open" if any_child_selected else ""
                subtable = (
                    f'<tr><td colspan="{len(headers)}" style="padding:0px 0px 4px 12px;">'
                    f"<details{open_attr}><summary style=\"font-size:11px; cursor:pointer; color:#888;\">"
                    f"Expanded instances ({len(parent.children)})</summary>"
                    f'<table style="width:100%; border-collapse:collapse;">{"".join(child_rows)}</table>'
                    f"</details></td></tr>"
                )
                rows_html.append(subtable)

        return (
            '<table style="width:100%; border-collapse:collapse;">'
            f"<thead><tr>{th_cells}</tr></thead>"
            f'<tbody>{"".join(rows_html)}</tbody>'
            "</table>"
        )

    def _refresh_tables(self):
        """Re-render both HTML tables to reflect current hover and selection state."""
        self.site_table_html.value = self._render_table_html(self.SITE_HEADERS, self.site_nodes)
        self.exchange_table_html.value = self._render_table_html(self.EXCHANGE_HEADERS, self.exchange_nodes)

    def _on_selector_changed(self, _change):
        if self._block_selection_events:
            return

        render_ids: list[int] = []
        for label in self.site_selector.value:
            render_ids.extend(self._site_option_to_ids.get(label, []))
        for label in self.exchange_selector.value:
            render_ids.extend(self._exchange_option_to_ids.get(label, []))

        # Deduplicate while preserving order
        self.current_selection = list(dict.fromkeys(render_ids))
        self._refresh_tables()

        for cb in self._selection_callbacks:
            cb()

    def _on_clear_selection_clicked(self, _btn):
        self.set_selection([])
        for cb in self._selection_callbacks:
            cb()

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

    def set_hover(self, render_ids: list[int]):
        """Set the hover state for specified ``render_ids`` and refresh table formatting.

        Parameters
        ----------
        render_ids : list of int
            Render IDs of currently hovered elements.
        """
        self.hover_ids = list(render_ids)
        self._refresh_tables()

    def set_selection(self, render_ids: list[int]):
        """Set the selection based on ``render_ids`` without emitting recursive selection signals.

        Parameters
        ----------
        render_ids : list of int
            Render IDs of elements to select in the tables.
        """
        self.current_selection = list(dict.fromkeys(render_ids))

        site_labels = [
            self._site_render_id_to_option[rid]
            for rid in self.current_selection
            if rid in self._site_render_id_to_option
        ]
        exchange_labels = [
            self._exchange_render_id_to_option[rid]
            for rid in self.current_selection
            if rid in self._exchange_render_id_to_option
        ]

        self._block_selection_events = True
        try:
            self.site_selector.value = tuple(site_labels)
            self.exchange_selector.value = tuple(exchange_labels)
        finally:
            self._block_selection_events = False

        self._refresh_tables()
