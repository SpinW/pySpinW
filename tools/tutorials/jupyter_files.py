import json

import base64

class JupyterFile:

    """ A Jupyter File """

    _metadata = {
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3"
            },
            "language_info": {
                "name": "python",
                "version": "3.12.0"
            }
        }
    }

    def __init__(self):
        self.notebook = {
            "cells": [],
            "metadata": self._metadata,
            "nbformat": 4,
            "nbformat_minor": 5
        }

    def add_code(self, lines: list[str] | str):
        """ Add a code cell"""

        if isinstance(lines, str):
            lines = [lines]

        # Remove front and back empty lines
        lines = [line.rstrip() for line in lines]

        ## Remove from front
        new_lines = []
        keep = False
        for line in lines:
            if line.strip() != "":
                keep = True

            if keep:
                new_lines.append(line)

        lines = new_lines

        ## Remover from back
        lines.reverse()

        new_lines = []
        keep = False
        for line in lines:
            if line.strip() != "":
                keep = True

            if keep:
                new_lines.append(line)

        lines = new_lines

        lines.reverse()

        # Don't make a cell if it is empty
        if not lines:
            return

        self.notebook["cells"].append(
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": "\n".join(lines)
            })

    def add_output_image_to_last(self, image_file):
        """ Add image output to last cell """

        with open(image_file, 'rb') as file:

            image_base64 = base64.b64encode(file.read()).decode("ascii")

            self.notebook["cells"][-1]["outputs"].append(
            {
                "output_type": "display_data",
                "data": {
                    "image/png": image_base64,
                    "text/plain": "Figure"
                },
                "metadata": {}
            })

    def add_text_output_to_last(self, lines):
        """ Add a text output to the last cell """

        self.notebook["cells"][-1]["outputs"].append(
            {
                "output_type": "display_data",
                "data": {
                    "text/plain": "".join(lines)
                },
                "metadata": {}
            })

    def add_text(self, lines: list[str] | str):
        """ Add a text cell """

        if isinstance(lines, str):
            lines = [lines]

        lines = [line for line in lines if line.strip() != ""]

        if not lines:
            return

        self.notebook["cells"].append(
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": lines
            })

    def write_notebook(self, filename):
        with open(filename, 'w') as file:
            json.dump(self.notebook, file)