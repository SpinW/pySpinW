""" Renderable base class """

from abc import ABC

import logging

logger = logging.getLogger("GL Subsystem")

class Renderable(ABC):
    """ Interface for everything that can be rendered with the OpenGL widget"""

    def render_wireframe(self):
        """ Override here with GL commands to render the solid object in a local frame"""
        logger.debug(f"{self.__class__} does not support wireframe rendering")


    def render_solid(self):
        """ Override here with GL commands to render the wireframe object in a local frame"""
        logger.debug(f"{self.__class__} does not support solid rendering")
