""" Metadata for exchanges """
from pyspinw.serialisation import SPWSerialisable, SPWSerialisationContext, SPWDeserialisationContext, expects_keys, \
    rgb_serialise, rgb_deserialise


class ExchangeMetadata(SPWSerialisable):
    """ Metadata for exchanges"""

    def __init__(self, color: tuple[float, float, float] | None = None):
        self.color = color

    def _serialise(self, context: SPWSerialisationContext) -> dict:
        return {
            "color": None if self.color is None else rgb_serialise(self.color)
        }

    @staticmethod
    @expects_keys("color")
    def _deserialise(json: dict, context: SPWDeserialisationContext):

        color =  None if json["color"] is None else rgb_deserialise(json["color"])

        return ExchangeMetadata(color=color)

    def copy(self):
        """ Return a copy of this exchange metadata """
        return ExchangeMetadata(self.color)