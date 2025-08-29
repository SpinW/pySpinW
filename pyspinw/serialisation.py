""" Serialisation mixin
"""

class SPWSerialisable:
    def serialise(self) -> dict:
        raise NotImplementedError("Serialisation not implemented")

    @staticmethod
    def deserialise(data: dict):
        raise NotImplementedError("Deserialisation not implemented")

