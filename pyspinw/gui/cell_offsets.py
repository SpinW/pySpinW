from pydantic import BaseModel


class CellOffset(BaseModel):
    i: int
    j: int
    k: int

    @property
    def as_tuple(self):
        return self.i, self.j, self.k