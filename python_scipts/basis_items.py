from dataclasses import dataclass
from enum import Enum

from python_scipts.input_items import InputMember

class InteractionType(Enum):
    RotorAtom = 1
    VibrotorAtom = 2
    RotorRotor = 3
    PlugIn = 9

@dataclass
class BasisItems:
    int_type: InteractionType
    levels: InputMember

    def to_string(self) -> str:
        return f"""&Basis
            ITYPE = {self.int_type.value},
            {self.levels.to_string()}"""

type AngL = int
type Vib = int

@dataclass
class VibrotorLevels:
    levels: list[float]
    states: list[tuple[AngL, Vib]]

    def to_string(self) -> str:
        return f"EUNITS = 7,\n NLEVEL = {len(self.levels)}, JLEVEL = " \
            + ", ".join(map(lambda x: f"{x[0]}, {x[1]}", self.states)) \
            + ",\nELEVEL = " + ", ".join(map(lambda x: str(x), self.levels))
