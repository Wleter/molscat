from dataclasses import dataclass
from enum import Enum, auto, nonmember
from typing import Optional, Protocol

from python_scipts.utils.types import Interval

def red_mass(mass1: float, mass2: float):
    return mass1 * mass2 / (mass1 + mass2)


class InputMember(Protocol):
    def to_string(self) -> str:
        ...

class PrintFlag(Enum):
    No = 0
    Default = 12

    def to_string(self) -> str:
        return f"IPRINT = {self.value}"

class DistanceType(Enum):
    ShortRange = auto()
    LongRange = auto()

class PropagatorType(Enum):
    Vogelaere = 2
    RMatrix = 3
    VariableStep = 4
    JohnsonLogDerivative = 5
    DiabaticLogDerivative = 6
    QDiabaticLogDerivative = 7
    SymplecticLogDerivative = 8
    Airy = 9
    WKB = -1

    def with_config(self, step_config: float, distances: DistanceType) -> "PropagatorType":
        self.step_config = step_config
        self.distances = distances

        return self

    step_config = nonmember(float)
    distances = nonmember(DistanceType)

    def distance_suffix(self) -> str:
        match self.distances:
            case DistanceType.ShortRange:
                return "S"
            case DistanceType.LongRange:
                return "L"
            case _:
                raise NameError

    def step_string(self) -> str:
        match self:
            case PropagatorType.Airy:
                return f"TOLHI{self.distance_suffix()}"
            case _:
                return f"DR{self.distance_suffix()}"

    def to_string(self) -> str:
        return f"IPROP{self.distance_suffix()} = {self.value}, \
            {self.step_string()} = {self.step_config}"

@dataclass
class Energies(InputMember):
    energies: list[float]

    def to_string(self) -> str:
        return f"EUNITS = 7, NNRG = {len(self.energies)}, Energy = " \
            + ", ".join(map(lambda x: str(x), self.energies))

class PropagationConfig(InputMember):
    r_min: float
    r_mid: float
    r_max: float
    near_config: PropagatorType
    far_config: PropagatorType

    def __init__(self, range: Interval, near_config: PropagatorType, far_config: PropagatorType):
        self.r_min = range.min_value
        self.r_max = range.max_value
        self.r_mid = range.mid_value
        self.near_config = near_config
        self.far_config = far_config

    def to_string(self) -> str:
        return f"""RMIN = {self.r_min}, RMID = {self.r_mid}, RMAX = {self.r_max}, IRMSET = 0,
        {self.near_config.to_string()},
        {self.far_config.to_string()}"""

class Symmetry(Enum):
    No = 0
    Odd = 1
    Even = 2

@dataclass
class InputItems(InputMember):
    label: str
    red_mass: float
    print_flag: PrintFlag
    propagation: PropagationConfig
    energies: Energies

    l_range: tuple[float, float] = (0, 0)
    symmetry: Symmetry = Symmetry.No
    optional: Optional[list[InputMember]] = None

    def to_string(self) -> str:
        input_string = f"""&INPUT
                LABEL = '{self.label}',
                URED = {self.red_mass},
                {self.print_flag.to_string()},
                {self.propagation.to_string()},
                JTOTL = {self.l_range[0]}, JTOTU = {self.l_range[1]}, IBFIX = {self.symmetry.value},
                {self.energies.to_string()}, RUNIT = 0.529177210903,
                LASTIN = 1"""

        if self.optional != None:
            input_string += ",\n".join(map(lambda x: x.to_string(), self.optional))

        return input_string
