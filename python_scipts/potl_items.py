from dataclasses import dataclass
from python_scipts.input_items import InputMember
from python_scipts.basis_items import AngL, Vib
from python_scipts.utils.units import CM_INV

@dataclass
class LennardJones:
    d6: float
    re: float

    def no_term(self) -> str:
        return "2"
    
    def powers(self) -> str:
        return "-12, -6"
    
    def coefs(self) -> str:
        c12 = self.d6 * self.re ** 12
        c6 = -2 * self.d6 * self.re ** 6
        return f"{c12}, {c6}"

class PotlSubItems(InputMember):
    pass

@dataclass
class CoupledLJPotentials(PotlSubItems):
    potentials: list[LennardJones]
    states: list[tuple[AngL, Vib, Vib]]

    def to_string(self) -> str:
        return f"MXLAM = {len(self.potentials)},\n LAMBDA = " \
            + ", ".join(map(lambda x: f"{x[0]}, {x[1]}, {x[2]}", self.states)) \
            + ",\n NTERM = " + ", ".join(map(lambda x: f"{x.no_term()}", self.potentials)) \
            + ",\n NPOWER = " + ", ".join(map(lambda x: f"{x.powers()}", self.potentials)) \
            + ",\n A = " + ", ".join(map(lambda x: f"{x.coefs()}", self.potentials))

@dataclass
class PotlItems(InputMember):
    potential: PotlSubItems
    
    def to_string(self) -> str:
        return f"""&POTL
            EPSIL = {1 / CM_INV},
            {self.potential.to_string()}"""