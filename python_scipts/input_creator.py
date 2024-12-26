from dataclasses import dataclass
from python_scipts.basis_items import BasisItems
from python_scipts.potl_items import PotlItems
from python_scipts.input_items import InputItems

@dataclass
class Input:
    input_items: InputItems
    basis_items: BasisItems
    potl_items: PotlItems

    def to_string(self) -> str:
        string = f"""{self.input_items.to_string()},
                /
                {self.basis_items.to_string()},
                /
                {self.potl_items.to_string()},
                /
                """

        trimmed_string = "\n".join(line.strip() for line in string.splitlines())

        return trimmed_string