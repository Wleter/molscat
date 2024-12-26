from dataclasses import dataclass

@dataclass
class Range:
    min_value: float
    max_value: float
    value_no: int

@dataclass
class Interval:
    min_value: float
    mid_value: float
    max_value: float