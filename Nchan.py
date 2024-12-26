import numpy as np
from python_scipts.basis_items import BasisItems, InteractionType, VibrotorLevels
from python_scipts.input_creator import Input
from python_scipts.input_items import DistanceType, Energies, InputItems, PrintFlag, PropagationConfig, PropagatorType, red_mass
from python_scipts.molscat_exec import CygwinExecutor
from python_scipts.potl_items import CoupledLJPotentials, LennardJones, PotlItems
from python_scipts.utils.types import Interval
from python_scipts.utils.units import ANGS, CM_INV, KELVIN

if __name__ == "__main__":
    int_type = InteractionType.VibrotorAtom
    
    N = 50

    wells = np.linspace(0.0019, 0.0022, N)

    levels = [(well / 0.0019 - 1.0) * KELVIN for well in wells]
    states = [(0, i + 1) for i in range(N)]
    vib_levels = VibrotorLevels(levels, states)

    potentials = []
    states = []
    for v in range(1, N + 1):
        states.append((0, v, v))
        potentials.append(LennardJones(wells[v - 1], 9.))
    for v in range(1, N):
        states.append((0, v, v + 1))
        potentials.append(LennardJones(wells[v] / 10, 6.))
    potential = CoupledLJPotentials(potentials, states)

    basis = BasisItems(int_type, vib_levels)

    interval = Interval(6.5, 10., 1e3)
    near_propagator = PropagatorType.DiabaticLogDerivative.with_config(1e-3, DistanceType.ShortRange)
    far_propagator = PropagatorType.Airy.with_config(1e-7, DistanceType.LongRange)

    propagators = PropagationConfig(interval, near_propagator, far_propagator)

    input = InputItems(
        f"{N}chan_scattering",
        red_mass(6.015122, 7.016004),
        PrintFlag.Default,
        propagators,
        Energies([1e-7 * KELVIN])
    )

    potential = PotlItems(potential)

    input = Input(input, basis, potential)

    CygwinExecutor().execute(input, f"{N}chan")