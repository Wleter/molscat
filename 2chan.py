from python_scipts.basis_items import BasisItems, InteractionType, VibrotorLevels
from python_scipts.input_creator import Input
from python_scipts.input_items import DistanceType, Energies, InputItems, PrintFlag, PropagationConfig, PropagatorType
from python_scipts.molscat_exec import CygwinExecutor
from python_scipts.potl_items import CoupledLJPotentials, LennardJones, PotlItems
from python_scipts.utils.types import Interval
from python_scipts.utils.units import ANGS, CM_INV, KELVIN

if __name__ == "__main__":
    int_type = InteractionType.VibrotorAtom

    levels = [0, 10 * CM_INV]
    states = [(0, 1), (0, 2)]
    vib_levels = VibrotorLevels(levels, states)

    basis = BasisItems(int_type, vib_levels)

    interval = Interval(0.8 * ANGS, 1.5 * ANGS, 6. * ANGS)
    near_propagator = PropagatorType.QDiabaticLogDerivative.with_config(0.001 * ANGS, DistanceType.ShortRange)
    far_propagator = PropagatorType.Airy.with_config(1e-5, DistanceType.LongRange)

    propagators = PropagationConfig(interval, near_propagator, far_propagator)

    input = InputItems(
        "2chan_scattering",
        20.,
        PrintFlag.Default,
        propagators,
        Energies([1e-7 * KELVIN])
    )

    potentials = []
    states = []
    N = 2
    for v in range(1, N + 1):
        states.append((0, v, v))
        potentials.append(LennardJones(400. * CM_INV, 5. * ANGS))
    for v in range(1, N):
        states.append((0, v, v + 1))
        potentials.append(LennardJones(40. * CM_INV, 5. * ANGS))

    potential = CoupledLJPotentials(potentials, states)
    potential = PotlItems(potential)

    input = Input(input, basis, potential)

    CygwinExecutor().execute(input, "2chan")