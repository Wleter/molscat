from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure

def plot() -> tuple[Figure, Axes]:
    fig, ax = plt.subplots()
    ax.grid()
    ax.tick_params(which='both', direction="in")

    return fig, ax