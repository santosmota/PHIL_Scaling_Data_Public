import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from matplotlib import cm
from matplotlib.ticker import LinearLocator

plt.rcParams.update(plt.rcParamsDefault)
plt.rcParams['mathtext.fontset'] = 'cm'  # 'cm' Computer modern # 'dejavuserif', 'dejavusans'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'cmr10'
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
# from matplotlib.ticker import FormatStrFormatter

###############################################################
# Annotations in 3D charts
# source: # https://gist.github.com/WetHat/1d6cd0f7309535311a539b42cccca89c
###############################################################
from mpl_toolkits.mplot3d.proj3d import proj_transform
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib.text import Annotation

class Annotation3D(Annotation):

    def __init__(self, text, xyz, *args, **kwargs):
        super().__init__(text, xy=(0, 0), *args, **kwargs)
        self._xyz = xyz

    def draw(self, renderer):
        x2, y2, z2 = proj_transform(*self._xyz, self.axes.M)
        self.xy = (x2, y2)
        super().draw(renderer)

def _annotate3D(ax, text, xyz, *args, **kwargs):
    '''Add anotation `text` to an `Axes3d` instance.'''

    annotation = Annotation3D(text, xyz, *args, **kwargs)
    ax.add_artist(annotation)

setattr(Axes3D, 'annotate3D', _annotate3D)


###############################################################
# Makes the surface curve of the error in l_{r}^{sd}
###############################################################
def superficies():
    ###############################
    # system data
    ###############################
    Fn = 50.0                           # frequency in Hz
    wn = Fn * 2 * np.pi                 # angular frequency in radians per second

    ###############################
    # fs = Full Scale converter
    Sbfs = 5e6                          # BESS rated power in VA
    Vbacfs = 690.0                      # BESS ac voltage
    Ibacfs = Sbfs / Vbacfs / 3**0.5     # BESS ac current

    Vbdcfs = 1100.0                     # BESS dc voltage
    Ibdcfs = Sbfs / Vbdcfs              # BESS dc current
    acdcfs = Vbacfs / Vbdcfs            # ac / dc voltage ratio

    Zbfs = Vbacfs / Ibacfs / 3**0.5     # BESS base resistance
    Lbfs = Zbfs / wn                    # BESS base inductance
    Cbfs = 1 / Zbfs / wn                # BESS base capacitance (ac side)

    # ltfs = 0.08                         # transformer pu values, not used
    # rtfs = 0.005

    Lrfs = 7.7465e-05                   # BESS reactor inductance
    lrfs = Lrfs / Lbfs                  # in pu

    # Cacfs = 0.0018386                   # capacitance, not used
    # cacfs = Cacfs / Cbfs

    ###############################
    # sd = Scaled Down converter
    Lrsd = 500e-6                       # SDC reactor
    # Ltsd = 316e-6                     # transformer
    # Cacsd = 50e-6                     # capacitance (ac)

    ################################
    # SDC Base voltages
    Vbacsd_start = 50.0
    Vbacsd_end = np.round(650*acdcfs, decimals=0)  # max limited by AC/DC ratio, and max 650Vdc

    # num=100,   # 100 samples makes a nice chart already
    Vbacsd = np.linspace(start=Vbacsd_start,
                         stop=Vbacsd_end,
                         num=100,  # num=int(Vbacsd_end - Vbacsd_start + 1),
                         endpoint=True)

    ################################
    # SDC Base currents
    Ibacsd_start = 25.0
    Ibacsd_end = 72.0

    # num=100,   # 100 samples makes a nice chart already
    Ibacsd = np.linspace(start=Ibacsd_start,
                         stop=Ibacsd_end,
                         num=100, # num=int(Ibacsd_end - Ibacsd_start + 1),
                         endpoint=True)

    ################################
    # Number of bases
    NIb = len(Ibacsd)
    NVb = len(Vbacsd)

    print('number of base voltages = ', NVb)
    print('number of base currents = ', NIb)

    ################################
    # Error plane
    err_limit = 5 * np.ones((NVb, NIb))

    ################################
    # Calculating the reactor inductance error
    x_I, y_V = np.meshgrid(Ibacsd, Vbacsd)          # creates the x and y matrices for the surface

    Zbsd = y_V / x_I / 3 ** 0.5                     # SDC base resistance (matrix)
    Lbsd = Zbsd / wn                                # SDC base inductance (matrix)

    lrsd = Lrsd / Lbsd                              # SDC l_{r}^{sd} (matrix)
    lrsd_err = 100.0 * np.abs(lrfs - lrsd) / lrfs   # SDC l_{r}^{sd} error (matrix)


    ###############################################################
    # size of the figure (in inches)
    ###############################################################
    figsizex = 5.5
    figsizey = 4.0 # 4.25

    fig, axes = plt.subplots(1, 1, sharex=True,
                             figsize=(figsizex, figsizey),
                             num='Err',
                             subplot_kw={"projection": "3d"})

    ###############################################################
    # Surface plots
    ###############################################################
    surf_err = axes.plot_surface(x_I, y_V, lrsd_err, cmap=cm.turbo, alpha=0.75,
                                 linewidth=0.0, antialiased=False)

    axes.plot_surface(x_I, y_V, err_limit, color='black', alpha=0.35, # cmap=cm.coolwarm,
                      linewidth=0.0, antialiased=False)

    ###############################################################
    # axis titles and annotations
    ###############################################################
    axes.set_ylabel(r'SDC base voltage (V)')
    axes.set_xlabel(r'SDC base current (A)')
    axes.set_zlabel(r'SDC reactor $l_r^{sd}$ error (%)')

    axes.annotate3D('81V, 72A, 5.4%', (72, 81, 5.37),
                  xytext=(30, 10),
                  textcoords='offset points',
                  bbox=dict(boxstyle="round", fc="lightyellow"),
                  arrowprops=dict(arrowstyle="->", color='black', lw=1.5))

    axes.annotate3D('74V, 72A, 5.0%', (72, 74, 5.0),
                    xytext=(-50, -20),
                    textcoords='offset points',
                    bbox=dict(boxstyle="round", fc="lightyellow"),
                    arrowprops=dict(arrowstyle="->", ec='black', lw=1.5))

    axes.annotate3D('Case 1', (72, 78, 4.0),
                    xytext=(30, -10),
                    textcoords='offset points',
                    bbox=dict(boxstyle="round", fc="white"),
                    # arrowprops=dict(arrowstyle="->", ec='black', lw=1.5)
                    )

    axes.annotate3D('363V, 72A, 78.8%', (72, 363, 80),
                    xytext=(-10, -30),
                    textcoords='offset points',
                    bbox=dict(boxstyle="round", fc="lightyellow"),
                    arrowprops=dict(arrowstyle="->", color='black', lw=1.5))

    axes.annotate3D('Case 2', (72, 390, 80),
                    xytext=(20, -10),
                    textcoords='offset points',
                    bbox=dict(boxstyle="round", fc="white"),
                    # arrowprops=dict(arrowstyle="->", color='black', lw=1.5)
                    )

    axes.annotate3D('5% error plane', (30, 350, 5),
                    xytext=(-80, -10),
                    textcoords='offset points',
                    bbox=dict(boxstyle="round", fc="lightyellow"),
                    # arrowprops=dict(arrowstyle="->", color='black', lw=1.5)
                    )

    ###############################################################
    # tight layout an show
    ###############################################################
    # axes.view_init(elev=10.0, azim=30.0, roll=0.0)
    axes.view_init(elev=22.0, azim=13.0, roll=0.0)

    clb = fig.colorbar(surf_err, shrink=0.6, aspect=20)
    clb.ax.set_title('$l_r^{sd}$ error (%)')

    fig.tight_layout()
    fig.savefig('3dplot.pdf', format='pdf')

    plt.show()


def main():
    superficies()


if __name__ == '__main__':
    main()
