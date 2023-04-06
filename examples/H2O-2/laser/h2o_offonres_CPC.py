#!/usr/bin/env python
from numpy import loadtxt as Read
import matplotlib as mpl
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvas
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"

Blue = '#0000ff'
Red = '#ff0000'

x11, y11 = Read('Offres-tdlda/pescel.H2O', skiprows=1, usecols=(0,2), unpack=True)
x12, y12 = Read('Offres-rta/pescel.H2O', skiprows=1,usecols=(0,2), unpack=True)
x21, y21 = Read('Onres-tdlda/pescel.H2O', skiprows=1, usecols=(0,2), unpack=True)
x22, y22 = Read('Onres-rta/pescel.H2O', skiprows=1, usecols=(0,2), unpack=True)

x31, y31 = Read('Offres-tdlda/penergies.H2O', skiprows=13, usecols=(0,10), unpack=True)
x32, y32 = Read('Offres-rta/penergies.H2O', skiprows=13, usecols=(0,10), unpack=True)
x41, y41 = Read('Onres-tdlda/penergies.H2O', skiprows=13, usecols=(0,10), unpack=True)
x42, y42 = Read('Onres-rta/penergies.H2O', skiprows=13, usecols=(0,10), unpack=True)
y31 = 13.605693123*y31 # Convert from Rydbergs to electronvolts
y32 = 13.605693123*y32 # Convert from Rydbergs to electronvolts
y41 = 13.605693123*y41 # Convert from Rydbergs to electronvolts
y42 = 13.605693123*y42 # Convert from Rydbergs to electronvolts

x51, y51 = Read('Offres-tdlda/pdip.H2O', skiprows=5, usecols=(0,3), unpack=True)
x52, y52 = Read('Offres-rta/pdip.H2O', skiprows=6, usecols=(0,3), unpack=True)
x61, y61 = Read('Onres-tdlda/pdip.H2O', skiprows=5, usecols=(0,3), unpack=True)
x62, y62 = Read('Onres-rta/pdip.H2O', skiprows=5, usecols=(0,3), unpack=True)

fig = Figure()
FigCanvas(fig)

Part11 = r'\begin{align*} '
Part12 = r'\omega_\mathrm{laser} &= 10.2~\mathrm{eV} \\'
Part13 = r'I_\mathrm{laser} &= 5.6\times 10^{13}~\mathrm{W/cm^2} \\'
Part14 = r'T_\mathrm{pulse} &= 36~\mathrm{fs}'
Part15 = r'\end{align*}'
LabelOnres = Part11 + Part12 + Part13 + Part14 + Part15
Part21 = r'\omega_\mathrm{laser} &= 11.4~\mathrm{eV} \\'
LabelOffres = Part11 + Part21 + Part13 + Part14 + Part15

ax1 = fig.add_subplot(321)
ax1.plot(x11, y11, label=r'TDLDA', color=Blue)
ax1.plot(x12, y12, label=r'RTA', color=Red)
ax1.set_ylabel(r'Ionization')
ax1.set_xlim(0, 73)
ax1.set_ylim(0, 0.058)
ax1.yaxis.set_ticks_position('both')
ax1.tick_params(labelbottom=False)
ax1.set_yticks([0, 0.01, 0.02, 0.03, 0.04, 0.05])
ax1.set_yticklabels(['0', '0.01', '0.02', '0.03', '0.04', '0.05'])
ax1.text(30, 0.00375, LabelOffres, fontsize=16)

ax2 = fig.add_subplot(322)
ax2.plot(x21, y21, label=r'TDLDA', color=Blue)
ax2.plot(x22, y22, label=r'RTA', color=Red)
ax2.set_xlim(0, 73)
ax2.set_ylim(0, 0.2)
ax2.yaxis.set_ticks_position('both')
ax2.tick_params(labelbottom=False)
ax2.tick_params(labelleft=False)
ax2.tick_params(labelright=True)
ax2.set_yticks([0.05, 0.1, 0.15])
ax2.set_yticklabels(['0.05', '0.10', '0.15'])
ax2.text(30, 0.015, LabelOnres, fontsize=16)

ax3 = fig.add_subplot(323, sharex=ax1)
ax3.plot(x31, y31, label=r'TDLDA', color=Blue)
ax3.plot(x32, y32, label=r'RTA', color=Red)
ax3.set_ylabel(r'$E_\mathrm{abs}$ (eV)')
ax3.set_xlim(0, 73)
ax3.set_ylim(0, 2)
ax3.yaxis.set_ticks_position('both')
ax3.tick_params(labelbottom=False)
ax3.set_yticks([0, 0.5, 1.0, 1.5])
ax3.set_yticklabels(['0', '0.5', '1.0', '1.5'])

ax4 = fig.add_subplot(324, sharex=ax2)
ax4.plot(x41, y41, label=r'TDLDA', color=Blue)
ax4.plot(x42, y42, label=r'RTA', color=Red)
ax4.set_xlim(0, 73)
ax4.set_ylim(0, 10)
ax4.yaxis.set_ticks_position('both')
ax4.tick_params(labelbottom=False)
ax4.tick_params(labelleft=False)
ax4.tick_params(labelright=True)
ax4.set_yticks([2, 4, 6, 8])
ax4.set_yticklabels(['2', '4', '6', '8'])

ax5 = fig.add_subplot(325, sharex=ax1)
ax5.plot(x51, y51, label=r'TDLDA', color=Blue)
ax5.plot(x52, y52, label=r'RTA', color=Red)
ax5.set_ylabel(r'$D_z$ ($a_0$)')
ax5.set_xlabel(r'Time (fs)')
ax5.set_xlim(0, 73)
ax5.set_ylim(-0.05, 0.05)
ax5.yaxis.set_ticks_position('both')
ax5.set_yticks([-0.025, 0, 0.025])
ax5.set_yticklabels(['-0.025', '0', '0.025'])
ax5.legend(loc='lower right')

ax6 = fig.add_subplot(326, sharex=ax2)
ax6.plot(x61, y61, label=r'TDLDA', color=Blue)
ax6.plot(x62, y62, label=r'RTA', color=Red)
ax6.set_xlabel(r'Time (fs)')
ax6.set_xlim(0, 73)
ax6.set_ylim(-0.1, 0.1)
ax6.yaxis.set_ticks_position('both')
ax6.tick_params(labelleft=False)
ax6.tick_params(labelright=True)
ax6.set_yticks([-0.05, 0, 0.05])
ax6.set_yticklabels(['-0.05', '0', '0.05'])

fig.subplots_adjust(hspace=0, wspace=0.1)
fig.savefig('h2o_offonres_CPC', bbox_inches='tight')
