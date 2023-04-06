#!/usr/bin/env python
import numpy as np
import matplotlib as mpl
from numpy import loadtxt as Read
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvas
from matplotlib.figure import Figure
from matplotlib import collections as Coll
mpl.rcParams['text.latex.preamble'] = [
r'\usepackage{mathrsfs}']

x11,y11=Read('pdip.Na41+',skiprows=5,usecols=(0,1),unpack=True)
x12,y12=Read('pdip.Na41+',skiprows=5,usecols=(0,3),unpack=True)
x21,y21=Read('pspectr.Na41+',skiprows=3,usecols=(1,3),unpack=True)
x22,y22=Read('pspectr.Na41+',skiprows=3,usecols=(1,9),unpack=True)

fig=Figure()
FigCanvas(fig)

ax1=fig.add_subplot(211)
# ax1.plot(x11,y11, color='#1f78b4', label=r'$\{x,y\}$-direction')
# ax1.plot(x12,y12, color='#e31a1c', label=r'$z$-direction')
ax1.plot(x11,y11, label=r'$\{x,y\}$-direction')
ax1.plot(x12,y12, label=r'$z$-direction')
ax1.set_xlabel(r'$\mathrm{Time}~t~/~[\mathrm{fs}]$')
ax1.set_ylabel(r'$D(t)~/~[a_0]$')
ax1.set_xlim(0,100)
ax1.set_ylim(-0.004,0.006)
ax1.set_yticks([-0.002, 0, 0.002, 0.004])
ax1.set_yticklabels(['-0.002', '0', '0.002', '0.004'])
ax1.xaxis.set_ticks_position('top')
ax1.xaxis.set_label_position('top')

ax2=fig.add_subplot(212)
# ax2.plot(x21,y21, color='#1f78b4', label=r'$\{x,y\}$-direction')
# ax2.plot(x22,y22, color='#e31a1c', label=r'$z$-direction')
ax2.plot(x21,y21, label=r'$\{x,y\}$-direction')
ax2.plot(x22,y22, label=r'$z$-direction')
ax2.axvline(x=5.29561, color='k', linestyle='--', linewidth=1)
ax2.set_xlabel(r'$\mathrm{Frequency}~\omega~/~[\mathrm{eV}]$')
ax2.set_ylabel(r'$\mathscr{F}[\widetilde{D}(t)](\omega)~/~[\mathrm{a.u.}]$')
ax2.set_xlim(2,6)
ax2.legend()
ax2.yaxis.set_ticks_position('none')
ax2.set_yticklabels([])

fig.subplots_adjust(hspace=0, wspace=0.1)
fig.savefig('Na41p', bbox_inches='tight')
