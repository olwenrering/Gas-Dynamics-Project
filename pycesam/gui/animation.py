#!/usr/bin/env python3
# coding: utf-8

# Code d'Evolution Stellaire, Adaptatif et Modulaire for the 2020 decade.
# Copyright (c) 1997-2023 The Cesam2k20 authors
# SPDX-License-Identifier : GPL-3.0-or-later

import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from pycesam import *
from tqdm import tqdm


class VarAnimation:
    """
    Class for showing an animation of the evolution of a given variable 
    """

    def __init__(self, mdl, n_ts, xx, yy, var=None, xvar=None, fig=None, ax=None, xmin=None, xmax=None, ymin=None, ymax=None, **kwargs):

        self.mdl = mdl
        self.n_ts = n_ts
        self.xx = xx
        self.yy = yy

        self.xmin = xmin*0.98 if xmin else np.nanmin(xx[-1])*0.98
        self.xmax = xmax*1.02 if xmax else np.nanmax(xx[-1])*1.02
        self.ymin = ymin*0.98 if ymin else np.nanmin(yy[-1])*0.98
        self.ymax = ymax*1.02 if ymax else np.nanmax(yy[-1])*1.02

        self.axl = ax[0]
        self.axr = ax[1]
        self.plvar, = self.axl.plot(xx[0], yy[0], marker="o", markersize=2, c="tab:blue")
        self.hr, = self.axr.plot(mdl.log_teff, mdl.log_l, c="tab:orange", linewidth=2, zorder=0)
        self.cur = self.axr.scatter(mdl.log_teff[0], mdl.log_l[0], s=25, c="tab:blue", zorder=5)
        self.fig = fig

        if xvar: self.axl.set_xlabel(xvar)
        if var: self.axl.set_ylabel(var)
        self.axr.set_xlabel(r"$\log_{10} T_\mathrm{eff}$")
        self.axr.set_ylabel(r"$\log_{10} L/L_\odot$")
        self.axr.invert_xaxis()

        self.animation = animation.FuncAnimation(
            self.fig, self.update, frames=self.n_ts, interval=50, blit=True, repeat=True, repeat_delay=500)
        self.paused = False

        fig.canvas.mpl_connect('button_press_event', self.toggle_pause)

    def toggle_pause(self, *args, **kwargs):
        if self.paused:
            self.animation.resume()
        else:
            self.animation.pause()
        self.paused = not self.paused

    def update(self, i):

        self.plvar.set_xdata(self.xx[i])
        self.plvar.set_ydata(self.yy[i])
        self.axl.set_xlim([self.xmin, self.xmax])
        self.axl.set_ylim([self.ymin, self.ymax])
        self.cur.set_offsets([self.mdl.log_teff[i], self.mdl.log_l[i]])
        self.fig.tight_layout()
        lg = self.axr.legend(handles=[self.hr], labels=[f"{self.mdl.age[i]:.2f} Myrs"])
        return self.plvar, self.cur, lg
    
    def save(self, name):

        writer = animation.PillowWriter(fps=15, metadata=dict(artist='Me'), bitrate=1800)
        self.animation.save(name+'.gif', progress_callback= tqdm( range( self.n_ts ), desc=f'Reading *{self.mdl.name}*.osc files' ), writer=writer)