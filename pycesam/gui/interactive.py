#!/usr/bin/env python3
# coding: utf-8

# Code d'Evolution Stellaire, Adaptatif et Modulaire for the 2020 decade.
# Copyright (c) 1997-2023 The Cesam2k20 authors
# SPDX-License-Identifier : GPL-3.0-or-later


from traits.api             import Bool, Float, Range, Int, Trait, HasTraits, Instance
from traitsui.api           import View, Item, HGroup, VGroup, Handler, RangeEditor
from traitsui.menu          import LiveButtons
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import subprocess

from pycesam.gui.mpl import MPLFigureEditor, MPLHandler


# -----------------------------------------------------------------------------

class Echelle(HasTraits):
    """
    Class for plotting frequency Echelle diagrams interactively.
    """

    dnu = Range(0.0, 1.0e6, 150.0, label=u'\u0394\u03BD')
    nu0 = Range(0.0, 1.0e6, 0.0, label=u'\u03BD\u2080')
    nu0_obs = Range(0.0, 1.0e6, 0.0, label=u'\u03BD\u2080 obs')
    figure = Instance(matplotlib.figure.Figure, ())
    dnu_min = Float(75)
    dnu_max = Float(300)
    obs = Bool(False)
    traits_view = View(Item('figure', editor=MPLFigureEditor(),
                            show_label=False),
                       VGroup(Item('dnu',
                                   editor=RangeEditor(low_name='dnu_min',
                                                      high_name='dnu_max')),
                              Item('nu0'),
                              Item('nu0_obs', visible_when='obs')),
                       width=660,
                       height=580,
                       resizable=True,
                       title='Echelle diagram',
                       handler=MPLHandler()
                       )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def __init__(self, ldeg, freqs, lsep=None, dnu=None, freq_obs=None,
                 with_lines=False, l_obs=None):

        matplotlib.use('Qt5Agg')  # wxPython does not work...

        self.lines = []
        self.freqs = []
        self.l = np.array(ldeg)
        if freq_obs is not None:
            self.obs = True
            self.lines_obs = []
            self.freqs_obs = []
            self.l_obs = l_obs

        self.lmax = max(self.l) + 1
        self.lmin = min(self.l)

        for l in range(self.lmin, self.lmax):
            self.freqs.append(freqs[self.l == l])
            if self.obs: self.freqs_obs.append(freq_obs[self.l_obs == l])

        self.figure = plt.figure()
        self.ax = self.figure.add_subplot(111)
        self.line_styles = ['bo','g^','rs','cD']

        if lsep is not None:
            nn = len(lsep[0])
            nmin = min(15, nn)
            nmax = min(25, nn)
            self.dnu = np.mean(lsep[0][15:25])
            print(f'dnu = {self.dnu} muHz')
        if dnu:
            self.dnu = dnu

        self.leg = []
        for i in range(self.lmin, self.lmax):
            if with_lines:
                line, = self.ax.plot([], [], self.line_styles[(i-self.lmin)%4]+'-', ms=8.0)
            else:
                line, = self.ax.plot([], [], self.line_styles[(i-self.lmin)%4], ms=8.0)

            self.lines.append(line)
            self.leg.append(f'$\\ell={i}$')

        if self.obs:
            for i in range(self.lmin, self.lmax):
                if with_lines:
                    line, = self.ax.plot([], [], self.line_styles[(i-self.lmin)%4]+':', ms=4.0)
                else:
                    line, = self.ax.plot([], [], self.line_styles[(i-self.lmin)%4], ms=4.0)

                self.lines_obs.append(line)


        self.ax.legend(self.leg)
        self.ax.set_xlabel(r'$\nu (\mu{\rm Hz})$')

        self.update_echelle(init=True)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def update_echelle(self, init=False):
        self.nus = []
        self.nub = []
        if self.obs:
            self.nus_obs = []
            self.nub_obs = []

        for i, j in enumerate(self.freqs):
            self.nus.append((j - self.nu0) % self.dnu)
            ind = np.array(map(int, j/self.dnu))
            self.nub.append(ind*self.dnu + self.nu0)
        if self.obs:
            for i, j in enumerate(self.freqs_obs):
                self.nus_obs.append((j - self.nu0_obs) % self.dnu)
                ind = np.array(map(int, j/self.dnu))
                self.nub_obs.append(ind*self.dnu + self.nu0_obs)

        for l in range(self.lmin, self.lmax):
            il = l - self.lmin
            self.lines[il].set_data(self.nus[il], self.nub[il])
            if self.obs:
                self.lines_obs[il].set_data(self.nus_obs[il], self.nub_obs[il])

        #        xmin = min(self.nus[0])
        if init:
            xmin = 0.0
            xmax = max(self.nus[0])
            ymin = min(self.nub[0])
            ymax = max(self.nub[0])
            self.ax.set_xlim(xmin, xmax)
            self.ax.set_ylim(ymin, ymax)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def _dnu_changed(self):
        self.dnu_min = 0.75*self.dnu
        self.dnu_max = 1.5*self.dnu
        self.update_echelle()
        self.figure.canvas.draw()

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def _nu0_changed(self):
        self.update_echelle()
        self.figure.canvas.draw()

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def _nu0_obs_changed(self):
        self.update_echelle()
        self.figure.canvas.draw()

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

class EchellePer(HasTraits):
    """
    Class for plotting period Echelle diagrams interactively.
    """

    dper = Range(0.0, 1.0e6, 150.0, label=u'\u0394 P')
    per_0 = Range(0.0, 1.0e6, 0.0, label=u'P\u2080')
    figure = Instance(matplotlib.figure.Figure, ())
    per_min = Float(0.1)
    per_max = Float(30000)
    traits_view = View(Item('figure', editor=MPLFigureEditor(),
                            show_label=False),
                       VGroup(Item('dper',
                                   editor=RangeEditor(low_name='per_min',
                                                      high_name='per_max')),
                              Item('per_0')),width=660,
                       height=580,
                       resizable=True,
                       title='Period Echelle diagram',
                       handler=MPLHandler()
                       )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def __init__(self, ldeg, pers, dper, lplot=None):
        matplotlib.use('Qt5Agg')  # wxPython does not work...

        self.lines = []
        self.pers = []
        self.l = np.array(ldeg)

        self.lmax = max(self.l) + 1
        self.lmin = min(self.l)

        for l in range(self.lmin, self.lmax):
            self.pers.append(pers[self.l == l])

        self.figure = plt.figure()
        self.ax = self.figure.add_subplot(111)
        self.line_styles = ['o','^','s','D']

        self.dper = dper

        self.leg = []
        for i in range(self.lmin, self.lmax):
            il = i - self.lmin
            line, = self.ax.plot([], [], self.line_styles[il%4])
            self.lines.append(line)
            self.leg.append(f'$\\ell={i}$')

        self.ax.legend(self.leg)
        self.ax.set_xlabel(r'$P ({\rm s})$')

        self.update_echelle()

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def update_echelle(self):
        self.ps = []
        self.pb = []
        for i, j in enumerate(self.pers):
            self.ps.append((j - self.per_0) % self.dper)
            #ind = arange(len(j))
            ind = np.array(map(int, j/self.dper))
            self.pb.append(1.0/(ind*self.dper + self.per_0)*1e6)
            # self.nub.append(j)

        for i in range(self.lmin, self.lmax):
            il = i - self.lmin
            self.lines[il].set_data(self.ps[il], self.pb[il])

#        xmin = min(self.nus[0])
        xmin = 0.0
        xmax = max(self.ps[0])
        ymin = min(self.pb[0])
        ymax = max(self.pb[0])
        self.ax.set_xlim(xmin, xmax)
        self.ax.set_ylim(ymin, ymax)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def _dper_changed(self):
        self.dper_min = 0.75*self.dper
        self.dper_max = 1.5*self.dper
        self.update_echelle()
        self.figure.canvas.draw()

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def _per_0_changed(self):
        self.update_echelle()
        self.figure.canvas.draw()


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

class PlotEigenFunction(HasTraits):
    """
    Class for plotting eigenfunctions interactively.
    """

    l = Range(0, 10000, 0)
    n = Range(-200, 200, 0)
    l_max = Int(10000)
    l_min = Int(0)
    n_max = Int(200)
    n_min = Int(0)
    figure = Instance(matplotlib.figure.Figure, ())
    hold = Bool(False)
    funct = Trait(0,
                  {u'y\u2081' : 0,
                   u'y\u2082' : 1,
                   u'y\u2083' : 2,
                   u'y\u2084' : 3,
                   u'z\u2081' : 4,
                   u'z\u2082' : 5}, label='Function to plot')

    traits_view = View(Item('figure', editor=MPLFigureEditor(),
                            show_label=False),
                       HGroup(Item('l',
                                   editor=RangeEditor(low_name='l_min',
                                                      high_name='l_max')),
                              Item('n',
                                   editor=RangeEditor(low_name='n_min',
                                                      high_name='n_max')),
                              Item('funct')),
                       Item('hold'),
                       width=660, height=580, resizable=True, handler=MPLHandler()
                       )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def __init__(self, name):
        matplotlib.use('Qt5Agg')  # wxPython does not work...

        self.name = name
        self.figure = plt.figure()
        self.ax = self.figure.add_subplot(111)
        self.line = self.ax.plot([], [])
        self.read_amde_all()

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def read_amde_all(self):
        """
        Reads an amde file. Stores the result in `self.xeig[npoints]` and
        `self.yeig[6,npoints]`.
         """


        cmd = f'{self.name}\n'

        print(f'Reading {self.name}')
        s = subprocess.Popen('convert-amde_all.x', stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE)
        res = s.communicate(input=cmd)[0].split('\n')

        if res:
            if res[0]:
                nline = 0
                imode = 0
                self.index = {}
                self.xeig = []
                self.yeig = []
                self.ls = []
                self.ns = []

                while True:
                    dat = map(int, res[nline].split())
                    ntot = dat[0]
                    self.ls.append(dat[1])
                    self.ns.append(dat[2])
                    print(f'Reading l = {dat[1]:d}, n = {dat[2]:d}')
                    self.index[(dat[1],dat[2])] = imode

                    nline += 1
                    self.xeig.append(np.zeros(ntot))
                    self.yeig.append(np.zeros((6, ntot)))
                    for i in range(ntot):
                        dat = map(float, res[nline].split())
                        self.xeig[-1][i] = dat[0]
                        self.yeig[-1][:,i] = dat[1:]
                        nline += 1

                    imode += 1
                    if nline+1 == len(res): break

                self.l_max = max(self.ls)
                self.l_min = min(self.ls)
                self.n_max = max(self.ns)
                self.n_min = min(self.ns)
                self.l = self.ls[0]
                self.n = self.ns[0]

                return

        print(f"{RED}Error:{NO_COLOR} mode not calculated.")


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def _l_changed(self):
        self.update_fig()


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def _n_changed(self):
        self.update_fig()

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def _funct_changed(self):
        self.update_fig()

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def update_fig(self):
        xmin = min(self.xeig[self.index[(self.l, self.n)]])
        xmax = max(self.xeig[self.index[(self.l, self.n)]])
        ymin = min(self.yeig[self.index[(self.l, self.n)]][self.funct_])
        ymax = max(self.yeig[self.index[(self.l, self.n)]][self.funct_])
        self.ax.set_xlim(xmin, xmax)
        self.ax.set_ylim(ymin, ymax)
        if self.hold:
            line, = self.ax.plot([], [])
            self.line.append(line)
        elif len(self.line) > 1:
            self.figure.clf()
            self.ax = self.figure.add_subplot(111)
            self.line = self.ax.plot([], [])

        self.line[-1].set_data(self.xeig[self.index[(self.l, self.n)]],
                               self.yeig[self.index[(self.l, self.n)]][self.funct_])
        self.figure.canvas.draw()

