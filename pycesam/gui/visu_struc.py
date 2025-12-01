#!/usr/bin/env python3
# coding: utf-8

# Code d'Evolution Stellaire, Adaptatif et Modulaire for the 2020 decade.
# Copyright (c) 1997-2023 The Cesam2k20 authors
# SPDX-License-Identifier : GPL-3.0-or-later


from pycesam.gui import *
from pycesam.gui.mpl import MPLFigureEditor, MPLHandler

from traits.api import List, Str
from traitsui.api import EnumEditor

from matplotlib.ticker import AutoMinorLocator


class PlotStruc(HasTraits):
    mod = Instance(CModelGUI)
    figure = Instance(plt.figure.Figure, ())
    xnames = List(comparison_mode=2)
    ynames = List(comparison_mode=2)
    ntot = Int(1)
    imod = Range(0, 1000000, 0, label='Model')
    modnames = List(comparison_mode=2)
    all_osc = Bool(False)
    age = Str()
    Xc = Str()
    xn = Str('r')
    yn = Str('m/M')
    x_log = Bool(False)
    x_inv = Bool(False)
    x_norm = Bool(False)
    y_log = Bool(False)
    y_inv = Bool(False)
    y_norm = Bool(False)

    traits_view = View(HGroup(Item('age', style='readonly'),
        Item('Xc', style='readonly')),
        Item('figure', editor=MPLFigureEditor(),
             show_label=False, width=1200,
         height=925,
         resizable=True),
        Item('imod', visible_when='all_osc',
             editor=RangeEditor(mode='spinner', high_name='ntot')),
        HGroup(
            HGroup(Item('xn', editor=EnumEditor(name='xnames')),
                   Item('x_log', label='log'),
                   Item('x_inv', label='inv')),
            HGroup(Item('yn', editor=EnumEditor(name='ynames')),
                   Item('y_log', label='log'),
                   Item('y_inv', label='inv')),
            ),
        width=1200, height=925,
        resizable=True, handler=MPLHandler()
        )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def __init__(self, mod, N2=False):

        self.mod = mod
        if '_plato' in self.mod.params.nom_output_:
            raise ValueError("Cesam2k20 plato file only contains accoustic quantities.\n \
                plot_struct cannot be used on them")
        self.xnames = mod.unom_vars + ['tau','r_acou']
        self.ynames = mod.unom_vars + ['tau','r_acou']
        self.xvar = self.xnames.index(self.xn)
        self.yvar = self.xnames.index(self.yn)
        self.all_osc = isinstance(mod.var, list)
        self.set_xlims()
        self.set_ylims()

        if self.all_osc:
            self.ntot = len(mod.var)-1
            self.modnames = list(map(str, range(self.ntot)))
            self.age = '%f Myr' % self.mod.age[self.imod]
            self.Xc = '%f' % self.mod.ab_c['H1'][self.imod]
            if N2:
                for i in range(self.ntot):
                    ra = self.mod.var[i][0][:]
                    ma = self.mod.var[i][1][:]
                    grav = self.ctes.ggrav*ma*self.mod.glob[i][0]/ra**3
                    self.mod.var[i][14][:] *= grav
        else:
            if N2:
                ra = self.mod.var[0][:]
                ma = self.mod.var[1][:]
                grav = self.ctes.ggrav*ma*self.mod.glob[0]/ra**3
                self.mod.var[14][:] *= grav

        if N2:
            self.mod.nom_vars[14] = r'$N^2$'
            self.mod.unom_vars[14] = u'N\u00b2'

        self.figure = plt.figure()
        self.figure.subplots_adjust(bottom=0.05,top=0.99,left=0.05,right=0.99,hspace=0.02,wspace=0.22)

        if self.all_osc:
            ################  HR ######################
            self.ax_hr = self.figure.add_subplot(431)
            self.format_ticks( self.ax_hr )
            xmin = min(self.mod.log_teff)-0.001*min(self.mod.log_teff)
            xmax = max(self.mod.log_teff)+0.001*max(self.mod.log_teff)
            ymin = min(self.mod.log_l)-0.01*min(self.mod.log_l)
            ymax = max(self.mod.log_l)+0.01*max(self.mod.log_l)
            self.ax_hr.set_xlim(xmax, xmin)
            self.ax_hr.set_ylim(ymin, ymax)
            self.ax_hr.set_xlabel(r'$\log T_{\rm eff}$')
            self.ax_hr.set_ylabel(r'$\log L/L_{\odot}$')
            self.line_hr, = self.ax_hr.plot(self.mod.log_teff,
                                            self.mod.log_l)

            self.pt_hr, = self.ax_hr.plot([self.mod.log_teff[self.imod]],
                                          [self.mod.log_l[self.imod]], 'o')
            ################  choice  ######################
            self.ax = self.figure.add_subplot(432)
            self.format_ticks( self.ax )
            self.line, = self.ax.plot(self.mod.var[self.imod][self.xvar],
                                      self.mod.var[self.imod][self.yvar])
            ################  grad T  ######################
            self.ax_gradT = self.figure.add_subplot(433)
            self.format_ticks( self.ax_gradT )
            self.line_gradT1, = self.ax_gradT.plot(self.mod.var[self.imod][self.xvar],
                                      self.mod.var[self.imod][10],label='$\\nabla_\mathrm{ad}$')
            self.line_gradT2, = self.ax_gradT.plot(self.mod.var[self.imod][self.xvar],
                                      self.mod.var[self.imod][21],label='$\\nabla_\mathrm{rad}$')
            self.line_gradT3, = self.ax_gradT.plot(self.mod.var[self.imod][self.xvar],
                                      self.mod.var[self.imod][5],label='$\\nabla_\mathrm{T}$')
            self.ax_gradT.set_xlabel(self.mod.nom_vars[self.xvar])
            self.ax_gradT.set_ylabel(r'$\nabla T$')
            self.ax_gradT.legend(loc=3)
            self.ax_gradT.set_ylim([0.1,0.5])
            ################  X, Y, Z  ##################
#            self.ax_XYZ = self.figure.add_subplot(456, sharex=self.ax_hr)
#            self.ax_XYZ.yaxis.set_minor_locator(AutoMinorLocator(5))
#            self.ax_XYZ.xaxis.set_minor_locator(AutoMinorLocator(5))
#            self.ax_XYZ.tick_params(which='major',bottom=True, top=True, left=True, right=True, direction='in')
#            self.ax_XYZ.tick_params(which='minor',bottom=True, top=True, left=True, right=True, direction='in')
#            self.line_XYZ1, = self.ax_XYZ.plot(self.mod.log_teff,
#                                      self.mod.ab_s[0],label='$X$')
#            self.line_XYZ2, = self.ax_XYZ.plot(self.mod.log_teff,
#                                      self.mod.ab_s[1]+self.mod.ab_s[2],label='$Y$')
#            self.line_XYZ3, = self.ax_XYZ.plot(self.mod.log_teff,
#                                      1.-self.mod.ab_s[0]-self.mod.ab_s[1]-self.mod.ab_s[2],label='$Z$')
#            self.pt_XYZ1, = self.ax_XYZ.plot([self.mod.log_teff[self.imod]],
#                                          self.mod.ab_s[0,self.imod], 'o')
#            self.pt_XYZ2, = self.ax_XYZ.plot([self.mod.log_teff[self.imod]],
#                                          self.mod.ab_s[1,self.imod]+self.mod.ab_s[2,self.imod], 'o')
#            self.pt_XYZ3, = self.ax_XYZ.plot([self.mod.log_teff[self.imod]],
#                                          1.-self.mod.ab_s[0,self.imod]-self.mod.ab_s[1,self.imod]-self.mod.ab_s[2,self.imod], 'o')
#            self.ax_XYZ.set_xlabel(r'$\log T_{\rm eff}$')
#            self.ax_XYZ.set_ylabel(r'$X, Y, Z$')
#            self.ax_XYZ.legend()
            ################  [Fe/H]  ##################
            self.FeHZ=np.log10((1.-self.mod.ab_s['H1']-self.mod.ab_s['He3']-self.mod.ab_s['He4'])/self.mod.ab_s['H1'])-np.log10(0.01781)-0.0045
            if 'Fe56' in self.mod.ab_s:
                self.FeHtrue=np.log10(self.mod.ab_s['Fe56']/self.mod.ab_s['H1'])-np.log10(56.)+4.55
#
            self.ax_FeH = self.figure.add_subplot(4,3,4, sharex=self.ax_hr)
            self.format_ticks( self.ax_FeH )
            self.line_FeH1, = self.ax_FeH.plot(self.mod.log_teff,
                                      self.FeHZ,label='from Z')
            if 'Fe56' in self.mod.ab_s:
                self.line_FeH2, = self.ax_FeH.plot(self.mod.log_teff,
                                      self.FeHtrue,label='true')
            self.pt_FeH1, = self.ax_FeH.plot([self.mod.log_teff[self.imod]],
                                          self.FeHZ[self.imod], 'o')
            if 'Fe56' in self.mod.ab_s:
                self.pt_FeH2, = self.ax_FeH.plot([self.mod.log_teff[self.imod]],
                                         self.FeHtrue[self.imod], 'o')
            self.ax_FeH.set_xlabel(r'$\log T_{\rm eff}$')
            self.ax_FeH.set_ylabel(r'[Fe/H]')
            self.ax_FeH.legend()

#            self.FeHZ=np.log10((1.-self.mod.ab_s[0]-self.mod.ab_s[1]-self.mod.ab_s[2])/self.mod.ab_s[0])-np.log10(0.01781)-0.0045
#            self.FeHtrue=np.log10(self.mod.ab_s[18]/self.mod.ab_s[0])-np.log10(56.)+4.55
#
#            self.ax_FeH = self.figure.add_subplot(4,3,4, sharex=self.ax_hr)
#            self.format_ticks( self.ax_FeH )
#            self.line_FeH1, = self.ax_FeH.plot(self.mod.log_teff,
#                                      self.FeHZ,label='from Z')
#            self.line_FeH2, = self.ax_FeH.plot(self.mod.log_teff,
#                                      self.FeHtrue,label='true')
#            self.pt_FeH1, = self.ax_FeH.plot([self.mod.log_teff[self.imod]],
#                                          self.FeHZ[self.imod], 'o')
#            self.pt_FeH2, = self.ax_FeH.plot([self.mod.log_teff[self.imod]],
#                                         self.FeHtrue[self.imod], 'o')
#            self.ax_FeH.set_xlabel(r'$\log T_{\rm eff}$')
#            self.ax_FeH.set_ylabel(r'[Fe/H]')
#            self.ax_FeH.legend()
            ################  Mcz  ##################
            self.Mcz=[]
            for k in range(len(self.mod.age)):
                  if len(self.mod.mcz[k])==1:
                        self.Mcz.append(self.mod.mcz[k][0])
                  if len(self.mod.mcz[k])>1:
                        self.Mcz.append(self.mod.mcz[k][len(self.mod.mcz[k])-1])
#
            self.ax_Mcz = self.figure.add_subplot(4,3,7, sharex=self.ax_hr)
            self.format_ticks( self.ax_Mcz )
            self.line_Mcz, = self.ax_Mcz.plot(self.mod.log_teff, self.Mcz)
            self.pt_Mcz, = self.ax_Mcz.plot([self.mod.log_teff[self.imod]], self.Mcz[self.imod], 'o')
            self.ax_Mcz.set_xlabel(r'$\log T_{\rm eff}$')
            self.ax_Mcz.set_ylabel(r'$M_{cz}/M^\ast$')
            ################  Rcz  ##################
            self.Rcz=[]
            for k in range(len(self.mod.age)):
                  if len(self.mod.rcz[k])==1:
                        self.Rcz.append(self.mod.rcz[k][0]/self.mod.rstar[k])
                  if len(self.mod.rcz[k])>1:
                        self.Rcz.append(self.mod.rcz[k][len(self.mod.rcz[k])-1]/self.mod.rstar[k])
#
            self.ax_Rcz = self.figure.add_subplot(4,3,10, sharex=self.ax_hr)
            self.format_ticks( self.ax_Rcz )
            self.line_Rcz, = self.ax_Rcz.plot(self.mod.log_teff,self.Rcz)
            self.pt_Rcz, = self.ax_Rcz.plot([self.mod.log_teff[self.imod]],self.Rcz[self.imod], 'o')
            self.ax_Rcz.set_xlabel(r'$\log T_{\rm eff}$')
            self.ax_Rcz.set_ylabel(r'$R_{cz}/R^\ast$')
            ################  sound speed  ##################
            self.cs=(self.mod.var[self.imod][9]*self.mod.var[self.imod][3]/self.mod.var[self.imod][4])**0.5 #Gamma1*P/rho
            ################  Tau (profondeur acoustique) ##################
            self.Tau=[]
            tau_temp=0.
            self.Tau.append(tau_temp)
            for k in range(len(self.mod.var[self.imod][0])-1):  
                  tau_temp=tau_temp+(1/self.cs[k+1]+1/self.cs[k])*0.5*(self.mod.var[self.imod][0][k]-self.mod.var[self.imod][0][k+1])
                  self.Tau.append(tau_temp) #profondeur acoustique
                  #if 1.-self.mod.var[self.imod][1][k]<=Mcz/self.mod.mstar[self.imod]:
                  #      Tau_cz=tau_temp #profondeur acoustique de la cz
                  #      ind_cz=k
            ################  racou (rayon acoustique)  ##################
            self.racou=[]
            for k in range(len(self.Tau)):
                  self.racou.append(self.Tau[len(self.Tau)-1]-self.Tau[k])
        #T_cz=Tau[len(Tau)-1]-Tau_cz #rayon acoustique
        #T_czsT0=T_cz/Tau[len(Tau)-1] #ratio rayon acoustique / rayon total
        #dcs2dr=derivative(self.mod.var[j][0],cs**2)
            ################  mu et gradmu  ##################
            mui = np.array([self.mod.var[self.imod][i] for i in [22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 37, 38, 39, 40]])
            denomi = np.array([1.0, 3.0, 4.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 20.0, 23.0, 24.0, 27.0, 28.0, 30.0, 32.0, 40.0, 56.0])
            self.muI=1.0 / np.sum( mui / denomi[:, None], axis=0)
            self.muE=2./(1.+self.mod.var[self.imod][22])
            self.mu=1./(1./self.muI+1./self.muE)

            self.gradmu=self.derivative(np.log(self.mod.var[self.imod][3]),np.log(self.mu))
#
            self.ax_mu = self.figure.add_subplot(4,3,5, sharex=self.ax_gradT)
            self.format_ticks( self.ax_mu )
            self.ax_mu.set_xlabel(self.mod.nom_vars[self.xvar])
            self.ax_mu.set_ylabel(r'$\mu$')
            self.line_mu, = self.ax_mu.plot(self.mod.var[self.imod][self.xvar],self.mu)
            self.ax_mu.set_ylim([0.5,0.8])
#
            self.ax_gradmu = self.figure.add_subplot(4,3,6, sharex=self.ax_gradT)
            self.format_ticks( self.ax_gradmu )
            self.ax_gradmu.set_xlabel(self.mod.nom_vars[self.xvar])
            self.ax_gradmu.set_ylabel(r'$\nabla \mu$')
            self.line_gradmu, = self.ax_gradmu.plot(self.mod.var[self.imod][self.xvar],self.gradmu)
            self.ax_gradmu.set_ylim([-0.05,0.25])
            ################  Gamma1  ######################
            self.ax_G1 = self.figure.add_subplot(4,3,8, sharex=self.ax_gradT)
            self.format_ticks( self.ax_G1 )
            self.ax_G1.set_xlabel(self.mod.nom_vars[self.xvar])
            self.ax_G1.set_ylabel(r'$\Gamma_1$')
            self.line_G1, = self.ax_G1.plot(self.mod.var[self.imod][self.xvar],
                                      self.mod.var[self.imod][9])
            self.ax_G1.set_ylim([1.5,1.7])
            ################  gradGamma1  ######################
            self.gradG1=self.derivative(np.log(self.mod.var[self.imod][3]),np.log(self.mod.var[self.imod][9]))
#
            self.ax_gradG1 = self.figure.add_subplot(4,3,9, sharex=self.ax_gradT)
            self.format_ticks( self.ax_gradG1 )
            self.ax_gradG1.set_xlabel(self.mod.nom_vars[self.xvar])
            self.ax_gradG1.set_ylabel(r'$\nabla \Gamma_1$')
            self.line_gradG1, = self.ax_gradG1.plot(self.mod.var[self.imod][self.xvar],self.gradG1)
            self.ax_gradG1.set_ylim([-0.01,0.02])
            ################  sound speed  ##################
            self.ax_cs = self.figure.add_subplot(4,3,11, sharex=self.ax_gradT)
            self.format_ticks( self.ax_cs )
            self.ax_cs.set_xlabel(self.mod.nom_vars[self.xvar])
            self.ax_cs.set_ylabel(r'$c_s$ (cm.s$^{-1})$')
            self.line_cs, = self.ax_cs.plot(self.mod.var[self.imod][self.xvar],self.cs)
            ################  dcs2dTau  ##################
            self.cs2=(self.mod.var[self.imod][9]*self.mod.var[self.imod][3]/self.mod.var[self.imod][4])
            self.dcs2dTau=self.derivative(self.Tau,self.cs2)
#
            self.ax_dcs2dTau = self.figure.add_subplot(4,3,12, sharex=self.ax_gradT)
            self.format_ticks( self.ax_dcs2dTau )
            self.ax_dcs2dTau.set_xlabel(self.mod.nom_vars[self.xvar])
            self.ax_dcs2dTau.set_ylabel(r'$dc_s^2/d\tau$')
            self.line_dcs2dTau, = self.ax_dcs2dTau.plot(self.mod.var[self.imod][self.xvar],self.dcs2dTau)
############################################################


        else:
            self.ax = self.figure.add_subplot(111)
            self.line, = self.ax.plot(self.mod.var[self.xvar],
                                      self.mod.var[self.yvar])

        self.ax.set_xlabel(self.mod.nom_vars[self.xvar])
        self.ax.set_ylabel(self.mod.nom_vars[self.yvar])
        self.list_ax = [self.ax, self.ax_gradT, self.ax_mu, self.ax_gradmu, self.ax_cs, self.ax_dcs2dTau, \
            self.ax_G1, self.ax_gradG1]


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def set_xlims(self):
        if self.all_osc:
          if self.xn == 'tau':
            self.xmin = min(self.Tau)
            self.xmax = max(self.Tau)
          elif self.xn == 'r_acou':
            self.xmin = min(self.racou)
            self.xmax = max(self.racou)
          else:
            self.xmin = min(self.mod.var[self.imod][self.xvar])
            self.xmax = max(self.mod.var[self.imod][self.xvar])
        else:
            self.xmin = min(self.mod.var[self.xvar])
            self.xmax = max(self.mod.var[self.xvar])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def set_ylims(self):
        if self.all_osc:
            self.ymin = min(self.mod.var[self.imod][self.yvar])
            self.ymax = max(self.mod.var[self.imod][self.yvar])
        else:
            self.ymin = min(self.mod.var[self.yvar])
            self.ymax = max(self.mod.var[self.yvar])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def format_ticks(self, ax):
        ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        ax.tick_params(which='major',bottom=True, top=True, left=True, right=True, direction='in')
        ax.tick_params(which='minor',bottom=True, top=True, left=True, right=True, direction='in')

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def _xn_changed(self):
        for ax in self.list_ax:
            if self.xmax == self.xmin:
                self.xmax * 1.05
                self.xmin * 0.95
            if self.x_inv:
                ax.set_xlim(self.xmax, self.xmin)
            else:
                ax.set_xlim(self.xmin, self.xmax)

            if self.xn == 'tau' :
                ax.set_xlabel(r'$\tau$')
            elif self.xn == 'r_acou' :
                ax.set_xlabel(r'$r_\mathrm{acou}$')
            else:
                self.xvar = self.xnames.index(self.xn)
                ax.set_xlabel(self.mod.nom_vars[self.xvar])
        self.set_xlims()
        self.update()

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def _yn_changed(self):
        self.yvar = self.xnames.index(self.yn)
        self.ax.set_ylabel(self.mod.nom_vars[self.yvar])
        self.set_ylims()
        if self.y_inv:
            self.ax.set_ylim(self.ymax, self.ymin)
        else:
            self.ax.set_ylim(self.ymin, self.ymax)
        self.update()

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def _x_log_changed(self, new):
        if new:
            sc = 'log'
        else:
            sc = 'linear'
        self.ax.set_xscale(sc)
        self.update()

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def _y_log_changed(self, new):
        if new:
            sc = 'log'
        else:
            sc = 'linear'
        self.ax.set_yscale(sc)
        self.update()

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def _x_inv_changed(self):
        if self.x_inv:
            self.ax.set_xlim(self.xmax, self.xmin)
        else:
            self.ax.set_xlim(self.xmin, self.xmax)
        self.update()

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def _y_inv_changed(self):
        if self.y_inv:
            self.ax.set_ylim(self.ymax, self.ymin)
        else:
            self.ax.set_ylim(self.ymin, self.ymax)
        self.update()


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def _imod_changed(self):
        self.age = '%f Myr' % self.mod.age[self.imod]
        self.Xc = '%f' % self.mod.ab_c['H1'][self.imod]
        self.set_xlims()
        self.set_ylims()
        if self.x_inv:
            self.ax.set_xlim(self.xmax, self.xmin)
            self.ax_gradT.set_xlim(self.xmax, self.xmin)
        else:
            self.ax.set_xlim(self.xmin, self.xmax)
            self.ax_gradT.set_xlim(self.xmin, self.xmax)

        if self.y_inv:
            self.ax.set_ylim(self.ymax, self.ymin)
        else:
            self.ax.set_ylim(self.ymin, self.ymax)

        self.update()

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def update(self):
        if self.all_osc:
            if self.xn in ['tau', 'r_acou']:
                ################  sound speed  ##################
                self.cs=(self.mod.var[self.imod][9]*self.mod.var[self.imod][3]/self.mod.var[self.imod][4])**0.5 #Gamma1*P/rho
                ################  Tau (profondeur acoustique) ##################
                self.Tau=[]
                tau_temp=0.
                self.Tau.append(tau_temp)
                for k in range(len(self.mod.var[self.imod][0])-1):  
                      tau_temp=tau_temp+(1/self.cs[k+1]+1/self.cs[k])*0.5*(self.mod.var[self.imod][0][k]-self.mod.var[self.imod][0][k+1])
                      self.Tau.append(tau_temp) #profondeur acoustique
                ################  racou (rayon acoustique)  ##################
                self.racou=[]
                for k in range(len(self.Tau)):
                      self.racou.append(self.Tau[len(self.Tau)-1]-self.Tau[k])
                ############################################################

                ################  HR ######################
                self.pt_hr.set_data([self.mod.log_teff[self.imod]],
                                    [self.mod.log_l[self.imod]])
                if self.xn == 'tau':
                    xx = self.Tau
                elif self.xn == 'r_acou' :
                    xx = self.racou

                ################  Choice  ######################
                self.line.set_data(xx, self.mod.var[self.imod][self.yvar])

                ################  gradT  ######################
                self.line_gradT1.set_data(xx, self.mod.var[self.imod][10])
                self.line_gradT2.set_data(xx, self.mod.var[self.imod][21])
                self.line_gradT3.set_data(xx, self.mod.var[self.imod][5])

                ################  [Fe/H]  ######################
                self.line_FeH1.set_data(self.mod.log_teff, self.FeHZ)
                self.line_FeH2.set_data(self.mod.log_teff, self.FeHtrue)
                self.pt_FeH1.set_data([self.mod.log_teff[self.imod]], self.FeHZ[self.imod])
                self.pt_FeH2.set_data([self.mod.log_teff[self.imod]], self.FeHtrue[self.imod])
                
                ################  Mcz  ##################
                self.line_Mcz.set_data(self.mod.log_teff, self.Mcz)
                self.pt_Mcz.set_data([self.mod.log_teff[self.imod]], self.Mcz[self.imod])
                
                ################  Rcz  ##################
                self.line_Rcz.set_data(self.mod.log_teff, self.Rcz)
                self.pt_Rcz.set_data([self.mod.log_teff[self.imod]], self.Rcz[self.imod])
                
                ################  mu et gradmu  ##################
                mui = np.array([self.mod.var[self.imod][i] for i in [22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 37, 38, 39, 40]])
                denomi = np.array([1.0, 3.0, 4.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 20.0, 23.0, 24.0, 27.0, 28.0, 30.0, 32.0, 40.0, 56.0])
                self.muI=1.0 / np.sum( mui / denomi[:, None], axis=0)
                self.muE=2./(1.+self.mod.var[self.imod][22])
                self.mu=1./(1./self.muI+1./self.muE)
                self.gradmu=self.derivative(np.log(self.mod.var[self.imod][3]),np.log(self.mu))

                self.line_mu.set_data(xx,self.mu)
    
                self.line_gradmu.set_data(xx,self.gradmu)
                ################  Gamma1  ##################
                self.line_G1.set_data(xx,self.mod.var[self.imod][9])
                ################  gradGamma1  ######################
                self.gradG1=self.derivative(np.log(self.mod.var[self.imod][3]),np.log(self.mod.var[self.imod][9]))
    
                self.line_gradG1.set_data(xx,self.gradG1)
                ################  sound speed  ##################
                self.line_cs.set_data(xx,self.cs)
                ################  dcs2dTau  ##################
                self.cs2=(self.mod.var[self.imod][9]*self.mod.var[self.imod][3]/self.mod.var[self.imod][4])
                self.dcs2dTau=self.derivative(self.Tau,self.cs2)
    
                self.line_dcs2dTau.set_data(self.Tau,self.dcs2dTau)

            else: 
                ################  HR ######################
                self.pt_hr.set_data([self.mod.log_teff[self.imod]],
                                    [self.mod.log_l[self.imod]])
                ################  Choice  ######################
                self.line.set_data(self.mod.var[self.imod][self.xvar],
                                   self.mod.var[self.imod][self.yvar])
                ################  gradT  ######################
                self.line_gradT1.set_data(self.mod.var[self.imod][self.xvar],
                                          self.mod.var[self.imod][10])
                self.line_gradT2.set_data(self.mod.var[self.imod][self.xvar],
                                          self.mod.var[self.imod][21])
                self.line_gradT3.set_data(self.mod.var[self.imod][self.xvar],
                                          self.mod.var[self.imod][5])
                ################  XYZ  ######################
                # self.line_XYZ1.set_data(self.mod.log_teff,self.mod.ab_s[0])
                # self.line_XYZ2.set_data(self.mod.log_teff,self.mod.ab_s[1]+self.mod.ab_s[2])
                # self.line_XYZ3.set_data(self.mod.log_teff,1.-self.mod.ab_s[0]-self.mod.ab_s[1]-self.mod.ab_s[2])
                # self.pt_XYZ1.set_data([self.mod.log_teff[self.imod]],self.mod.ab_s[0,self.imod])
                # self.pt_XYZ2.set_data([self.mod.log_teff[self.imod]],self.mod.ab_s[1,self.imod]+self.mod.ab_s[2,self.imod])
                # self.pt_XYZ3.set_data([self.mod.log_teff[self.imod]],1.-self.mod.ab_s[0,self.imod]-self.mod.ab_s[1,self.imod]-self.mod.ab_s[2,self.imod])
                ################  [Fe/H]  ######################
                self.line_FeH1.set_data(self.mod.log_teff,self.FeHZ)
                self.line_FeH2.set_data(self.mod.log_teff,self.FeHtrue)
                self.pt_FeH1.set_data([self.mod.log_teff[self.imod]],self.FeHZ[self.imod])
                self.pt_FeH2.set_data([self.mod.log_teff[self.imod]],self.FeHtrue[self.imod])
                ################  Mcz  ##################
                self.line_Mcz.set_data(self.mod.log_teff,self.Mcz)
                self.pt_Mcz.set_data([self.mod.log_teff[self.imod]],self.Mcz[self.imod])
                ################  Rcz  ##################
                self.line_Rcz.set_data(self.mod.log_teff,self.Rcz)
                self.pt_Rcz.set_data([self.mod.log_teff[self.imod]],self.Rcz[self.imod])
                ################  mu et gradmu  ##################
                mui = np.array([self.mod.var[self.imod][i] for i in [22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 37, 38, 39, 40]])
                denomi = np.array([1.0, 3.0, 4.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 20.0, 23.0, 24.0, 27.0, 28.0, 30.0, 32.0, 40.0, 56.0])
                self.muI=1.0 / np.sum( mui / denomi[:, None], axis=0)
                self.muE=2./(1.+self.mod.var[self.imod][22])
                self.mu=1./(1./self.muI+1./self.muE)
                self.gradmu=self.derivative(np.log(self.mod.var[self.imod][3]),np.log(self.mu))
    
                self.line_mu.set_data(self.mod.var[self.imod][self.xvar],self.mu)
    
                self.line_gradmu.set_data(self.mod.var[self.imod][self.xvar],self.gradmu)
                ################  Gamma1  ##################
                self.line_G1.set_data(self.mod.var[self.imod][self.xvar],
                                          self.mod.var[self.imod][9])
                ################  gradGamma1  ######################
                self.gradG1=self.derivative(np.log(self.mod.var[self.imod][3]),np.log(self.mod.var[self.imod][9]))
    
                self.line_gradG1.set_data(self.mod.var[self.imod][self.xvar],self.gradG1)
                ################  sound speed  ##################
                self.cs=(self.mod.var[self.imod][9]*self.mod.var[self.imod][3]/self.mod.var[self.imod][4])**0.5 #Gamma1*P/rho
    
                self.line_cs.set_data(self.mod.var[self.imod][self.xvar],self.cs)
                ################  dcs2dTau  ##################
                self.cs2=(self.mod.var[self.imod][9]*self.mod.var[self.imod][3]/self.mod.var[self.imod][4])
                self.dcs2dTau=self.derivative(self.Tau,self.cs2)
    
                self.line_dcs2dTau.set_data(self.mod.var[self.imod][self.xvar],self.dcs2dTau)
                ##################################
        else:
            self.line.set_data(self.mod.var[self.xvar], self.mod.var[self.yvar])
        self.figure.canvas.draw()

    def derivative(self,x,y):
        dydx=[]
        for i in range(len(x)):
            if i==0:
                dydx.append((y[i+1]-y[i])/(x[i+1]-x[i]))
            elif i==len(x)-1:
                dydx.append((y[i]-y[i-1])/(x[i]-x[i-1]))
            else:
                dydx.append((y[i+1]-y[i-1])/(x[i+1]-x[i-1]))
        return dydx


def plot_struc(mod, N2=False):
    pl = PlotStruc(mod, N2=N2)
    pl.configure_traits()
