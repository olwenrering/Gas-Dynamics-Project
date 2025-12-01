#!/usr/bin/env python3
# coding: utf-8

# Code d'Evolution Stellaire, Adaptatif et Modulaire for the 2020 decade.
# Copyright (c) 1997-2023 The Cesam2k20 authors
# SPDX-License-Identifier : GPL-3.0-or-later


from pycesam import *
from pycesam.gui.cparams    import CParametersGUI#, CParametersGUI_Tk
from pycesam.gui.cesam_run  import CRunGUI
from pycesam.gui.freqs      import FreqsGUI
from pycesam.gui.CSetupGUI  import CSetupGUI
from pycesam.gui.animation  import VarAnimation
from traits.api             import Bool, DelegatesTo, Float, Range, Int, Trait
from traitsui.api           import View, Item, Tabbed, HGroup, VGroup, Group, Handler, RangeEditor
from traitsui.menu          import LiveButtons

import matplotlib as mplt
mplt.use('Qt5Agg')
import matplotlib.pyplot as plt
plt.ion()
import pycesam.gui.interactive as itr

#paramgui = CParametersGUI



class CModelGUI(CModel):
    """
    <CModelGUI.__init__(name, reinit=False, read=True, obs=False, freq=False, create=False)>

    Class that represent a Cesam2k20 evolutionary track with graphic interface.

    :param name: Name of the model
    :type name: string

    :kparam reinit: If True, deletes all calculated files.
    :ktype reinit: boolean

    :kparam read: If True, reads files .HR and .osc
    :ktype read: boolean

    :kparam obs: If True, compute observable quantities such as magnitude and colors.
    :ktype obs: boolean

    :kparam freq: If True, GUI is preset with computation of frequency and no computation of model.
    :ktype freq: boolean

    :kparam create: If True, the don file is created but Cesam2k20 is not ran.
    :ktype create: boolean

    :kparam job: Some options we want to enforce.
    :ktype job: dict

    :raise NameError: if name length is 0 or > 31 characters.
    """

    calc_mod    = Bool(True, label='Calculate model')
    calc_frq    = Bool(False, label='Calculate frequencies')

    params      = Instance(CParametersGUI)
    run         = Instance(CRunGUI)
    fparams     = Instance(FreqsGUI)

    traits_view = View(Item('calc_mod'), Item('calc_frq'),
                       Tabbed(Item('params',  label='Model Parameters',
                                   style='custom', show_label=False),
                              Item('run',     label='Run options',
                                   style='custom', show_label=False),
                              Item('fparams', label='Oscillation options',
                                   style='custom', show_label=False)),
                       buttons = LiveButtons, kind='live',
                       title='Cesam2k20 parameters',height=700, resizable=True )

    def __init__( self, name, reinit=False, read=True, obs=False, model=True, freq=False, create=False,
        job=None ):
        """
        <CModelGUI.__init__(name, reinit=False, read=True, obs=False, freq=False, create=False)>

        Class that represent a Cesam2k20 evolutionary track and allow graphical representation
        in a Jupyter Notebook.

        :param name: Name of the model.
        :type name: string

        :kparam reinit: If True, deletes all calculated files.
        :ktype reinit: boolean

        :kparam read: If True, reads files .HR and .osc
        :ktype read: boolean

        :kparam obs: If True, compute observable quantities such as magnitude and colors.
        :ktype obs: boolean

        :kparam freq: If True, GUI is preset with computation of frequency and no computation of model.
        :ktype freq: boolean

        :kparam create: If True, the don file is created but Cesam2k20 is not ran.
        :ktype create: boolean

        :kparam job: Some options we want to enforce.
        :ktype job: dict

        :raise NameError: if name length is 0 or > 31 characters.
        """

        self.calc_mod = model and not create
        self.calc_frq = freq

        super().__init__( name, reinit=reinit, read=read, job=job, fromGUI=True )

        self.__init( )

    def __str__(self):
        """
        Gives the string representation of a CModelGUI instance.

        :return: Mass, X0, Y0 and Z0 of the model if has been correctly calculated. Otherwise, only mass.
        :rtype:  string
        """
        if self.finished:
            return (f"Cesam2k20 model (GUI) \'{self.name}\'\n\tmodel calculated\n\tParameters:\n"  \
                    f"\t\tMass = {self.params.mtot:.3f}\n\t\tX0 =  {self.params.x0:.5f}\n\t\t" \
                        f"Y0 = {self.params.y0:.5f}\n\t\tZ0 =  {self.params.z0:.5f}")
        else:
            return f"Cesam2k20 model (GUI) \'{self.name}\'\n\tmodel not calculated\n\t" \
                    f"Mass = {self.params.mtot:.3f}"


    def __call__(self, mkdon=True, debug=False, edit=True, log=False, devnull=False, progress_bar=True, executable=None, **kwargs):
        """
        <CModelGUI.__call__( mkdon=True, debug=False, edit=True, log=False )>

        Calculates models or frequencies.

        :kparam mkdon: Creates (or re-creates) a .don file if `mkdon==True`; uses existing
            .don file otherwise. Default: True.
        :ktype mkdon: boolean, optional

        :kparam debug: If True, run the debug Cesam2k20 executable (`cesam2k20_dbg.x`). Default: False.
        :ktype debug: boolean, optional

        :kparam edit: If True, allow to edit .don file or oscillation parameter file.
        :ktype edit: boolean, optional

        :kparam log: If True, direct standard output to .log file. Default: False.
        :ktype log: boolean, optional

        :kparam devnull: If True, redirect standard output to /dev/null.
        :ktype devnull: boolean

        :kparam valgrind: If True, calls Cesam2k20 with valgrind analysis software.
        :ktype valgrind: boolean, optional

        :kparam coverage: If True, evaluate test coverage by calling `cesam2k20_cov.x` with gcov.
        :ktype coverage: boolean, optional
        """

        self.fparams.write_fsettings( )
        if self.run.job_ == 'freqs':
            self.calc_freqsGUI(edit=edit)

        else:
            if edit:
                ok = self.run.configure_traits()
                if not ok: return
            else:
                ok = True

            self.finished = False
            if mkdon: self.params.mkdon()
            if self.run.job_ != 'rep':
                if os.path.exists(self.hr):    os.remove(self.hr)
                if os.path.exists(self.hrnew): os.remove(self.hrnew)

            self.run.write_rsettings()

            self.result = None

            valgrind = kwargs.get( 'valgrind', False )
            coverage = kwargs.get( 'coverage', False )

            if self.run.pre_pms: self._CModel__precompute_pms(valgrind=valgrind, coverage=coverage, debug=False, log=log, devnull=devnull)

            print("\n-----------------------------------------")
            print(f"Model: {self.name}")
            print(f"    Started at {time.asctime()} ")
            print(f"Computing evolution of model {self.name}...")
            t1 = time.time()

            # runs Cesam2k20
            self.parent, self.child = mulp.Pipe()
            self.process = mulp.Process(
                target=self.run_cesam,
                kwargs={'pipe':self.child, 'debug':debug, 'log':log, 'valgrind':valgrind, 'coverage':coverage, 'devnull':devnull, 'executable':executable} )
            self.process.start()

            if progress_bar: self._CModel__display_progress_bar( )

            sortie, err = self.parent.recv( )

            self._CModel__process_out( sortie, err, t1, True, debug, log, valgrind )

    def __init( self ):
        """
        <CModel.__init( )>

        Method used to initialize all subclasses used by CModelGUI and read files if neaded.
        """

        self.params   = CParametersGUI( self.name )
        self.run      = CRunGUI( self.name, job=self.job )
        self.ctes     = CConstants( )

        self.set_osc()
        self.set_osc2d()
        self.set_rep()
        self.init_fparamsNB('p')
        # read .frun file
        if os.path.exists(self.fparams.frun_file):
            self.fparams.read_fsettings( )

        self.ctes.set( self.params )

        self.is_finished( )

        self.hr_exists = False
        if os.path.exists(self.hr) or os.path.exists(self.hrnew):
            self.hr_exists = True


        if self.reinit and self.finished:
            self._reinit()
        elif self.read:
            if self.finished: self.read_osc()
            if self.hr_exists:
                self.read_hr()

    def init_fparamsNB( self, p_or_g, keep=False ):

        """
        <CModelGUI.init_fparamsNB(self, p_or_g, keep=False)>
        Initializes parameters for frequency computations. See FreqsGUI for documentation.

        :param p_or_g: Chooses whether to calculate p or g modes.
        :type p_or_g: string

        :kparam keep: If True, keeps some options if `CModelGUI.init_fparamsNB` called again.
        :ktype keep: boolean
        """

        if not keep:
            self.fparams = FreqsGUI( self.name, self.all_osc )
            self.fparams.modes = p_or_g + '-modes'

        # adipls
        self.amdl = self.name + '.amdl'
        self.agsm = self.name + '.agsm'
        self.ssm  = self.name + '.ssm'
        self.amde = self.name + '.amde'
        self.rkr  = self.name + '.rkr'
        self.gkr  = self.name + '.gkr'

    def hr_axis(self, ax=None, nologx=False, nology=False, **kwargs):
        """
        <CModelGUI.hr_axis( ax=None )>

        Sets the axes for plotting tracks on the HRD.

        :kparam ax : axis instance, optional
        :ktype ax: matplotlib.axes._subplots.AxesSubplot

        :kparam nologx: True if the linear value instead of the log value of the x
            quantity must be plotted. Default: False.
        :ktype nologx: boolean

        :kparam nology: True if the linear value instead of the log value of the y
            quantity must be plotted. Default: False.
        :ktype nology: boolean


        :kwarg fontsize: Font size used for witting axe labels. Default value:
            `plt.rcParams['font.size']`.
        :kwtype fontsize: float

        :kwarg show_ax_label: Wether we display the ax labels or not. Default: True.
        :kwtype show_ax_label: boolean
        """

        func = lambda x, nolog: 10**x if nolog else x
        x_label = r'$T_{\rm eff}~[{\rm K}]$' if nologx else r'$\log T_{\rm eff}~[{\rm K}]$'
        y_label = r'$L / L_\odot$' if nology else r'$\log L / L_\odot$'

        imin      = kwargs.get( "imin", 0 )
        ima       = len( self.log_teff )
        imax      = min( ima, kwargs.get( "imax", ima ) )

        xx = func( self.log_teff[imin:imax], nologx )
        yy = func( self.log_l[imin:imax],    nology )

        mil, mal = min(yy), max(yy)
        add_l    = 0.03 * (mal - mil)
        add_teff = 0.03 * (max(xx) - min(xx))

        if nology:
            lims = (max(xx)+add_teff, min(xx)-add_teff,
                    mil - add_l,  mal + add_l)
        else:
            lims = (max(xx)+add_teff, min(xx)-add_teff,
                    mil - add_l,  mal + np.sign(mal) * add_l)
        fs = kwargs.get( 'fontsize', plt.rcParams['font.size'] )
        if ax:
            if ax.lines:
                xlim = ax.get_xlim()
                ylim = ax.get_ylim()
                lims = (max(lims[0], xlim[0]), min(lims[1], xlim[1]),
                        min(lims[2], ylim[0]), max(lims[3], ylim[1]))
            ax.set_xlim(lims[0], lims[1])
            ax.set_ylim(lims[2], lims[3])
            if kwargs.get( 'show_ax_label', True ):
                ax.set_xlabel(x_label, fontsize=fs)
                ax.set_ylabel(y_label, fontsize=fs)
        else:
            plt.axis(lims)
            if kwargs.get( 'show_ax_label', True ):
                plt.xlabel(x_label, fontsize=fs)
                plt.ylabel(y_label, fontsize=fs)

    def kiel_axis(self, ax=None, nologx=False, nology=False, **kwargs):
        """
        <CModelGUI.kiel_axis(ax=None)>

        Sets the axes for plotting tracks on the Kiel Diagram.

        :kparam ax: axis instance, optional
        :ktype ax: `matplotlib.axes._subplots.AxesSubplot`

        :kparam nologx: True if the linear value instead of the log value of the x
            quantity must be plotted. Default: False.
        :ktype nologx: boolean

        :kparam nology: True if the linear value instead of the log value of the y
            quantity must be plotted. Default: False.
        :ktype nology: boolean

        :kwarg fontsize: Font size used for witting axe labels. Default value:
            `plt.rcParams['font.size']`.
        :kwtype fontsize: float

        :kwarg show_ax_label: Wether we display the ax labels or not. Default: True.
        :kwtype show_ax_label: boolean
        """

        func = lambda x, nolog: 10**x if nolog else x
        xlabel = r'$T_{\rm eff}~[{\rm K}]$' if nologx else r'$\log T_{\rm eff}~[{\rm K}]$'
        ylabel = r'$g~[{\rm cm}\,{\rm s}^{-2}]$' if nology else r'$\log g~[{\rm cm}\,{\rm s}^{-2}]$'

        xx = func( self.log_teff, nologx )
        yy = func( self.log_g,    nology )

        mig, mag = min(yy), max(yy)
        add_g    = 0.03 * (mag - mig)
        add_teff = 0.03 * (max(xx) - min(xx))

        lims = (max(xx)+add_teff, min(xx)-add_teff,
                mag + np.sign(mag) * add_g, mig + np.sign(mig) * add_g)
        fs = kwargs.get( 'fontsize', plt.rcParams['font.size'] )
        if ax:
            if ax.lines:
                xlim = ax.get_xlim()
                ylim = ax.get_ylim()
                lims = (max(lims[0], xlim[0]), min(lims[1], xlim[1]),
                        max(lims[2], ylim[0]), min(lims[3], ylim[1]))
            ax.set_xlim(lims[0], lims[1])
            ax.set_ylim(lims[2], lims[3])
            if kwargs.get( 'show_ax_label', True ):
                ax.set_xlabel( xlabel, fontsize=fs )
                ax.set_ylabel( ylabel, fontsize=fs )
        else:
            plt.axis(lims)
            if kwargs.get( 'show_ax_label', True ):
                plt.xlabel( xlabel, fontsize=fs )
                plt.ylabel( ylabel, fontsize=fs )


    def hr_box(self, t, dt, l, dl, logt=True, logl=True, ls='solid'):
        """
        <CModelGUI.hr_box(t, dt, l, dl, logt=True, logl=True, linestyle='solid')>

        Draws an error box on the HRD.

        :param t:     Effective temperature (or log Teff if `logt=True`)
        :type t:      float

        :param d :    Error in Teff (or in log Teff if `logt=True`)
        :type dt:     float

        :param l:     Luminosity (or log L if `logl=True`)
        :type l:      float

        :param d :    Error in L (or in log L if `logl=True`)
        :type dl:     float

        :kparam logt : If True, temperature in log (default)
        :ktype logt:   bool, optional

        :kparam logl : If True, luminosity in log (default)
        :ktype logl:   bool, optional

        :kparam ls:    Linestyle to use in plot. Options: `'solid'`, `'dashed'`, `'dashdot'`,
            `'dotted'`. Default: `'solid'`.
        :ktype ls:     string
        """

        if not logt:
            tmin = np.log10(t - dt)
            tmax = np.log10(t + dt)
        else:
            tmin = t - dt
            tmax = t + dt

        if not logl:
            lmin = np.log10(l - dl)
            lmax = np.log10(l + dl)
        else:
            lmin = l - dl
            lmax = l + dl

        ax = [tmax, tmax, tmin, tmin]
        ay = [lmax, lmin, lmin, lmax]

        plt.fill(ax, ay, fill=False, linestyle=ls)

    def get_teff_omega_model( self ):
        r        = 10**self.log_r * self.ctes.rsun
        req      = 10**self.log_req * self.ctes.rsun
        rpol     = 10**self.log_rpol * self.ctes.rsun
        omegak   = np.sqrt( self.ctes.ggrav * self.mstar * self.ctes.msun / req**3 )
        omega    = self.omega_s / omegak
        mred     = self.mstar*self.ctes.msun * (1 - self.omega_s**2/(8*np.pi/3 * self.ctes.ggrav * self.mstar*self.ctes.msun / r**3))

        Fpol     = np.exp( 2 * omega**2 * (rpol/req)**3 / 3 )
        Feq      = (1 - omega**2)**(-2/3)

        lsgm     = 10**self.log_l*self.ctes.lsun / (4 * np.pi * sigma * self.ctes.ggrav * mred)
        teff_eq  = (lsgm * Feq  * self.geffeq)**0.25
        teff_pol = (lsgm * Fpol * self.geffpol)**0.25
        return np.log10( teff_eq ), np.log10( teff_pol )




    def plot_hr(self, ls='-', mods=None, ax=None, nologx=False, nology=False, update=False, **kwargs):
        """
        <CModelGUI.plot_hr(ls=None, mods=None, ax=None, log=True, **kwargs)

        Plot Hertsprung-Russel diagramm.

        :kparam ls:      Matplotlib supported line style. Default: `'-'`.
        :ktype ls:       string

        :kparam mods:    Draw logteff, logl of specific times steps. Default: None
        :ktype mods:     list

        :kparam ax:      Instance of matplotlib.axes.axes. Default: Current axis
        :ktype ax:       matplotlib.axes.axes

        :kparam nologx: True if the linear value instead of the log value of the x
            quantity must be plotted. Default: False.
        :ktype nologx: boolean

        :kparam nology: True if the linear value instead of the log value of the y
            quantity must be plotted. Default: False.
        :ktype nology: boolean

        :kparam update: True if the user wants to update an existing HR diagram. Default: False.
        :ktype update: boolean

        :kwarg label: Label ateached to the evolutionnary track. You need to call
            `plt.legend()` to display it. Default: None
        :kwtype label: str

        :kwarg linewidth: Set width of line. Default: Use the one currently set.
        :kwtype linewidth: float

        :kwarg imin: Index of first time step to be plotted. Defaut: 0
        :kwtype imin: int

        :kwarg imax: Index of last time step to be plotted. Defaut: number of
            total time step in the model
        :kwtype imax: int

        :kwarg color: Color of the line. Default: Next color in the cycle.
        :kwtype color: str

        :kwarg highlight_mods: Time steps highlighted by a dot
        :kwtype highlight_mods: array of int

        :kwarg show_ax_label: Wether we display the ax labels or not. Default: True.
        :kwtype show_ax_label: boolean

        :kwarg markersize: Marker size if marker are used
        :kwtype markersize: int

        :kwarg teff_theta: Shade area for Teff varying with latitude.
        :kwtype teff_theta: bool
        """

        if not self.hr_exists:
            print('Warning: HR file does not exist.')
            return
        if self.hr_empty:
            print('Warning: HR file exists but is empty.')
            return
        if (not self.log_l.any()) or (not self.log_teff.any()):
            print('Warning: HR file empty.')
            return

        func = lambda x, nolog: 10**x if nolog else x
        xx = func( self.log_teff, nologx )
        yy = func( self.log_l,    nology )

        imin      = kwargs.get( "imin", 0 )
        ima       = len( self.log_teff )
        imax      = min( ima, kwargs.get( "imax", ima ) )
        ls        = ls if ls is not None else '-'
        ax        = plt.gca() if ax is None else ax

        kwargs['label']      = kwargs.get( "label", None )
        kwargs['linewidth']  = kwargs.get( "linewidth", plt.rcParams['lines.linewidth']  )
        kwargs['marker']     = kwargs.get( "marker", '' )
        kwargs['color']      = kwargs.get( "color", ax._get_lines.get_next_color())
        kwargs['markersize'] = kwargs.get( "markersize", plt.rcParams['lines.markersize'] )
        kwargs['teff_theta'] = kwargs.get( "teff_theta", False )

        self.hr_axis( ax, nologx=nologx, nology=nology, **kwargs)
        rmkeys = ['imin', 'imax', 'show_ax_label']
        for k in set(kwargs.keys()).intersection(rmkeys): del kwargs[k]

        # span a region for teff depending on colatitude
        if kwargs['teff_theta'] and hasattr(self, "omega_s"):
            teffeq, teffpol = self.get_teff_omega_model( )
            teffpol = func( teffpol, nologx )
            teffeq = func( teffeq, nologx )
            ax.fill_betweenx( yy[imin:imax], teffpol[imin:imax], teffeq[imin:imax], facecolor=kwargs['color'],
                alpha=0.5 )

        # Do we only update an existing HRD ?
        if update and hasattr(self, "hrline"):
            self.hrline.set_data( xx[imin:imax], yy[imin:imax] )
        # or do we create it ?
        else:
            self.hrline, = ax.plot( xx[imin:imax], yy[imin:imax], ls=ls, label=kwargs['label'],
                linewidth=kwargs['linewidth'], color=kwargs['color'], marker=kwargs['marker'] )

        if mods is not None:
            ax.plot(xx[mods], yy[mods], 'o', **kwargs)

        highlight_mods = kwargs.get("highlight_mods", None )
        rmkeys = ['highlight_mods', 'label']
        for k in set(kwargs.keys()).intersection(rmkeys): del kwargs[k]

        if highlight_mods is not None:
            for hm in highlight_mods:
                ax.plot([xx[hm]], [yy[hm]], 'o', color=kwargs['color'],
                    markersize=kwargs['markersize'] )

        return ax


    def plot_kiel(self, ls=None, mods=None, ax=None, nologx=False, nology=False, update=False,
            **kwargs):
        """
        <CModelGUI.plot_kiel(ls=None, mods=None, ax=None, **kwargs)>

        Plot Kiel diagramm (Teff-logg).

        :kparam ls: Matplotlib supported line style. Default: None.
        :ktype ls: str

        :kparam mods: Draw logteff, logg of specific times steps. Default: None
        :ktype mods: list

        :kparam ax: Instance of matplotlib.axes.axes. Default: None
        :ktype ax: matplotlib.axes.axes

        :kparam nologx: True if the linear value instead of the log value of the x
            quantity must be plotted. Default: False.
        :ktype nologx: boolean

        :kparam nology: True if the linear value instead of the log value of the y
            quantity must be plotted. Default: False.
        :ktype nology: boolean

        :kparam update: True if the user wants to update an existing Kiel diagram. Default: False.
        :ktype update: boolean

        :kwarg label: Give label to curve. Default: ''
        :kwtype label:  str

        :kwarg linewidth: Line width of tracks. Default: 2.
        :kwtype linewidth: float

        :kwarg imin: index of first model in track. Default: 0.
        :kwtype imin: int

        :kwarg iax: index of last model in track. Default: -1.
        :kwtype iax: int

        :kwarg color: Color of the track. Default: next color in color cycle.
        :kwtype color: str

        :kwarg teff_theta: Shade area for Teff varying with latitude.
        :kwtype teff_theta: bool
        """
        if not self.hr_exists:
            print('Warning: HR file does not exist.')
            return
        if self.hr_empty:
            print('Warning: HR file exists but is empty.')
            return
        if (not self.log_l.any()) or (not self.log_teff.any()):
            print('Warning: HR file empty.')
            return

        func = lambda x, nolog: 10**x if nolog else x
        xx = func( self.log_teff, nologx )
        yy = func( self.log_g,    nology )

        imin      = kwargs.get( "imin", 0 )
        ima       = len( self.log_teff )
        imax      = min( ima, kwargs.get( "imax", ima ) )
        ls        = ls if ls is not None else '-'
        ax        = plt.gca() if ax is None else ax

        kwargs['label']      = kwargs.get( "label", None )
        kwargs['linewidth']  = kwargs.get( "linewidth", plt.rcParams['lines.linewidth']  )
        kwargs['marker']     = kwargs.get( "marker", '' )
        kwargs['color']      = kwargs.get( "color", ax._get_lines.get_next_color())
        kwargs['markersize'] = kwargs.get( "markersize", plt.rcParams['lines.markersize'] )
        kwargs['teff_theta'] = kwargs.get( "teff_theta", False )


        self.kiel_axis( ax, nologx=nologx, nology=nology, **kwargs)
        rmkeys = ['imin', 'imax', 'show_ax_label']
        for k in set(kwargs.keys()).intersection(rmkeys): del kwargs[k]

        if kwargs['teff_theta'] and hasattr(self, "omega_s"):
            teffeq, teffpol = self.get_teff_omega_model( )
            teffpol = func( teffpol, nologx )
            teffeq = func( teffeq, nologx )
            ax.fill_betweenx( yy[imin:imax], teffpol[imin:imax], teffeq[imin:imax], facecolor=kwargs['color'],
                alpha=0.5 )

        # Do we only update an existing Kiel diagram ?
        if update:
            if not hasattr(self, "kielline"):
                raise CESAMError( "You called plot_hr with the argument update=True, but your instance of CModelGUI"
                    "does not have an kielline attribute yet. You must call plot_hr with argument update=False first.")

            self.kielline.set_data( xx[imin:imax], yy[imin:imax] )
        # or do we create it ?
        else:
            self.kielline, = ax.plot( xx[imin:imax], yy[imin:imax], ls=ls, label=kwargs['label'],
                linewidth=kwargs['linewidth'], color=kwargs['color'], marker=kwargs['marker'] )

        if mods is not None:
            ax.plot(xx[mods], yy[mods], 'o', **kwargs)

        highlight_mods = kwargs.get("highlight_mods", None )
        rmkeys = ['highlight_mods', 'label']
        for k in set(kwargs.keys()).intersection(rmkeys): del kwargs[k]

        if highlight_mods is not None:
            for hm in highlight_mods:
                ax.plot([xx[hm]], [yy[hm]], 'o', color=kwargs['color'],
                    markersize=kwargs['markersize'] )

        return ax


    def plot_cmd( self, ls='-' ):
        """
        <CModelGUI.plot_cmd( ls=None )>

        Plots a color-magnitude diagram. If B-V and $M_V$ have not been previously
        computed for model, computes then by calling `self.read_hr(obs=True)`

        :kparam ls: Line style. Default: '-'.
        :ktype ls: str
        """

        try:
            axs = (min(self.bv),max(self.bv),
                    max(self.mv),min(self.mv))
        except AttributeError:
            self.read_hr(obs=True)
            axs = (min(self.bv),max(self.bv),
                    max(self.mv),min(self.mv))

        plt.axis(axs)
        plt.xlabel(r'$B-V$')
        plt.ylabel(r'$M_V$')

        plt.plot(self.bv, self.mv, ls=ls)


    def plot_cz( self, r_or_m='m', hatch=False, ax=None, colors=['k'], **kwargs ):
        """
        <CModelGUI.plot_cz( r_or_m='m', hatch=False )>
        Plots a Kippenhahn Diagram of the model's convection zones.

        :kparam r_or_m: If 'r', eulerian plot If 'm', Lagragean plot. Default: 'm'.
        :ktype obs: str

        :kparam hatch: Whether or not we draw hatches on different zones.
            Default: False.
        :ktype hatch: boolean

        :kparam ax:      Instance of matplotlib.axes.axes. Default: Current axis
        :ktype ax:       matplotlib.axes.axes

        :kparam colors:  List of colors to be cycled through when plotting zones limits. Default: ['k']
        :ktype colors:   list or numpy.ndarray

        :kwarg imin: index of first model in track. Default: 0.
        :kwtype imin: int

        :kwarg iax: index of last model in track. Default: -1.
        :kwtype iax: int

        :raises ValueError: `r_or_m` is not 'm' or 'r'.
        """

        nzones_tot = len(self.iconv_start)
        npts       = [self.iconv_end[i] - self.iconv_start[i] + 1 for i in range(nzones_tot)]

        if not isinstance(colors, (np.ndarray, list) ):
            raise TypeError( "colors should be a list, or a numpy.ndarray.")
        nc         = len(colors)

        if (r_or_m == 'm'):
            yy          = self.mstar
            yconv_start = self.mconv_start
            yconv_end   = self.mconv_end
            ymax        = max(self.mstar)
            ylabel      = r"$m/M_\odot$"
        elif r_or_m == 'r':
            yy          = 10**self.log_r
            yconv_start = self.rconv_start
            yconv_end   = self.rconv_end
            ymax        = max(10**self.log_r)
            ylabel      = r"$r/R_\odot$"
        else:
            raise ValueError( f"r_or_m must be 'm' or 'r'. Received {r_or_m}." )

        imin      = kwargs.get( "imin", 0 )
        ima       = len( self.age ) -1
        imax      = min( ima, kwargs.get( "imax", ima ) )

        y_pts = [np.zeros(2*npts[i]+1) for i in range(nzones_tot)]
        x_pts = [np.zeros(2*npts[i]+1) for i in range(nzones_tot)]


        for i in range(nzones_tot):
            age = self.age[self.iconv_start[i]:self.iconv_end[i]+1]

            y_pts[i][:npts[i]]   = yconv_start[i]
            y_pts[i][npts[i]:-1] = yconv_end[i][::-1]
            y_pts[i][-1]         = yconv_start[i][0]

            x_pts[i][:npts[i]]   = age
            x_pts[i][npts[i]:-1] = age[::-1]
            x_pts[i][-1]         = age[0]

        y_border = np.append(0.0, yy)
        y_border = np.append(y_border, 0.0)
        y_border = np.append(y_border, 0.0)

        x_border = np.append(self.age[0], self.age)
        x_border = np.append(x_border, self.age[-1])
        x_border = np.append(x_border, self.age[0])

        border = np.array([[x, y] for x, y in zip(x_border, y_border)] )

        ax = plt.subplot(111) if ax is None else ax

        if hatch:
            hatch1, hatch2, fc1, fc2 = None, '//', 'none', 'none'
            alpha = 0.99
        else:
            hatch1, hatch2, fc1, fc2  = None, None, '0.75', '0.35'
            alpha = 1.0

        pc = plt.Polygon( border, hatch=hatch1, fc=fc1, alpha=alpha )
        ax.add_patch(pc)

        for i in range(nzones_tot):
            pts = np.array([[x, y] for x, y in zip(x_pts[i], y_pts[i])] )
            pc = plt.Polygon( pts, hatch=hatch2, fc=fc2, alpha=alpha )
            ax.add_patch(pc)
            ax.plot(x_pts[i], y_pts[i], color=colors[i%nc])

        ax.plot(self.age, yy, 'k')

        ax.set_xlim( self.age[imin], self.age[imax] )
        ax.set_ylim( 0.0, ymax*1.05 )
        ax.set_xlabel( r"Age [Myrs]" )
        ax.set_ylabel( ylabel )

    def plot_burn( self, r_or_m='m', hatch=False, color=True, contours=False, with_cz=True, fig=None, ax=None ):
        """
        <CModelGUI.plot_burn( r_or_m='m', hatch=False )>
        Plots a Kippenhahn Diagram of the model's convection zones.

        :kparam r_or_m: If 'r', eulerian plot If 'm', Lagragean plot. Default: 'm'.
        :ktype obs: str

        :kparam hatch: Whether or not we draw hatches on different zones.
            Default: False.
        :ktype hatch: boolean

        :kparam color: Whether or not we use a color scale (reds).
            Default: True.
        :ktype hatch: boolean

        :kparam contours: Whether or not we draw contours.
            Default: False.
        :ktype hatch: boolean

        :kparam with_cz: Whether or not we draw convection zones (hatched).
            Default: True.
        :ktype hatch: boolean

        :kparam fig: Instance of `matplotlib.figure.Figure`. Default: Current figure.
        :ktype fig: matplotlib.figure.Figure

        :kparam ax: Instance of `matplotlib.axes.axes`. Default: Current axis
        :ktype ax: matplotlib.axes.axes

        :raises ValueError: `r_or_m` is not 'm' or 'r'.
        """
        fig = plt.gcf( ) if fig is None else fig

        if hasattr(self, "iburn_start"):
            if color:
                red4 = (1.0, 0.0, 0)
                red3 = (1.0, 0.5, 0.5)
                red2 = (1.0, 0.67, 0.67)
                red1 = (1.0, 0.83, 0.83)
                colors = []
                for i in range(self.n_contours_burn):
                    k = (self.n_contours_burn - i - 1) / self.n_contours_burn
                    colors.append( (1.0, k, k) )

                color_back = (0.8, 0.9, 0.9)
            else:
                colors = ['0.5', '0.4', '0.3', '0.2']
                color_back = '0.8'

            for contour in range(4):

                nzones_tot = len(self.iburn_start[contour])
                npts       = [self.iburn_end[contour][i] - self.iburn_start[contour][i] + \
                              1 for i in range(nzones_tot)]

                if (r_or_m == 'm'):
                    yy          = self.mstar
                    yburn_start = self.mburn_start[contour]
                    yburn_end   = self.mburn_end[contour]
                    ymax        = max(self.mstar)
                    ylabel      = r"$m/M_\odot$"
                elif r_or_m == 'r':
                    yy          = 10**self.log_r
                    yburn_start = self.rburn_start[contour]
                    yburn_end   = self.rburn_end[contour]
                    ymax        = max(10**self.log_r)
                    ylabel      = r"$r/R_\odot$"
                else:
                    raise ValueError( f"r_or_m must be 'm' or 'r'. Received {r_or_m}." )

                y_pts = [np.zeros(2*npts[i]+1) for i in range(nzones_tot)]
                x_pts = [np.zeros(2*npts[i]+1) for i in range(nzones_tot)]

                for i in range(nzones_tot):
                    age = self.age[self.iburn_start[contour][i]:self.iburn_end[contour][i]+1]

                    y_pts[i][:npts[i]]   = yburn_start[i]
                    y_pts[i][npts[i]:-1] = yburn_end[i][::-1]
                    y_pts[i][-1]         = yburn_start[i][0]

                    x_pts[i][:npts[i]]   = age
                    x_pts[i][npts[i]:-1] = age[::-1]
                    x_pts[i][-1]         = age[0]

                y_border = np.append(0.0, yy)
                y_border = np.append(y_border, 0.0)
                y_border = np.append(y_border, 0.0)

                x_border = np.append(self.age[0], self.age)
                x_border = np.append(x_border, self.age[-1])
                x_border = np.append(x_border, self.age[0])

                border = np.array([[x, y] for x, y in zip(x_border, y_border)] )

                ax = plt.subplot(111) if ax is None else ax

                # ax = plt.subplot(111)

                if hatch:
                    hatch1, hatch2, fc1, fc2 = None, '//', None, None
                else:
                    hatch1, hatch2, fc1, fc2  = None, None, color_back, colors[contour]

                if contour == 0:
                    pc = plt.Polygon(border, hatch=hatch1, fc=fc1)
                    ax.add_patch(pc)

                for i in range(nzones_tot):
                    pts = np.array([[x, y] for x, y in zip(x_pts[i], y_pts[i])] )
                    pc = plt.Polygon(pts, hatch=hatch2, fc=fc2)
                    ax.add_patch(pc)
                    if contours: ax.plot(x_pts[i], y_pts[i], 'k')

                ax.plot(self.age, yy, 'k')

            ax.set_xlim( self.age[0], self.age[-1] )
            ax.set_ylim( 0.0, ymax*1.05 )
            ax.set_xlabel( r"Age [Myrs]" )
            ax.set_ylabel( ylabel )

            if with_cz:
                self.plot_cz(r_or_m=r_or_m, hatch=True)

            cmap = mplt.colors.ListedColormap([color_back]+colors)
            try:
                bounds = np.log10(self.burn_contours[:self.n_contours_burn])
            except AttributeError:
                bounds = [0.0, 1.0, 4.0, 6.0]
            norm = mplt.colors.BoundaryNorm(bounds, cmap.N, extend='both')
            fig.colorbar(
                    mplt.cm.ScalarMappable(cmap=cmap, norm=norm),
                    location='right',
                    ax=ax,
                    extend='both',
                    ticks=bounds,
                    spacing='proportional',
                    label=r'$\log \varepsilon\,({\rm erg} \, {\rm g}^{-1} \,{\rm s}^{-1})$')
        else:
            print('HR file without burn zones.')

    def plot_central(self, limits=False, phases=False):
        """
        Plots Tc as a function of rhoc.

        :param limits: Sets whether to plot limits of degenerate zones, radiation pressure, etc.
        :type limits: bool

        :param phases: Sets whether to show evolutionary phase.
        :type limits: bool
        """


        if phases:
            names = ['PMS', 'MS', 'RGB', 'Core He burn', 'AGB']
            try:
                inds = [self.pms, self.ms, self.rgb, self.cohe, self.agb]
                legend = []

                for ind, name in zip (inds, names):
                    if len(ind) > 0:
                        plt.plot(self.rhoc[ind], self.Tc[ind])
                        legend.append(name)

                plt.legend(legend)
            except AttributeError:
                print('No evolutionary phases available, will plot without.')
                plt.plot(self.rhoc, self.Tc)
        else:
            plt.plot(self.rhoc, self.Tc)

        plt.loglog()

        xmin, xmax = plt.xlim()
        ymin, ymax = plt.ylim()

        if limits:

            ne = 0.5*self.rhoc/amu
            pF = (3.0*ne/8/np.pi)**(1.0/3)*hpl

            # find degenerate zones
            Npsi = 4  # find where EF = 4*kB*T
            Tnorel = 0.5/Npsi/mel/kbol*pF**2
            Trel   = clight/Npsi/kbol*pF
            T = np.where(Trel < Tnorel, Trel, Tnorel)
            plt.plot(self.rhoc, T, 'k--')

            # find Prad = Pgas (using a perfect gas with mu = 1)
            mu = 1.0  # typical, not to be taken too seriously
            T_rad = (3*rgas/mu/aradia*self.rhoc)**(1.0/3.0)
            plt.plot(self.rhoc, T_rad, 'k--')


            plt.xlim((xmin, xmax))
            plt.ylim((ymin, ymax))


        plt.xlabel(r'$\rho_{\rm c}\,({\rm g} \,{\rm cm}^{-3})$')
        plt.ylabel(r'$T_{\rm c}\,({\rm K})$')


    def plot_osc(self, ix, iy, ls='-'):
        """
        Plot a column of the .osc file as a function of another one.

        :param ix: Index of column represented in abscissa.
        :type ix: int

        :param iy: Index of column represented in ordinates.
        :type iy: int

        :kparam ls: Line style. Default: '-'.
        :ktype ls: str
        """
        xxlabel = self.nom_vars[ix]
        yylabel = self.nom_vars[iy]

        plt.plot(self.var[ix], self.var[iy], ls=ls)

        plt.xlabel(xxlabel)
        plt.ylabel(yylabel)


    def plot_all_eps(self, plot_log=False, age_lims=None, npnt=1000,
        ncontours=500, zmax=None, zmin=None, cz=False):

        if not self.all_osc:
            print(f"{RED}Error:{NO_COLOR} files .osc not available...")
            return

        if age_lims:
            age_i = age_lims[0]
            age_f = age_lims[1]
        else:
            age_i = min(self.age)
            age_f = max(self.age)


        mstar = max(self.mstar)
        x = np.array([])
        y = np.linspace(0.0, mstar, npnt)
        eps, eps_g = [], []

        for i in range(1, self.nmod):
            if age_i <= self.age[i] <= age_f:
                ma = self.var[i][1][::-1]*self.mstar[i]
                ma1 = self.var[i-1][1][::-1]*self.mstar[i-1]
                dt = (self.age[i] - self.age[i-1])*1e6*365.2425*24*3600
                f = interpolate.interp1d(ma1, self.var[i-1][2][::-1])
                t1 = f(y)
                f = interpolate.interp1d(ma, self.var[i][2][::-1])
                t2 = f(y)
                dtdt = (t2 - t1)/dt
                f = interpolate.interp1d(ma1, self.var[i-1][3][::-1])
                p1 = f(y)
                f = interpolate.interp1d(ma, self.var[i][3][::-1])
                p2 = f(y)
                dpdt = (p2 - p1)/dt

                f = interpolate.interp1d(ma, self.var[i][12][::-1])
                cp = f(y)
                f = interpolate.interp1d(ma, self.var[i][4][::-1])
                rho = f(y)
                f = interpolate.interp1d(ma, self.var[i][11][::-1])
                delta = f(y)

                eg = cp*dtdt - delta/rho*dpdt
                eps_g.append(eg)

                f = interpolate.interp1d(ma, self.var[i][8][::-1])
                eps.append(f(y))
                x = np.append(x, self.age[i])

        eps = np.array(eps).T
        eps_g = -np.array(eps_g).T
        eps = eps - eps_g
        if plot_log:
            if zmax is not None:
                eps = np.where(eps > zmax, np.log10(zmax), np.log10(eps))
                eps_g = np.where(eps_g > zmax, np.log10(zmax), np.log10(eps_g))
                if zmin is not None:
                    eps = np.where(eps < zmin, zmin, eps)
                    eps_g = np.where(eps_g < zmin, zmin, eps_g)
            elif zmin is not None:
                eps = np.where(eps < zmin, np.log10(zmin), np.log10(eps))
                eps_g = np.where(eps_g < zmin, np.log10(zmin), np.log10(eps_g))
            else:
                eps = np.log10(eps)
                eps_g = np.log10(eps_g)
        else:
            if zmax is not None:
                eps = np.where(eps > zmax, zmax, eps)
                eps_g = np.where(eps_g > zmax, zmax, eps_g)
            if zmin is not None:
                eps = np.where(eps < zmin, zmin, eps)
                eps_g = np.where(eps_g < zmin, zmin, eps_g)

        figs = []
        ax = []

        zz = [eps, eps_g]
        legs = [r'$\varepsilon$', r'$\varepsilon_g$']
        for i in range(2):
            figs.append(plt.figure())
            ax.append(figs[-1].add_subplot(111))

            plt.contour(x, y, zz[i], ncontours)
            plt.contourf(x, y, zz[i], ncontours)
            plt.xlabel(r'Age [Myrs]')
            plt.ylabel(r'$m/M_{\odot}$')

            cb = plt.colorbar()
            if plot_log:
                cb.set_label('$' + r'\log ' + legs[i][1:])
            else:
                cb.set_label(legs[i])

            if cz: self.plot_cz('m', ls='w', lw=1.0, hatch='/', fill=False)


    def plot_all_osc(self, ivar, plot_log=False, age_lims=None, npnt=1000,
                     ncontours=500, zmax=None, zmin=None, cz=False):

        if not self.all_osc:
            print(f"{RED}Error:{NO_COLOR} files .osc not available...")
            return

        age_i = age_lims[0] if age_lims else min( self.age )
        age_f = age_lims[1] if age_lims else max( self.age )


        mstar = max(self.mstar)
        x = np.array([])
        y = np.linspace(0.0, mstar, npnt)
        z = []

        for i in range(self.nmod):
            if age_i <= self.age[i] <= age_f:
                ma = self.var[i][1][::-1]*self.mstar[i]
                z1 = self.var[i][ivar][::-1]
                f  = interpolate.interp1d(ma, z1)
                z.append(f(y))
                x  = np.append(x, self.age[i])

        z = np.array(z).T
        if plot_log:
            if zmax is not None:
                z = np.where(z > zmax, np.log10(zmax), np.log10(z))
                if zmin is not None:
                    z = np.where(z < zmin, zmin, z)
            elif zmin is not None:
                z = np.where(z < zmin, np.log10(zmin), np.log10(z))
            else:
                z = np.log10(z)
        else:
            if zmax is not None:
                z = np.where(z > zmax, zmax, z)
            if zmin is not None:
                z = np.where(z < zmin, zmin, z)

        plt.contour(x, y, z, ncontours)
        plt.contourf(x, y, z, ncontours)
        plt.xlabel(r'${\rm age}\, {\rm (Myrs)}$')
        plt.ylabel(r'$m/M_{\odot}$')

        cb = plt.colorbar()
        if plot_log:
            cb.set_label('$' + r'\log ' + self.nom_vars[ivar][1:])
        else:
            cb.set_label(self.nom_vars[ivar])

        if cz: self.plot_cz('m', ls='w', lw=1.0, hatch='/', fill=False)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def plot_all_osc_r(self, ivar, plot_log=False, age_lims=None, npnt=2000,
                       ncontours=500, zmax=None, zmin=None, cz=False, rmax=None,
                       fill=True,label=False):

        if not self.all_osc:
            print(f"{RED}Error:{NO_COLOR} files .osc not available...")
            return

        if age_lims:
            age_i = age_lims[0]
            age_f = age_lims[1]
        else:
            age_i = min(self.age)
            age_f = max(self.age)


        rstar = max(self.ctes.rsun*10**self.log_r)
        if rmax is not None:
            rstar = min(rmax, rstar)

        x = np.array([])
        y = np.linspace(0.0, rstar, npnt)
        z = []

        for i in range(self.nmod):
            if age_i <= self.age[i] <= age_f:
                ra = self.var[i][0][::-1]
                z1 = self.var[i][ivar][::-1]
                f = interpolate.interp1d(ra, z1)
                mask = y < max(ra)
                zz1 = np.zeros(npnt)
                zz1[mask] = f(y[mask])
                zz1[~mask] = np.nan
                z.append(zz1)
                x = np.append(x, self.age[i])

        z = np.array(z).T
        if label: levels = np.linspace(z[~np.isnan(z)].min(), z[~np.isnan(z)].max(), 10)
        if plot_log:
            if zmax is not None:
                z = np.where(z > zmax, np.log10(zmax), np.log10(z))
                if zmin is not None:
                    z = np.where(z < zmin, zmin, z)
            elif zmin is not None:
                z = np.where(z < zmin, np.log10(zmin), np.log10(z))
            else:
                z = np.log10(z)
        else:
            if zmax is not None:
                z = np.where(z > zmax, zmax, z)
            if zmin is not None:
                z = np.where(z < zmin, zmin, z)

        if fill:
            plt.contour(x, y/self.ctes.rsun, z, ncontours)
            plt.contourf(x, y/self.ctes.rsun, z, ncontours)
        else:
            plt.contour(x, y/self.ctes.rsun, z, ncontours, colors='k')
            ctb = plt.contour(x, y/self.ctes.rsun, z, levels, colors='r')
            if label: plt.clabel(ctb, levels)

        plt.xlabel(r'${\rm age}\, {\rm (Myrs)}$')
        plt.ylabel(r'$r/R_{\odot}$')

        if fill:
            cb = plt.colorbar()
            if plot_log:
                cb.set_label('$' + r'\log ' + self.nom_vars[ivar][1:])
            else:
                cb.set_label(self.nom_vars[ivar])

        if cz:
            self.plot_cz('r', ls='w' if fill else 'b', lw=1.0, hatch='/', fill=False)

    def plot_ani(self, var=None, xvar=None, filetype="conv.dat", fig=None, ax=None, nologx=False, nology=False, save=False, save_name=None, **kwargs):
        """
        <CModelGUI.plot_ani(var=None, xvar=None, filetype="conv.dat", fig=None, ax=None, logx=True, logy=False, save=False, save_name=None, **kwargs)

        Plot evolution of a given variable with HR diagram on the side.

        :kparam var:        Name of the variable of the y-axis. Default: None
        :ktype var:         string

        :kparam xvar:       Name of the variable of the x-axis. Default: None
        :ktype mods:        string

        :kparam filetype:   Filetype that stores the variables. "conv.dat" or "osc"
        :ktype mods:        string

        :kparam fig:        Instance of matplotlib.figure Default: Current figure
        :ktype fig:         matplotlib.figure

        :kparam ax:         Instance of matplotlib.axes.axes. Default: Current axis
        :ktype ax:          matplotlib.axes.axes

        :kparam nologx:       True if linear scale should be used for the x axis.  Default: False.
        :ktype nologx:        boolean

        :kparam nology:       True if linear scale should be used for the y axis. Default: False.
        :ktype nology:        boolean

        :kparam save:       True if the generate animation should be saved. Default: False.
        :ktype update:      boolean

        :kparam save_name:  Name of .gif file. If not given. the name of model is used. Default: None
        :ktype update:      string
        """

        if not self.hr_exists:
            print('Warning: HR file does not exist.')
            return
        if self.hr_empty:
            print('Warning: HR file exists but is empty.')
            return
        if (not self.log_l.any()) or (not self.log_teff.any()):
            print('Warning: HR file empty.')
            return

        ivar, varstr = self.get_var_info(var, filetype) if var is not None else (19, r'$\alpha_mathrm{MLT}$')

        ixvar, xvarstr = self.get_var_info(xvar, filetype) if xvar is not None else (5, r'$\rho$')

        if fig is None:
            fig, ax = plt.subplots(1, 2, figsize=(9, 4))
        elif ax is None:
            ax = fig.add_subplot(1, 2)

        if filetype == "conv.dat":

            n_ts = sum(1 for f in os.listdir(os.getcwd()) if f.endswith(self.name+"-conv.dat"))
            print(f"{n_ts} {self.name}-conv.dat files found.")

            func = lambda x, nolog: x if nolog else np.log10(x)

            xx = []
            yy = []
            xmin = np.nan
            xmax = np.nan
            ymin = np.nan
            ymax = np.nan
            for i in tqdm( range( n_ts ), desc=f'Reading *{self.name}-conv.dat files' ):
                xx.append(func(np.loadtxt(f"{i:05d}-{self.name}-conv.dat", skiprows=2, usecols=(ixvar,)), nologx))
                xmin = np.nanmin([np.min(xx[-1]), xmin])
                xmax = np.nanmax([np.max(xx[-1]), xmax])
                yy.append(func(np.loadtxt(f"{i:05d}-{self.name}-conv.dat", skiprows=2, usecols=(ivar,)), nology))
                ymin = np.nanmin([np.min(yy[-1]), ymin])
                ymax = np.nanmax([np.max(yy[-1]), ymax])

        elif filetype == "osc":
            self.read_osc()

            if len(self.osc) < 2:
                print("1 osc file is not enough for animation")
                return

            xx = self.var[ixvar, :]
            yy = self.var[ivar, :]

        else:
            print("Neither .osc files nor *conv.dat found")
            return

        varstr = varstr if nology else r"$\log_{10}$"+varstr
        xvarstr = xvarstr if nologx else r"$\log_{10}$"+xvarstr

        ani = VarAnimation(self, n_ts=n_ts, xx=xx, yy=yy, var=varstr, xvar=xvarstr, xmin=xmin*0.99, xmax=xmax*1.01,
                           ymin=ymin*0.99, ymax=ymax*1.01, fig=fig, ax=ax)

        if save:
            ani.save(save_name) if save_name else ani.save(self.name)

        plt.show()

    def get_var_info(self, var, filetype=None):
        # to be updated later

        osc_list = {'r': (r'$r$', 'radius'),
                    'm':  (r'$m/M_{\star}$', 'mass'),
                    'T': (r'$T$', 'temperature'),
                    'Ptot': (r'$P_{\rm tot}$', 'pressure'),
                    'rho': (r'$\rho$', 'density'),
                    'nabla': (r'$\nabla$', 'temperature gradient'),
                    'L': (r'$L$', 'luminosity'),
                    'kappa': (r'$\kappa$', 'opacity'),
                    'vareps': (r'$\varepsilon$', 'varespsilon'),
                    'gamma1': (r'$\Gamma_1$', 'adiabatic index'),
                    'nabla_ad': (r'$\nabla_{\rm ad}$', 'adiabatic temperature gradient'),
                    'delta': (r'$\delta$', 'isobaric expansion coefficient'),
                    'cp': (r'$c_p$', 'specific heat capacity'),
                    'mue': (r'$1/\mu_e$', 'mue'),
                    'A': (r'$\mathcal{A}$', 'calculation of the Brunt-Väisäla frequency'),
                    'omega': (r'$\Omega$', 'angular velocity'),
                    'kappaT': (r'$\kappa_T$', 'kappaT'),
                    'kapparho': (r'$\kappa_{\rho}$', 'kapparho'),
                    'varepsT': (r'$\varepsilon_T$', 'varespT'),
                    'varepsrho': (r'$\varepsilon_{\rho}$', 'varesprho'),
                    'Ptot2gas': (r'$P_{\rm tot}/P_{\rm gas}$', 'pressure'),
                    'nablarad': (r'$\nabla_{\rm rad}$', 'radiative temperature gradient')}
        conv_dat_list = {'r': (r'$r$', 'radius'),
                         'm23':  (r'$m^{2/3}$', 'mass coordinate'),
                         'm/Msun':  (r'$m/M_\odot$', 'mass'),
                         'T': (r'$T$', 'temperature'),
                         'Ptot': (r'$P_{\rm tot}$', 'pressure'),
                         'rho': (r'$\rho$', 'density'),
                         'cp': (r'$c_p$', 'specific heat capacity'),
                         'delta': (r'$\delta$', 'isobaric expansion coefficient'),
                         'Fconv': (r'$F_\mathrm{conv}$','convective flux'),
                         'g': (r'$g$', 'gravity'),
                         'kap': (r'$\kappa$', 'kappa'),
                         'Hp': (r'$H_p$', 'pressure scale height'),
                         's': (r'$s$', 'entropy'),
                         'nabconv': (r'$\nabla_\mathrm{conv}$', 'actual temperature gradient'),
                         'nabrad': (r'$\nabla_\mathrm{rad}$', 'radiative temperature gradient'),
                         'nabad': (r'$\nabla_\mathrm{ad}$', 'adiabatic temperature gradient'),
                         'vconv': (r'$v_\mathrm{conv}$', 'convective velocity'),
                         'mmw': (r'$\mu$', 'mean molecular weight'),
                         'kappa': (r'$\kappa$', 'optical depth'),
                         'amlt': (r'$\alpha_\mathrm{MLT}$', 'mixing length parameter'),
                         'deltab': (r'$\delta_\mathbf{B}$', 'magnetic inhibition parameter'),
                         'br': (r'$B_r$', 'radial component of the magnetic field')}

        if filetype == 'osc':
            try:
                ivar = list(osc_list.keys()).index(var)
                varstr = osc_list[var][0]
            except ValueError:
                print("Variable not available. Supported variables in osc files:")
                for variable, explanation in osc_list.items():
                    print(f"{variable:>10}: "+ explanation[0] +f", {explanation[1]}")
                print('Please enter a valid variable name from the list.')
                raise
        elif filetype == 'conv.dat':
            try:
                ivar = list(conv_dat_list.keys()).index(var)
                varstr = conv_dat_list[var][0]
            except ValueError:
                print("Variable not available. Supported variables in conv.dat files:")
                for variable, explanation in conv_dat_list.items():
                    print(f"{variable:>10}: "+ explanation[0] +f", {explanation[1]}")
                print('Please enter a valid variable name from the list.')
                raise
        else:
            raise CESAMError('Only osc or conv.dat filetypes are supported.')

        return ivar, varstr


    def edit_params(self):
        ok = self.params.configure_traits()
        if ok:
            self.params.mkdon()


    def calc_freqsGUI(self, mods=None, edit=True):

        self.fparams.write_fsettings( )

        if not self.finished and \
            not (self.fparams.adipls_new_amdl or self.fparams.osc_code_name_ == 'acor'):
            print('Model not yet calculated...')
            return 0

        if edit:
            self.init_fparamsNB(self.fparams.modes_, keep=True)
            ok = self.fparams.edit_traits()
        else:
            ok = True

        if self.fparams.reduce_ != 'none':
            self.read_hr( )

        if ok:
            if self.fparams.osc_code_name_=='adipls':
                self.calc_freqs_adipls( mods=mods, from_gui=True )
            elif self.fparams.osc_code_name_=='acor':
                self.calc_freqs_acor( mods=mods )


    def redistrbGUI(self, amdl, gui=True):

        ok = True
        if gui:
            ok = self.fparams.configure_traits(view='remesh_view')

        if ok:
            self.redistrb( amdl )


    def remesh(self, rep=None, don=None):

        ok = self.fparams.configure_traits(view='remesh_view')

        if ok:

            if not rep:
                rep = self.file_dat
            if not don:
                don = self.name

            if self.fparams.remesh_ == 'n':
                cmd = f'{rep}\n{don}\nn\ny\n'
            elif self.fparams.remesh_ in ('p', 'g', 'r', 'o'):
                cmd = f'{rep}\n{don}\ny\no\n{self.fparams.cts[0]}, {self.fparams.cts[1]}, {self.fparams.cts[2]}\n{self.fparams.npoints:d}\ny\n'
            else:
                cmd = f'{rep}\n{don}\ny\na\n{self.params.npoints}\ny\n'

            print("Re-meshing model and converting to amdl...", end='')
            s = subp.Popen('rep-osc_v3.x', stdin=subp.PIPE,
                                 stdout=subp.PIPE, stderr=subp.PIPE)
            res = s.communicate(input=cmd)[0]
            print(f"{GREEN}[Done]{NO_COLOR}\n")
            messg = res.split('\n')[-2]
            print(f"{GREEN_U}{messg}{NO_COLOR}")


    def plot_ulines(self, x, xstep=10, ny=50, mag=False, nc=20,
                    colors=None, smooth=False, arrows=False,
                    u_csi=None, quadrant=False):

        ntheta = 1000
        #nr = 200
        n = len(x)
        nq = 10

        if quadrant:
            theta = np.linspace(0.0, 0.5*np.pi, ntheta)
        else:
            theta = np.linspace(0.0, 2*np.pi, ntheta)

        if not isinstance(u_csi, np.ndarray):
            if smooth:
                y = self.U_csis
            else:
                y = self.U_csi
        else:
            y = u_csi

        r,t = np.meshgrid(x, theta)
        rr,tt = np.meshgrid(y, theta)

        Z = (np.cos(t)**3 - np.cos(t))*rr
        Z = np.log(np.abs(Z))
        X = r*np.sin(tt)
        Y = r*np.cos(tt)
        saved = plt.rcParams['contour.negative_linestyle']
        plt.rcParams['contour.negative_linestyle'] = 'solid'
        plt.contour(X, Y, np.abs(Z), nc, colors=colors)
        plt.rcParams['contour.negative_linestyle'] = saved

        if arrows:
            theta = [0.0, 0.5*np.pi, np.pi, 1.5*np.pi]
            xx = x[::n//nq]
            r,t = np.meshgrid(xx, theta)

            U2 = self.U2s[::n//nq]

            ur,tt = np.meshgrid(U2, theta)

            # U = self.lp2cos(t)*ur
            U = np.sign(ur)*np.cos(2*t)

            X = r*np.sin(t)
            Y = r*np.cos(t)

            Ux = U*np.sin(t)
            Uy = U*np.cos(t)

            plt.quiver(X, Y, Ux, Uy, pivot='tip')



    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def plot_uvec(self, x, xstep=10, ny=50, mag=False, normalize=False,
                  smooth=False, maxlength=None, key=None, qx=0.5, qy=0.93,
                  u=None, v=None, quadrant=False):
        dtheta = np.pi/ny
        ntheta = ny

        if quadrant:
            theta = np.arange(0.0,0.5*np.pi+dtheta,dtheta)
        else:
            theta = np.linspace(0.0, 2*np.pi, ntheta)

        xx = x[::xstep]

        if not isinstance(u, np.ndarray):
            if smooth:
                u = self.U2s
            else:
                u = self.U2
        if not isinstance(v, np.ndarray):
            if smooth:
                v = self.V2s
            else:
                v = self.V2

        r,t = np.meshgrid(xx, theta)

        ur,tt = np.meshgrid(u[::xstep], theta)
        vr,tt = np.meshgrid(v[::xstep], theta)

        U = lp2cos(t)*ur
        V = dlp2cos(t)*vr

        X = r*np.sin(t)
        Y = r*np.cos(t)

        if mag and not normalize:
            Z = np.sqrt(Ux**2 + Uy**2)
            plt.pcolor(X, Y, Z)
            plt.colorbar()

        mask = U != 0.0
        if normalize:
            Ux = (U*np.sin(t) + V*np.cos(t))/np.sqrt(U**2 + V**2)
            Uy = (U*np.cos(t) - V*np.sin(t))/np.sqrt(U**2 + V**2)
            q = plt.quiver(X[mask], Y[mask], Ux[mask], Uy[mask], scale=40.0,
                       angles='xy')
        else:
            if maxlength:
                length = np.sqrt(U**2 + V**2)
                lims = length > maxlength
                U[lims] *= maxlength/length[lims]
                V[lims] *= maxlength/length[lims]
            else:
                maxlength = max(np.abs(u))
            Ux = (U*np.sin(t) + V*np.cos(t))/maxlength
            Uy = (U*np.cos(t) - V*np.sin(t))/maxlength
            q = plt.quiver(X[mask], Y[mask], Ux[mask], Uy[mask], scale=40.0,
                       angles='xy')

        if key and not normalize:
            leng = 1.4*maxlength
            a2 = int(np.floor(np.log10(np.abs(leng))))
            a1 = leng/10**a2

            strg = f"${a1:.2f} \\times 10^{{{a1}}}\, {{\rm {key}}}$"
            plt.quiverkey(q, qx, qy, 1.4, strg, coordinates='figure')

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def streamplot_u(self, x, smooth=False, density=5, color=None,
                     thickness=None, u=None, v=None, npt=500):

        Y, X = np.mgrid[min(x):max(x):500j, min(x):max(x):500j]


        if not isinstance(u, np.ndarray):
            if smooth:
                u = self.U2s
            else:
                u = self.U2
        if not isinstance(v, np.ndarray):
            if smooth:
                v = self.V2s
            else:
                v = self.V2

        R = np.sqrt(X**2 + Y**2)
        theta = np.arctan(X/Y)

        fU = interpolate.interp1d(x, u)
        fV = interpolate.interp1d(x, v)

        U2 = np.zeros((npt, npt))
        V2 = np.zeros((npt, npt))

        msk = (R < max(x))

        U2[msk] = fU(R[msk])
        V2[msk] = fV(R[msk])

        Ur = 0.5*U2*(3.0*np.cos(theta)**2 - 1.0)
        Utheta = -3.0*V2*np.cos(theta)*np.sin(theta)

        Ux = Ur*np.sin(theta) + Utheta*np.cos(theta)
        Uy = Ur*np.cos(theta) - Utheta*np.sin(theta)

        U = np.sqrt(Ux**2 + Uy**2)
        msk = [U > 0.0]

        th = np.zeros(U.shape)
        th_min = -np.log10(min(U[msk]))
        th[msk] = np.log10(U[msk]) + th_min

        if not color:
            plt.streamplot(X, Y, Ux, Uy, density=density, color=U)
            plt.colorbar()
        elif not thickness:
            plt.streamplot(X, Y, Ux, Uy, density=density, color=color)
        else:
            plt.streamplot(X, Y, Ux, Uy, density=density, color=color,
                       linewidth=th)


    def plot_echelle(self, lstyles='', dnu=None, interactive=False, nu0=0.0,
                     freq_obs=None, l_obs=None, nu0_obs=0.0, legn=True,
                     with_lines=False, freq_errors=None):
        if interactive:
            ec = itr.Echelle(self.l, self.freq, lsep=self.ls, dnu=dnu, freq_obs=freq_obs,
                         l_obs=l_obs, with_lines=with_lines)
            ec.configure_traits()
        else:
            freqs = []

            lmax = max(self.l) + 1
            lmin = min(self.l)

            for l in range(lmin, lmax):
                freqs.append(self.freq[self.l == l])

            if len(lstyles)==0: lstyles = ['bo','g^','rs','cD']
            nstyles = len(lstyles)

            if dnu is None:
                nn = len(self.ls[0])
                nmin = min(15, nn)
                nmax = min(25, nn)
                dnu = np.mean(self.ls[0][nmin:nmax])
                print(f'dnu = {dnu} muHz')

            nus = []
            nub = []
            for i, j in enumerate(freqs):
                nus.append((j - nu0) % dnu)
                ind = np.array(map(int, j/dnu))
                nub.append(ind*dnu + nu0)


            xmin = 0.0
            xmax = max(nus[0])
            ymin = min(nub[0])
            ymax = max(nub[0])
            plt.xlim(xmin, xmax)
            plt.ylim(ymin, ymax)
            leg = []
            for l in range(lmin, lmax):
                il = l - lmin
                if with_lines:
                    plt.plot(nus[il], nub[il], lstyles[il%nstyles]+'-', ms=8)
                else:
                    plt.plot(nus[il], nub[il], lstyles[il%nstyles], ms=8)

                leg.append(f'$\\ell={l:d}$')

            if legn: plt.legend(leg)
            plt.xlabel(r'$\nu (\mu{\rm Hz})$')

            if freq_obs is not None:
                freqs = []

                for l in range(lmin, lmax):
                    freqs.append(freq_obs[l_obs == l])
                if freq_errors is not None:
                    freq_err = []
                    for l in range(lmin, lmax):
                        freq_err.append(freq_errors[l_obs == l])

                nus = []
                nub = []
                for i, j in enumerate(freqs):
                    nus.append((j - nu0_obs) % dnu)
                    ind = np.array(map(int, j/dnu))
                    nub.append(ind*dnu + nu0_obs)


                for l in range(lmin, lmax):
                    il = l - lmin
                    if with_lines:
                        plt.plot(nus[il], nub[il], lstyles[il%nstyles]+':', ms=4.0)
                    else:
                        plt.plot(nus[il], nub[il], lstyles[il%nstyles], ms=4.0)
                    if freq_errors is not None:
                        plt.errorbar(nus[il], nub[il], xerr=freq_err[il],
                                 fmt=lstyles[il%nstyles], ms=4)



    def plot_echelle_p(self, dper):
        ec = itr.EchellePer(self.l, self.per, dper)
        ec.configure_traits()

    def plot_amde(self, filename=None):
        if not filename: filename = self.amde
        pef = itr.PlotEigenFunction(filename)
        pef.configure_traits()

