#!/usr/bin/env python3
# coding: utf-8

from pycesam import tools
from pycesam.gui import *
from pycesam.gui.mpl import MPLFigureEditor
from traits.api import List, Str, Button, observe
from traitsui.api import EnumEditor, ButtonEditor, View, Item, HGroup, VGroup
import matplotlib


class MPLHandler(Handler):

    def close(self, info, is_ok):
        plt.close(info.object.figure)

        return True


class GetNumModel:

    def __init__(self, mdl, ax):
        self.mdl = mdl
        self.ax = ax
        self.mdl.hr_axis(ax=self.ax)
        self.line, = self.ax.plot(mdl.log_teff, mdl.log_l, picker=2)
        self.pt_hr, = self.ax.plot([], [], 'o')
        self.pt = 0

    def onpick(self, event):
        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind
        self.pt = sum(ind)/len(ind)
        self.pt_hr.set_data([xdata[self.pt]], [ydata[self.pt]])
        self.ax.figure.canvas.draw()
        print('Model #%d, age = %.3e Myr' % (self.pt, self.mdl.age[self.pt]))

    def connect(self):
        self.ax.figure.canvas.mpl_connect('pick_event', self.onpick)


class PlotOscGUI(HasTraits):

    __dict_val = {'select_model' :
    {'Model index' : 'index',
     'Age' : 'age',
     'Fraction of remaining central H1' : 'xc',
     'Central density' : 'rhoc',
     'Central temperature' : 'Tc'}}
    models       = List(Instance(CModelGUI))
    figure       = Instance(matplotlib.figure.Figure, ())
    varnames     = List(comparison_mode=2)
    nmod         = Int(1)
    imod         = Range(0,   1000000, 0,   label='Model index')
    tage         = Range(0.0, 13800.0, 0.0, label='Model age')
    txc          = Range(0.0, 2.0,     1.0, label='Fraction of remaining central H1')
    trhoc        = Range(0.0, 1e6,     100, label='Central density')
    tTc          = Range(0.0, 1e10,    1e7, label='Central temperature')
    select_model = Trait('Model index', __dict_val['select_model'],
        label='Select model according to:')
    modnames     = List(comparison_mode=2)
    xn           = Str('r')
    yn           = Str('m/M')
    x_log        = Bool(False)
    x_inv        = Bool(False)
    y_log        = Bool(False)
    y_inv        = Bool(False)
    switch       = Button(u'\u21C5')

    traits_view = \
        View(
            Item('figure', editor=MPLFigureEditor(), show_label=False),
            Item( 'select_model' ),
            Item('imod',  visible_when="select_model_ == 'index'",
                editor=RangeEditor(mode='spinner', high_name='nmod', low=0)),
            Item('tage',  visible_when="select_model_ == 'age'",
                editor=RangeEditor(mode='xslider', high=13800, is_float=True)),
            Item('txc',   visible_when="select_model_ == 'xc'",
                editor=RangeEditor(mode='slider', is_float=True)),
            Item('trhoc', visible_when="select_model_ == 'rhoc'",
                editor=RangeEditor(mode='xslider', low=0.0, high=1e6, is_float=True)),
            Item('tTc',   visible_when="select_model_ == 'Tc'",
                editor=RangeEditor(mode='xslider', low=0.0, high=1e10, is_float=True)),
        VGroup(
            HGroup(Item('xn', editor=EnumEditor(name='varnames')),
                   Item('x_log', label='log'),
                   Item('x_inv', label='inv')),
            HGroup(
                Item( '30' ),
                Item( 'switch', show_label=False )),
            HGroup(Item('yn', editor=EnumEditor(name='varnames')),
                   Item('y_log', label='log'),
                   Item('y_inv', label='inv'))),
        width=800, height=600,
        resizable=True, handler=MPLHandler()
        )

    def __init__( self, models ):
        """
        Controls the GUI that displays structure of models.

        :param models: list of Cesam2k20 models.
        :type models: list of <CModelGUI>

        :member xn: Name of the variable on x-axis.
        :mtype xn: str

        :member yn: Name of the variable on y-axis.
        :mtype yn: str

        :member x_log: True if x-axis is in log scale.
        :mtype x_log: bool

        :member x_inv: True if x-axis is inverted.
        :mtype x_inv: bool

        :member y_log: True if y-axis is in log scale.
        :mtype y_log: bool

        :member y_inv: True if y-axis is inverted.
        :mtype y_inv: bool

        """

        self.models  = models
        self.Nmodels = len( self.models )
        self.imods   = [self.imod for _ in range(self.Nmodels)]
        self.get_varnames( )
        self.all_osc = np.array( [m.params.nom_output_[:3] for m in self.models] ) == 'all'
        self.ages    = [0.0 if self.all_osc[i] else self.models[i].age[-1] for i in range(self.Nmodels)]
        self.set_nmod( )
        self.build_dictionnary( )
        self.set_composite_vars( )
        self.initial_plot()

    def get_varnames( self ):
        """
        Retrieves the unicode name of variables avalable in the structure file.

        :raises PYCESAMGUIError: If all osc files are not in the same format ('nadia', 'adia', etc.)
        """
        output_short = [m.params.nom_output_[4:] for m in self.models]
        if len( set(output_short) ) > 1:
            raise PYCESAMGUIError( "all outputs should be of the same type ('nadia', "
                "'adia', 'plato' or 'invers'). You can still mix output for all or only "
                "last model." )

        self.varnames = self.models[-1].unom_vars

    def build_dictionnary( self ):
        """
        From the variables available for each model, builds a dictionnary with mode names
        and unicode name of each variables as key.
        """
        self.var_dict = {}
        for i, m in enumerate( self.models ):
            self.var_dict[m.name] = {}
            for j, key in enumerate( self.varnames ):
                if self.all_osc[i]:
                    self.var_dict[m.name][key] = []
                    for k in range( m.nmod ):
                        self.var_dict[m.name][key].append( m.var[k][j] )
                else:
                    self.var_dict[m.name][key] = [m.var[j]]

    def set_nmod( self ):
        """
        Get maximum number of vailable structures for all models.
        """
        self.nmod = self.models[0].nmod
        if self.Nmodels > 1:
            for i, m in enumerate( self.models[1:], start=1 ):
                if self.all_osc[i]:
                    self.nmod = max( self.nmod, m.nmod )

        self.nmod -= 1

    def __get_N2( self, model, i ):
        """
        Computes Brunt-Vaissala frequency.

        :param model: Cesam2k20 model.
        :type model: instance of <CModelGUI>

        :param i: Index of the model in the list.
        :type i: int

        :return: Brunt-Vaissala frequency.
        :rtype: list of np.ndarray
        """
        if self.all_osc[i]:
            return model.get_N2( )
        else:
            return [model.get_N2( )]

    def __get_rR( self, model, i):
        """
        Computes r / Rstar.

        :param model: Cesam2k20 model.
        :type model: instance of <CModelGUI>

        :param i: Index of the model in the list.
        :type i: int

        :return: r/Rstar.
        :rtype: list of np.ndarray
        """
        if self.all_osc[i]:
            var = []
            for k in range( model.nmod ):
                var.append( model.var[k][0] / (model.rstar[k]*model.ctes.rsun) )
        else:
            var = [model.var[0] / (model.rstar[-1]*model.ctes.rsun)]
        return var

    def __get_rRsun( self, model, i):
        """
        Computes r/Rsun.

        :param model: Cesam2k20 model.
        :type model: instance of <CModelGUI>

        :param i: Index of the model in the list.
        :type i: int

        :return: r/Rsun.
        :rtype: list of np.ndarray
        """
        if self.all_osc[i]:
            var = []
            for k in range( model.nmod ):
                var.append( model.var[k][0] / model.ctes.rsun )
        else:
            var = [model.var[0] / model.ctes.rsun]
        return var

    def __get_mMsun( self, model, i):
        """
        Computes m/Msun.

        :param model: Cesam2k20 model.
        :type model: instance of <CModelGUI>

        :param i: Index of the model in the list.
        :type i: int

        :return: m/Msun.
        :rtype: list of np.ndarray
        """
        if self.all_osc[i]:
            var = []
            for k in range( model.nmod ):
                var.append( model.var[k][1] * model.mstar[k] )
        else:
            var = [model.var[1] * model.mstar[-1]]
        return var

    def __set_new_var( self, func, rname, uname ):
        """
        Defines a new variable visible through the GUI.

        :param func: Function that will be applied to compute the new quantity.
            Should be of the form func( model, i ), with model a CModelGUI
            instance and i an integer.
        :type func: function

        :param rname: Latex name of the variable.
        :type rname: raw-string

        :param uname: Unicode name of the variable.
        :type uname: unicode-string
        """
        self.varnames.append( uname )
        for i, m in enumerate( self.models ):
            new_var = func( m, i )
            m.nom_vars.append( rname )
            m.unom_vars.append( uname )
            self.var_dict[m.name][uname] = new_var

    def set_composite_vars( self ):
        """
        Defines all the new varaibles you need.
        """
        self.__set_new_var( self.__get_N2, r'$N^2$', u'N\u00b2' )
        self.__set_new_var( self.__get_rR, r'$r/R_\star$', u'r/R' )
        self.__set_new_var( self.__get_rRsun, r'$r/R_\odot$', u'r/Rsun' )
        self.__set_new_var( self.__get_mMsun, r'$m/M_\odot$', u'm/Msun' )

    def initial_plot( self ):
        """
        Creates the initial plots with default settings.
        """
        self.ix    = self.varnames.index(self.xn)
        self.iy    = self.varnames.index(self.yn)
        self.lines = []
        self.pts_hr = []
        self.xsc = 'linear'
        self.ysc = 'linear'

        any_all_osc = any( self.all_osc )

        if any_all_osc:
            self.figure, ((self.ax_osc, self.ax_hr)) = plt.subplots(1,2)
        else:
            self.figure, ((self.ax_osc)) = plt.subplots(1,1)
            self.ax_hr = None

        self.colors     = list( plt.rcParams['axes.prop_cycle'] )

        for i, m in enumerate( self.models ):
            color = self.colors[i]['color']
            imod  = int( self.imods[i] )
            imod_hr  = self.imods[i] if self.all_osc[i] else -1
            if any_all_osc:

                # plot osc variable
                line, = self.ax_osc.plot( self.var_dict[m.name][self.xn][imod],
                    self.var_dict[m.name][self.yn][imod], color=color, label=m.name )
                self.lines.append( line )

                # plot HR diagram
                m.plot_hr( ax=self.ax_hr, color=color )

                pt_hr, = self.ax_hr.plot([m.log_teff[imod_hr]], [m.log_l[imod_hr]], 'o', color=color,
                    label=self.ages[i])
                self.pts_hr.append( pt_hr )
            else:
                line, = self.ax_osc.plot(self.var_dict[m.name][self.xn][imod],
                                    self.var_dict[m.name][self.yn][imod], color=color)
                self.lines.append( line )

        self.ax_osc.set_xlabel(m.nom_vars[self.ix])
        self.ax_osc.set_ylabel(m.nom_vars[self.iy])
        self.ax_osc.legend( loc='best', fontsize=6 ).get_frame().set_linewidth(0.0)
        if self.ax_hr is not None:
            self.ax_hr.legend(  loc='best', fontsize=6 ).get_frame().set_linewidth(0.0)


    @observe('xn')
    def _xn_changed(self):
        """
        Changes variale on x-axis. Called whenever `self.xn` is changed.
        """
        for i, m in enumerate( self.models ):
            self.ix = self.varnames.index(self.xn)
            self.ax_osc.set_xlabel(m.nom_vars[self.ix])

            imod = int( self.imods[i] if self.all_osc[i] else 0 )
            self.update( i, imod )
        self.__redraw( )

    @observe( 'yn' )
    def _yn_changed(self):
        """
        Changes variable on y-axis. Called whenever `self.yn` is changed.
        """
        for i, m in enumerate( self.models ):
            self.iy = self.varnames.index(self.yn)
            self.ax_osc.set_ylabel(m.nom_vars[self.iy])
            imod = int( self.imods[i] if self.all_osc[i] else 0 )
            self.update( i, imod )
        self.__redraw( )

    @observe( 'x_log' )
    def _x_log_changed( self ):
        """
        Changes, for x-axis, log scale to linear and vice-versa.
        Called whenever `self.x_log` is changed.
        """
        self.xsc = 'log' if self.x_log else 'linear'
        self.ax_osc.set_xscale(self.xsc)
        self.figure.canvas.draw()

    @observe( 'y_log' )
    def _y_log_changed( self ):
        """
        Changes, for y-axis, log scale to linear and vice-versa.
        Called whenever `self.y_log` is changed.
        """
        self.ysc = 'log' if self.y_log else 'linear'
        self.ax_osc.set_yscale(self.ysc)
        self.figure.canvas.draw()

    @observe( 'x_inv' )
    def _x_inv_changed( self ):
        """
        Inverses x-axis. Called whenever `self.x_inv` is changed.
        """
        self.ax_osc.invert_xaxis()
        self.figure.canvas.draw()

    @observe( 'y_inv' )
    def _y_inv_changed( self ):
        """
        Inverses y-axis. Called whenever `self.y_inv` is changed.
        """
        self.ax_osc.invert_yaxis()
        self.figure.canvas.draw()

    def _switch_fired( self ):
        """
        Switches variables on x and y axes. Called whenever `self.switch` is changed.
        """
        self.xn, self.yn = swap( self.xn, self.yn )
        self.x_inv, self.y_inv = swap( self.x_inv, self.y_inv )
        self.x_log, self.y_log = swap( self.x_log, self.y_log )

    def _tage_changed(self):
        """
        Finds index of a model's time step with age closest to a target age.
        Then redraws figure.
        """
        for i, m in enumerate( self.models ):
            if self.all_osc[i]:
                self.imods[i] = np.argmin( np.abs( m.age - self.tage) )
                self.__update_age( i )

                self.update( i, self.imods[i] )

        self.__redraw( )

    def _txc_changed(self):
        """
        Finds index of a model's time step with core hydrogen to initial core
        hydrogen ratio closest to a target ratio.
        Then redraws figure.
        """
        for i, m in enumerate( self.models ):
            if self.all_osc[i]:
                self.imods[i] = m.nmod - np.argmin( np.abs( m.ab_c['H1'][::-1] - self.txc*m.ab_c['H1'][0]) ) - 1
                self.__update_age( i )

                self.update( i, self.imods[i] )

        self.__redraw( )

    def _trhoc_changed(self):
        """
        Finds index of a model's time step with core density closest to a target core density.
        Then redraws figure.
        """
        for i, m in enumerate( self.models ):
            if self.all_osc[i]:
                rhoc = np.array([m.var[k][4][-1] for k in range( m.nmod )])
                self.imods[i] = np.argmin( np.abs( rhoc - self.trhoc ) )
                self.__update_age( i )

                self.update( i, self.imods[i] )

        self.__redraw( )

    def _tTc_changed(self):
        """
        Find index of a model's time step with core temperature closest to a target core temperature.
        Then redraws figure.
        """
        for i, m in enumerate( self.models ):
            if self.all_osc[i]:
                Tc = np.array([m.var[k][2][-1] for k in range( m.nmod )])
                self.imods[i] = np.argmin( np.abs( Tc - self.tTc ) )
                self.__update_age( i )

                self.update( i, self.imods[i] )

        self.__redraw( )

    def _imod_changed(self):
        """
        Redraws figure with structure corresponding to desired `self.imod`.
        """
        for i, m in enumerate( self.models ):
            self.imods[i] = int( min( self.imod, m.nmod-1) if self.all_osc[i] else 0 )
            self.__update_age( i )

            self.update( i, self.imods[i] )

        self.__redraw( )

    def __update_age( self, i ):
        """
        Gets new value of model's age every time imods changes.

        :param i: Index of the mdoel.
        :type i: int
        """
        iage = self.imods[i] if self.all_osc[i] else -1
        self.ages[i] = 'Age = {:8.2f}'.format(self.models[i].age[iage])


    def __redraw( self ):
        """
        Redraws the figure with new settings.
        """
        self.ax_osc.relim()
        self.ax_osc.autoscale_view()
        if self.ax_hr is not None:
            self.ax_hr.legend(  loc='best', fontsize=6 ).get_frame().set_linewidth(0.0)
        self.figure.canvas.draw()


    def update( self, i, imod ):
        """
        Update line data of given model.

        :param i: Index of the model in the list.
        :type i: int

        :param imod: Index of the time step whose values should be displayed.
        :param imod: int
        """
        m = self.models[i]
        self.lines[i].set_data(self.var_dict[m.name][self.xn][imod],
                               self.var_dict[m.name][self.yn][imod])
        if self.all_osc[i]:
            self.pts_hr[i].set_data([m.log_teff[imod]], [m.log_l[imod]])
            self.pts_hr[i].set_label(self.ages[i])


def polar2cartesian(r, t, grid, x, y, order=3):
    X, Y = np.meshgrid(x, y)

    new_r = np.sqrt(X*X+Y*Y)
    new_t = np.arctan2(X, Y)

    ir = interpolate.interp1d(r, np.arange(len(r)), bounds_error=False)
    it = interpolate.interp1d(t, np.arange(len(t)), bounds_error=False)

    new_ir = ir(new_r.ravel())
    new_it = it(new_t.ravel())

    new_ir[new_r.ravel() > r.max()] = len(r) - 1
    new_ir[new_r.ravel() < r.min()] = 0

    return plt.map_coordinates(grid, np.array([new_ir, new_it]),
                           order=order).reshape(new_r.shape)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def movie(files, output="output.avi", fps=25, erase=True):
    """
    <movie(files, output="output.avi", fps=25, erase=True)>

    Makes a movie from image files.

    Parameters
    ----------
    files : string
        A string that contains the specifications describing the
        image files. Example: '*.png'
    output : string
        Name of the output video file. Default: `output.avi`
    fps : integer
        Frames per second, default: 25
    erase : boolean
        If True, deletes image files after creating movie.
        Default: delete.
    """

    command = 'mencoder mf://%s -mf type=png:w=800:h=600:fps=%d -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o %s' % (files, fps, output)
    os.system(command)
    if erase: os.system("rm -f %s" % files)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def lp2cos(x):
    return 0.5*(3.0*np.cos(x)**2 - 1.0)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def dlp2cos(x):
    return -3.0*np.cos(x)*np.sin(x)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def plotpolar(x, y, label=None, quadrant=False, vmin=None, vmax=None, cmap=None):
    """
    <plotpolar(x, y, label=None, quadrant=False, vmin=None, vmax=None)>
    """

    ntheta = 200
    nr = 200
    n = len(x)

    if quadrant:
        theta = np.linspace(0.0, 0.5*np.pi, ntheta)
    else:
        theta = np.linspace(0.0, 2*np.pi, ntheta)

    x1 = x[::n/nr]
    x1 = np.append(x1, x[-1])
    y1 = y[::n/nr]
    y1 = np.append(y1, y[-1])

    r,t = np.meshgrid(x1, theta)
    rr,tt = np.meshgrid(y1, theta)
    Z = lp2cos(t)*rr

    X = r*np.sin(tt)
    Y = r*np.cos(tt)
    if cmap:
        plt.pcolor(X, Y, Z, vmin=vmin, vmax=vmax, cmap=cmap)
    else:
        plt.pcolor(X, Y, Z, vmin=vmin, vmax=vmax)
    cb = plt.colorbar()
    if label: cb.set_label(label)


def get_num_model(mdl, mods=None, mod_nums=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    if mods:
        if not isinstance(mods, list): mods = [mods]
        if mod_nums:
            if not isinstance(mod_nums, list): mod_nums = [mod_nums]
            if len(mod_nums) != len(mods):
                try:
                    raise IndexError
                except IndexError:
                    print('mods and mod_nums must have the same size.')
            else:
                for md, mn in zip(mods, mod_nums):
                    md.plot_hr(ax=ax,mods=[mn])
        else:
            for md in mods:
                md.plot_hr(ax=ax)

    num_mod = GetNumModel(mdl, ax)
    num_mod.connect()
    plt.show()

    return num_mod.pt

def set_plot_page(profile, dpi=100.0):
    inch = 25.4
    # plt.rcParams['font.serif'] = 'Times New Roman'
    plt.rcParams['figure.dpi'] = dpi
    if profile == 'aa_onecolumn':
        aspect = 4.0/3.0
        plt.rcParams['font.size'] = 12.0
        w = 180.0/inch
        plt.rcParams['figure.figsize'] = [w, w/aspect]
    elif profile == 'aa_twocolumn':
        aspect = 4.0/3.0
        plt.rcParams['font.size'] = 12.0
        w = 88.0/inch
        plt.rcParams['figure.figsize'] = [w, w/aspect]
    elif profile == 'presentation_full':
        aspect = 4.0/3.0
        plt.rcParams['font.size'] = 24.0
        w = 280.0/inch
        plt.rcParams['figure.figsize'] = [w, w/aspect]
    elif profile == 'presentation_half':
        aspect = 4.0/3.0
        plt.rcParams['font.size'] = 16.0
        w = 120.0/inch
        plt.rcParams['figure.figsize'] = [w, w/aspect]
    else:
        print('Not a valid option, will do nothing.')

