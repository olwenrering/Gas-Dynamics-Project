#!/usr/bin/env python3
# coding: utf-8

# Code d'Evolution Stellaire, Adaptatif et Modulaire for the 2020 decade.
# Copyright (c) 1997-2023 The Cesam2k20 authors
# SPDX-License-Identifier : GPL-3.0-or-later


# try to import the new traits library
try:
    from traits.api           import HasTraits, Trait, Float, Int, Bool, File, \
        List, Range, String
# look for the old traits library if the new one is not found
except ImportError:
    from enthought.traits.api import HasTraits, Trait, Float, Int, Bool, File, \
        List, Range, String
import os
import glob
import numpy      as np
import subprocess as subp
from threading    import Thread

from pycesam.constants import FPARAMS_DICT
from pycesam.tools import str2type

GREEN     = "\033[32;01m"
NO_COLOR  = "\033[0m"

class Freqs(HasTraits):

    # Dictionaries for mapped lists
    _dict_val        = FPARAMS_DICT

    osc_code_name    = Trait('ADIPLS', _dict_val['osc_code_name'], label='Oscillation code')
    to_zip           = Bool(False, label="Store agsm files in a zip.")
    ##################
    # ADIPLS options
    # Type of modes ('p' or 'g')
    adipls_new_amdl  = Bool(False, label="Use another amdl file ?")
    adipls_amdl_name = String('', label='Name of amdl file')
    modes            = Trait('p-modes', _dict_val['modes'], label='Kind of modes')
    reduce           = Trait('None', _dict_val['reduce'], label='Reduce to some modes ?')
    nred             = Range(1, 10000, 5, label='Number of width of the Gaussian before/after nu_max')
    # Correction for radial order for l=1. Default: 0.
    irsord          = Int(20, label='Setting for order of modes')
    # Method for integrating eqs. Default: 1 (shooting)
    mdintg          = Trait('Shooting method', _dict_val['mdintg'], label='Method of integration')
    # Number of points to use. Default: 2400
    npoints         = Int(2400, label='Number of points')
    iekinr          = Trait('With total photospheric displacement', _dict_val['iekinr'], label='Normalization of energy')
    icaswn          = 10010
    # If True, calculates rotational kernel
    rotkr           = Bool(False, label='Rotational kernel')
    # If True, calculates Gamma_1 kernel
    gm1kr           = Bool(False, label='Gamma1 kernel')
    # If True, outputs eigenfunctions
    amde            = Bool(False, label='Output eigenfunctions')
    iscan           = Int(250, label='Number of steps')

    lmin            = Range(0, 10000, 0, label='Minimum l')
    lmax            = Range(0, 10000, 3, label='Maximum l')
    dels            = Range(0, 10000, 1, label='Step in l')
    nsel            = Int(4, label='Number of ls')
    step            = Range(1, 10000, 1, label='Step')

    ivarf           = Int(1)
    all_osc         = Bool(False)
    nsig            = Trait('Frequency', _dict_val['nsig'], label='Step in')
    sig1            = Float(4.0, label=u'\u03C3\u2081\u00B2')
    sig2            = Float(3000.0, label=u'\u03C3\u2082\u00B2')
    remesh          = Trait('Just increase mesh density', _dict_val['remesh'], label='Re-meshing options')

    cts             = List([10.0, 1.0e-2, 1.5e-2], label='Re-meshing parameters')

    iriche          = Bool(True, label='Richardson extrapolation')
    istsbc          = Trait('Exponentially decaying solution', _dict_val['istsbc'], label='Surface condition')
    fctsbc          = Trait(u'\u03B4p = 0', _dict_val['fctsbc'], label='Function fsbc')
    fsbc            = Float( 0.0 )

    ##################
    # ACOR options
    input_file      = Trait('Compute frequencies from 1D output files', _dict_val['input_file'], label='Input file')
    step1d          = Range(1, 10000, 1, label='Step 1D')
    step2d          = Range(1, 10000, 1, label='Step 2D')
    only_last       = Bool(False, label="Only last model")
    mod_flag        = Trait('Model name', _dict_val['mod_flag'], label='ACOR frequency file name')
    mod_flag_custom = String('', label='Custom model flag')
    mod_suff        = String('', label='Model suffixe (added after model flag)')
    acor_star_type  = Trait('Low mass star', _dict_val['acor_star_type'], label='Type of star')
    try:
        acor_exe_path = String(os.environ['acordir'] + '/bin/', label='Path to ACOR executable files')
    except KeyError:
        try:
            acor_exe_path = String(os.environ['cesacordir'] + '/bin/', label='Path to ACOR executable files')
        except KeyError:
            acor_exe_path = String(os.environ['CESDIR'] + '/ACOR/bin/', label='Path to ACOR executable files')

    out_freq_path       = String(os.getcwd(), label='Path location for the storage of output frequencies')
    out_eigen_path      = String(os.getcwd(), label='Path location for the storage of eigenfunctions')
    acor_rot            = Bool(True, label='Rotation ?')
    acor_rot_grad       = Bool(True, label='Rotation gradients ?')
    acor_rot_toro       = Bool(True, label='Toroidal part ?')
    acor_rot_prof       = Trait('Differential', _dict_val['acor_rot_prof'], label='Type of rotation profile')
    acor_profile_ext    = Bool(False, label='Rotation profile in an external file ?')
    profile_filename    = String('', label='Profile file name')
    unif_om             = Float(1.0, label=u'Uniform \u03a9 (muHz)')
    diff_om             = Float(1.0, label=u'Central \u03a9 (muHz)')
    biz_rad             = Float(1.0, label=u'Bizonal: transition radius (r/R)')
    biz_oms             = Float(1.0, label=u'Bizonal: surface \u03a9 (rad/s)')
    biz_omc             = Float(1.0, label=u'Bizonal: central \u03a9 (rad/s)')
    lat_omeq            = Float(1.0, label=u'Latitudinal: Equatorial \u03a9 (muHz)')
    lat_dom             = Float(1.0, label=u'Latitudinal:  \u0394\u03a9 (muHz)')
    acor_freq_target    = Trait('Scanning', _dict_val['acor_freq_target'], label='Frequency targets')
    acor_freq_units     = Trait('microHz', _dict_val['acor_freq_units'], label='Frequency units')
    acor_freq_start     = Float(5.0, label=u'Starting Frequency')
    acor_freq_step      = Float(5.0, label=u'Step in Frequency')
    acor_nsteps         = Float(100.0, label=u'Number of steps')
    acor_freq_file_type = Trait('ADIPLS', _dict_val['acor_freq_file_type'], label='Type of frequency file')
    acor_params_select  = Trait('lmax and m', _dict_val['acor_params_select'], label='Select mode parameters')
    acor_freq_filename  = String('', label='Frequency file name')
    acor_m              = Range(-1000, 1000, 0, label='Azimutal order m')
    acor_parity         = Bool(True, label='l of same parity as m ?')
    acor_nsh            = Range(0, 10000, 1, label='Number of spherical harmonics')
    acor_rad_res        = Range(1, 10000, 1, label='Radial resolution (1:the model one, 2:half of it, etc..)')
    acor_newt_it        = Range(1, 10000, 15, label='Number of Newton iterations')
    acor_follow         = Bool(False, label='Follow mode with respect of rotation ?')
    acor_eigen          = Bool(False, label='Print out eigenfunctions ?')
    acor_rm_input       = Bool( True, label='Remove ACOR_input file at the end?' )


    def __init__( self, name, all_osc, verbose=True ):
        """
        Class that represents ADIPLS or ACOR input parameters

        :param name: Generic name of the model
        :type name: str

        :param all_osc: True if all osc files are written.
        :type all_osc: boolean

        :kparam verbose: If True, prints additional information.
        :ktype verbose: dict
        """
        self.frun_file = name + '.frun'
        self.verbose = verbose

        self.all_osc = all_osc
        self.on_trait_change(self.manage_changes, 'modes_')
        self.on_trait_change(self.l_changes, 'lmin, lmax, dels')
        self.on_trait_change(self.set_cts, 'remesh')
        self.l_changes( )

    def manage_changes(self):
        """
        Adapts settings when the value of `self.modes` changes.
        """
        if self.modes_ == 'p':
            self.ivarf = 1
            self.nsig = 'Frequency'
            self.sig1 = 1.0
            self.sig2 = 6000.0
        elif self.modes == 'g':
            self.ivarf = 2
            self.nsig = '1/Frequency'
            self.sig1 = 0.009
            self.sig2 = 33.0

    def l_changes(self):
        """
        Adapts settings when the value of `self.lmax`, `self.lmin` or `self.dels` changes.
        """
        self.nsel = int( (self.lmax - self.lmin + 1)/self.dels )

    def set_cts(self):
        """
        Adapts remesh settings when `self.remesh` changes
        """
        if self.remesh_ == 'p':
            self.cts = [10.0, 1.0e-2, 1.5e-2]
        elif self.remesh_ == 'g':
            self.cts = [2.5e-2, 1.0e-1, 1.0e-4]
        elif self.remesh_ == 'r':
            self.cts = [1.0e4, 1.0e-2, 1.0e-2]
            self.nsig = 'Red giants'
            self.iscan = 10
            self.mdintg = '4th-order shooting method'
            self.modes = 'all'

    def _fctsbc_changed( self ):
        """
        Method automatically called whenever value of `self.fctsbc` is changed
        """
        if self.fctsbc_ == 0:
            self.fsbc = 0.0
        elif self.fctsbc_ == 1:
            self.fsbc = 1.0

    def _istsbc_changed( self ):
        """
        Method automatically called whenever value of `self.istsbc` is changed
        """
        if self.istsbc_ == 1:
            self.fctsbc_ = 0
            self.fsbc = 0.0
        elif self.istsbc_ == 0:
            if self.fctsbc_ == 0:
                self.fsbc = 0.0
            elif self.fctsbc_ == 1:
                self.fsbc = 1.0

    def write_fsettings(self):
        """
        Writes frequency settings to a `.frun` file.
        """

        adipls_vars = ['adipls_new_amdl', 'adipls_amdl_name', 'modes', 'reduce',
            'nred', 'step', 'lmin', 'lmax', 'dels', 'nsel', 'sig1',
            'sig2', 'iscan', 'nsig', 'remesh', 'cts', 'npoints', 'irsord',
            'mdintg', 'iriche', 'iekinr', 'rotkr', 'gm1kr', 'amde', 'istsbc',
            'fctsbc', 'fsbc', 'to_zip']
        acor_vars = ['input_file', 'step1d', 'step2d', 'only_last',
            'acor_star_type', 'mod_flag', 'mod_flag_custom', 'mod_suff',
            'acor_exe_path', 'out_freq_path', 'acor_eigen', 'out_eigen_path',
            'acor_rot', 'acor_rot_grad', 'acor_rot_toro', 'acor_rot_prof',
            'acor_profile_ext', 'profile_filename', 'unif_om', 'diff_om',
            'biz_rad', 'biz_oms', 'biz_omc', 'lat_omeq', 'lat_dom',
            'acor_freq_target', 'acor_freq_units', 'acor_freq_start',
            'acor_freq_step', 'acor_nsteps', 'acor_freq_file_type',
            'acor_freq_filename', 'acor_params_select', 'acor_nsh', 'acor_m',
            'acor_parity', 'lmax', 'acor_rad_res', 'acor_newt_it',
            'acor_follow', 'acor_rm_input']
        remesh_vars = ['remesh', 'cts', 'npoints']

        fvars = adipls_vars if self.osc_code_name == 'ADIPLS' else acor_vars

        lines = []
        lines.append('{:25s}{:s}\n'.format('osc_code_name', str(self.__getattribute__( 'osc_code_name' ))))
        for var in fvars:
            lines.append('{:25s}{:s}\n'.format(var, str(self.__getattribute__(var))))

        for var in remesh_vars:
            lines.append('{:25s}{:s}\n'.format(var, str(self.__getattribute__(var))))

        with open(self.frun_file, 'w') as f:
            f.writelines(lines)


    def read_fsettings(self):
        """
        Reads frequency settings from a `.frun` file.
        """

        if self.verbose: print(f"Reading {self.frun_file}...", end='')

        with open(self.frun_file, 'r') as f:
            lines = f.readlines()

        isplit = 25

        for line in lines:
            key = line[:isplit].strip()
            val = line[isplit:].strip()

            self.__setattr__( key, str2type( self.__getattribute__( key ), val) )

        if self.verbose: print(f"{GREEN}[Done]{NO_COLOR}")




class CFacor:

    def __init__(self, filename=None):
        """
        Class that defines a data structure for ACOR's oscillation outputs.

        :kparam filename: Name of the file in which oscillation results are stored.
        :ktype filename: str
        """
        self.filename   = filename
        self.number, self.m, self.l, self.nt, self.npt, self.ngt, self.n, self.np, \
            self.ng, self.index = [np.array([], dtype='int') for _ in range(10)]
        self.Sigma_in, self.Sigma_out, self.Sigma_out2, self.freq, self.ekin, self.I, \
            self.Knl, self.Beta, self.I_core, self.Knl_core, self.Beta_core, \
            self.Omk = [np.array([]) for _ in range(12)]


        self.listef_names = {0: r"numero", 1: r"m", 2: r"l", 3: r"n", 4: r"np", 5: r"n",
            6: r"$\sigma_{\rm in}$", 7: r"$\sigma_{\rm out}$", 8: r"$\sigma^2_{\rm out}$",
            9: r"$\sigma_{\rm out}$ (muHz)", 10: r"Ec (mode)", 11: r"$I$ (mode)",
            12: r"$K_{n\ell}$ (mode)", 13: r"$\beta$ (mode)", 14: r"$I_{\rm core}$ (mode)",
            15: r"$K_{n\ell, \rm core}$ (mode)", 16: r"$\beta_{\rm core}$ (mode)",
            17: r"$\Omega_{\rm k} \times 10^{-6} / 2 \np.pi"}

    def merge( self, other_cfacor ):
        """
        Merges a secondary `CFacor` instance to current instance.
        In practice: merges the two attributes `self.__dict__`.

        :param other_cfacor: Instance of `CFacor` whose attributes must be injected into
            current instance.
        :type other_cfacor: <CFacor>
        """
        this_item  = self.__dict__
        other_item = other_cfacor.__dict__

        items = list(this_item.keys())
        items.remove( 'filename' )
        items.remove( 'listef_names' )
        for item in items:
            this_item[item] = np.append( this_item[item], other_item[item] )



class CRunAcor( Thread ):

    def __init__( self, fparams, evol_code="CESAM2k20", nmod=None, d2=False,
        aparam_name='ACOR_parameters', write=True, quiet=False ):
        """
        Class that allows to run two instances of ACOR in parallel. Useful when we compute
        frequencies with odd/even parities.

        :param fparams: Instance of Freqs class that represents frequency parameters.
        :type fparams: <Freqs>

        :param evol_code: Name of the evolution code that computed the structure model. Optional.
        :type evol_code: stringnselsse

        :param nmod: Name of the model that contains the structure. Optional.
        :type nmod: string

        :param d2: If True, the structure model is a 2D model. Optional.
        :type d2: boolean

        :param aparam_name: Name of the ACOR parameter file. Default: `'ACOR_parameters'`.
        :type aparam_name: string

        :param write: If True, the ACOR parameter file is written when instance is created. Optional.
        :type write: boolean

        :param quiet: If True, standard output will not be printed when frequencies are computed. Optional.
        :type quiet: boolean
        """
        self.fparams     = fparams
        self.evol_code   = evol_code
        self.nmod        = nmod
        self.d2          = d2
        self.aparam_name = aparam_name
        self.write       = write
        self.quiet       = quiet
        self.cfacor      = CFacor( )
        self.fname       = self.fparams.mod_flag_custom
        self.modname     = self.fparams.mod_flag_custom
        if self.fparams.input_file_[:2]=='d1':
            if self.fparams.acor_rot:
                self.fname += "_"+self.fparams.mod_suff
            else:
                self.fname += "_norot"

        self.Mmpar_name = "_M%i_m%i_par%i" % (self.fparams.acor_nsh, self.fparams.acor_m,
            not self.fparams.acor_parity)
        self.fname += self.Mmpar_name

        if self.write:
            self.__write_aparam( )

        Thread.__init__(self)

    def __write_aparam( self ):
        """
        Private method that writes the file with ACOR parameters.
        """

        model_path = os.getcwd( )
        with open(self.aparam_name, 'w') as f:
            f.write("### paths and stellar model  ;\n")
            f.write("# directory that contains the model (path);\n")
            f.write(f"{model_path.replace('//', '/')};\n")
            f.write("# evolution code used for the stellar model?;\n")
            f.write(f"{self.evol_code};\n")
            f.write("# file containing the stellar model;\n")
            f.write(self.nmod+";\n")
            f.write("# stellar model flag;\n")
            f.write(self.fparams.mod_flag_custom+"\n")
            f.write("# 2D structural model? 1->yes 0->no;\n")
            f.write("   %i;\n" % (self.d2))
            f.write("# type of star? 1->low mass 2->massive 3 -> Gamma Dors;\n")
            f.write("   %i;\n" % (self.fparams.acor_star_type_))
            f.write("# path location of the Interface and the ACOR program;\n")
            f.write(self.fparams.acor_exe_path.replace('//', '/')+";\n")
            f.write("# path location for storage of the output frequencies;\n")
            f.write(self.fparams.out_freq_path.replace('//', '/')+";\n")
            f.write("# path location for storage of the eigenfunctions ;\n")
            f.write(self.fparams.out_eigen_path.replace('//', '/')+";\n")
            f.write("# Following of the modes with respect to rotation? 1->yes 2->no;\n")
            f.write("   %i;\n" % (self.fparams.acor_follow))
            f.write("# Print out the eigenfunctions? 1->yes 0->no;\n")
            f.write("   %i;\n" % (self.fparams.acor_eigen))
            f.write(";\n")
            f.write("### rotation caracteristics  ;\n")
            f.write("# ACOR with or without rotation: 0 -> without , 1 -> with ;\n")
            f.write("   %i;\n" % (self.fparams.acor_rot))
            f.write("# Computation with or without accounting for rotation gradient?   0-> without,    1-> with;\n")
            f.write("   %i;\n" % (self.fparams.acor_rot_grad))
            f.write("# Computation with or without toroidal part? 0-> without,    1-> with;\n")
            f.write("   %i;\n" % (self.fparams.acor_rot_toro))
            f.write("# Type of rotation profile? 1->uniform   2->differential  3->bizonal  4 -> bizonal with smooth gradient   5 -> latitudinal diff: Om = Om_eq - DOm cos^2 theta;\n")
            f.write("   %i;\n" % (self.fparams.acor_rot_prof_))
            f.write("# Rotation profile given in the model of structure or in an external file? 1 -> int  2 -> ext;\n")
            f.write("   %i;\n" % (self.fparams.acor_profile_ext+1))
            f.write("# Name of the file? (format: fractionnal radius, Omega(rad/s), no header);\n")
            if self.fparams.acor_profile_ext:
                f.write("bla;\n")
            else:
                f.write(self.fparams.profile_filename+";\n")
            f.write("# Rotation suffix for the files name?;\n")
            f.write(self.fparams.mod_suff+";\n")
            f.write("# If uniform, value of Om (muHz)? if differential, value of the central Om (muHz)?;\n")
            if self.fparams.acor_rot_prof_ == 1:
                f.write("   %f;\n" % self.fparams.unif_om)
            else:
                f.write("   %f;\n" % self.fparams.diff_om)
            if self.fparams.acor_rot_prof_ == 3 or self.fparams.acor_rot_prof_ == 4:
                f.write("# If bizonal, surface rotation rate? (rad/s) If latitudinal, equatorial value (muHz)?;\n")
                f.write("   %f;\n" % self.fparams.biz_oms)
                f.write("# If bizonal, central rotation rate? (rad/s) If latitudinal, Dom value (muHz)?;\n")
                f.write("   %f;\n;\n" % self.fparams.biz_omc)
            else:
                f.write("# If bizonal, surface rotation rate? (rad/s) If latitudinal, equatorial value (muHz)?;\n")
                f.write("   %f;\n" % self.fparams.lat_omeq)
                f.write("# If bizonal, central rotation rate? (rad/s) If latitudinal, Dom value (muHz)?;\n")
                f.write("   %f;\n;\n" % self.fparams.lat_dom)
            f.write("### frequency target;\n")
            f.write("# Choice for the frequency targets: 1->scanning 2->file 3->delta Pi  4->Delta Nu ;\n")
            f.write("   %i;\n" % (self.fparams.acor_freq_target_))
            f.write("# if scanning answer the four next questions;\n")
            f.write("# input frequency unit: 1->microHz, 2->Omk;\n")
            f.write("   %i;\n" % (self.fparams.acor_freq_units_))
            f.write("# starting frequency:;\n")
            f.write("   %f;\n" % (self.fparams.acor_freq_start))
            f.write("# step in frequency: ;\n")
            f.write("   %f;\n" % (self.fparams.acor_freq_step))
            f.write("# number of modes, or number of large separations:;\n")
            f.write("   %i;\n" % (self.fparams.acor_nsteps))
            f.write("# if frequencies in a file: kind of file?;\n")
            f.write("# 1->ADIPLS  2->ACOR  3->Marie-Jo  4->SebS?;\n")
            f.write("   %i;\n" % (self.fparams.acor_freq_file_type_))
            f.write("# name of the file?;\n")
            f.write(self.fparams.acor_freq_filename + ";\n")
            f.write(";\n")
            f.write("### Mode caracteristics:  ;\n")
            f.write("# azimutal order m: ;\n")
            f.write("   %i;\n" % (self.fparams.acor_m))
            f.write("# parity: 0 -> S. H. of l same parity as m 1 -> S. H. of opposite parity as m;\n")
            f.write("   %i;\n" % (not self.fparams.acor_parity))
            f.write(";\n")
            f.write("### computations parameters:;\n")
            f.write("# number of spherical harmonics:;\n")
            f.write("   %i;\n" % (self.fparams.acor_nsh))
            f.write("# radial resolution: 1->the model one 2->half of it, etc.. ;\n")
            f.write("   %i;\n" % (self.fparams.acor_rad_res))
            f.write("# number of Newton iterations:;\n")
            f.write("   %i;\n" % (self.fparams.acor_newt_it))


    def run( self ):
        """
        Run `script_run_ACOR_unified` executable. `CRunAcor.run` is only executed by the Thread class
        """
        # If a list_frequences file already exists with same name, we remove it
        os.system(f"rm -f liste_frequences_{self.fname} acor_freq_{self.fname}")

        command = f'{self.fparams.acor_exe_path}/script_run_ACOR_unified {self.aparam_name}'
        if self.quiet:
            res = subp.Popen( command.split( ), stdin=subp.PIPE, stdout=subp.PIPE,
                shell=False ).communicate()[0]
        else:
            res = subp.Popen( command.split( ) ).communicate()[0]

        self.clean( )
        self.read_acor_liste_freq( )

    def clean( self ):
        """
        Clean directory of not needed files.
        """
        # In the future, ACOR_input should not be removed here.
        if self.fparams.acor_rm_input:
            os.system( "rm -f ACOR_input" )
        os.system( "rm -f Propagation_ACOR Propagation_cavities" )
        os.system( "rm -f fort.*" )

    def read_acor_liste_freq( self ):
        """
        Reads a `acor_liste_freq` file.
        """
        suff_w = '' if self.fparams.acor_rot_toro else '_notor'

        if self.fparams.acor_rot:
            if self.fparams.acor_rot_prof_ == 1:
                suff_rot ='unif*muHz'+suff_w
            elif self.fparams.acor_rot_prof_ == 2:
                suff_rot = '*'
                suff_rot += '' if self.fparams.acor_rot_grad else '_nograd'
                suff_rot += suff_w
            elif self.fparams.acor_rot_prof_ == 3:
                suff_rot ='bizone_CoS*'
                suff_rot += '_S*'
            elif self.fparams.acor_rot_prof_ == 4:
                suff_rot = 'bizone_CoS*'
                suff_rot += '_C*muHz'
            elif self.fparams.acor_rot_prof_ == 5 or self.fparams.acor_rot_prof_ == 6:
                suff_rot = self.fparams.mod_suff.strip()
        else:
            suff_rot='_norot'

        files = glob.glob( f"liste_frequences_{self.modname}_{suff_rot}{self.Mmpar_name}" )#lfreq_name+suff_rot )
        lfreq_name = files[0]

        data = np.loadtxt( lfreq_name, skiprows=1)
        if data.shape[1] == 21:
            self.cfacor.number, self.cfacor.m, self.cfacor.l, self.cfacor.nt, self.cfacor.npt, \
                self.cfacor.ngt, self.cfacor.n, self.cfacor.np, \
                self.cfacor.ng = np.transpose( data[:,0:9].astype(int) )
            self.cfacor.Sigma_in, self.cfacor.Sigma_out, self.cfacor.Sigma_out2, \
                self.cfacor.freq, self.cfacor.ekin, self.cfacor.I, self.cfacor.Knl, \
                self.cfacor.Beta, self.cfacor.I_core, self.cfacor.Knl_core, \
                self.cfacor.Beta_core, self.cfacor.Omk = np.transpose( data[:,9:21] )
        else:
            self.cfacor.number, self.cfacor.m, self.cfacor.l, self.cfacor.n, self.cfacor.np, self.cfacor.ng = np.transpose( data[:,0:6].astype(int) )

            self.cfacor.Sigma_in, self.cfacor.Sigma_out, self.cfacor.Sigma_out2, \
                self.cfacor.freq, self.cfacor.ekin, self.cfacor.I, self.cfacor.Knl, \
                self.cfacor.Beta, self.cfacor.I_core, self.cfacor.Knl_core, \
                self.cfacor.Beta_core, self.cfacor.Omk = np.transpose( data[:,6:18] )
