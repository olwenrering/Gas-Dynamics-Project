#!/usr/bin/env python3
# coding: utf-8

# Code d'Evolution Stellaire, Adaptatif et Modulaire for the 2020 decade.
# Copyright (c) 1997-2023 The Cesam2k20 authors
# SPDX-License-Identifier : GPL-3.0-or-later


# try to import the new traits library
try:
    from traits.api import HasTraits, Trait, Float, File, Bool
# look for the old traits library if the new one is not found
except ImportError:
    from enthought.traits.api import HasTraits, Trait, Float, File, Bool
import os

from pycesam.constants import RUN_DICT
from pycesam.tools import isbool, str2bool

NO_COLOR  = "\033[0m"
DGRAY     = "\033[30;01m"
RED       = "\033[31;01m"
GREEN     = "\033[32;01m"
YELLOW    = "\033[33;01m"
BLUE      = "\033[34;01m"
MAGENTA   = "\033[35;01m"
CYAN      = "\033[36;01m"
LGRAY     = "\033[37;01m"
DGRAY_U   = "\033[40;01m"
RED_U     = "\033[41;01m"
GREEN_U   = "\033[42;01m"
YELLOW_U  = "\033[43;01m"
BLUE_U    = "\033[44;01m"
MAGENTA_U = "\033[45;01m"
CYAN_U    = "\033[46;01m"
LGRAY_U   = "\033[47;01m"


class CRun(HasTraits):

    _dict_val = RUN_DICT
    job       = Trait( 'From PMS', _dict_val['job'], label='Evolution' )

    c_iben    = Float( 5.0e-6, label="Iben's constant" )
    type_file = Trait( 'ASCII', 'Binary' )

    mod_init  = File( os.environ['CESDIR'] + '/SUN_STAR_DATA/INIT/pms/8d-5.pms' )
    dt0       = Float( 0.0, label='Initial time step (Myrs)' )
    pre_pms   = Bool( False, label='Precompute PMS model' )
    nohup     = Bool( False, label="Detach from shell" )
    li_zams   = Bool( True, label='Force the computation with a low Li when starting from ZAMS' )


    def __init__( self, name, job=None, verbose=True ):
        """
        <CRun.__init__( name, job=None, verbose=True )>
        Initialise options to control Cesam2k20's run.

        :param name: Name of the model.
        :type name: str

        :kparam job: OPtion to set the starting point from elsewhere than the .run file.
        :ktype job: dict

        :kparam verbose: If True, prints additional information.
        :ktype verbose: dict
        """
        self.data_dir = os.environ['CESDIR'] + '/SUN_STAR_DATA/'
        self.verbose  = verbose
        self.name     = name
        self.run_file = name + '.run'
        self.pms      = name + '_B.pms'
        self.zams     = name + '_B.zams'
        self.rep      = name + '_B.rep'

        # True if mod_init has been changed by hand and must not be touched again
        self.mod_init_changed = False

        if os.path.exists(self.run_file): # read .run file
            self.read_rsettings()

        if job is not None:
            self.job     = job['job']
            self.manage_type_file_changes( )
            if 'mod_init' in job.keys() and job['mod_init']:
                self.mod_init = job['mod_init']
            if 'type_file' in job.keys():
                self.type_file = job['type_file']
            if 'c_iben' in job.keys():
                self.c_iben = job['c_iben']
            if 'dt0' in job.keys():
                self.dt0 = job['dt0']

        self.on_trait_change(self.manage_type_file_changes, 'job_, type_file')
        self.on_trait_change(self.manage_mod_init_changes, 'mod_init')

    def manage_mod_init_changes( self ):
        """
        <CRun.manage_mod_init_changes( )>
        Automatically called by a decorator of `traits` module, any time `mod_init` is changed.
        """
        self.mod_init_changed = True

    def manage_type_file_changes( self ):
        """
        <CRun.manage_type_file_changes( )>
        Automatically called by a decorator of `traits` module, any time `type_file` is changed.
        """
        curr_dir = os.getcwd()
        if self.job_ == 'rep':
            self.type_file = 'Binary'
            if not self.mod_init_changed and not self.mod_init.endswith(('.rep', '.reph5')):
                self.mod_init = self.rep
        elif self.job_ == 'pms':
            self.fltr = ['*.pms']
            if not self.mod_init_changed:
                if self.type_file == 'Binary':
                    self.mod_init =  self.pms
                else:
                    self.mod_init = self.data_dir + '/INIT/pms/8d-5.pms'
        elif self.job_ == 'zams' and not self.mod_init_changed:
            if self.type_file == 'Binary':
                self.mod_init = self.zams
            else:
                self.mod_init = self.data_dir + '/INIT/zams/m010.zams'

    def write_rsettings( self ):
        """
        <CRun.write_rsettings( )>
        Write run settings to a .run file.
        """
        vars = ['job', 'c_iben', 'type_file', 'mod_init', 'dt0', 'pre_pms', 'li_zams']

        lines = []
        for var in vars:
            lines.append('{:8s}     {:16s}\n'.format(var, str(self.__getattribute__(var))))

        with open(self.run_file, 'w') as f:
            f.writelines(lines)

    def copy( self, other ):
        """
        <CRun.copy( other )>
        Copy run settings to another CRun instance, possibly from another model.

        :param other: Other instance of CRun, which have same parameters as `self`.
        :type other: <pycesam.CRun>
        """
        self.job       = other.job
        self.c_iben    = other.c_iben
        self.type_file = other.type_file
        self.mod_init  = other.mod_init
        self.dt0       = other.dt0

    def read_rsettings(self):
        """
        <CRun.read_rsettings( )>
        Read run settings in a .run file.
        """
        if self.verbose: print(f"Reading {self.run_file}...", end='')

        with open(self.run_file, 'r') as f:
            lines = f.readlines()

        for line in lines:
            words  = line.split()
            key    = words[0]
            isplit = len(key) + 1
            val    = line[isplit:].strip()
            if key in ['c_iben', 'dt0'] : val = float(val)
            if key in ['pre_pms', 'li_zams']:
                if val in ['True', 'False']: val = eval(val)

            if isbool(val):
                self.__setattr__(key, str2bool(val) )
            else:
                self.__setattr__(key, val)

        if self.verbose: print(f"{GREEN}[Done]{NO_COLOR}")


