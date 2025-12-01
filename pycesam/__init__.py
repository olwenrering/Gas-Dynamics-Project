#!/usr/bin/env python3
# coding: utf-8

# Code d'Evolution Stellaire, Adaptatif et Modulaire for the 2020 decade.
# Copyright (c) 1997-2023 The Cesam2k20 authors
# SPDX-License-Identifier : GPL-3.0-or-later


import os, glob
import subprocess as subp
import multiprocessing as mulp
import numpy as np
import re
import time
import h5py
import sys
import zipfile
from contextlib import nullcontext
from copy import deepcopy
from scipy import interpolate
from tqdm.auto import tqdm
# from tqdm.notebook import tqdm
# from tqdm.rich import tqdm
from time import sleep

from traits.api import HasTraits, Instance

from pycesam.tools import *
from pycesam.constants import *
from pycesam.cparams import CParameters
from pycesam.CConstants import CConstants
from pycesam.cesam_run import CRun
from pycesam.freqs import Freqs, CFacor, CRunAcor
from pycesam.FortranIO import *

try:
    sys.path.insert(0, os.environ['CESDIR']+'/python3')
    from cesam2k20_fortran import cesam2k20_fortran
    cesfort_available = True
except:
    cesfort_available = False

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



class CModel(HasTraits):
    """
    Class that represents a Cesam2k20 evolution model.
    """
    params  = Instance(CParameters)
    run     = Instance(CRun)
    fparams = Instance(Freqs)

    def __init__( self, name, reinit=False, read=True, job=None, fromGUI=False, verbose=True ):
        """
        <CModel.__init__(name, reinit=False, read=True, job=None)>

        Class that represents a Cesam2k20 evolution model.

        :param name: Name of the model.
        :type name: string

        :kparam reinit: If True, deletes all calculated files. Default: False.
        :ktype reinit: boolean

        :kparam read: If True, reads files .HR and .osc. Default: True.
        :ktype read: boolean

        :kparam job: Some options we want to enforce.
        :ktype job: dict

        :kparam fromGUI: True if called from `super()` function in CModelGUI. Default: False
        :ktype fromGUI: boolean

        :kparam verbose: If True, information is written for the user. Default: True.
        :ktype verbose: boolean

        :raises NameError: if name length is $0$ or $> 31$ characters.
        :raises CESAMError: If name of model starts with a hyphen `'-'`.
        """

        if not name:
            raise CESAMError( 'Name of the model must not be an empty string.' )
        elif name[0] == '-':
            raise CESAMError("pycesam will never allow you to create a model with name beginning with '-'."
                "By all means necessary. Never. This is for you own safety.")
        elif len(name) > 49:
            raise CESAMError('Name too long! Must be < 31 characters.')
        elif len(name) == 0:
            raise CESAMError('Name must not be empty.')

        self.name      = name
        self.don       = name + '.don'
        self.hr        = name + '.HR'
        self.hrnew     = name + '.HRnew'
        self.lis       = name + '.lis'
        self.file_atm  = name + '_B.atm'
        self.file_pms  = name + '_B.pms'
        self.file_zams = name + '_B.zams'
        self.file_dat  = name + '_B.dat'
        self.err_file  = name + '.err'
        self.job       = job
        self.read      = read
        self.reinit    = reinit
        self.verbose   = verbose

        self.HRoffset  = 0
        self.nmod      = 0

        if not fromGUI: self.__init( )


    def __str__(self):
        """
        Gives the string representation of a CModel instance.

        :return: Mass, X0, Y0 and Z0 of the model if has been correctly calculated. Otherwise, only mass.
        :rtype:  string
        """
        if self.finished:
            return ("Cesam2k20 model \'%s\'\n\tmodel calculated\n\tParameters:\n" + \
                   "\t\tMass = %.3f\n\t\tX0 =  %.5f\n\t\tY0 =  %.5f\n\t\tZ0 =  %.5f") % \
                   (self.name, self.params.mtot, self.params.x0, self.params.y0, self.params.z0)
        else:
            return "Cesam2k20 model \'%s\'\n\tmodel not calculated\n\tMass = %.3f" % (self.name, self.params.mtot)


    def __repr__(self):
        """
        Gives the representation of a CModel instance (see `self.__str__`).

        :return: Mass, X0, Y0 and Z0 of the model if has been correctly calculated. Otherwise, only mass.
        :rtype:  string
        """
        message = self.__str__()

        return message


    def __getitem__(self, key):
        """
        Return value of a given attributename.

        :raises KeyError: if attribute does not exist in class.
        """
        if key in self.params.class_trait_names():
            return self.params.__getattribute__(key)
        elif key in self.fparams.class_trait_names():
            return self.fparams.__getattribute__(key)
        else:
            raise KeyError


    def __setitem__(self, key, value):
        """
        Set value to a given attributename.

        :raises KeyError: if attribute does not exist in class.
        """
        if key in self.params.class_trait_names():
            self.params.__setattr__(key, value)
        elif key in self.fparams.class_trait_names():
            self.fparams.__setattr__(key, value)
        else:
            raise KeyError

    def __eq__(self, other):
        """
        Checks whether model parameters are equal to parameters of another model.

        :param other: Instance of class CModel.
        :type other: <CModel>
        """
        return self.params.__getstate__() == other.params.__getstate__()


    def __len__(self):
        """
        Returns number of time steps in the model
        """
        return self.nmod


    def __ne__(self, other):
        """
        Checks whether model parameters are not equal to parameters of another model.

        :param other: Instance of class CModel.
        :type other: <CModel>
        """
        return not self.__eq__(other)

    def __call__( self, mkdon=True, mkeos=True, debug=False, log=False, devnull=False, verbose=None,
        progress_bar=True, executable=None, **kwargs ):
        """
        <CModel.__call__( mkdon=True, debug=False, log=False, verbose=True )>

        Calculates models or frequencies.

        :kparam mkdon: Creates (or re-creates) a .don file if `mkdon==True`; uses existing
            .don file otherwise. Default: True.
        :ktype mkdon: boolean

        :kparam mkeos: Creates (or re-creates) the EoS file if `mkeos==True`; uses existing
            EoS file otherwise. Default: True.
        :ktype mkeos: boolean

        :kparam debug: If True, run the debug Cesam2k20 executable (`cesam2k20_dbg.x`). Default: False.
        :ktype debug: boolean

        :kparam log: If True, direct standard output to .log file. Default: False.
        :ktype log: boolean

        :kparam devnull: If True, redirect standard output to /dev/null.
        :ktype devnull: boolean

        :kparam verbose: If True, information is written for the user. Default: True.
        :ktype verbose: boolean

        :kwarg valgrind: If True, calls Cesam2k20 with valgrind analysis software.
        :kwtype valgrind: boolean

        :kwarg coverage: If True, evaluate test coverage by calling `cesam2k20_cov.x` with gcov.
        :kwtype coverage: boolean
        """


        verbose = self.verbose if verbose is None else verbose

        self.fparams.write_fsettings( )
        if self.run.job_ == 'freqs':
            self.calc_freqs()

        else:
            self.finished = False
            self.params.mkdon( overwrite=mkdon, replace_eos=mkeos )
            if self.run.job_ != 'rep':
                if os.path.exists(self.hr):
                    os.remove(self.hr)
                if os.path.exists(self.hrnew):
                    os.remove(self.hrnew)

            self.run.write_rsettings()

            self.result = None

            valgrind = kwargs.get( 'valgrind', False )
            coverage = kwargs.get( 'coverage', False )

            if self.run.pre_pms: self.__precompute_pms(valgrind=valgrind, coverage=coverage, debug=debug, log=log, devnull=devnull )

            if verbose:
                print("\n-----------------------------------------")
                print(f"Model: {self.name}")
                print(f"    Started at {time.asctime()} ")
                print(f"Computing evolution of model {self.name}..." )
            t1 = time.time()

            # runs Cesam2k20
            self.parent, self.child = mulp.Pipe()
            self.process = mulp.Process(
                target=self.run_cesam,
                kwargs={'pipe':self.child, 'debug':debug, 'log':log, 'valgrind':valgrind, 'coverage':coverage, 'devnull':devnull, 'executable':executable} )
            self.process.start()

            progress_bar = progress_bar and (not self.__is_backgrounded())
            if progress_bar: self.__display_progress_bar( )

            sortie, err = self.parent.recv( )
            self.__process_out( sortie, err, t1, verbose, debug, log, valgrind )

            #if coverage:
            #    cwd = os.getcwd()
            #    cesdir = os.environ['CESDIR']
            #    os.system( f'mv {cesdir}/src/*gcda {cwd}')
            #    os.system( f'mv {cesdir}/src/*gcno {cwd}')
            #    s = subp.Popen(['gcov', '-b', '-c', 'cesam2k20_cov.x'], stdin=subp.PIPE, stdout=subp.PIPE)

    def __display_progress_bar( self ):
        bar_format  = '{desc:<40}: {percentage:3.0f}%|{bar:20}| {n:5.2e}/{total:5.2e} {unit:>10}'
        barx_format = '{desc:<40}: {percentage:3.0f}%|{bar:20}| {n:.2f}/{total:.2f}   {unit:>16}'
        positions = iter(range(4,0,-1))
        with tqdm(desc="Progress on age", unit='Myrs', total=self.params.agemax, colour=u'#ffbf00', 
                bar_format=bar_format, leave=False, position=next(positions) ) as pbar_age, \
            tqdm(desc="Progress on central density", unit='g/cm3', total=self.params.rhoc, colour=u'#006a4e', 
                bar_format=bar_format, leave=False, position=next(positions) ) if self.params.rhoc > -1.0 else nullcontext() as pbar_rhoc, \
            tqdm(desc="Progress on central temperature", unit='K', total=self.params.t_stop, colour=u'#8a2be2', 
                bar_format=bar_format, leave=False, position=next(positions) ) if self.params.t_stop > -1.0 else nullcontext() as pbar_tc, \
            tqdm(desc=f"Decrease of H1 mass fraction (Xc0={self.params.x0:.2f})", unit='mass fraction', total=self.params.x0-self.params.x_stop, colour=u'#00cc99', 
                bar_format=barx_format, leave=False, position=next(positions) ) if self.params.x_stop > -1.0 else nullcontext() as pbar_xc:

            size = 0
            while self.process.is_alive():
                hr = self.hr if os.path.exists( self.hr ) else self.hrnew
                size2 = os.path.getsize( hr ) if os.path.exists( hr ) else 0
                if size != size2:
                    self.read_hr( update=True, verbose=False )
                    if self.age.size > 0:
                        pbar_age.update(self.age[-1]-pbar_age.n)
                        # pbar_nmod.update(self.nmod-pbar_nmod.n)
                        if pbar_rhoc is not None:
                            pbar_rhoc.update(self.rhoc[-1]-pbar_rhoc.n)
                            pbar_rhoc.refresh()
                        if pbar_tc is not None:
                            pbar_tc.update(self.Tc[-1]-pbar_tc.n)
                            pbar_tc.refresh()
                        if pbar_xc is not None:
                            pbar_xc.n = self.params.x0-self.ab_c['H1'][-1]
                            pbar_xc.refresh()
                    size = size2
                sleep(1)
        return

    def __is_backgrounded( self ):
        try:
            # Get foreground process group of the terminal
            fg_pgrp = os.tcgetpgrp(sys.stdout.fileno())
            # Compare to our process group
            our_pgrp = os.getpgrp()
            return fg_pgrp != our_pgrp
        except OSError:
            # No controlling terminal (e.g., nohup, systemd)
            return False

    def __precompute_pms(self, valgrind=False, coverage=False, debug=False, log=False, devnull=False):
        m_max = 4.0
        m_min = 0.5

        if self.params.mtot > m_max or self.params.mtot < m_min:
            if self.params.mtot > m_max:
                step = 0.5
                m_init = m_max
            else:
                step = 0.1
                m_init = m_min

            m_goal = self.params.mtot
            agemax = self.params.agemax

            self.params.mtot   = m_init
            self.params.agemax = 0.0
            self.params.mkdon()

            print(f'Calculating PMS model for mass = {self.params.mtot}')
            sortie = self.run_cesam( debug=debug, log=log, valgrind=valgrind, coverage=coverage, devnull=devnull )

            nsteps = np.abs(m_init - m_goal)/step + 1

            nsteps = int(round(nsteps))
            print(f'nsteps = {nsteps}')

            self.run.type_file = 'Binary'
            self.run.mod_init  = self.file_pms
            for i in range(1,nsteps-1):
                mass = m_init - i/(nsteps - 1)*(m_init - m_goal)
                self.params.mtot = mass

                self.params.mkdon()
                print(f'Calculating PMS model for mass = {self.params.mtot}')
                sortie = self.run_cesam( debug=debug, log=log, valgrind=valgrind, coverage=coverage, devnull=devnull )

            print(f"Done. Now compute model as usual.")

            self.run.type_file = 'Binary'
            self.params.mtot = m_goal
            self.params.agemax = agemax
            self.params.mkdon()

        else:
            print(f'0.5 < self.params.mtot = {self.params.mtot} < 5.0: No need to step in mass.')


    def __init( self ):
        """
        <CModel.__init( )>

        Method used to initialize all subclasses used by CModel and read files if neaded.
        """
        self.params   = CParameters( self.name, verbose=self.verbose )
        self.run      = CRun( self.name, job=self.job, verbose=self.verbose )
        self.ctes     = CConstants( )


        self.set_osc()
        self.set_osc2d()
        self.set_rep()
        self.init_fparams('p')
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



    def __process_out( self, out, err, tstart, verbose, debug, log, valgrind ):
        """
        <CModel.__process_out( out, err, tstart, verbose, debug, log, valgrind )>

        Method that process the output of Cesam2k20 (Was the evolution a success or was it aborted?),
        and read the output files if needed

        :param out: Contains all messages printed to standard output by Cesam2k20.
        :type out: array of str

        :param err: Contains all messages printed to standard error by Cesam2k20.
        :type err: array of str

        :param tsart: Time at which the computation started.
        :type tsart: float

        :param verbose: If True, information is written for the user. Default: True.
        :type verbose: boolean

        :param debug: If True, run the debug Cesam2k20 executable (`cesam2k20_dbg.x`). Default: False.
        :type debug: boolean

        :param log: If True, direct standard output to .log file. Default: False.
        :type log: boolean

        :param valgrind: If True, calls Cesam2k20 with valgrind analysis software.
        :type valgrind: boolean
        """
        err = err.split('\n')

        for line in err:
            if 'Succes' in line:
                self.result   = line
                self.error    = None
                self.finished = True

                t2 = time.time()
                hr, mn, sc = time.gmtime(t2 - tstart)[3:6]
                if verbose:
                    print(f"\n{GREEN}[Done]{NO_COLOR}")
                    print(f"   Finished at {time.asctime()}")
                    if hr == 0:
                        if mn == 0:
                            print(f"   Time: {sc:02.2f}s")
                        else:
                            print(f"   Time: {mn:02d}m {sc:02.2f}s")
                    else:
                        print(f"   Time: {hr:d}h {mn:02d}m {sc:02.2f}s")
                    print(f"{GREEN_U}---> {self.result}{NO_COLOR}")

                # If 'deformed' was in success message, then, we don't need to try anything
                if 'deformed' in line: return

                self.read_hr( verbose=verbose )
                if 'manu' not in line:
                    self.set_osc()
                    self.set_osc2d()
                    self.set_rep()
                    if self.params.nom_output_.startswith('osc_') or self.params.nom_output_.startswith('all_'):
                        self.read_osc( read_vc=False, verbose=verbose )

                if verbose and not self.hr_empty:
                    print(f"log Teff = {self.log_teff[-1]} ; Teff = {10.0**self.log_teff[-1]} K")
                    print(f"log L/L0 = {self.log_l[-1]} ; L/L0 = {10.0**self.log_l[-1]}")
                    print(f"log R/R0 = {self.log_r[-1]} ; R/R0 = {10.0**self.log_r[-1]}")
                    print(f"log g = {self.log_g[-1]}")
                    print(f"age = {self.age[-1]} Myrs")
                    if ('manu' not in line):
                        if self.params.nom_output_.startswith('osc_') and not self.params.nom_output_.endswith('h5'):
                            print(f"Tc   = {self.var[2][-1]:e} K")
                            print(f"rhoc = {self.var[4][-1]:e} g/cm^3")
                        if self.params.nom_output_.startswith('all_'):
                            print(f"Tc   = {self.var[-1][2][-1]:e} K")
                            print(f"rhoc = {self.var[-1][4][-1]:e} g/cm^3")
                elif self.hr_empty:
                    print("It appears that the computation of the model finished correctly\n",
                        "but no HR or HRnew file was generated.")
                return

        # only goes here if the calculation was interrupted tries to find and print error messages
        if out is None:
            self.error = ["NO ERROR MESSAGE: STDOUT WAS REDIRECTED TO /dev/null"]
        else:
            out = out.split('\n')
            if(len(out) >= 10):
                self.error = out[-10:]
            else:
                self.error = out


        if verbose:
            print(f"{RED}[Error]{NO_COLOR}")
            print(f"{RED_U}Not finished, error message: {NO_COLOR}")
            print(f"{RED}#########################################{NO_COLOR}")
            for line in err:
                print(line)
            for line in self.error:
                print(line)

            print(f"{RED}#########################################{NO_COLOR}")

            self.read_hr( verbose=verbose )
            # if not self.err_hr and not self.hr_exists:
            if self.hr_exists:
                print(f"  age = {self.age[-1]} Myrs")
            else:
                print(f"{RED}Didn\'t even start...{NO_COLOR}")

            if len( self.error ):
                for i in range( min( len(self.error), 10 ) ):
                    self.error[i] += '\n'

        with open(self.err_file, 'w') as f:
            f.writelines(self.error)

        #print( "\nYou found a bug! Send the .don file to troubleshootingcesam2k20@disroot.org ?" )
        #go = nonBlockingInput( "[N]/y (10 secs to choose): ", 10, default='N' )
#
        #if go.lower() == 'y':
        #    email_adress = input("Enter your email adress:" )
#
        #    text = 'Cesam2k20 bug report\n'
        #    self.read_hr( verbose=verbose )
        #    if not self.err_hr:
        #        text += f"  age = {self.age[-1]} Myrs\n"
        #    else:
        #        text += f"Didn\'t even start\n"
        #    text += f"Commit hash: {get_hash()}"
        #    try:
        #        email( self, email_adress, text )
        #    except ConnectionRefusedError:
        #        print( "The email could not be sent because of a ConnectionRefusedError.",
        #            "You can still send your .don file, the error message, and ideally the",
        #            f"commit hash ({get_hash()}) of your version of Cesam2k20 to the adress",
        #            "troubleshootingcesam2k20@disroot.org" )


        self.set_osc()
        self.set_osc2d()
        self.set_rep()

    def __tail( self, file, n=1 ):
        lines = ''
        with open( file, 'rb' ) as f:
            try:  # catch OSError in case of a one line file
                f.seek( -n, os.SEEK_END )
                while f.read(1) != b'\n':
                    f.seek( -n, os.SEEK_CUR )
            except OSError:
                f.seek(0)
            for i in range( n ):
                lines += f.readline().decode('latin-1')

        return lines


    def run_cesam( self, pipe=None, debug=False, log=False, devnull=False, executable=None, **kwargs ):
        """
        <CModel.run_cesam( debug=False, log=False )>

        Private method that runs one of the Cesam2k20 executable. If needed, it stores the
        standard output to file.

        :kparam debug: If True, run in debug mode (runs `'cesam2k20_dbg.x'` instead of `'cesam2k20.x'`)
        :ktype debug: boolean

        :kparam log: If True, stores standard output to file.
        :ktype log: boolean

        :kparam devnull: If True, redirect standard output to /dev/null.
        :ktype devnull: boolean

        :kwargs valgrind: If True, calls Cesam2k20 with valgrind analysis software.
        :kwtype valgrind: boolean

        :kwargs coverage: If True, evaluate test coverage by calling `cesam2k20_cov.x` with gcov.
        :kwtype coverage: boolean
        """

        bbin = self.run.type_file == 'Binary'
        if self.run.job_ == 'pms':
            cmd = ["3", f"{'o' if bbin else 'n'}", f"{self.run.mod_init}", f"{self.name}", f"{self.run.c_iben}"]
        elif self.run.job_ == 'rep':
            cmd = ["1", "o", f"{self.run.mod_init}", f"{self.name}", f"{'o' if bbin else 'n'}", f"{self.run.dt0}"]
        elif self.run.job_ == 'zams':
            cmd = ["2", f"{'o' if bbin else 'n'}", f"{self.run.mod_init}", f"{self.name}", f"{'o' if self.run.li_zams else 'n'}"]
        else:
            print(f'{RED}Error:{NO_COLOR} must start from PMS, ZAMS or previous model')
            return

        if debug:
            stdout = open( self.name + '.log', 'w+', encoding='latin-1' )
            s = subp.Popen( ['cesam2k20_dbg.x' if executable is None else executable, *cmd], stdout=stdout, stderr=subp.PIPE )
        elif kwargs.get( 'valgrind', False ):
            command = f'valgrind -v --leak-check=full --show-leak-kinds=all --track-origins=yes \
                --read-inline-info=yes --xml=yes --xml-file={self.name}-valgrind.xml cesam2k20_dbg.x {cmd}'
            stdout = open( self.name + '.profile_log', 'w+', encoding='latin-1' )
            s = subp.Popen( [*command.split(), *cmd], stdout=stdout, stderr=subp.PIPE )
        elif kwargs.get( 'coverage', False ):
            stdout = open( self.name + '.cov_log', 'w+', encoding='latin-1' )
            s = subp.Popen( ['cesam2k20_cov.x' if executable is None else executable, *cmd], stdout=stdout, stderr=subp.PIPE )
        else:
            if log:
                stdout = open( self.name + '.log', 'w+', encoding='latin-1' )
            elif devnull:
                stdout = subp.DEVNULL
            else:
                stdout = subp.PIPE
            s = subp.Popen( ['cesam2k20.x' if executable is None else executable, *cmd], stdout=stdout, stderr=subp.PIPE )
        sortie, err = s.communicate( )#input=cmd.encode() )

        if debug or log:
            sortie = self.__tail( self.name + '.log', n=10 )
        elif kwargs.get( 'valgrind', False ):
            sortie = self.__tail( self.name + '.profile_log', n=10 )
        elif kwargs.get( 'coverage', False ):
            sortie = self.__tail( self.name + '.cov_log', n=10 )
        else:
            sortie = sortie if sortie is None else sortie.decode('latin-1')
            if sortie is not None:
                sortie = sortie.split('\n')
                if(len(sortie) >= 10):
                    sortie = sortie[-10:]
                sortie = '\n'.join( sortie )

        if pipe is None:
            return sortie, err.decode('latin-1')
        else:
            pipe.send( (sortie, err.decode('latin-1')) )
            pipe.close()


    def is_finished( self ):
        """
        Checks if the computation of the evolutionary track has been carried to the end.
        Modifies boolean attribute self.finished.
        """
        if self.all_osc:
            if self.osc:
                self.finished = len(self.osc) > 0
            elif self.osch5:
                self.finished = len([self.osch5]) > 0
            elif self.osc2d:
                self.finished = len(self.osc2d) > 0
        elif self.osc or self.osch5  or self.osc2d:
            if self.osc:
                self.finished = os.path.exists(self.osc)
            elif self.osch5:
                self.finished = os.path.exists(self.osch5)
            elif self.osc2d:
                self.finished = os.path.exists(self.osc2d)
        else:
            self.finished = False


    def _reinit(self):
        """
        Removes all calculated files. Does not remove file .don.
        """
        files = [self.hr, self.hrnew, self.lis, self.file_atm, self.file_pms, self.file_zams]
        for i in files:
            if os.path.exists(i): os.remove(i)
        self.finished = False


    def set_osc(self):
        """
        Find the names of .osc files; stores names in `self.osc`.
        """
        self.osc     = None
        self.osch5   = None
        self.all_osc = False
        save_mods    = True

        if self.params.nom_output_.endswith('h5'):
            if '_adia' in self.params.nom_output_:
                self.suf = '-ad.osch5'
            elif '_nadia' in self.params.nom_output_:
                self.suf = '-nad.osch5'
            elif '_invers' in self.params.nom_output_:
                self.suf = '-inv.osch5'
            elif '_plato' in self.params.nom_output_:
                self.suf = '-plato.osch5'

            if ((self.params.nom_output_.startswith('all_')) and save_mods): # saves all models as .osch5
                self.all_osc = True
                if self.params.nb_max_modeles > 0: # saves all _B.rep
                    self.params.nb_max_modeles *= -1

            self.osch5 = self.name + self.suf

        else:
            if '_adia' in self.params.nom_output_:
                self.suf = '-ad.osc'
            elif '_nadia' in self.params.nom_output_:
                self.suf = '-nad.osc'
            elif '_invers' in self.params.nom_output_:
                self.suf = '-inv.osc'
            elif '_plato' in self.params.nom_output_:
                self.suf = '-plato.osc'

            if ((self.params.nom_output_ == 'no_output') or self.params.nom_output_.endswith('ascii')) :
                self.osc = []
            else:
                if ((self.params.nom_output_.startswith('all_')) and save_mods): # saves all models
                                                               # as .osc
                    self.all_osc = True
                    self.osc = [i for i in os.listdir('.') if i[6:] == self.name + self.suf
                                and i[:5].isdigit()]
                    self.osc.sort()

                    if self.params.nb_max_modeles > 0:
                        self.params.nb_max_modeles *= -1 # saves all _B.rep

                elif save_mods: # no model saved
                    self.osc = self.name + self.suf
                else:
                    self.osc = []

            self.nosc = len( self.osc )

    def set_osc2d(self):
        """
        Find the names of .osc2d files; stores names in `self.osc2d`.
        """
        self.osc2d = None
        self.all_osc2d = False
        save_mods = True

        suf = '.osc2d'

        if '2d' in self.params.nom_output_:
            if ((self.params.nom_output_.startswith('all_')) and save_mods): # saves all models
                                                           # as .osc
                self.all_osc2d = True
                self.osc2d = [i for i in os.listdir('.') if i[6:] == self.name + suf
                            and i[:5].isdigit()]
                self.osc2d.sort()

                if self.params.nb_max_modeles > 0:
                    self.params.nb_max_modeles *= -1 # saves all _B.rep

            elif ((self.params.nom_output_ != 'no_output') and save_mods): # no model saved
                self.osc2d = self.name + suf
            else:
                self.osc2d = []
        else:
            self.osc2d = []


    def set_rep(self):
        """
        Find the names of .rep files; stores names in `self.file_rep`.
        """
        self.file_rep = None
        self.file_dat = None
        if self.params.all_rep:
            self.file_rep = sorted( glob.glob( self.name + '*_B.rep' ) )
            self.file_dat = sorted( glob.glob( self.name + '*_B.dat' ) )
            if len(self.file_dat) == 1: self.file_dat = self.file_dat[0]

        else:
            self.file_rep = self.name + '_B.rep'
            self.file_dat = self.name + '_B.dat'

    def set_agsm( self ):
        """
        Find the names of .agsm files; stores names in `self.agsm`.
        """
        self.agsm = [i for i in os.listdir('.') if i[6:6+len(self.name)] == self.name and
                         i[:5].isdigit() and i[-5:] == '.agsm']
        self.agsm.sort()

    def set_acor( self ):
        """
        Find the names of .acor files; stores names in `self.acor`.
        """
        self.acor = [i for i in os.listdir('.') if i[6:6+len(self.name)] == self.name and
                         i[:5].isdigit() and i[-5:] == '.acor']
        self.acor.sort()

    def set_amdls( self ):
        """
        Find the names of .amdl files; stores names in self.amdls.
        """

        self.amdls = None
        save_mods = True

        if ((self.params.nom_output_.startswith('all_')) and save_mods): # saves all models
                                                       # as .osc
            self.amdls = [i for i in os.listdir('.') if i[6:] == self.name + '.amdl'
                        and i[:5].isdigit()]
            self.amdls.sort()

            if self.params.nb_max_modeles > 0:
                self.params.nb_max_modeles *= -1 # saves all _B.rep

        elif save_mods: # no model saved
            self.amdls = self.name + '.amdl'
        else:
            self.amdls = []


    def init_fparams(self, p_or_g, keep=False):
        """
        Initializes parameters for frequency calculations.

        :param p_or_g: Can be either `'p'` or `'g'`.
        :type p_or_g: character

        :kparam keep: If True : keeps default parameters.
        :ktype keep: bool
        """

        if not keep:
            self.fparams = Freqs( self.name, self.all_osc, verbose=self.verbose )
            self.fparams.modes = p_or_g + '-modes'

        # adipls
        self.amdl = self.name + '.amdl'
        self.agsm = []
        self.ssm  = []
        self.amde = self.name + '.amde'
        self.rkr  = self.name + '.rkr'
        self.gkr  = self.name + '.gkr'



    def isoc(self, age, x):
        """
        Returns a variable x interpolated in age.

        :param age: Desired age.
        :type age: float

        :param x: Table of values to interpolate. Must have the same size as self.age.
        :type x: array or list

        :return: Value of quantity x at desired age.
        :rtype: float

        :example: To get `log_teff` at t = 20.0 Myrs for model m:

        >>> log_teff = m.isoc(20.0, m.log_teff)
        """
        n = len(self.age)
        for i in range(n):
            if self.age[i] > age:
                xx = x[i-1] + (x[i] - x[i-1])/(self.age[i] - self.age[i-1])*\
                     (age - self.age[i-1])
                return xx

        print(f'age = {age} greater than agemax = {self.age[-1]}')


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def read_hr(self, obs=False, zsx_sol=None, verbose=None, update=False, process_cz=True):
        """
        :kparam verbose: If True, information is written for the user. Default: True.
        :ktype verbose: boolean
        """
        verbose = self.verbose if verbose is None else verbose
        self.err_hr = False
        if os.path.exists(self.hr):
            self.hr_exists = True
            self.__read_hr_old( obs=obs, zsx_sol=zsx_sol, verbose=verbose )
        elif os.path.exists(self.hrnew):
            self.hr_exists = True
            self.__read_hr_new( obs=obs, zsx_sol=zsx_sol, verbose=verbose, update=update, process_cz=process_cz)
        else:
            self.hr_empty = True
            self.hr_exists = False
            self.err_hr = True
            print(f"{RED}Error:{NO_COLOR} files {self.hr} or {self.hrnew} do not exist.")
            return

    def __read_hr_new(self, obs=False, zsx_sol=None, verbose=True, update=False, process_cz=True):
        """
        Reads the model's HR file. Computes colors and visual magnitudes if
        required. Teff-color transformations and bolometric corrections from
        VandenBerg & Clem, 2003, AJ, 126, 778 (uses the program described by
        VandenBerg et al., 2006, ApJS, 162, 375).

        :kparam obs: If True, computes colors and visual magnitudes. Default: False
        :ktype obs: bool

        :kparam zsx_sol: Solar Z/X. If None, will use solar Z/X correponding to model parameter 'Initial abundancies'.
            If unable, will use default value of 0.0181 (Asplund et al. 2009).
        :ktype zsx_sol: float

        :kparam verbose: If True, print some details.
        :ktype verbose: bool

        :member age: Age of each models.
        :mtype age: np.ndarray

        :member log_teff: $\log T_{\rm eff}$ of each models.
        :mtype log_teff: np.ndarray

        :member log_l: $\log T_{\rm eff}$ of each models.
        :mtype log_l: np.ndarray

        :member log_r: $\log T_{\rm eff}$ of each models.
        :mtype log_r: np.ndarray

        :member log_g: $\log T_{\rm eff}$ of each models.
        :mtype log_g: np.ndarray

        :member FeH: $[\rm Fe/H]$ of each models.
        :mtype FeH: np.ndarray

        :member z: Metallicity of each models.
        :mtype z: np.ndarray

        :member mstar: $M_\star$ of each models.
        :mtype mstar: np.ndarray

        :member omega_s: Surface angular velocity of each models.
        :mtype omega_s: np.ndarray

        :member vrot: Rotation velocity of each models.
        :mtype vrot: np.ndarray

        :member alpha: $\alpha$ parameter of each models.
        :mtype alpha: np.ndarray

        :member senv: Entropy of adiabat of each models.
        :mtype senv: np.ndarray

        :member mu: Mean molecular weight at bottom of CZ of each models.
        :mtype mu: np.ndarray
        """

        __dict_vars = {
            "                # model" : "imod",
            "             age (Myrs)" : "age",
            "                  phase" : "evol_phase",
            "               log Teff" : "log_teff",
            "             log L/Lsun" : "log_l",
            "             log R/Rsun" : "log_r",
            "           log Req/Rsun" : "log_req",
            "          log Rpol/Rsun" : "log_rpol",
            "             Mstar/Msun" : "mstar",
            "               G_eff,eq" : "geffeq",
            "              G_eff,pol" : "geffpol",
            "         Mdot (Msun/yr)" : "mdot",
            "           Mdot routine" : "mdot_prescription",
            "                 Tc (K)" : "Tc",
            "          rhoc (g/cm^3)" : "rhoc",
            "          Pc (dyn/cm^2)" : "Pc",
            "            L_pp (Lsun)" : "L_pp",
            "           L_cno (Lsun)" : "L_cno",
            "           L_3al (Lsun)" : "L_3al",
            "            L_gr (Lsun)" : "L_grav",
            "            L_nu (Lsun)" : "L_neut",
            "           neutrinos PP" : "nu_pp",
            "          neutrinos Be7" : "nu_be7",
            "           neutrinos B8" : "nu_b8",
            "          neutrinos N13" : "nu_n13",
            "          neutrinos O15" : "nu_015",
            "          neutrinos F17" : "nu_f17",
            "        tau_conv (days)" : "tau_conv",
            "   Omega_center (rad/s)" : "omega_c",
            "     Omega_surf (rad/s)" : "omega_s",
            "                  alpha" : "alpha",
            "    s_ad [1e-9 erg/g/K]" : "senv",
            "           enveloppe mu" : "mu"}

        str_vars = ["evol_phase", "mdot_prescription"]  # list of variables that are strings


        self.err_hr   = False
        self.hr_empty = os.path.getsize( self.hrnew ) < 1
        if self.hr_empty:
            self.err_hr = True
            print(f"{RED}Error:{NO_COLOR} file {self.hrnew} exists but is empty.")
            return

        if zsx_sol is None:
            try:
                zsx_sol = zsx_sun[self.params.nom_abon_]
            except KeyError:
                zsx_sol = 0.0181
                print(f'No solar metalicity found corresponding to {self.params.nom_abon}.')
                print(f'Using Z/X = {zsx_sol} to compute [Fe/H].')

        const1 = np.log10(self.ctes.gmsun/self.ctes.rsun**2)

        pwd = os.getcwd()
        vicdir = f"{os.environ.get('CESDIR')}/vicmodels"

        if verbose: print(f"Reading {self.hrnew}...", end='')

        read_from_start = not update or not self.nmod

        # Read header only once
        if read_from_start:
            head = []
            with open(self.hrnew, 'r', encoding='latin-1') as f:
                head.append(f.readline())
                head.append(f.readline())
                head_offset = f.seek( 0, 1 )
            # Byte offset used when reading HR file
            self.HRoffset  = head_offset

            self.nchim         = int(head[0][23:26])
            self.__hrnvars     = int(head[0][7:15])
            self.__hrnvars_tot = self.__hrnvars + 2*self.nchim
            self.nom_elem      = []
            col = 36
            for i in range(self.nchim):
                self.nom_elem.append(head[0][col:col+4].strip())
                col += 4

            col += 21
            try:
                self.n_contours_burn = int(head[0][col:col+4].strip())
                col += 22
                length = 23
                self.burn_contours = np.zeros(self.n_contours_burn)
                for i in range(self.n_contours_burn):
                    self.burn_contours[i] = float(head[0][col:col+length])
                    col += length
            except ValueError:
                print("\nHR file without burn zones.")
                self.n_contours_burn = 0

            self.__vars_float = []
            self.__vars_str   = []
            col        = 0
            length     = len(list(__dict_vars.keys())[0])
            self.__ifloat     = []
            self.__istr       = []

            # get types of each columns and number of each type
            for i in range(self.__hrnvars):
                names = head[1][col:col+length]
                col += length
                if __dict_vars[names] in str_vars:
                    self.__vars_str.append(__dict_vars[names])
                    self.__istr.append(i)
                else:
                    self.__ifloat.append(i)
                    self.__vars_float.append(__dict_vars[names])

            for i in range(self.__hrnvars,self.__hrnvars_tot):
                self.__ifloat.append(self.__ifloat[-1]+1)

            self.__hrnvars_str   = len(self.__vars_str)
            self.__hrnvars_float = self.__hrnvars_tot - self.__hrnvars_str

            self.rcz, self.mcz, self.rov, self.lconv = [[] for _ in range(4)]
            self.ab_c, self.ab_s, self.abn_s = [{} for _ in range(3)]

            if self.n_contours_burn > 0:
                self.n_burn_zones = np.zeros( (self.n_contours_burn, 0), dtype=int )
                self.m_burn, self.r_burn = [[] for _ in range(2)]

            self.nmod = 0

        with open(self.hrnew, 'r', encoding='latin-1') as f:
            f.seek( self.HRoffset )
            lines         = f.readlines()
            self.HRoffset = f.seek( 0, 1 )

        nlines     = len(lines)
        old_nmod = self.nmod
        self.nmod += nlines

        data              = np.zeros( (self.__hrnvars_float, nlines ) )
        data_str          = [[] for _ in range(self.__hrnvars_str)]
        if read_from_start:
            if self.n_contours_burn > 0:
                self.n_burn_zones = np.zeros( (self.n_contours_burn, self.nmod), dtype=int )
                self.m_burn, self.r_burn = [[] for _ in range(2)]
        else:
            if self.n_contours_burn > 0:
                self.n_burn_zones = np.append( self.n_burn_zones, np.zeros( (self.n_contours_burn, nlines), dtype=int ),
                    axis=1 )

        max_ifloat = np.max( self.__ifloat )
        for i in range(nlines):
            line = lines[i].split()
            if len( line ) < max_ifloat:
                raise CESAMError(f"Line of file {self.hrnew} starting at line {old_nmod+i+3} is probably truncated.")
            for j in range(self.__hrnvars_float):
                try:
                    test = float(line[self.__ifloat[j]])
                    data[j, i] = float(line[self.__ifloat[j]])
                except ValueError:
                    data[j, i] = 0.0

            for j in range(self.__hrnvars_str):
                data_str[j].append(line[self.__istr[j]].strip())

            nconv = int(line[self.__hrnvars_tot])
            ivar  = self.__hrnvars_tot + 1
            lconvd, mczd, rczd, rovd, mst, rst = [[] for _ in range(6)]
            for k in range(nconv):
                lconvd.append(line[ivar])
                ivar += 1


            for k in range(nconv):
                mczd.append(float(line[ivar]))
                rczd.append(float(line[ivar+1]))
                try:
                    rovd.append(float(line[ivar+2]))
                except ValueError:
                    rovd.append(-1.0)

                ivar += 3

            self.lconv.append(lconvd)
            self.mcz.append(mczd)
            self.rcz.append(rczd)
            self.rov.append(rovd)

            if self.n_contours_burn > 0:
                for contour in range(self.n_contours_burn):
                    self.n_burn_zones[contour, self.nmod-nlines+i] = int(line[ivar])
                    mstd, rstd = [[] for _ in range(2)]
                    ivar += 1
                    for zone in range(self.n_burn_zones[contour, self.nmod-nlines+i]):
                        mstd.append([float(line[ivar]), float(line[ivar+1])])
                        rstd.append([float(line[ivar+2]), float(line[ivar+3])])
                        ivar += 4

                    mst.append(mstd)
                    rst.append(rstd)

                self.m_burn.append(mst)
                self.r_burn.append(rst)

        for i in range(self.__hrnvars-self.__hrnvars_str):
            if not read_from_start:
                self.__setattr__( self.__vars_float[i], np.append(
                    getattr( self, self.__vars_float[i] ), data[i,:], axis=0 ) )
            else:
                self.__setattr__( self.__vars_float[i], data[i,:] )

        for i in range(self.__hrnvars_str):
            if not read_from_start:
                self.__setattr__( self.__vars_str[i], np.append(
                    getattr( self, self.__vars_str[i] ), data_str[i], axis=0 ) )
            else:
                self.__setattr__( self.__vars_str[i], data_str[i] )

        if hasattr( self, 'imod' ):
            self.imod = self.imod.astype(int)

        self.log_g = const1 + np.log10(self.mstar) - 2.0*self.log_r

        for i in range(self.nchim):
            key = self.nom_elem[i]
            if key in self.ab_c.keys():
                self.ab_c[key]  = np.append( self.ab_c[key], data[self.__hrnvars-self.__hrnvars_str+2*i,:] )
                self.ab_s[key]  = np.append( self.ab_s[key], data[self.__hrnvars-self.__hrnvars_str+2*i+1,:] )
            else:
                self.ab_c[key]  = data[self.__hrnvars-self.__hrnvars_str+2*i,:]
                self.ab_s[key]  = data[self.__hrnvars-self.__hrnvars_str+2*i+1,:]

            self.abn_s[key] = self.ab_s[key]

            if key in nucleo.keys():
                self.abn_s[key] = self.abn_s[key]/nucleo[key]
            else:
                self.abn_s[key] = self.abn_s[key]/float(key[-2:])

        # abundances per mole (H = 12 dex)
        sum_mass = np.zeros(self.nmod)

        for key in self.ab_s:
            sum_mass += np.array( self.abn_s[key] )

        for key in self.ab_s:
            self.abn_s[key] /= sum_mass
            self.abn_s[key] = np.log10(self.abn_s[key]/self.abn_s['H1']) + 12


        try:
            self.L_nuc = self.L_pp + self.L_cno + self.L_3al
        except AttributeError:
            print('\nNo L_nuc in HR file. ')

        try:
            self.pms  = [i for i in range(self.nmod) if self.evol_phase[i] == 'PMS']
            self.ms   = [i for i in range(self.nmod) if self.evol_phase[i] == 'MS']
            self.rgb  = [i for i in range(self.nmod) if self.evol_phase[i] == 'RGB']
            self.cohe = [i for i in range(self.nmod) if self.evol_phase[i] == 'CoHe']
            self.agb  = [i for i in range(self.nmod) if self.evol_phase[i] == 'AGB']
        except AttributeError:
            print('No evolutionary phases found in HRnew file.')

        if process_cz:
            self.__process_cz()

            if self.n_contours_burn > 0:
                self.__process_bz( )
            else:
                print('\nNo burn zones in HR file.')

        x = sum([self.ab_s[i][:] for i in self.ab_s if i in ['H1', 'H2']])
        y = sum([self.ab_s[i][:] for i in self.ab_s if i in ['He3', 'He4']])
        self.z = 1.0 - x - y

        self.FeH  = np.log10(self.z/x) - np.log10(zsx_sol)


        self.rstar  = 10.0**self.log_r

        if verbose: print(f"{GREEN}[Done]{NO_COLOR}")

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def __read_hr_old(self, obs=False, zsx_sol=None, verbose=True):
        """
        Reads the model's HR file. Computes colors and visual magnitudes if
        required. Teff-color transformations and bolometric corrections from
        VandenBerg & Clem, 2003, AJ, 126, 778 (uses the program described by
        VandenBerg et al., 2006, ApJS, 162, 375).

        :kparam obs: If True, computes colors and visual magnitudes. Default: False
        :ktype obs: bool

        :kparam zsx_sol: Solar Z/X. If None, will use solar Z/X correponding to model parameter 'Initial abundancies'.
            If unable, will use default value of 0.0181 (Asplund et al. 2009).
        :ktype zsx_sol: float

        :kparam verbose: If True, print some details.
        :ktype verbose: bool

        :member age: Age of each models.
        :mtype age: np.ndarray

        :member log_teff: $\log T_{\rm eff}$ of each models.
        :mtype log_teff: np.ndarray

        :member log_l: $\log T_{\rm eff}$ of each models.
        :mtype log_l: np.ndarray

        :member log_r: $\log T_{\rm eff}$ of each models.
        :mtype log_r: np.ndarray

        :member log_g: $\log T_{\rm eff}$ of each models.
        :mtype log_g: np.ndarray

        :member FeH: $[\rm Fe/H]$ of each models.
        :mtype FeH: np.ndarray

        :member z: Metallicity of each models.
        :mtype z: np.ndarray

        :member mstar: $M_\star$ of each models.
        :mtype mstar: np.ndarray

        :member omega_s: Surface angular velocity of each models.
        :mtype omega_s: np.ndarray

        :member vrot: Rotation velocity of each models.
        :mtype vrot: np.ndarray

        :member alpha: $\alpha$ parameter of each models.
        :mtype alpha: np.ndarray

        :member senv: Entropy of adiabat of each models.
        :mtype senv: np.ndarray

        :member mu: Mean molecular weight at bottom of CZ of each models.
        :mtype mu: np.ndarray
        """

        self.err_hr = False
        self.hr_empty  = False
        zsx_sol = zsx_sun[self.params.nom_abon_] if zsx_sol is None else zsx_sol
        old = 0  # reads old format

        self.hr_empty = os.path.getsize( self.hr ) < 1
        if self.hr_empty:
            self.err_hr = True
            print(f"{RED}Error:{NO_COLOR} file {self.hr} exists but is empty.")
            return

        const1 = np.log10(self.ctes.gmsun/self.ctes.rsun**2)

        pwd = os.getcwd()
        vicdir = f"{os.environ.get('CESDIR')}/vicmodels"

        self.imod     = np.array([], dtype=int)
        self.age, self.log_teff, self.log_l, self.log_r, self.log_g, \
                  self.FeH,      self.z,     self.mstar, self.omega_s, \
                  self.vrot,     self.alpha, self.senv,  self.mu \
                  = [np.array([]) for _ in range(13)]
        self.rcz, self.mcz, self.rov, self.lconv, self.nom_elem = [[] for _ in range(5)]
        self.ab_c, self.ab_s, self.abn_s = [{} for _ in range(3)]

        if verbose: print(f"Reading {self.hr}...", end='')
        with open(self.hr, 'r', encoding='latin-1') as f:
            cont = f.readlines()
        size = len(cont)

        first = True

        i = 0
        if size == 0:
            if verbose:
                print(f"{YELLOW}[Empty]{NO_COLOR} ", end='')
            self.__close_hr(obs=obs, err=True)
            self.__process_cz()

            return

        while i < size:
            try:
                if cont[i][0:4] == 'Wrot':
                    i += 1
            except IndexError:
                print(f"{RED} Problem at line: {NO_COLOR}{i}. ", end='')
                self.__close_hr(obs=obs, err=True)
                self.__process_cz()

                return

            try:
                self.age = np.append(self.age, float(cont[i][1:22]))
                self.nchim = int(cont[i][23:25])
                krot = int(cont[i][27])
                nconv = int(cont[i][30])
                self.imod = np.append( self.imod, int( cont[i][32:37] ) )
                lconvd = []
                for k in range(nconv):
                    if old == 0:
                        old = False
                        try:
                            lconvd.append(cont[i][38+2*k])
                        except IndexError:
                            old = True
                            lconvd.append(cont[i][36+2*k])
                    else:
                        if old:
                            lconvd.append(cont[i][36+2*k])
                        else:
                            lconvd.append(cont[i][38+2*k])

                self.lconv.append(lconvd)
            except (IndexError, ValueError):
                print(f"{RED} Problem at line: {NO_COLOR}{i}. ", end='')
                self.__close_hr(obs=obs, err=True)
                self.__process_cz()

                return

            i += 1

            try:
                step = 14
                if old:
                    step -= 1

                self.log_teff = np.append(self.log_teff, float(cont[i][     0:  step]))
                self.log_l    = np.append(self.log_l,    float(cont[i][  step:2*step]))
                self.log_r    = np.append(self.log_r,    float(cont[i][2*step:3*step]))
                self.mstar    = np.append(self.mstar,    float(cont[i][3*step:4*step]))

                const2 = const1 + np.log10(self.mstar[-1]) - 2.0*self.log_r[-1]
                self.log_g = np.append(self.log_g, const2)
            except (IndexError, ValueError):
                print(f"{RED} Problem at line: {NO_COLOR}{i}. ", end='')
                self.__close_hr(obs=obs, err=True)
                self.__process_cz()

                return

            mczd = []
            rczd = []
            rovd = []

            # convective zones

            for j in range(nconv):
                if j % 2 == 0:
                    try:
                        if old:
                            mczd.append(float(cont[i][52:65]))
                            rczd.append(float(cont[i][65:78]))
                        else:
                            mczd.append(float(cont[i][57:70]))
                            rczd.append(float(cont[i][71:84]))
                    except IndexError:
                        self.__close_hr(obs=obs, err=True)
                        self.__process_cz()

                        return

                    i += 1   # next line

                    try:
                        if old:
                            rovd.append(float(cont[i][0:13]))
                        else:
                            rovd.append(float(cont[i][0:14]))
                    except (ValueError, IndexError):
                        rovd.append(0.0)
                else:
                    try:
                        if old:
                            mczd.append(float(cont[i][13:26]))
                            rczd.append(float(cont[i][26:39]))
                            try:
                                rovd.append(float(cont[i][39:52]))
                            except ValueError:
                                rovd.append(0.0)
                        else:
                            mczd.append(float(cont[i][15:28]))
                            rczd.append(float(cont[i][29:42]))
                            try:
                                rovd.append(float(cont[i][43:56]))
                            except ValueError:
                                rovd.append(0.0)
                    except IndexError:
                        self.__close_hr(obs=obs, err=True)
                        self.__process_cz()

                        return
                    except ValueError:
                        print(f"{RED} Problem at line: {NO_COLOR}{i}. ", end='')

            self.mcz.append(mczd)
            self.rcz.append(rczd)
            try:
                self.rov.append(rovd)
            except ValueError:
                print(f"{RED} Problem at line: {NO_COLOR}{i}. ", end='')
                rovd.append(0.0)

            try:
                self.alpha = np.append( self.alpha, float(cont[i+1][0:14]) )
                self.senv  = np.append( self.senv,  float(cont[i+1][14:28]) )
                self.mu    = np.append( self.mu,    float(cont[i+1][28:42]) )
                i += 1
            except ValueError: # Old format where alpha, s_env and mu were not specified
                pass
            except IndexError: # Reach end of file
                print(f"{RED} Problem at line: {NO_COLOR}{i}. ", end='')

            i+=1 # advance line
            for j in range(self.nchim):
                try:
                    nom_elem = cont[i][0:4].strip()
                    try:
                        ab_c = float(cont[i][4:16])   # central abundances
                    except ValueError:
                        ab_c = 0.0
                    try:
                        ab_s = float(cont[i][16:28])  # surface abundances
                    except ValueError:
                        ab_s = 0.0

                    if nom_elem in self.ab_c:
                        self.ab_c[nom_elem] = np.append(self.ab_c[nom_elem], ab_c)
                        self.ab_s[nom_elem] = np.append(self.ab_s[nom_elem], ab_s)
                    else:
                        self.ab_c[nom_elem] = np.array([ab_c])
                        self.ab_s[nom_elem] = np.array([ab_s])

                except IndexError:
                    print(f"{RED} Problem at line: {NO_COLOR}{i}. ", end='')
                    self.__close_hr(obs=obs, err=True)
                    self.__process_cz()

                    return

                i += 1  # advance line

            if first:
                self.nom_elem = list( self.ab_c.keys() )
                first = False


            # rotation
            if krot > 1:
                try:
                    self.vrot = np.append(self.vrot, float(cont[i][4:16]))
                    self.omega_s = np.append(self.omega_s, float(cont[i][16:28]))
                except IndexError:
                    print(f"{RED} Problem at line: {NO_COLOR}{i}. ", end='')
                    self.__close_hr(obs=obs, err=True)
                    self.__process_cz()

                    return

                i += 1  # advance line


            x = 0.0
            z = 0.0
            for key in self.ab_s:
                if key in ['H1', 'H2']:
                    x += self.ab_s[key][-1]
                elif not key in ['He3', 'He4']:
                    z += self.ab_s[key][-1]

            self.z = np.append(self.z, z)
            if z > 0.0:
                logz = np.log10(z/x) - np.log10(zsx_sol)
                self.FeH = np.append(self.FeH,logz)

        self.__close_hr(obs=obs, verbose=verbose, old=old)
        self.__process_cz()


    def __close_hr(self, err=False, obs=False, verbose=True, old=False):
        """
        Properly close a HR file and remove data that have been partially read.

        :kparam err: If True, an error was found and partially read data must be removed.
        :ktype err: boolean

        :kparam obs: If True, some observable parameters are computed such as magnitude and colors.
        :ktype obs: boolean

        :kparam verbose: If False, do not print details, except if err is True.
        :ktype verbose: boolean

        :kparam old: If True, HR file is in olf format
        :ktype old: boolean
        """

        self.err_hr = err
        if verbose and not err:
            print(f"{GREEN}[Done]{NO_COLOR}")
            if old: print('HR file in old format.')
        elif err:
            print(f"{YELLOW}[Closing HR]{NO_COLOR}")
            if old: print('HR file in old format.')

        self.nmod = len(self.age)
        if err:
            self.nmod -= 1

            self.age      = self.age[:self.nmod]
            self.mcz      = self.mcz[:self.nmod]
            self.rcz      = self.rcz[:self.nmod]
            self.rov      = self.rov[:self.nmod]
            self.mstar    = self.mstar[:self.nmod]
            self.log_g    = self.log_g[:self.nmod]
            self.log_teff = self.log_teff[:self.nmod]
            self.log_l    = self.log_l[:self.nmod]
            self.log_r    = self.log_r[:self.nmod]
            self.FeH      = self.FeH[:self.nmod]
            self.z        = self.z[:self.nmod]
            try:
                self.alpha = self.alpha[:self.nmod]
                self.senv  = self.senv[:self.nmod]
                self.mu    = self.mu[:self.nmod]
            except NameError:
                pass
            for key in self.ab_s:
                self.ab_c[key] = self.ab_c[key][:self.nmod]
                self.ab_s[key] = self.ab_s[key][:self.nmod]
            return

        self.rstar  = 10.0**self.log_r
        self.nu_max = 3050.0*self.mstar/(self.rstar**2)/np.sqrt((10.0**self.log_teff)/5777.0)
        self.Dnu    = 134.7*np.sqrt(self.mstar/(self.rstar**3))

        # abundances per mole (H = 12 dex)
        sum_mass = np.zeros(self.nmod)
        for key in self.ab_s:
            try:
                self.abn_s[key] = self.ab_s[key]/nucleo[key]
            except KeyError:
                self.abn_s[key] = self.ab_s[key]/float(key[-2:])

        for key in self.ab_s:
            sum_mass += self.abn_s[key]

        for key in self.ab_s:
            self.abn_s[key] /= sum_mass
            self.abn_s[key] = np.log10(self.abn_s[key]/self.abn_s['H1']) + 12

        try:
            self.p_rot = 2.0*np.pi/self.omega_s/3600.0/24.0
        except:
            return

        if obs:
            # calculate observational CMD
            self.mv = np.array([])
            self.bv = np.array([])
            self.vr = np.array([])
            self.vi = np.array([])

            cmd = ''
            for i in range(self.nmod):
                cmd += f'{self.FeH[i]}, {self.log_g[i]}, {self.log_teff[i]}\n'

            cmd += '-20.0, -20.0, -20.0\n'

            os.chdir(vicdir)
            s = subp.Popen('./teff-color.x', stdin=subp.PIPE, stdout=subp.PIPE)

            sortie, err = s.communicate(input=cmd.encode())
            sortie = sortie.decode().split('\n')
            os.chdir(pwd)

            i = 0
            for line in sortie[2:-1]:
                bv, vr, vi, corr = line.split()

                self.bv = np.append(self.bv, float(bv))
                self.vr = np.append(self.vr, float(vr))
                self.vi = np.append(self.vi, float(vi))

                Mbol = mbol_sun - 2.5*self.log_l[i] - float(corr) # in fact it is Mv!!!
                self.mv = np.append(self.mv, Mbol)  # <-----------------+
                i += 1

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def __process_cz(self):
        """
        Processe limits of convection zones, turning them into easier to plot quantities
        :member m_conv: $m/M_{\odot}$ of each limit for each time step.
        :mtype m_conv: list

        :member r_conv: $m/M_{\odot}$ of each limit for each time step.
        :mtype r_conv: list

        :member iconv_start: index of the beginning of CZ.
        :mtype iconv_start: list

        :member iconv_end: index of the end of CZ.
        :mtype iconv_end: list

        :member mconv_start: lower mass of the CZ
        :mtype mconv_start: list

        :member rconv_start: lower radius of the CZ
        :mtype rconv_start: list

        :member mconv_end: upper mass of the CZ
        :mtype mconv_end: list

        :member rconv_end: upper radius of the CZ
        :mtype rconv_end: list
        """

        if self.nmod < 1:
            return
        mconv_temp = []
        rconv_temp = []
        ncz   = []
        for i in range(self.nmod):
            # determine number of cz
            nlims = len(self.mcz[i])

            mconv_temp.append([])
            rconv_temp.append([])
            mconv_current = []
            rconv_current = []
            ncz_current   = 0

            if nlims > 0:
                if self.lconv[i][0] == 'F': # 1st limit is an end of a CZ -> central CZ
                    mconv_current = [0.0]   # start of CZ is at center
                    rconv_current = [0.0]
                    ncz_current   = 1

            for j in range(nlims):
                if self.lconv[i][j] == 'F': # end of a CZ
                    mconv_current.append(self.mstar[i] - self.mcz[i][j])
                    mconv_temp[-1].append(mconv_current)
                    mconv_current = []

                    rconv_current.append(self.rcz[i][j])
                    rconv_temp[-1].append(rconv_current)
                    rconv_current = []
                else:  # start of a new CZ. Increares nzc,
                    ncz_current += 1
                    mconv_current.append(self.mstar[i] - self.mcz[i][j])
                    rconv_current.append(self.rcz[i][j])
            if len(mconv_current) == 1:
                mconv_current.append(self.mstar[i])
                mconv_temp[-1].append(mconv_current)

                rconv_current.append(10**self.log_r[i])
                rconv_temp[-1].append(rconv_current)
            ncz.append(ncz_current)


        tol = 0.15

        nzones_tot = len(mconv_temp[0])
        alive = [[True for i in range(nzones_tot)]]
        nzones = [nzones_tot]
        self.mconv = [mconv_temp[0]]
        self.rconv = [rconv_temp[0]]
        self.iconv_start = []
        self.iconv_end   = []
        if nzones_tot > 0:
            self.iconv_start = [0]*nzones_tot
            self.iconv_end   = [0]*nzones_tot

        for i in range(1,self.nmod):
            nzones.append(len(mconv_temp[i]))
            self.mconv.append([None for k in range(nzones_tot)])
            self.rconv.append([None for k in range(nzones_tot)])
            diff = nzones[i] - nzones[i-1]

            if diff < 0: # CZ disappear, nzone_tot does not change
                alive_current = [False for k in range(nzones_tot)]
                checked_j = [False for j in range(nzones[i])]
                for k in range(nzones_tot):
                    if alive[i-1][k]:

                        for j in range(nzones[i]): # check which ones are dead
                            if not checked_j[j]:
                                diff_s = np.abs(mconv_temp[i][j][0]  - self.mconv[i-1][k][0])/self.mstar[i]
                                diff_e = np.abs(mconv_temp[i][j][1]  - self.mconv[i-1][k][1])/self.mstar[i]
                                if np.abs(self.mconv[i-1][k][1] - self.mstar[i])/self.mstar[i] < 1.0e-8: #surface
                                    if np.abs(mconv_temp[i][j][1] - self.mstar[i])/self.mstar[i] < 1.0e-8:
                                        alive_current[k] = True
                                elif self.mconv[i-1][k][0]/self.mstar[i] < 1.0e-8:
                                    alive_current[k] = diff_s < tol
                                else:
                                    alive_current[k] = (diff_s < tol) & (diff_e < tol)

                                if alive_current[k]:
                                    checked_j[j] = True
                                    self.mconv[i][k] = mconv_temp[i][j]
                                    self.rconv[i][k] = rconv_temp[i][j]
                                    break  # k is alive, no need to check other j's

                alive.append(alive_current)
                for zone in range(nzones_tot):
                    if alive[i-1][zone] and not alive[i][zone]:
                        self.iconv_end[zone] = i - 1

            else: # CZ appear, nzone_tot increases
                alive.append(alive[i-1])
                if diff == 0 and nzones_tot == 1:
                    self.mconv[i][0] = mconv_temp[i][0]
                    self.rconv[i][0] = rconv_temp[i][0]
                else:
                    alive_current = [False for k in range(nzones_tot)]

                    for j in range(nzones[i]): # check which ones appear
                        for k in range(nzones_tot):
                            found = False

                            if not alive[i][k]: continue #CZs do not ressuscitate
                            if alive_current[k]: continue # already found

                            diff_s = np.abs(mconv_temp[i][j][0]  - self.mconv[i-1][k][0])/self.mstar[i]
                            diff_e = np.abs(mconv_temp[i][j][1]  - self.mconv[i-1][k][1])/self.mstar[i]
                            if np.abs(self.mconv[i-1][k][1] - self.mstar[i])/self.mstar[i] < 1.0e-8: #surrface
                                if np.abs(mconv_temp[i][j][1] - self.mstar[i])/self.mstar[i] < 1.0e-8:
                                    alive_current[k] = True
                                else:
                                    alive_current[k] = False
                            elif self.mconv[i-1][k][0]/self.mstar[i] < 1.0e-6:
                                alive_current[k] = diff_s < tol
                            else:
                                alive_current[k] = (diff_s < tol) & (diff_e < tol)
                            if alive_current[k]:
                                self.mconv[i][k] = mconv_temp[i][j]
                                self.rconv[i][k] = rconv_temp[i][j]
                                found = True
                                break  # found k corresponding to j is alive, no need to check other k's


                        if not found:  # new CZ
                            for l in range(i):
                                self.mconv[l].append(None)  # add a new empty CZ
                                self.rconv[l].append(None)  # to previous mconv
                            nzones_tot += 1
                            self.mconv[i].append(mconv_temp[i][j])
                            self.rconv[i].append(rconv_temp[i][j])
                            alive[i].append(True)
                            self.iconv_start.append(i)
                            self.iconv_end.append(0)
                            alive_current.append(True)

                    alive[i] = alive_current
                    for k in range(nzones_tot):
                        if not alive[i][k] and alive[i-1][k]:
                            self.iconv_end[k] = i - 1


        for i in range(nzones_tot):
            if self.iconv_end[i] == 0: self.iconv_end[i] = self.nmod - 1

        self.mconv_start = [[] for i in range(nzones_tot)]
        self.mconv_end = [[] for i in range(nzones_tot)]

        self.rconv_start = [[] for i in range(nzones_tot)]
        self.rconv_end = [[] for i in range(nzones_tot)]

        for i in range(nzones_tot):
            start = self.iconv_start[i] - 1 if self.iconv_start[i] > 0 else 0
            end   = self.iconv_end[i]   + 1 if self.iconv_end[i]   < self.nmod - 1 else self.nmod - 1
            for j in range(self.iconv_start[i], self.iconv_end[i]+1):
                if j == self.iconv_start[i] and j > 0:
                    if np.abs(self.mconv[j][i][0])/self.mstar[j] < 1e-8: # central CZ
                        self.mconv[j-1][i] = [self.mconv[j][i][0], self.mconv[j][i][0]]
                        self.mconv_start[i].append(self.mconv[j][i][0]) # start with mstart = mend = 0
                        self.mconv_end[i].append(self.mconv[j][i][0])

                        self.rconv[j-1][i] = [self.rconv[j][i][0], self.rconv[j][i][0]]
                        self.rconv_start[i].append(self.rconv[j][i][0])
                        self.rconv_end[i].append(self.rconv[j][i][0])

                    elif np.abs(self.mconv[j][i][1] - self.mstar[j])/self.mstar[j] < 1.0e-6: # surface CZ
                        self.mconv[j-1][i] = [self.mconv[j][i][1], self.mconv[j][i][1]]
                        self.mconv_start[i].append(self.mconv[j][i][1]) # start with mstart = mend = mstar
                        self.mconv_end[i].append(self.mconv[j][i][1])

                        self.rconv[j-1][i] = [self.rconv[j][i][1], self.rconv[j][i][1]]
                        self.rconv_start[i].append(self.rconv[j][i][1])
                        self.rconv_end[i].append(self.rconv[j][i][1])
                    else:  # start with average between start and end
                        m_ave = 0.5*(self.mconv[j][i][0] + self.mconv[j][i][1])
                        self.mconv[j-1][i] = [m_ave, m_ave]
                        self.mconv_start[i].append(m_ave)
                        self.mconv_end[i].append(m_ave)

                        r_ave = 0.5*(self.rconv[j][i][0] + self.rconv[j][i][1])
                        self.mconv[j-1][i] = [r_ave, r_ave]
                        self.rconv_start[i].append(r_ave)
                        self.rconv_end[i].append(r_ave)

                self.mconv_start[i].append( self.mconv[j][i][0])
                self.mconv_end[i].append(   self.mconv[j][i][1])
                self.rconv_start[i].append( self.rconv[j][i][0])
                self.rconv_end[i].append(   self.rconv[j][i][1])


            if self.iconv_end[i] < self.nmod - 1:
                if np.abs(self.mconv[self.iconv_end[i]][i][0])/self.mstar[self.iconv_end[i]] < 1e-6: # central CZ
                    self.mconv[j+1][i] = [self.mconv[j][i][0], self.mconv[j][i][0]]
                    self.mconv_start[i].append(self.mconv[self.iconv_end[i]][i][0])
                    self.mconv_end[i].append(self.mconv[self.iconv_end[i]][i][0])

                    self.rconv[j+1][i] = [self.rconv[j][i][0], self.rconv[j][i][0]]
                    self.rconv_start[i].append(self.rconv[self.iconv_end[i]][i][0])
                    self.rconv_end[i].append(self.rconv[self.iconv_end[i]][i][0])
                elif np.abs(self.mconv[self.iconv_end[i]][i][1] - self.mstar[self.iconv_end[i]])/self.mstar[self.iconv_end[i]] < 1.0e-6: # surface CZ
                    self.mconv[j+1][i] = [self.mconv[j][i][1], self.mconv[j][i][1]]
                    self.mconv_start[i].append(self.mconv[self.iconv_end[i]][i][1])
                    self.mconv_end[i].append(self.mconv[self.iconv_end[i]][i][1])

                    self.rconv[j+1][i] = [self.rconv[j][i][1], self.rconv[j][i][1]]
                    self.rconv_start[i].append(self.rconv[self.iconv_end[i]][i][1])
                    self.rconv_end[i].append(self.rconv[self.iconv_end[i]][i][1])
                else:
                    m_ave = 0.5*(self.mconv[self.iconv_end[i]][i][0] + self.mconv[self.iconv_end[i]][i][1])
                    self.mconv[j+1][i] = [m_ave, m_ave]
                    self.mconv_start[i].append(m_ave)
                    self.mconv_end[i].append(m_ave)

                    r_ave = 0.5*(self.rconv[self.iconv_end[i]][i][0] + self.rconv[self.iconv_end[i]][i][1])
                    self.rconv[j+1][i] = [r_ave, r_ave]
                    self.rconv_start[i].append(r_ave)
                    self.rconv_end[i].append(r_ave)

            self.iconv_start[i] = start
            self.iconv_end[i]   = end

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def __process_bz( self ):
        """
        Processe limits of burning zones, turning them into easier to plot quantities
        """
        self.mburn = [[], [], [], []]
        self.rburn = [[], [], [], []]
        self.iburn_start = [[], [], [], []]
        self.iburn_end   = [[], [], [], []]
        self.mburn_start = [[], [], [], []]
        self.mburn_end   = [[], [], [], []]
        self.rburn_start = [[], [], [], []]
        self.rburn_end   = [[], [], [], []]

        tol = 0.05
        for contour in range(4):
            nzones = []
            alive = []
            start = -1
            for i in range(self.nmod):
                nzones_tot = self.n_burn_zones[contour,i]
                alive.append([True for i in range(nzones_tot)])
                nzones.append(nzones_tot)
                if nzones_tot > 0:
                    self.mburn[contour].append(self.m_burn[i][contour])
                    self.rburn[contour].append(self.r_burn[i][contour])
                    self.iburn_start[contour] = [i]*nzones_tot
                    self.iburn_end[contour]   = [0]*nzones_tot

                    start = i
                    break
                else:
                    self.mburn[contour].append([None])
                    self.rburn[contour].append([None])

            if start < 0: continue

            for i in range(start+1,self.nmod):
                nzones.append(self.n_burn_zones[contour,i])
                self.mburn[contour].append([None for k in range(nzones_tot)])
                self.rburn[contour].append([None for k in range(nzones_tot)])
                diff = nzones[i] - nzones[i-1]

                if diff < 0: # CZ disappear, nzone_tot does not change
                    alive_current = [False for k in range(nzones_tot)]
                    checked_j = [False for j in range(nzones[i])]
                    for k in range(nzones_tot):
                        if alive[i-1][k]:

                            for j in range(nzones[i]): # check which ones are dead
                                if not checked_j[j]:
                                    diff_s = np.abs(self.m_burn[i][contour][j][0]  - self.mburn[contour][i-1][k][0])/self.mstar[i]
                                    diff_e = np.abs(self.m_burn[i][contour][j][1]  - self.mburn[contour][i-1][k][1])/self.mstar[i]
                                    if np.abs(self.mburn[contour][i-1][k][1] - self.mstar[i])/self.mstar[i] < 1.0e-8: #surface
                                        if np.abs(self.m_burn[i][contour][j][1] - self.mstar[i])/self.mstar[i] < 1.0e-8:
                                            alive_current[k] = True
                                    elif self.mburn[contour][i-1][k][0]/self.mstar[i] < 1.0e-8:
                                        alive_current[k] = diff_s < tol
                                    else:
                                        alive_current[k] = (diff_s < tol) & (diff_e < tol)

                                    if alive_current[k]:
                                        checked_j[j] = True
                                        self.mburn[contour][i][k] = self.m_burn[i][contour][j]
                                        self.rburn[contour][i][k] = self.r_burn[i][contour][j]
                                        break  # k is alive, no need to check other j's

                    alive.append(alive_current)
                    for zone in range(nzones_tot):
                        if alive[i-1][zone] and not alive[i][zone]:
                            self.iburn_end[contour][zone] = i - 1

                else: # CZ appear, nzone_tot increases
                    alive.append(alive[i-1])
                    if (nzones[i] == 0):
                        for j in range(nzones_tot):
                            alive[i][j] = False
                            self.mburn[contour][i][j] = None
                            self.rburn[contour][i][j] = None
                    if diff == 0 and nzones_tot == 1 and nzones[i] > 0:
                        self.mburn[contour][i][0] = self.m_burn[i][contour][0]
                        self.rburn[contour][i][0] = self.r_burn[i][contour][0]
                    else:
                        alive_current = [False for k in range(nzones_tot)]

                        for j in range(nzones[i]): # check which ones appear
                            found = False
                            for k in range(nzones_tot):
                                found = False

                                if not alive[i][k]: continue #CZs do not ressuscitate
                                if alive_current[k]: continue # already found

                                diff_s = np.abs(self.m_burn[i][contour][j][0]  - self.mburn[contour][i-1][k][0])/self.mstar[i]
                                diff_e = np.abs(self.m_burn[i][contour][j][1]  - self.mburn[contour][i-1][k][1])/self.mstar[i]
                                if np.abs(self.mburn[contour][i-1][k][1] - self.mstar[i])/self.mstar[i] < 1.0e-8: #surrface
                                    if np.abs(self.m_burn[i][contour][j][1] - self.mstar[i])/self.mstar[i] < 1.0e-8:
                                        alive_current[k] = True
                                    else:
                                        alive_current[k] = False
                                elif self.mburn[contour][i-1][k][0]/self.mstar[i] < 1.0e-6:
                                    alive_current[k] = diff_s < tol
                                else:
                                    alive_current[k] = (diff_s < tol) & (diff_e < tol)
                                if alive_current[k]:
                                    self.mburn[contour][i][k] = self.m_burn[i][contour][j]
                                    self.rburn[contour][i][k] = self.r_burn[i][contour][j]
                                    found = True
                                    break  # found k corresponding to j is alive, no need to check other k's


                            if not found:  # new CZ
                                for l in range(i):
                                    self.mburn[contour][l].append(None)  # add a new empty CZ
                                    self.rburn[contour][l].append(None)  # to previous mburn[contour]
                                nzones_tot += 1
                                self.mburn[contour][i].append(self.m_burn[i][contour][j])
                                self.rburn[contour][i].append(self.r_burn[i][contour][j])
                                alive[i].append(True)
                                self.iburn_start[contour].append(i)
                                self.iburn_end[contour].append(0)
                                alive_current.append(True)

                        alive[i] = alive_current
                        for k in range(nzones_tot):
                            if not alive[i][k] and alive[i-1][k]:
                                self.iburn_end[contour][k] = i - 1



            for i in range(nzones_tot):
                if self.iburn_end[contour][i] == 0: self.iburn_end[contour][i] = self.nmod - 1


            for i in range(nzones_tot):
                self.mburn_start[contour].append([])
                self.mburn_end[contour].append([])

                self.rburn_start[contour].append([])
                self.rburn_end[contour].append([])
                start = self.iburn_start[contour][i] - 1 if self.iburn_start[contour][i] > start else start
                end   = self.iburn_end[contour][i]   + 1 if self.iburn_end[contour][i]   < self.nmod - 1 else self.nmod - 1
                for j in range(self.iburn_start[contour][i], self.iburn_end[contour][i]+1):
                    if j == self.iburn_start[contour][i] and j > start:
                        if np.abs(self.mburn[contour][j][i][0])/self.mstar[j] < 1e-8: # central CZ
                            self.mburn[contour][j-1][i] = [self.mburn[contour][j][i][0], self.mburn[contour][j][i][0]]
                            self.mburn_start[contour][i].append(self.mburn[contour][j][i][0]) # start with mstart = mend = 0
                            self.mburn_end[contour][i].append(self.mburn[contour][j][i][0])

                            self.rburn[contour][j-1][i] = [self.rburn[contour][j][i][0], self.rburn[contour][j][i][0]]
                            self.rburn_start[contour][i].append(self.rburn[contour][j][i][0])
                            self.rburn_end[contour][i].append(self.rburn[contour][j][i][0])

                        elif np.abs(self.mburn[contour][j][i][1] - self.mstar[j])/self.mstar[j] < 1.0e-6: # surface CZ
                            self.mburn[contour][j-1][i] = [self.mburn[contour][j][i][1], self.mburn[contour][j][i][1]]
                            self.mburn_start[contour][i].append(self.mburn[contour][j][i][1]) # start with mstart = mend = mstar
                            self.mburn_end[contour][i].append(self.mburn[contour][j][i][1])

                            self.rburn[contour][j-1][i] = [self.rburn[contour][j][i][1], self.rburn[contour][j][i][1]]
                            self.rburn_start[contour][i].append(self.rburn[contour][j][i][1])
                            self.rburn_end[contour][i].append(self.rburn[contour][j][i][1])
                        else:  # start with average between start and end
                            m_ave = 0.5*(self.mburn[contour][j][i][0] + self.mburn[contour][j][i][1])
                            self.mburn[contour][j-1][i] = [m_ave, m_ave]
                            self.mburn_start[contour][i].append(m_ave)
                            self.mburn_end[contour][i].append(m_ave)

                            r_ave = 0.5*(self.rburn[contour][j][i][0] + self.rburn[contour][j][i][1])
                            self.mburn[contour][j-1][i] = [r_ave, r_ave]
                            self.rburn_start[contour][i].append(r_ave)
                            self.rburn_end[contour][i].append(r_ave)

                    if self.mburn[contour][j][i] is not None:
                        self.mburn_start[contour][i].append(self.mburn[contour][j][i][0])
                        self.mburn_end[contour][i].append(self.mburn[contour][j][i][1])
                        self.rburn_start[contour][i].append(self.rburn[contour][j][i][0])
                        self.rburn_end[contour][i].append(self.rburn[contour][j][i][1])
                    else:
                        self.mburn_start[contour][i].append(None)
                        self.mburn_end[contour][i].append(None)
                        self.rburn_start[contour][i].append(None)
                        self.rburn_end[contour][i].append(None)


                if self.iburn_end[contour][i] < self.nmod - 1:
                    if np.abs(self.mburn[contour][self.iburn_end[contour][i]][i][0])/self.mstar[self.iburn_end[contour][i]] < 1e-6: # central CZ
                        self.mburn[contour][j+1][i] = [self.mburn[contour][j][i][0], self.mburn[contour][j][i][0]]
                        self.mburn_start[contour][i].append(self.mburn[contour][self.iburn_end[contour][i]][i][0])
                        self.mburn_end[contour][i].append(self.mburn[contour][self.iburn_end[contour][i]][i][0])

                        self.rburn[contour][j+1][i] = [self.rburn[contour][j][i][0], self.rburn[contour][j][i][0]]
                        self.rburn_start[contour][i].append(self.rburn[contour][self.iburn_end[contour][i]][i][0])
                        self.rburn_end[contour][i].append(self.rburn[contour][self.iburn_end[contour][i]][i][0])
                    elif np.abs(self.mburn[contour][self.iburn_end[contour][i]][i][1] - self.mstar[self.iburn_end[contour][i]])/self.mstar[self.iburn_end[contour][i]] < 1.0e-6: # surface CZ
                        self.mburn[contour][j+1][i] = [self.mburn[contour][j][i][1], self.mburn[contour][j][i][1]]
                        self.mburn_start[contour][i].append(self.mburn[contour][self.iburn_end[contour][i]][i][1])
                        self.mburn_end[contour][i].append(self.mburn[contour][self.iburn_end[contour][i]][i][1])

                        self.rburn[contour][j+1][i] = [self.rburn[contour][j][i][1], self.rburn[contour][j][i][1]]
                        self.rburn_start[contour][i].append(self.rburn[contour][self.iburn_end[contour][i]][i][1])
                        self.rburn_end[contour][i].append(self.rburn[contour][self.iburn_end[contour][i]][i][1])
                    else:
                        m_ave = 0.5*(self.mburn[contour][self.iburn_end[contour][i]][i][0] + self.mburn[contour][self.iburn_end[contour][i]][i][1])
                        self.mburn[contour][j+1][i] = [m_ave, m_ave]
                        self.mburn_start[contour][i].append(m_ave)
                        self.mburn_end[contour][i].append(m_ave)

                        r_ave = 0.5*(self.rburn[contour][self.iburn_end[contour][i]][i][0] + self.rburn[contour][self.iburn_end[contour][i]][i][1])
                        self.rburn[contour][j+1][i] = [r_ave, r_ave]
                        self.rburn_start[contour][i].append(r_ave)
                        self.rburn_end[contour][i].append(r_ave)

                self.iburn_start[contour][i] = start
                self.iburn_end[contour][i]   = end

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def read_allrot(self, smooth=False):
        """
        <CModel.read_allrot(smooth=False)>

        Read all files `*_coeff_rota.dat`.

        :kparam smooth: If True, `U2_t`, `V2_t` and `U_csi_t` are smoothed.
        :ktype smooth: bool

        :member r_t: Radius.
        :mtype r_t: list of lists

        :member m_t: Mass.
        :mtype m_t: list of lists

        :member Omega_t: Angular velocity.
        :mtype Omega_t: list of lists

        :member U2_t: Vertical meridional cirulation.
        :mtype U2_t: list of lists

        :member V2_t: Horizontal meridional circulation.
        :mtype V2_t: list of lists

        :member Lambda_t: Fluctuation of mean molecular weight.
        :mtype Lambda_t: list of lists

        :member Psi_t: Fluactions of temperature.
        :mtype Psi_t: list of lists

        :member Theta_t: Fluctuation of density.
        :mtype Theta_t: list of lists

        :member T_t: Temperature.
        :mtype T_t: list of lists

        :member U_csi_t: $\frac{1}{2}\rho r^2 U_2$
        :mtype U_csi_t: list of lists
        """
        n = [name for name in os.listdir('.')
             if f"{self.name}_coeff_rota.dat" in name]
        n.sort()
        self.r_t = []
        self.m_t = []
        self.Omega_t = []
        self.U2_t = []
        self.V2_t = []
        self.Lambda_t = []
        self.Psi_t = []
        self.Theta_t = []
        self.T_t = []
        self.U_csi_t = []

        for i in n:
            self.read_rot(i)
            self.r_t.append(self.r)
            self.m_t.append(self.m)
            self.Omega_t.append(self.Omega)
            self.Lambda_t.append(self.Lambda)
            self.Psi_t.append(self.Psi)
            self.Theta_t.append(self.Theta)
            self.T_t.append(self.T)
            if smooth:
                self.U2_t.append(self.U2s)
                self.V2_t.append(self.V2s)
                self.U_csi_t.append(self.U_csis)
            else:
                self.U2_t.append(self.U2)
                self.V2_t.append(self.V2)
                self.U_csi_t.append(self.U_csi)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def read_rot2d(self, filename=None, old=False):
        """
        DEPRECIATED
        Reads the rorational coefficients file filename. If `filename=None`
        (default), reads `self.name+'_coeff_rota.dat2d'`; if filename is an
        integer, reads file `number+'-*_coeff_rota.dat2d'`.

        :raises TypeError: if type is not int or str.

        Sets:
              self.r
              self.m
              self.Omega
              self.U2
              self.V2
              self.Theta
              self.Upsi
              self.Lambda
              self.Deff
              self.Dh
              self.Dv
              self.T
              self.rho
              self.mu
              self.Psi

        """

        if not filename:
            filename = f'{self.name}_coeff_rota.dat2d'
        elif isinstance(filename, int):
            filename = f"{filename:05d}-{self.name}_coeff_rota.dat2d"
        elif not isinstance(filename, str):
            raise TypeError(f"Argument must be int or str, not {type(filename)}.")

        print(f"Reading file {filename}... ")
        with open(filename, 'r', encoding='latin-1') as f:
            cont = f.readlines()

        i = 5
        if cont[i][1:3] == 'ar':
            i +=1
        nchim = int(cont[i][1:3])
        i += 1
        n_rot = int(cont[i][6:11])
        nb_var = int(cont[i][18:20])
        ncoeff = int(cont[i][28:30])

        self.r = np.zeros((n_rot))
        self.m = np.zeros((n_rot))
        m_l = (nb_var - 7) // 6
        self.U      = np.zeros((m_l, n_rot))
        self.Ups    = np.zeros((m_l, n_rot))
        self.Lambda = np.zeros((m_l, n_rot))
        self.Omega  = np.zeros((m_l, n_rot))
        self.Theta  = np.zeros((m_l, n_rot))
        self.Psi    = np.zeros((m_l, n_rot))

        self.Deff = np.zeros((n_rot))
        self.Dh = np.zeros((n_rot))
        self.Dv = np.zeros((n_rot))
        self.dOdm = np.zeros((n_rot))
        self.theta_m = np.zeros((n_rot))
        i+=1

        for l in range(n_rot):
            data = cont[i+l].split()

            self.r[l] = float( data[0] )
            self.m[l] = float( data[1] )
            for ll in range(m_l):
                self.U[     ll, l] = float( data[2+ll] )
                self.Ups[   ll, l] = float( data[2+m_l+ll] )
                self.Lambda[ll, l] = float( data[2+2*m_l+ll] )
                self.Omega[ ll, l] = float( data[2+3*m_l+ll] )
                self.Theta[ ll, l] = float( data[2+4*m_l+ll] )
                self.Psi[   ll, l] = float( data[2+5*m_l+ll] )
            self.Deff[l]    = float( data[2+6*m_l] )
            self.Dh[l]      = float( data[2+6*m_l+1] )
            self.Dv[l]      = float( data[2+6*m_l+2] )
            self.dOdm[l]     = float( data[2+6*m_l+3] )
            self.theta_m[l] = float( data[2+6*m_l+4] )

        print(f"{GREEN}[Done]{NO_COLOR}")


    def read_rot(self, filename=None, target_age=None, old=False):
        """
        Reads the rorational coefficients file filename.

        :kparam filename: Either a string corresponding to the name of the file, or
            an integer corresponding to the number of the time step (reads `number+'-*_coeff_rota.dat'`).
            If filename=None (default) and `target_age=None`, reads `self.name+'_coeff_rota.dat'`.
        :ktype filename: string or int

        :kparam target_age: look of the model with closest age to `target_age`.
        :ktype target_age: float

        :kparam old: If True, the file is written in an old format.
        :ktype old: boolean

        :member r: Radius.
        :mtype r: list

        :member m: Mass.
        :mtype m: list

        :member Omega: Angular velocity.
        :mtype Omega: list

        :member U2: Vertical meridional cirulation.
        :mtype U2: list

        :member V2: Horizontal meridional circulation.
        :mtype V2: list

        :member Lambda: Fluctuation of mean molecular weight.
        :mtype Lambda: list

        :member Upsi: Value of Upsilon, see Cesam2k20's doc.
        :mtype Upsi: list

        :member Psi: Fluactions of temperature.
        :mtype Psi: list

        :member Theta: Fluctuation of density.
        :mtype Theta: list

        :member Deff: Effective turbulent diffusion coefficient.
        :mtype Deff: list

        :member Dh: Horizontal turbulent diffusion coefficient.
        :mtype Dh: list

        :member Dv: Vertical turbulent diffusion coefficient.
        :mtype Dv: list

        :member T: Temperature.
        :mtype T: list

        :member rho: Density.
        :mtype rho: list

        :member mu: Mean molecular weight.
        :mtype mu: list

        :member U2s: Smoothed vertical meridional cirulation.
        :mtype U2s: list

        :member V2s: Smoothed horizontal meridional circulation.
        :mtype V2s: list

        :member U_csi: $\frac{1}{2}\rho r^2 U_2$.
        :mtype U_csi: list

        :member U_csis: Smoothed `U_csi`.
        :mtype U_csis: list


        :raises TypeError(1): If filename is not int, str or float

        :raises TypeError(2): If filename and `target_age` are None at the same time.

        :raises ValueError: If targeted age is higher than maximum age of the model
        """

        if filename is None and target_age is None:
            filename = f'{self.name}_coeff_rota.dat'

        elif filename is not None and target_age is not None:
            raise TypeError("filename and target_age can not be None at the same time.")

        elif target_age is not None:
            if target_age > self.age[-1]:
                raise ValueError(f"Age = {target_age} Myrs greater than agemax = {self.age[-1]}Myrs.")
            i1 = np.argmin( np.abs( self.age - target_age ) )
            filename = f"{i1:05d}-{self.name}_coeff_rota.dat"

        elif isinstance(filename, int):
            filename = f"{filename:05d}-{self.name}_coeff_rota.dat"

        elif not isinstance(filename, str):
            raise TypeError(f"Argument must be int or str, not {type(filename)}.")

        print(f"Reading file {filename}... ")
        with open(filename, 'r', encoding='latin-1') as f:
            cont = f.readlines()

        i = 5
        if cont[i][1:3] == 'ar':
            i +=1
        nchim = int(cont[i][1:3])
        for j in range(nchim):
            if j > 13:
                i1 = 1
                k = 5*(j - 14) + 1
            else:
                k = 5*j + 4
                i1 = 0
            pass

        i += 1+i1
        n_rot = int(cont[i][6:11])
        nb_var = int(cont[i][18:20])
        ncoeff = int(cont[i][28:30])

        ecrit = np.zeros((nb_var, n_rot))
        self.X = np.zeros((nchim, n_rot))

        for l in range(n_rot):
            for j in range(nb_var):
                if j%5 == 0:
                    i += 1
                    k = 0

                if old:
                    ecrit[j, l] = float(cont[i][k:k+19])
                    k += 19
                else:
                    ecrit[j, l] = float(cont[i][k:k+22])
                    k += 22

        self.r = ecrit[0, :]
        self.m = ecrit[1, :]
        self.Omega = ecrit[2, :]

        if nb_var > 4:
            self.U2     = ecrit[3, :]
            self.Theta  = ecrit[4, :]
            self.Upsi   = ecrit[5, :]
            self.Lambda = ecrit[6, :]
            self.Deff   = ecrit[8, :]
            self.Dh     = ecrit[9, :]
            self.Dv     = ecrit[10, :]
            self.T      = ecrit[11, :]
            self.rho    = ecrit[12, :]
            self.gradmu = ecrit[13, :]
            self.mu     = ecrit[29,:]
            self.n2t    = ecrit[27,:]
            self.n2mu   = ecrit[28,:]
            self.K      = ecrit[22,:]    # diffusivity
            self.chi    = ecrit[19,:]  # conductivity
            if ncoeff > 30:
                self.eta    = ecrit[30,:]
                self.nu     = ecrit[31,:]
                self.Br     = ecrit[32,:]
                self.Bphi   = ecrit[33,:]
                self.kappa  = ecrit[34,:]
                self.pr     = ecrit[35,:]
                self.delta  = ecrit[36,:]
                self.grad   = ecrit[37,:]
                self.gradad = ecrit[38,:]
                self.cp     = ecrit[39,:]
                self.gamma1 = ecrit[40,:]
                self.lu     = ecrit[41,:]
                self.hp     = ecrit[42,:]
                self.X      = ecrit[43:43+nchim,:]


            self.Theta *=  self.Omega

            phi             = 1.0
            delta           = 1.0
            self.Psi        = (phi*self.Lambda - self.Theta)/delta

            self.V2         = np.zeros(n_rot) # init.

            index           = self.U2 != 0
            U2rad           = self.U2[index]
            U2rads          = savitzky_golay(U2rad,kernel=61)
            self.U2s        = deepcopy(self.U2)
            self.U2s[index] = U2rads
            self.V2[0]      = 0.0

            for l in range(1,n_rot-1):
                self.V2[l]  = self.rho[l+1]*(self.r[l+1]**2)*self.U2[l+1]
                self.V2[l] -= self.rho[l-1]*(self.r[l-1]**2)*self.U2[l-1]
                self.V2[l] /= (self.r[l+1] - self.r[l-1])*\
                              6.0*self.r[l]*self.rho[l]

            l = n_rot-1
            self.V2[l]  = self.rho[l]*(self.r[l]**2)*self.U2[l]
            self.V2[l] -= self.rho[l-1]*(self.r[l-1]**2)*self.U2[l-1]
            self.V2[l] /= (self.r[l] - self.r[l-1])*6.0*self.r[l]*self.rho[l]

            self.V2s = savitzky_golay(self.V2,kernel=61)

            self.U_csi = 0.5*self.rho*self.r**2*self.U2
            self.U_csis = 0.5*self.rho*self.r**2*self.U2s

            ra   = self.r*self.ctes.rsun
            ma   = self.m*self.ctes.msun
            grav = self.ctes.ggrav*ma/ra**2

            alpha        = 1.0 + 0.75*grav*self.Theta/ra/self.Omega**2
            self.Omega2  = 0.2*ra/self.Dh*(2.0*self.V2 - alpha*self.U2)*self.Omega
            self.Omega2s = 0.2*ra/self.Dh*(2.0*self.V2s - alpha*self.U2s)*self.Omega

            self.Krad    = 4.0*aradia*clight*self.T**3/3.0/self.kappa/self.rho**2/self.cp

            self.n2omega = 4.0*self.Omega**2 + 3.0*grav*self.Theta/ra
            self.U2s[self.U2 == 0.0] = 0.0
            self.V2s[self.V2 == 0.0] = 0.0

        else:
            self.rho = ecrit[3, :]


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def j_ang( self, infile='rota', n_ini=None, n_end=None ):
        """
        <CModel.j_ang(self, infile='rota', nn=None)>
        Calculates the evolution of total angular momentum.

        :kparam infile: Kind of infile file. Either `'rota'` (default): `_coeff_rota.dat` file; or
            `'osc'`: osc file.
        :ktype infile: string

        :kparam n_ini: Starts computation from `n_ini`. Default: None.
        :ktype n_ini: int

        :kparam n_end: Stops computation at `n_end`. Default: None.
        :ktype n_end: int
        """
        if infile == 'rota':
            n = sorted( [int(name[:5]) for name in os.listdir('.') if self.name + "_coeff_rota.dat" in name] )
            ni, ne = n[0], n[-1]

        elif infile == 'osc':
            ni, ne = 0, self.nosc-1
        else:
            raise ValueError("infile should be either 'rota' or 'osc'." )

        if n_ini is None:
            n_ini = ni
        if n_end is None:
            n_end = ne

        self.j = np.zeros(n_end - n_ini + 1)
        for i in range(n_ini, n_end+1):
            if infile == 'rota':
                self.read_rot(i)
                r     = self.r
                omega = self.Omega
                m     = self.m

            elif infile == 'osc':
                r     = self.var[i][0][::-1] / self.ctes.rsun
                omega = self.var[i][15][::-1]
                m     = self.var[i][1][::-1]

            ji = 0.5 * r[1:-1]**2 * omega[1:-1] * (m[2:] - m[0:-2])
            self.j[i-n_ini] = ji.sum()
            self.j[i-n_ini] += 0.5 * (r[0]**2*omega[0]*(m[1] - m[0]) + r[-1]**2*omega[-1]*(m[-1] - m[-2]))


        self.j *= 2.0/3.0


    def read_osc_glob( self, filename ):

        """
        Reads the glob part of an .osc file. Can be used to retrieve the value of the gravitational
        constant used to compute the model.

        :param filename: Name of files to be read.
        :type filename: string or list of string

        :return: Global values
        :rtype: array of floats

        :member glob[0]:  $M_\star M_\odot$
        :member glob[1]:  $R_{\rm tot} R_\odot$
        :member glob[2]:  $L_{\rm tot} L_\odot$
        :member glob[3]:  $Z0$
        :member glob[4]:  $X0$
        :member glob[5]:  $\alpha$
        :member glob[6]:  $X$ in CZ
        :member glob[7]:  $Y$ in CZ
        :member glob[8]:  $\partial^2 p / \partial r^2$
        :member glob[9]:  $\partial^2 \rho / \partial r^2$
        :member glob[10]: Age
        :member glob[11]: $w_{\rm rot}$ initial (global rotation velocity)
        :member glob[12]: $w_{\rm rot}$ initial
        :member glob[13]: Gravitational constant
        :member glob[14]: Solar mass
        :member glob[15]: Solar radius
        :member glob[16]: Solar luminosity
        :member glob[17]: Entropy of the adiabat
        :member glob[18]: Mean molecular weight at bottom of the adiabat
        """
        if not os.path.exists(filename):
            raise CESAMError(f"osc file does not exist.")

        if type( filename ) != str and type( filename ) != np.str_:
            raise CESAMError(f"No output specified.")

        with open(filename, 'r', encoding='latin-1') as f:
            cont = f.readlines()

        i          = 1
        osc_header = cont[0:4]
        revision   = osc_header[1][16:31]

        nchim      = int(cont[4][:3])

        # name of chemical elements
        for j in range(nchim):
            if j > 13:
                i1 = 1
                k = 5*(j - 14) + 1
            else:
                k = 5*j + 4
                i1 = 0
            pass

        nglob = int(cont[5+i1][10:20])
        krot  = int(cont[5+i1][40:50])

        globi = np.zeros(nglob)

        i     = 5 + i1
        # global variables
        shift = 19
        j     = 0
        k     = 0
        if (float(cont[i+1][0:0+shift]) < 1e29):
            shift = 24

        for j in range(nglob):
            if j%5 == 0:
                i += 1
                k = 0
            globi[j] = float(cont[i][k:k+shift])
            k += shift
        return globi

    def read_osc( self, filename=None, mods=None, read_vc=True, verbose=None,
        only_glob=False, fast=True, safe=False ):
        """
        Wrapper that reads the structure files. If `nom_output` is set to an HDF5 output, `read_osc` calls
        `read_osc_hdf5`, otherwise, it calls `read_osc_ascii` (except for `-plato.osc` files, it calls
        `read_osc_ascii_plato`).

        :kparam filename: Name of files to be read.
        :ktype filename: string or list of string

        :kparam mods: Index of the models to be read.
        :ktype mods: list of integer

        :kparam read_vc: If True, computes convective quantities (`vconv`, `Fconv`, `Pturb`...).
        :ktype  read_vc: boolean

        :kparam verbose: If True, print details.
        :ktype verbose: boolean

        :kparam only_glob: If True, only reads global parameters, not structure.
        :ktype only_glob: boolean

        :kparam safe: If True, try to read HDF5 records until something is corrupted. Then returns the index of the
            last model that was read properly.
        :ktype safe: boolean

        :member glob: array or list of array
            Global parameters for a model or for all models. Dimension: `nglob (, nmod)`

        :member var: array or list of array
            Structure variables for a model or for all models. Dimension: `nvar, ntot (, nmod)`
            nvar = 22 for adiabatic oscillations

        :member glob[0]:  $M_\star M_\odot$
        :member glob[1]:  $R_{\rm tot} R_\odot$
        :member glob[2]:  $L_{\rm tot} L_\odot$
        :member glob[3]:  $Z0$
        :member glob[4]:  $X0$
        :member glob[5]:  $\alpha$
        :member glob[6]:  $X$ in CZ
        :member glob[7]:  $Y$ in CZ
        :member glob[8]:  $\partial^2 p / \partial r^2$
        :member glob[9]:  $\partial^2 \rho / \partial r^2$
        :member glob[10]: Age
        :member glob[11]: $w_{\rm rot}$ initial (global rotation velocity)
        :member glob[12]: $w_{\rm rot}$ initial
        :member glob[13]: Gravitational constant
        :member glob[14]: Solar mass
        :member glob[15]: Solar radius
        :member glob[16]: Solar luminosity
        :member glob[17]: Entropy of the adiabat
        :member glob[18]: Mean molecular weight at bottom of the adiabat

        :member var[0,i]:   $r R_\odot$
        :member var[1,i]:   $m/M_\star$
        :member var[2,i]:   $T$
        :member var[3,i]:   $P_{\rm tot}$
        :member var[4,i]:   $\rho$
        :member var[5,i]:   Actual gradient $d ln T / d ln P$
        :member var[6,i]:   $L$
        :member var[7,i]:   $\kappa$
        :member var[8,i]:   Thermal + gravitational energy
        :member var[9,i]:   Adiabatic index $\Gamma_1$
        :member var[10,i]:  Adiabatic gradient
        :member var[11,i]:  $\delta$
        :member var[12,i]:  $c_p$
        :member var[13,i]:  $1/mu_{\rm elec.}$
        :member var[14,i]:  Brunt-Vaissala frequency, (0 at centre)
        :member var[15,i]:  Angular velocity [radian/sec]
        :member var[16,i]:  $d \ln \kappa / d \ln T$
        :member var[17,i]:  $d \ln \kappa / d \ln ro$
        :member var[18,i]:  $d \epsilon_{\rm nuc} / d \ln T$
        :member var[19,i]:  $d \epsilon_{\rm nuc} / d \ln ro$
        :member var[20,i]:  $P_{\rm tot}/P_{\rm gas}$
        :member var[21,i]:  Radiative gradient

        :member var[22+j,i]: `xchim[j]*nucleo[j], j=0,nchim-1`
        """
        verbose = self.verbose if verbose is None else verbose
        h5_file = False
        nstop   = None
        if filename is not None:
            h5_file = filename.endswith( 'h5' )
        if self.params.nom_output_.endswith('h5') or h5_file:
            nstop = self.__read_osc_hdf5(filename=filename, mods=mods, read_vc=read_vc, verbose=verbose,
                only_glob=only_glob, safe=safe)
        else:
            if self.params.nom_output_.endswith('_plato'):
                self.__read_osc_ascii_plato(filename=filename, mods=mods, read_vc=read_vc, verbose=verbose,
                    only_glob=only_glob)
            else:
                if fast and cesfort_available:
                    self.__read_osc_ascii_fast(filename=filename, mods=mods, read_vc=read_vc, verbose=verbose,
                        only_glob=only_glob)
                else:
                    self.__read_osc_ascii(filename=filename, mods=mods, read_vc=read_vc, verbose=verbose,
                        only_glob=only_glob)

        return nstop


    def __read_osc_ascii(self, filename=None, mods=None, read_vc=True, verbose=True,
        only_glob=False):

        """
        Reads the file `*.osc`

        :kparam filename: Name of files to be read.
        :ktype filename: string or list of string

        :kparam mods: Index of the models to be read.
        :ktype mods: list of integer

        :kparam read_vc: If True, computes convective quantities (`vconv`, `Fconv`, `Pturb`...).
        :ktype  read_vc: boolean

        :kparam verbose: If True, print details.
        :ktype verbose: boolean

        :kparam only_glob: If True, only reads global parameters, not structure.
        :ktype only_glob: boolean

        :member glob: See `self.read_osc`.
        :mtype glob: array or list of array

        :member var: See `self.read_osc`.
        :mtype var: array or list of array
        """
        if not self.osc:
            if verbose:
                print(f"{RED}Error:{NO_COLOR} no output specified.")
            return

        self.nom_elem = []

        if filename:
            if type(filename) == str:
                files = [filename]
            else:
                files = filename
        elif mods:
            if self.params.nom_output_.endswith(('_adia', '_adia2d', '_adia2dh5')):
                suf = '-ad'
            elif self.params.nom_output_.endswith(('_nadia', '_nadia2d', '_nadia2dh5')):
                suf = '-nad'
            elif self.params.nom_output_.endswith(('_invers', '_invers2d', '_invers2dh5')):
                suf = '-inv'
            else:
                raise CESAMError( "Output not supported: " + self.params.nom_output_ + '.' )
            files = [f'{i:05d}-{self.name}{suf}.osc' for i in mods]
        else:
            files = self.osc if self.all_osc else [self.osc]

        nfile = 0
        self.ntot = []
        self.glob = []
        self.var = []

        nfiles = len(files)
        for i in tqdm( range( nfiles ), desc=f'Reading *{self.name}*.osc files' ):
            fich = files[i]

            with open(fich, 'r', encoding='latin-1') as f:
                cont = f.readlines()
            nfile += 1

            i = 1

            self.osc_header = cont[0:4]
            self.revision = self.osc_header[1][19:26]

            if nfile == 1:
                self.nchim = int(cont[4][:3])

                # name of chemical elements
                for j in range(self.nchim):
                    if j > 13:
                        i1 = 1
                        k = 5*(j - 14) + 1
                    else:
                        k = 5*j + 4
                        i1 = 0
                    self.nom_elem.append(cont[4+i1][k:k+4].strip())

                self.nglob = int(cont[5+i1][10:20])
                self.nvar = int(cont[5+i1][20:30])
                self.krot = int(cont[5+i1][40:50])

            ntoti = int(cont[5+i1][:10])
            globi = np.zeros(self.nglob)
            vari = np.zeros((self.nvar + self.nchim, ntoti))

            i = 5 + i1

            # global variables
            shift = 19
            j=0
            k = 0
            if (float(cont[i+1][0:0+shift]) < 1e29):
                shift = 24

            for j in range(self.nglob):
                if j%5 == 0:
                    i += 1
                    k = 0
                globi[j] = float(cont[i][k:k+shift])
                k += shift

            if not only_glob:

                # variables
                for l in range(ntoti):
                    for j in range(self.nvar + self.nchim):
                        if j%5 == 0:
                            i += 1
                            k = 0

                        try:
                            vari[j, l] = float(cont[i][k:k+shift])
                        except ValueError:
                            vari[j, l] = 0.0

                        k += shift

                # Convert log(m/mstar) -> m/mstar
                index = vari[1, :] > -1.e10
                vari[1, index] = np.exp(vari[1, index])
                index = vari[1, :] < -1.e10
                vari[1, index] = 0.0

            if self.params.nom_output_.startswith('all_'):#len(files) > 1:
                self.ntot.append(ntoti)
                self.glob.append(globi)
                if not only_glob:
                    self.var.append(vari)
            else:
                self.ntot = ntoti
                self.glob = globi
                if not only_glob:
                    self.var = vari
                    if read_vc: self.calc_vconv('osc', var=vari, glob=globi)

        #-------------- finished looping over files

        self.nom_glob = [r'$M_{\star}$', r'$R_{\star}$', r'$L_{\star}$',
                         r'$Z_0$', r'$X_0$', r'$\alpha$', 'X in CZ',
                         'Y in CZ', 'd2p', 'd2ro', 'age',
                         r'$w_{\rm rot}$ initial (global rotation velocity)',
                         r'$w_{\rm rot}$ initial', 'Gravitional constant',
                         r'$M_{\odot}$', r'$R_{\odot}$', r'$L_{\odot}$',
                         r'$S_{\rm ad}$', r'$\mu_{\rm ad}$']

        if not only_glob:
            # Names of Vars...
            self.nom_vars = [r'$r$', r'$m/M_{\star}$', r'$T$',
                             r'$P_{\rm tot}$', r'$\rho$', r'$\nabla$', r'$L$',
                             r'$\kappa$', r'$\varepsilon$', r'$\Gamma_1$',
                             r'$\nabla_{\rm ad}$', r'$\delta$', r'$c_p$',
                             r'$1/\mu_e$', r'$\mathcal{A}$', r'$\Omega$',
                             r'$\kappa_T$', r'$\kappa_{\rho}$',
                             r'$\varepsilon_T$', r'$\varepsilon_{\rho}$']
            self.unom_vars = ['r',                  'm/M',             'T',             'P',
                             u'\u03c1',            u'\u2207',          'L',            u'\u03ba',
                             u'\u03b5',            u'\u0393\u2081',   u'\u2207\u2090', u'\u03b4',
                              'Cp',                u'1/\u03bc\u2091', u'A', u'\u03a9', u'\u03baT',
                             u'\u03ba\u03c1',      u'\u03b5T',        u'\u03b5\u03c1']

            if self.params.cpturb != 0.0:
                self.nom_vars += [r'$P_{\rm tot}/P_{\rm gas}$']
                self.unom_vars += ['P/Pgas']
            else:
                self.nom_vars += [r'$\nabla_{\rm \mu}$']
                self.unom_vars += [u'\u2207\u03bc']

            self.nom_vars += [r'$\nabla_{\rm rad}$']
            self.unom_vars += [u'\u2207\u1d63\u2090']


            if '_invers' in self.params.nom_output_:
                self.nom_vars += [r'$\partial \Gamma_1/\partial \ln P$',
                                  r'$\partial \Gamma_1/\partial \ln T$',
                                  r'$\partial \Gamma_1/\partial Y$']

                self.unom_vars += [u'\u2202\u0393\u2081/\u2202lnP',
                                   u'\u2202\u0393\u2081/\u2202lnT',
                                   u'\u2202\u0393\u2081/\u2202Y']

            elif '_nadia' in self.params.nom_output_:

                self.nom_vars += [r'$\partial \Gamma_1/\partial \ln P$',
                                  r'$\partial \Gamma_1/\partial \ln T$',
                                  r'$\partial \Gamma_1/\partial Y$',
                                  r'$(\partial P/\partial \rho)_{T,X}$',
                                  r'$(\partial P/\partial T)_{\rho,X}$',
                                  r'$(\partial P/\partial X)_{T,\rho}$',
                                  r'$(\partial u/\partial \rho)_{T,X}$',
                                  r'$(\partial u/\partial T)_{\rho,X}$',
                                  r'$(\partial u/\partial X)_{T,\rho}$',
                                  r'$u$',
                                  r'$\left(\partial^2 P/\partial \rho^2\right)_{T,X}$',
                                  r'$\left(\partial^2 P/\partial \rho \partial T \right)_X$',
                                  r'$\left(\partial^2 P/\partial T^2\right)_{\rho,X}$',
                                  r'$\left(\partial^2 u/\partial \rho^2\right)_{T,X}$',
                                  r'$\left(\partial^2 u/\partial \rho \partial T \right)_X$',
                                  r'$\left(\partial^2 u/\partial T^2\right)_{\rho,X}$',
                                  r'${\rm d} K/{\rm d}X$',
                                  r'${\rm d}^2 K/{\rm d}T^2$',
                                  r'${\rm d}\varepsilon /{\rm d}X$',
                                  r'${\rm d} X/{\rm d}r$',
                                  r'$J-B$',r'$\Gamma$']


                self.unom_vars += [u'\u2202\u0393\u2081/\u2202lnP',
                                   u'\u2202\u0393\u2081/\u2202lnT',
                                   u'\u2202\u0393\u2081/\u2202Y',
                                   u'(\u2202P/\u2202\u03c1)T,X',
                                   u'(\u2202P/\u2202T)\u03c1,X',
                                   u'(\u2202P/\u2202X)T,\u03c1',
                                   u'(\u2202u/\u2202\u03c1)T,X',
                                   u'(\u2202u/\u2202T)\u03c1,X',
                                   u'(\u2202u/\u2202X)T,\u03c1', 'u',
                                   u'(\u2202\u00b2P/\u2202\u03c1\u00b2)T,X',
                                   u'(\u2202\u00b2P/\u2202\u03c1\u2202T)X',
                                   u'(\u2202\u00b2P/\u2202T\u00b2)\u03c1,X',
                                   u'(\u2202\u00b2u/\u2202\u03c1\u00b2)T,X',
                                   u'(\u2202\u00b2u/\u2202\u03c1\u2202T)X',
                                   u'(\u2202\u00b2u/\u2202T\u00b2)\u03c1,X',
                                   'dK/dX', u'd\u00b2K/dT\u00b2',
                                   u'd\u03b5/dX', 'dX/dr$', 'J-B',u'\u0393']
                if self.nvar > 44:
                  self.nom_vars += [r'$\varepsilon_g$']
                  self.unom_vars += [u'-T dS/dt']

            nel_for = []

            for i in self.nom_elem:
                nel_for.append(r'$X_{\rm %s}$' % i)

            self.nom_vars += nel_for
            self.unom_vars += self.nom_elem


    def __read_osc_ascii_fast(self, filename=None, mods=None, read_vc=True, verbose=True,
        only_glob=False):

        """
        Reads the file *.osc

        :kparam filename: Name of files to be read.
        :ktype filename: string or list of string

        :kparam mods: Index of the models to be read.
        :ktype mods: list of integer

        :kparam read_vc: If True, computes convective quantities (`vconv`, `Fconv`, `Pturb`...).
        :ktype  read_vc: boolean

        :kparam verbose: If True, print details.
        :ktype verbose: boolean

        :kparam only_glob: If True, only reads global parameters, not structure.
        :ktype only_glob: boolean

        :member glob: See `self.read_osc`
        :mtype glob: array or list of array

        :member var: See `self.read_osc`
        :mtype var: array or list of array
        """
        if not self.osc:
            if verbose:
                print(f"{RED}Error:{NO_COLOR} no output specified.")
            return

        if filename:
            if type(filename) == str:
                files = [filename]
            else:
                files = filename
        elif type( mods ) == list or type( mods ) == np.ndarray:
            if self.params.nom_output_.endswith(('_adia', '_adia2d', '_adia2dh5')):
                suf = '-ad'
            elif self.params.nom_output_.endswith(('_nadia', '_nadia2d', '_nadia2dh5')):
                suf = '-nad'
            elif self.params.nom_output_.endswith(('_invers', '_invers2d', '_invers2dh5')):
                suf = '-inv'
            else:
                raise CESAMError( "Output not supported: " + self.params.nom_output_ + '.' )
            files = [f'{i:05d}-{self.name}{suf}.osc' for i in mods]
        else:
            files = self.osc if self.all_osc else [self.osc]

        self.ntot       = []
        self.glob       = []
        self.var        = []
        self.osc_header = []

        nfiles = len(files)
        for ifile in tqdm( range( nfiles ), desc=f'Reading *{self.name}*.osc files' ):
            file = files[ifile]

            if ifile == 0:
                with open(file, 'r', encoding='latin-1') as f:
                    # read 4 lines of 95 characters (= 380 chars)
                    for i in range(4):
                        self.osc_header.append( f.readline() )
                    self.revision = self.osc_header[1][19:26]

                    line = f.readline()
                    self.nchim = int( line[:3] )
                    self.nom_elem = list( map( str.strip, line.split()[1:] ) )

                    if (self.nchim-1) // 14 > 0:
                        for i in range( (self.nchim-1) // 14 ):
                            self.nom_elem.extend( list( map( str.strip, f.readline().split() ) ) )

                    words = f.readline().split()
                    self.nglob = int( words[1] )
                    self.nvar  = int( words[2] )
                    self.krot  = int( words[3] )
                    ntoti = int( words[0] )
                    size = (len( f.readline() )-1) // 5 + 1

            cesam2k20_fortran.read_osc( file, only_glob, size-1 )
            ntoti = int( cesam2k20_fortran.itot )
            if not only_glob:
                vari  = np.copy( cesam2k20_fortran.var )
                vari[1,:] = np.where( vari[1,:] > -1.e10, np.exp(vari[1,:]), 0.0)

            if self.params.nom_output_.startswith('all_'): #nfiles > 1:
                self.ntot.append(ntoti)
                self.glob.append( np.copy( cesam2k20_fortran.glob ) )
                if not only_glob:
                    self.var.append( vari )
            else:
                self.ntot = ntoti
                self.glob = np.copy( cesam2k20_fortran.glob )
                if not only_glob:
                    self.var = vari
                    if read_vc: self.calc_vconv('osc', var=self.var, glob=self.glob)

        #-------------- finished looping over files

        self.nom_glob = [r'$M_{\star}$', r'$R_{\star}$', r'$L_{\star}$',
                         r'$Z_0$', r'$X_0$', r'$\alpha$', 'X in CZ',
                         'Y in CZ', 'd2p', 'd2ro', 'age',
                         'wrot initial (global rotation velocity)',
                         'w_rot initial']

        if not only_glob:
            # Names of Vars...
            self.nom_vars = [r'$r$', r'$m/M_{\star}$', r'$T$',
                             r'$P_{\rm tot}$', r'$\rho$', r'$\nabla$', r'$L$',
                             r'$\kappa$', r'$\varepsilon$', r'$\Gamma_1$',
                             r'$\nabla_{\rm ad}$', r'$\delta$', r'$c_p$',
                             r'$1/\mu_e$', r'$\mathcal{A}$', r'$\Omega$',
                             r'$\kappa_T$', r'$\kappa_{\rho}$',
                             r'$\varepsilon_T$', r'$\varepsilon_{\rho}$']
            self.unom_vars = ['r',                  'm/M',             'T',             'P',
                             u'\u03c1',            u'\u2207',          'L',            u'\u03ba',
                             u'\u03b5',            u'\u0393\u2081',   u'\u2207\u2090', u'\u03b4',
                              'Cp',                u'1/\u03bc\u2091', u'A', u'\u03a9', u'\u03baT',
                             u'\u03ba\u03c1',      u'\u03b5T',        u'\u03b5\u03c1']

            if self.params.cpturb != 0.0:
                self.nom_vars += [r'$P_{\rm tot}/P_{\rm gas}$']
                self.unom_vars += ['P/Pgas']
            else:
                self.nom_vars += [r'$\nabla_{\rm \mu}$']
                self.unom_vars += [u'\u2207\u03bc']

            self.nom_vars += [r'$\nabla_{\rm rad}$']
            self.unom_vars += [u'\u2207\u1d63\u2090']

            if '_invers' in self.params.nom_output_:
                self.nom_vars += [r'$\partial \Gamma_1/\partial \ln P$',
                                  r'$\partial \Gamma_1/\partial \ln T$',
                                  r'$\partial \Gamma_1/\partial Y$']

                self.unom_vars += [u'\u2202\u0393\u2081/\u2202lnP',
                                   u'\u2202\u0393\u2081/\u2202lnT',
                                   u'\u2202\u0393\u2081/\u2202Y']

            elif '_nadia' in self.params.nom_output_:

                self.nom_vars += [r'$\partial \Gamma_1/\partial \ln P$',
                                  r'$\partial \Gamma_1/\partial \ln T$',
                                  r'$\partial \Gamma_1/\partial Y$',
                                  r'$(\partial P/\partial \rho)_{T,X}$',
                                  r'$(\partial P/\partial T)_{\rho,X}$',
                                  r'$(\partial P/\partial X)_{T,\rho}$',
                                  r'$(\partial u/\partial \rho)_{T,X}$',
                                  r'$(\partial u/\partial T)_{\rho,X}$',
                                  r'$(\partial u/\partial X)_{T,\rho}$',
                                  r'$u$',
                                  r'$\left(\partial^2 P/\partial \rho^2\right)_{T,X}$',
                                  r'$\left(\partial^2 P/\partial \rho \partial T \right)_X$',
                                  r'$\left(\partial^2 P/\partial T^2\right)_{\rho,X}$',
                                  r'$\left(\partial^2 u/\partial \rho^2\right)_{T,X}$',
                                  r'$\left(\partial^2 u/\partial \rho \partial T \right)_X$',
                                  r'$\left(\partial^2 u/\partial T^2\right)_{\rho,X}$',
                                  r'${\rm d} K/{\rm d}X$',
                                  r'${\rm d}^2 K/{\rm d}T^2$',
                                  r'${\rm d}\varepsilon /{\rm d}X$',
                                  r'${\rm d} X/{\rm d}r$',
                                  r'$J-B$',r'$\Gamma$']


                self.unom_vars += [u'\u2202\u0393\u2081/\u2202lnP',
                                   u'\u2202\u0393\u2081/\u2202lnT',
                                   u'\u2202\u0393\u2081/\u2202Y',
                                   u'(\u2202P/\u2202\u03c1)T,X',
                                   u'(\u2202P/\u2202T)\u03c1,X',
                                   u'(\u2202P/\u2202X)T,\u03c1',
                                   u'(\u2202u/\u2202\u03c1)T,X',
                                   u'(\u2202u/\u2202T)\u03c1,X',
                                   u'(\u2202u/\u2202X)T,\u03c1', 'u',
                                   u'(\u2202\u00b2P/\u2202\u03c1\u00b2)T,X',
                                   u'(\u2202\u00b2P/\u2202\u03c1\u2202T)X',
                                   u'(\u2202\u00b2P/\u2202T\u00b2)\u03c1,X',
                                   u'(\u2202\u00b2u/\u2202\u03c1\u00b2)T,X',
                                   u'(\u2202\u00b2u/\u2202\u03c1\u2202T)X',
                                   u'(\u2202\u00b2u/\u2202T\u00b2)\u03c1,X',
                                   'dK/dX', u'd\u00b2K/dT\u00b2',
                                   u'd\u03b5/dX', 'dX/dr$', 'J-B',u'\u0393']
                if self.nvar > 44:
                  self.nom_vars += [r'$\varepsilon_g$']
                  self.unom_vars += [u'-T dS/dt']

            nel_for = []

            for i in self.nom_elem:
                nel_for.append(r'$X_{\rm %s}$' % i)

            self.nom_vars += nel_for
            self.unom_vars += self.nom_elem


    def __read_osc_ascii_plato(self, filename=None, mods=None, read_vc=True, verbose=True,
        only_glob=False):
        """
        Reads the file *.osc

        :kparam filename: Name of files to be read.
        :ktype filename: string or list of string

        :kparam mods: Index of the models to be read.
        :ktype mods: list of integer

        :kparam read_vc: If True, computes convective quantities (`vconv`, `Fconv`, `Pturb`...).
        :ktype  read_vc: boolean

        :kparam verbose: If True, print details.
        :ktype verbose: boolean

        :kparam only_glob: If True, only reads global parameters, not structure.
        :ktype only_glob: boolean

        :member glob: See `self.read_osc`
        :mtype glob: array or list of array

        :member var: Structure variables for a model or for all models.
        :mtype var: array or list of array

        :member var[0,i]:   $r R_\odot$
        :member var[1,i]:   $m/M_\star$
        :member var[2,i]:   $P_{\rm tot}$
        :member var[3,i]:   $\rho$
        :member var[4,i]:   Adiabatic index $\Gamma_1$
        :member var[5,i]:  Brunt-Vaissala frequency, (0 at centre)
        :member var[6,i]:  $Y$
        :member var[7,i]:  $Z$
        """
        if not self.osc:
            if verbose:
                print(f"{RED}Error:{NO_COLOR} no output specified.")
            return

        self.nom_elem = []

        if filename:
            if type(filename) == str:
                files = [filename]
            else:
                files = filename
        elif mods:
            suf = '-plato'
            files = [f'{i:05d}-{self.name}{suf}.osc' for i in mods]
        else:
            files = self.osc if self.all_osc else [self.osc]

        nfile = 0
        self.ntot = []
        self.glob = []
        self.var  = []

        nfiles = len(files)
        for i in tqdm( range( nfiles ), desc=f'Reading *{self.name}*.osc files' ):
            fich = files[i]

            with open(fich, 'r', encoding='latin-1') as f:
                cont = f.readlines()
            nfile += 1

            i = 1

            self.osc_header = cont[0:4]
            self.revision = self.osc_header[1][19:26]

            if nfile == 1:
                self.nchim = int(cont[4][:3])

                # name of chemical elements
                for j in range(self.nchim):
                    if j > 13:
                        i1 = 1
                        k = 5*(j - 14)
                    else:
                        k = 5*j + 4
                        i1 = 0
                    self.nom_elem.append(cont[4+i1][k:k+4].strip())

                self.nglob = int(cont[5+i1][10:20])
                self.nvar = int(cont[5+i1][20:30])
                self.krot = int(cont[5+i1][40:50])

            ntoti = int(cont[5+i1][:10])
            globi = np.zeros(self.nglob)
            vari = np.zeros((self.nvar + self.nchim, ntoti))

            i = 5 + i1

            # global variables
            shift = 19
            j=0
            k = 0
            if (float(cont[i+1][0:0+shift]) < 1e33):
                shift = 24

            for j in range(self.nglob):
                if j%5 == 0:
                    i += 1
                    k = 0
                globi[j] = float(cont[i][k:k+shift])
                k += shift

            # variables
            for l in range(ntoti):
                for j in range(self.nvar):
                    if j%5 == 0:
                        i += 1
                        k = 0

                    try:
                        vari[j, l] = float(cont[i][k:k+shift])
                    except ValueError:
                        vari[j, l] = 0.0

                    k += shift

            # Convert log(m/mstar) -> m/mstar
            index = vari[1, :] > -1.e10
            vari[1, index] = np.exp(vari[1, index])
            index = vari[1, :] < -1.e10
            vari[1, index] = 0.0

            if self.params.nom_output_.startswith('all_'): #len(files) > 1:
                self.ntot.append(ntoti)
                self.glob.append(globi)
                self.var.append(vari)
            else:
                self.ntot = ntoti
                self.glob = globi
                self.var = vari


        #-------------- finished looping over files

        # Names of Vars.
        self.nom_vars = [r'$r$', r'$m/M_{\star}$', r'$P_{\rm tot}$', r'$\rho$', r'$\Gamma_1$', r'$\mathcal{A}$', '$Y$', '$Z$']
        self.unom_vars = ['r', 'm/M', 'P', u'\u03c1', u'\u0393\u2081', u'A', 'Y', 'Z']

    def __read_osc_hdf5(self, filename=None, mods=None, read_vc=True, verbose=True,
        only_glob=False, safe=False ):
        """
        Reads the file *.osch5.

        :kparam filename: Name of files to be read.
        :ktype filename: string or list of string

        :kparam mods: Index of the models to be read.
        :ktype mods: list of integer

        :kparam read_vc: If True, computes convective quantities (`vconv`, `Fconv`, `Pturb`...).
        :ktype  read_vc: boolean

        :kparam verbose: If True, print details.
        :ktype verbose: boolean

        :kparam only_glob: If True, only reads global parameters, not structure.
        :ktype only_glob: boolean

        :member glob: array or list of array
            See `self.read_osc`

        :member var: array or list of array
            See `self.read_osc`
        """
        steps    = []

        if filename is not None:
            file = filename
        elif mods is not None:
            if self.params.nom_output_.endswith('_adiah5'):
                suf = '-ad'
            elif self.params.nom_output_.endswith('_nadiah5'):
                suf = '-nad'
            elif self.params.nom_output_.endswith('_platoh5'):
                suf = '-plato'
            else:
                suf = '-inv'
            # steps = [f'/{i:05d}-{self.name}{suf}.osch5' for i in mods]
            file = self.osch5
        else:
            file = self.osch5

        self.ntot  = []
        self.glob  = []
        self.var   = []
        if verbose: print(f'Reading {file}...', end='')


        with h5py.File(file, 'r') as f:
            # if not steps:
            steps = []
            for s in list(f.keys()):
                if mods is None or int( s[:5] ) in mods:
                    steps.append( s )

            nstop = -1
            for i, step in enumerate( steps ):
                try:
                    self.revision = f[f'{step}/header/header'][1][19:26]
                    self.nchim = f[f'{step}/header'].attrs['nchim']
                    self.nom_elem = np.array(f[f'{step}/header/elem_name'][:])

                    self.nglob = f[f'{step}/glob'].attrs['nglob']
                    self.nvar  = f[f'{step}/var'].attrs['nvar']
                    totvar  = f[f'{step}/var'].attrs['totvar']
                    self.krot  = f[f'{step}/glob'].attrs['Krot']

                    ntoti = f[f'{step}/var'].attrs['itot']
                    globi = f[f'{step}/glob/glob'][:]
                    vari  = f[f'{step}/var/var'][:]

                    # Convert log(m/mstar) -> m/mstar
                    tmp_var1 = vari[1, :]
                    vari[1, :] = np.where( tmp_var1 > -1.e10, np.exp(tmp_var1), 0.0 )


                    self.ntot.append( ntoti )
                    self.glob.append( globi )
                    self.var.append(  vari )

                    if read_vc and not self.params.nom_output_.endswith('_platoh5'): self.calc_vconv('osc', var=vari, glob=globi)
                except:
                    break

                nstop = i-1

            self.nom_vars  = [s.decode().replace('\\\\', '\\') for s in f.attrs['latex_var_name']]
            self.nom_glob  = [s.decode().replace('\\\\', '\\') for s in f.attrs['latex_glob_name']]
            self.unom_vars = [s.decode().replace('\\\\', '\\') for s in f.attrs['unicode_var_name']]

        if self.params.nom_output_.startswith('osc_'): #len( steps ) == 1:
            self.var  = self.var[0]
            self.glob = self.glob[0]
            self.ntot = self.ntot[0]
        self.nom_elem  = [s.decode() for s in self.nom_elem]
        if ' n' in self.nom_elem: self.nom_elem.remove( ' n')

        self.ntot  = np.array( self.ntot )
        self.glob  = np.array( self.glob )

        nel_for = []
        for elem in self.nom_elem:
                nel_for.append(r'$X_{\rm %s}$' % (elem))

        if not self.params.nom_output_.endswith('_platoh5'):
            self.nom_vars += nel_for
            self.unom_vars += self.nom_elem

        if verbose:
            print(f"{GREEN}[Done]{NO_COLOR}")

        if safe:
            return nstop
        else:
            return None



    def __read_record( self, f_osc):
        """
        Private method that reads a record from a fortran binary file.

        :param f_osc: Instance of an opened fortran binary file.
        :type f_osc: `_io.TextIOWrapper`
        """
        h = np.fromfile(f_osc, dtype='i4', count=2)[1]
        return np.fromfile(f_osc, dtype='float64', count=int(h / 8))

    def read_osc2d( self, filename, verbose=None ):
        """
        <CModel.read_osc2d( filename )>

        Reads a osc2d file. For most 2D arrays, first dimension is the angular sector, 2nd dimension is the radius.
        For the angle, index goes from `0` to `Nt` (pole to equator).

        :param filename   : Name of files to be read.
        :type filename    : string

        :kparam verbose   : If True, print out some comments
        :ktype verbose    : bool

        :member nph       : Index of layer corresponding to the surface
        :mtype nph        : int

        :member Nr        : Number of layers in the total structure (interior + atmosphere)
        :mtype Nr         : int

        :member Nt        : number of angular sectors
        :mtype Nt         : int

        :member ml        : number of legendre poly in the decomposition
        :mtype ml         : int

        :member ith       : index of characteristic angle
        :mtype ith        : int

        :member cte_p     : normalization cste for pressure ($p_{\rm file} = p_{\rm mod} cte_p$)
        :mtype cte_p      : float

        :member cte_rho   : normalization cste for density ($\rho_{\rm file} = \rho_{\rm mod} cte_\rho$)
        :mtype cte_rho    : float

        :member omega_k   : break up velocity, normalization cste for omega
            ($omega_{\rm file} = omega_{\rm mod} / omega_{\rm k}$)
        :mtype omega_k    : float

        :member Rph       : radius at the surface (at equator)
        :mtype Rph        : float

        :member Req       : radius with atm (at equator)
        :mtype Req        : float

        :member M         : total mass
        :mtype M          : float

        :member G         : Universal gravity cste
        :mtype G          : float


        :member theta     : angles of sectors
        :mtype theta      : np.ndarray

        :member weight    : weight for legendre quadrature
        :mtype weight     : np.ndarray

        :member cos_theta : cos of angles
        :mtype cos_theta  : np.ndarray


        :member r2d       : radius of characteristics (spherical near center) (shape: `(Nt, Nr)`). In unit of `Req`
        :mtype r2d        : np.ndarray

        :member pr2d      : pressure on characteristics (shape: `(Nt, Nr)`)
        :mtype pr2d       : np.ndarray

        :member rho2d     : density on characteristics (shape: `(Nt, Nr)`)
        :mtype rho2d      : np.ndarray

        :member omega2d   : rotation on characteristics (shape: `(Nt, Nr)`)
        :mtype omega2d    : np.ndarray

        :member gamma12d  : $Gamma_1$ on characteristics (shape: `(Nt, Nr)`)
        :mtype gamma12d   : np.ndarray


        :member drdz      : $d r2d / d \zeta$ (zeta is radius along characteristic angles) (shape: `(Nt, Nr)`)
        :mtype drdz       : np.ndarray

        :member dpdr      : $dp / dr$  (shape: `(Nt, Nr)`)
        :mtype dpdr       : np.ndarray

        :member drhodr    : $d\rho / dr$ (shape: `(Nt, Nr)`)
        :mtype drhodr     : np.ndarray

        :member domegadr  : $d\Omega / dr$ (shape: `(Nt, Nr)`)
        :mtype domegadr   : np.ndarray


        :member drdcos    : $dr / d \cos \theta$
        :mtype drdcos     : np.ndarray

        :member dpdcos    : $dp / d \cos \theta$
        :mtype dpdcos     : np.ndarray

        :member drhodcos  : $drho / d \cos \theta$
        :mtype drhodcos   : np.ndarray

        :member domegadcos: $d\Omega / d \cos \theta$
        :mtype domegadcos : np.ndarray


        :member d2rdz2    : $d^2 r / d \zeta^2$
        :mtype d2rdz2     : np.ndarray

        :member d2rdcosdz : $d^2 r / d \cos \theta / d \zeta$
        :mtype d2rdcosdz  : np.ndarray


        :member d2rdcos2  : $d^2 r / d \cos^2 \theta$
        :mtype d2rdcos2   : np.ndarray

        :member d2pdcos2  : $d^2 p / d \cos^2 \theta$
        :mtype d2pdcos2   : np.ndarray

        :member d2rhodcos2: $d^2 \rho / d \cos^2 \theta$
        :mtype d2rhodcos2 : np.ndarray


        :member Ar        : $(dp/dr) / \Gamma_{1, \rm 2D} / pr2d - (d\rho/dr) / rho2d$
        :mtype Ar         : np.ndarray

        :member m         : Mass inside characteristic shells
        :mtype m          : np.ndarray


        :member rho_l        : component of the decomposition of rho over ISOBARS (only if option `extended_osc2d` was activated)
        :mtype rho_l         : np.ndarray

        :member fp           : factor $f_p$ (only if option `extended_osc2d` was activated)
        :mtype fp            : np.ndarray

        :member ft           : factor $f_t$ (only if option `extended_osc2d` was activated)
        :mtype ft            : np.ndarray

        :member fd           : factor $f_d$ (only if option `extended_osc2d` was activated)
        :mtype fd            : np.ndarray

        :member Sp           : Surface of isobars (only if option `extended_osc2d` was activated)
        :mtype Sp            : np.ndarray

        :member geffav[0, :] : $<g_{\rm eff}>$ (only if option `extended_osc2d` was activated)
        :mtype geffav[0, :]  : np.ndarray

        :member geffav[1, :] : $<g_{\rm eff}^{-1}>$ (only if option `extended_osc2d` was activated)
        :mtype geffav[1, :]  : np.ndarray

        :member geffav[2, :] : $< g_{\rm eff}^{-1} r^2 \sin^2 \theta>$ (only if option `extended_osc2d` was activated)
        :mtype geffav[2, :]  : np.ndarray

        :member geff2d[0, :, :]: $|g_{\rm eff}|$ (only if option `extended_osc2d` was activated)
        :mtype geff2d[0, :, :] : np.ndarray

        :member geff2d[1, :, :]: $g_{\rm eff,r}$ (only if option `extended_osc2d` was activated)
        :mtype geff2d[1, :, :] : np.ndarray

        :member geff2d[2, :, :]: $g_{\rm eff, \theta}$ (only if option `extended_osc2d` was activated)
        :mtype geff2d[2, :, :] : np.ndarray

        :member theta_m      : $\theta_m$ as a function of mass
        :mtype theta_m       : np.ndarray


        :member U_l     : Legendre components of the vertical merid. circu. (only if transport by Mathis & Zahn 2004) in 2D is used)
        :mtype U_l      : nd.array

        :member V_l     : Legendre components of the horizontal merid. circu (only if transport by Mathis & Zahn 2004) in 2D is used)
        :mtype V_l      : nd.array

        :member Ups_l   : Legendre components of Upsilon (see doc.) (only if transport by Mathis & Zahn 2004) in 2D is used)
        :mtype Ups_l    : nd.array

        :member psi_l   : Legendre components of Psi ($T_l / T$) (only if transport by Mathis & Zahn 2004) in 2D is used)
        :mtype psi_l    : nd.array

        :member lambda_l: Legendre components of Lambda ($\mu_l / \mu$) (only if transport by Mathis & Zahn 2004) in 2D is used)
        :mtype lambda_l : nd.array

        :member deff    : Effective diff. coeff for chemicals (e.g. Maeder 2003) (only if transport by Mathis & Zahn 2004) in 2D is used)
        :mtype deff     : nd.array

        :member dh      : horizontal shear-induced diffusion coefficient (only if transport by Mathis & Zahn 2004) in 2D is used)
        :mtype dh       : nd.array

        :member dv      : vertical shear-induced diffusion coefficient (only if transport by Mathis & Zahn 2004) in 2D is used)
        :mtype dv       : nd.array


        :member omega_c: cut-off frequency (in units of `omega_k`)
        :member omega_c_Hz : cut-off frequency in microHz
        """

        verbose = self.verbose if verbose is None else verbose

        if verbose: print(f'Reading {filename}...', end='')
        with open(filename, 'rb') as f_osc:
            header = np.fromfile(f_osc, dtype='i4', count=1)
            self.nph, self.Nr, self.Nt, self.ml, self.ith, self.Krot, self.isph = np.fromfile(f_osc, dtype='i4', count=7)
            header = np.fromfile(f_osc, dtype='i4', count=2)
            self.cte_p, self.cte_rho, self.omega_k, self.Rph, self.Req, self.M, self.G = np.fromfile(f_osc, dtype='float64', count=7)

            self.theta      = self.__read_record( f_osc )
            self.weight     = self.__read_record( f_osc )
            self.cos_theta  = self.__read_record( f_osc )

            self.r2d        = np.reshape( self.__read_record( f_osc ), (self.Nt, self.Nr) )
            self.pr2d       = np.reshape( self.__read_record( f_osc ), (self.Nt, self.Nr) )
            self.rho2d      = np.reshape( self.__read_record( f_osc ), (self.Nt, self.Nr) )
            self.omega2d    = np.reshape( self.__read_record( f_osc ), (self.Nt, self.Nr) )
            self.gamma12d   = np.reshape( self.__read_record( f_osc ), (self.Nt, self.Nr) )

            self.drdz       = np.reshape( self.__read_record( f_osc ), (self.Nt, self.Nr) )
            self.dpdr       = np.reshape( self.__read_record( f_osc ), (self.Nt, self.Nr) )
            self.drhodr     = np.reshape( self.__read_record( f_osc ), (self.Nt, self.Nr) )
            self.domegadr   = np.reshape( self.__read_record( f_osc ), (self.Nt, self.Nr) )

            self.drdcos     = np.reshape( self.__read_record( f_osc ), (self.Nt, self.Nr) )
            self.dpdcos     = np.reshape( self.__read_record( f_osc ), (self.Nt, self.Nr) )
            self.drhodcos   = np.reshape( self.__read_record( f_osc ), (self.Nt, self.Nr) )
            self.domegadcos = np.reshape( self.__read_record( f_osc ), (self.Nt, self.Nr) )

            self.d2rdz2     = np.reshape( self.__read_record( f_osc ), (self.Nt, self.Nr) )
            self.d2rdcosdz  = np.reshape( self.__read_record( f_osc ), (self.Nt, self.Nr) )

            self.d2rdcos2   = np.reshape( self.__read_record( f_osc ), (self.Nt, self.Nr) )
            self.d2pdcos2   = np.reshape( self.__read_record( f_osc ), (self.Nt, self.Nr) )
            self.d2rhodcos2 = np.reshape( self.__read_record( f_osc ), (self.Nt, self.Nr) )

            self.Ar         = np.reshape( self.__read_record( f_osc ), (self.Nt, self.Nr) )
            self.m          = self.__read_record( f_osc )

            header = np.fromfile(f_osc, dtype='i4', count=2)
            iext, iext_rota = np.fromfile(f_osc, dtype='i4', count=2)

            if iext:
                self.rho_l         = np.reshape( self.__read_record( f_osc ), (self.ml, self.Nr) )
                self.fp            = self.__read_record( f_osc )
                self.ft            = self.__read_record( f_osc )
                self.fd            = self.__read_record( f_osc )
                self.Sp            = self.__read_record( f_osc )

                self.geffav        = np.zeros((3, self.Nr))
                self.geffav[0, :]  = self.__read_record( f_osc )
                self.geffav[1, :]  = self.__read_record( f_osc )
                self.geffav[2, :]  = self.__read_record( f_osc )

                self.geff2d        = np.zeros((3, self.Nt, self.Nr))
                self.geff2d[0,:,:] = np.reshape( self.__read_record( f_osc ), (self.Nt, self.Nr) )
                self.geff2d[1,:,:] = np.reshape( self.__read_record( f_osc ), (self.Nt, self.Nr) )
                self.geff2d[2,:,:] = np.reshape( self.__read_record( f_osc ), (self.Nt, self.Nr) )
                self.theta_m       = self.__read_record( f_osc )

            if iext_rota:
                if self.Krot == 6:
                    ml = self.ml
                else:
                    ml = 1
                self.U_l      = np.reshape( self.__read_record( f_osc ), (ml, self.Nr) )
                self.V_l      = np.reshape( self.__read_record( f_osc ), (ml, self.Nr) )
                self.Ups_l    = np.reshape( self.__read_record( f_osc ), (ml, self.Nr) )
                self.psi_l    = np.reshape( self.__read_record( f_osc ), (ml, self.Nr) )
                self.lambda_l = np.reshape( self.__read_record( f_osc ), (ml, self.Nr) )
                self.deff     = self.__read_record( f_osc )
                self.dh       = self.__read_record( f_osc )
                self.dv       = self.__read_record( f_osc )


        # we set right value to ith
        ind_sorted     = self.theta.argsort()
        self.ith       = np.where(ind_sorted == 0)[0][0]
        self.theta     = self.theta[ind_sorted]
        self.weight    = self.weight[ind_sorted]
        self.cos_theta = self.cos_theta[ind_sorted]

        for i in range(self.Nr):
            self.r2d[:,         i] = self.r2d[       :, i][ind_sorted]
            self.pr2d[:,        i] = self.pr2d[      :, i][ind_sorted]
            self.rho2d[:,       i] = self.rho2d[     :, i][ind_sorted]
            self.omega2d[:,     i] = self.omega2d[   :, i][ind_sorted]
            self.gamma12d[:,    i] = self.gamma12d[  :, i][ind_sorted]
            self.drdz[:,        i] = self.drdz[      :, i][ind_sorted]
            self.dpdr[:,        i] = self.dpdr[      :, i][ind_sorted]
            self.drhodr[:,      i] = self.drhodr[    :, i][ind_sorted]
            self.domegadr[:,    i] = self.domegadr[  :, i][ind_sorted]
            self.drdcos[:,      i] = self.drdcos[    :, i][ind_sorted]
            self.dpdcos[:,      i] = self.dpdcos[    :, i][ind_sorted]
            self.drhodcos[:,    i] = self.drhodcos[  :, i][ind_sorted]
            self.domegadcos[:,  i] = self.domegadcos[:, i][ind_sorted]
            self.d2rdz2[:,      i] = self.d2rdz2[    :, i][ind_sorted]
            self.d2rdcosdz[:,   i] = self.d2rdcosdz[ :, i][ind_sorted]
            self.d2rhodcos2[:,  i] = self.d2rhodcos2[:, i][ind_sorted]
            self.d2pdcos2[:,    i] = self.d2pdcos2[  :, i][ind_sorted]
            self.d2rdcos2[:,    i] = self.d2rdcos2[:,   i][ind_sorted]
            self.Ar[:,          i] = self.Ar[        :, i][ind_sorted]
            if iext:
                self.geff2d[0, :, i] = self.geff2d[0, :, i][ind_sorted]
                self.geff2d[1, :, i] = self.geff2d[1, :, i][ind_sorted]
                self.geff2d[2, :, i] = self.geff2d[2, :, i][ind_sorted]

        self.omega_c = - 0.5*np.sqrt( self.gamma12d[self.Nt-1, self.Nr-1] / self.rho2d[self.Nt-1, self.Nr-1] / self.pr2d[self.Nt-1, self.Nr-1] ) \
            * self.dpdr[self.Nt-1, self.Nr-1]*self.drdz[self.Nt-1, self.Nr-1]
        self.omega_c_Hz = self.omega_c * self.omega_k*1e6/(2*np.pi)

        if verbose: print(f"{GREEN}[Done]{NO_COLOR}")

    def calc_vconv(self, origin, l_hp=True, K0=1.7, Gphi=2, beta=0.05, i=None, var=None,
                   glob=None):
        """
        <CModel.calc_vconv( origin, l_hp=True, K0=1.7, Gphi=2, beta=0.05, i=None, var=None,
            glob=None )>

        Computes the convective velocity, flux and the entropy fluctuations
        according to a given formulation of convection.

        :param origin: Origin of the data. Either from osc file or rot file.
        :type origin: str

        :kparam l_hp: if True, mixing $\textrm{length} = \alpha H_p$.
        :ktype l_hp: float

        :kparam Gphi: Gough's anisotropy factor (default `Gphi=2`).
        :ktype Gphi: float

        :kparam K0: coefficient K0 (see H02), (default `K0=1.5`).
        :ktype K0: float

        :kparam i: Index of the time step for which we want to compute convective quantities.
        :ktype i: integer

        :kparam var: structure data of a model (single time step).
        :ktype var: 2d array

        :kparam glob: global data of a model (single time step).
        :ktype glob: 2d array

        :member Fconv: Convectif flux
        :mtype Fconv:  np.ndarray

        :member Frad:  Radiatif flux
        :mtype Frad:   np.ndarray

        :member lmix:  Mixing-length ($\textrm{length} = \alpha H_p$)
        :mtype lmix:   np.ndarray

        :member Pturb: Turbulent pressure
        :mtype Pturb:  np.ndarray

        :member vconv: Convective velocity
        :mtype vconv:  np.ndarray

        :example: To compute convective quantities of model number 100, from data
        stored in an osc file:
        >>> mdl = CModel( 'model' )
        >>> mdl.calc_vconv( 'osc', i=100 )
        >>> mdl.calc_vconv( 'osc', i=-1 ) # for only last model

        :reference H02:     Heiter et al, A&A, 2002
        :reference CGM:     Canuto Goldman Mazitelli ApJ 473, 550-559, 1996
        :reference CM:      Canuto Mazitelli ApJ 370, 295, 1991
        :reference PaperI:  Samadi & Goupil, 2001

        :history 28.04.03: Reza Samadi
        :history 30.04.04: Updated
        :history 11.02.11: Adapted JP

        :raises ValueError: You must either `i` or `var` and `glob`.
        :raises TypeError: Keyword i must be an integer.
        """

        if origin == 'osc':
            if self.all_osc and ((var is None or glob is None) and i is None):
                raise ValueError( "You must either `i` or `var` and `glob`. See documentation." )

            if var is None:
                if type( self.var ) == list:
                    if type( i ) != int: raise TypeError( "Keyword i must be an integer." )

                    var = self.var[i]
                else:
                    var = self.var


            if glob is None:
                if isinstance( self.glob, (np.ndarray, list) ):
                    if type( i ) != int: raise TypeError( "Keyword i must be an integer." )

                    glob = self.glob[i]
                else:
                    glob = self.glob

            pr     = var[3]
            ma     = var[1]*glob[0]
            rho    = var[4]
            ra     = var[0]
            te     = var[2]
            kappa  = var[7]
            delta  = var[11]
            grad   = var[5]
            gradad = var[10]
            cp     = var[12]
            gamma1 = var[9]
            lum    = var[6]
        elif origin == 'rot':
            pr      = self.pr
            ma      = self.m*self.ctes.msun
            rho     = self.rho
            r       = self.r*self.ctes.rsun
            T       = self.T
            kappa   = self.kappa
            delta   = self.delta
            gradT   = self.grad
            gradTad = self.gradad
            cp      = self.cp
            gamma1  = self.gamma1
            lum     = self.lu*self.ctes.lsun

        grav = self.ctes.ggrav*ma[:-1]/ra[:-1]**2
        hp   = pr[:-1]/rho[:-1]/grav     # echelle des hauteurs
        grav = np.append(grav, 0.0)
        hp   = np.append(hp, 0.0)


        alpha = glob[5]
        if l_hp:
            lmix = alpha*hp                   # mixing-length
        else:
            hp    = interpolate.interp1d(r[-1::-1], hp[-1::-1], kind='linear')
            lmix  = []
            nconv = len(self.lconv[-1])
            rstar = 10.0**self.log_r[-1]
            ricz  = []
            for i in range(nconv):
                if self.lconv[-1][i] == 'F':
                    if i == 0:
                        ricz.append([0.0, self.rcz[-1][i]*self.ctes.rsun])
                    else:
                        ricz.append([self.rcz[-1][i-1]*self.ctes.rsun,
                                     self.rcz[-1][i]*self.ctes.rsun])
                else:
                    if i == nconv-1:
                        ricz.append([self.rcz[-1][i]*self.ctes.rsun, rstar*self.ctes.rsun])
                    else:
                        ricz.append([self.rcz[-1][i]*self.ctes.rsun,
                                     self.rcz[-1][i+1]*self.ctes.rsun])
            lmix = []
            for i in range(self.ntot):
                lcz = False
                for rcz in ricz:
                    if rcz[0] < r[i] < rcz[1]:
                        rcz[1] = min(rcz[1], self.var[0][0])
                        hp1 = hp(rcz[0])
                        hp2 = hp(rcz[1])
                        z1 = r[i] - rcz[0] + beta*hp1
                        z2 = rcz[1] - r[i] + beta*hp2
                        lmix.append(z1*z2/(z1 + z2))
                        lcz = True
                        break
                if not lcz:
                    lmix.append(0.0)
            lmix = np.array(lmix)

        alpha_s = gradad*te*rho*gamma1 #  (dP/ds)_rho

        #  K -> conductibilite thermique (sigma = a * c /4)
        krad = 4.0/3.0*aradia*clight*te**3/kappa/rho

        # conductivite thermometrique:
        ki = krad/cp/rho
        # le 4a**2 de H02 (Eq. 6) ou Eq. 6 CM:
        aa = (2.0*(lmix**2)/9.0/ki*grav)**2*rho/2.0/pr*delta
        dgrad = grad - gradad


        if(self.params.nom_conv_ == 'conv_jmj'):
            # CESAM : V/A=2/9 * l
            V_Al   = 2.0/9.0
            # xi = 1/72(3*V/A/l) ^2
            xi0    = 1.0/72.0*(3.0*V_Al)**2
            # CESAM: phi0= 9/8
            phi0   = 1.0/4.0/V_Al
            #* 0.5/2.4
            tauedd = kappa*rho*lum
            with np.errstate(divide='ignore'):
                dvst = (1.0 + 2.0/3.0/V_Al/tauedd**2)
                xi   = xi0*dvst**2
                with np.errstate(invalid='ignore'):
                    BB         = xi*(lmix**4)*(rho*cp)**2*delta*grav*krad**(-2.0)/hp
                    Gamma      = (np.sqrt(1.0 + 4.0*BB*dgrad) - 1.0)/2.0
                    # vertical component of velocity
                    self.vconv = np.sqrt( Gamma / (Gamma + 1) / 8 * delta * alpha**2 * hp * grav * dgrad )
                self.Pturb = rho*self.vconv**2
                # Eq.56
                R2  = (6.0*gamma1/alpha)**2
                # Eq.(6) H02
                sig = aa*(grad - gradad)
                S   = 81.0/2.0*sig
                with np.errstate(invalid='ignore'):
                    # Eq. 4 CM
                    phi = phi0/dvst*(np.sqrt(1.0 + sig) - 1.0)**3/sig
                # Flux convectif Eq. 2 H02 ou Eq.3 CM
                    self.Fconv = krad*te/hp*dgrad*phi

        elif (self.params.nom_conv_ == 'conv_cgm_reza'):
            # formulation CGM ApJ 473:550 , 1996

            # Eq.(19) H02 :
            ca   = 10.8654
            cb   = 0.00489073
            ck   = 0.149888
            cm   = 0.189238
            cn   = 1.85011

            # Eq.(21) H02 :
            cc  = 0.0108071
            cd  = 0.00301208
            ce  = 0.000334441
            cf  = 0.000125
            cpp = 0.72
            cq  = 0.92
            cr  = 1.2
            ct  = 1.5

            # Eq.(6) H02
            S          = 81.0/2.0*aa*(grad - gradad)
            N          = len(S)
            mask       = S < 0.0
            F3, ss, F4, v2, self.vconv, F5, F1, F2, self.Pturb, phi, self.Fconv \
                       = [np.zeros( N ) for _ in range(11)]
            S[mask]    = 0.0

            with np.errstate(divide='ignore'):
                # Eq. 89 CGM
                F3         = 0.00101392*S**2/(np.sqrt(0.000017848*S**2 + 1.0) + 1.0)*(K0/1.5)**3
                # Eq. 90 CGM
                ss         = 0.000777055*S**0.868589
                F4         = 2.256815*(ss - 1.0) / (ss + 1.0) + 6.39899
                # mean square value of the convective velocity
                v2         = (ki/lmix)**2*F3*F4
                # according to Eq. 88 CGM
                # vertical component of velocity
                self.vconv = np.sqrt(v2/Gphi)
                # Eq. 92 CGM
                ss         = 0.00111378*S**0.868589
                F5         = 1.49168 + 0.45185*(-1.0 + ss) / (1.0 + ss)
                # Eq. 91 CGM
                self.Pturb = rho*(ki/lmix)**2*F3*F5
                # Eq. 18 H02
                F1         = (K0/1.5)**3*ca*(S**ck)*((1.0 + cb*S)**cm -1.0)**cn
                # Eq. 20 H02
                F2         = 1.0 + cc*S**cpp/(1.0 + cd*S**cq) + ce*S**cr/(1.0 + cf*S**ct)
                # Eq.17 H02
                phi        = F1*F2

                # Flux convectif Eq. 2 H02 ou Eq.3 CM
                self.Fconv = krad*te/hp*dgrad*phi

        with np.errstate(divide='ignore'):
            self.Frad = krad*te*grad/hp
        self.lmix = lmix


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def write_osc(self, filename):
        """
        <CModel.write_osc( filename )>

        Write an osc file.

        :param filename: Name of the osc file to be written.
        :type filename: string
        """

        cont = self.osc_header
        line = f"{self.nchim:3d} " + ' '.join(self.nom_elem)

        cont.append(line + '\n')

        line = "%10i%10i%10i%10i%10i" % (self.ntot, self.nglob, self.nvar,
                                         self.nchim, self.krot)

        cont.append(line + '\n')

        i = 0
        # global variables
        while True:
            line = ''
            if i+5 > self.nglob:
                for k in self.glob[i:]:    line += "%19.12E" % k
                cont.append(line + '\n')
                break
            else:
                for k in self.glob[i:i+5]: line += "%19.12E" % k
                cont.append(line + '\n')
            i += 5

        self.vars_out        = deepcopy(self.var)
        self.vars_out[1,:-1] = np.log(self.var[1,:-1])
        self.vars_out[1,-1]  = -1.0e38
        # global variables
        nw = self.nvar+self.nchim
        for j in range(self.ntot):
            i = 0
            while True:
                line = ''
                if i+5 > nw:
                    for k in self.vars_out[i:,j]: line += "%19.12E" % k
                    cont.append(line + '\n')
                    break
                else:
                    for k in self.vars_out[i:i+5,j]: line += "%19.12E" % k
                    cont.append(line + '\n')
                i += 5

        with open(filename, 'w', encoding='latin-1') as f:
            f.writelines(cont)



    def calc_freqs( self, oscs=None, amdls=None, mods=None, write_amdl=True,
        quiet=False, read=True, from_gui=False, out_name=None, verbose=True, nmax=-1 ):
        """
        <CModel.calc_freqs( oscs=None, mods=None, write_amdl=False, quiet=False, read=True )>

        Wrapper for the computation of frequencies that will either call ADIPLS or ACOR.

        :kparam oscs: .osc files for which we want to compute the frequencies.
            oscs and amdls cannot be not None at the same time.
            If oscs and amdls are None, will compute frequencies for all .osc files.
        :ktype oscs: list of strings

        :kparam amdls: .amdl files for which we want to compute the frequencies.
            oscs and amdls cannot be not None at the same time.
            If oscs and amdls are None, will compute frequencies for all .osc files.
        :ktype amdls: list of strings

        :kparam mods: Indices of the time steps for which frequencies must be computed.
        :ktype mods: list of integers

        :kparam write_amdl: If True, .amdl have already been computed and we do not need
            to compute them again.
        :ktype write_amdl: boolean

        :kparam quiet: If True, standard output is not printed to screen.
        :ktype quiet: boolean

        :kparam read: If True, frequencies are read after computation.
        :ktype read: boolean

        :kparam from_gui: If True, `calc_freqs` was called from GUI. Default: False.
        :ktype  from_gui: boolean

        :kparam out_name: Name of the output .agsm file (useful when you compute only one file).
            If not given, the filename is constructed using the osc file. Default: None
        :ktype  out_name: str

        :kparam verbose: If True, information is written for the user. Default: True.
        :ktype verbose: boolean
        """
        if not write_amdl and not self.fparams.adipls_new_amdl:
            if not self.finished:
                print('Model not yet calculated...')
                return 0

        if self.fparams.reduce_ != 'none':
            self.read_hr( )

        self.fparams.write_fsettings( )

        self.init_fparams( self.fparams.modes_, keep=True )
        if self.fparams.osc_code_name_=='adipls':
            self.calc_freqs_adipls( oscs=oscs, amdls=amdls, mods=mods, write_amdl=write_amdl,
                read=read, from_gui=from_gui, out_name=out_name, verbose=verbose, nmax=nmax )
        elif self.fparams.osc_code_name_=='acor':
            self.calc_freqs_acor( mods=mods, quiet=quiet )

    def calc_freqs_adipls( self, oscs=None, amdls=None, mods=None, write_amdl=True,
        read=True, from_gui=False, out_name=None, verbose=True, nmax=None ):
        """
        <CModel.calc_freqs_adipls( oscs=None, mods=None, write_amdl=False, read=True )>
        Computes frequencies using the ADIPLS oscillation code
        (J. Christensen-Dalsgaard et al., 2008).

        :kparam oscs: .osc files for which we want to compute the frequencies.
            oscs and amdls cannot be not None at the same time.
            If oscs and amdls are None, will compute frequencies for all .osc files.
        :ktype oscs: list of strings

        :kparam amdls: .amdl files for which we want to compute the frequencies.
            oscs and amdls cannot be not None at the same time.
            If oscs and amdls are None, will compute frequencies for all .osc files.
        :ktype amdls: list of strings

        :kparam mods: Indices of the time steps for which frequencies must be computed.
        :ktype mods: list of integers

        :kparam write_amdl: If True, .amdl have already been computed and we do not need
            to compute them again.
        :ktype write_amdl: boolean

        :kparam read: If True, frequencies are read after computation.
        :ktype read: boolean

        :kparam from_gui: If True, `calc_freqs` was called from GUI. Default: False.
        :ktype from_gui: boolean

        :kparam out_name: Name of the output .agsm file (useful when you compute only one file).
            If not given, the filename is constructed using the osc file. Default: None
        :ktype  out_name: str

        :kparam verbose: If True, information is written for the user. Default: True.
        :ktype verbose: boolean
        """

        gsm = []
        ssm = []

        if write_amdl and amdls is not None:
            raise CESAMError( "amdls has been given a not None value. write_amdl should be False." )

        amdl_given = False


        if self.fparams.to_zip:
            if os.path.exists( f'{self.name}.agsm.zip' ):       os.remove( f'{self.name}.agsm.zip' )
            if os.path.exists( f'{self.name}.ssm.zip' ):        os.remove( f'{self.name}.ssm.zip' )
            if os.path.exists( f'{self.name}.amdl.zip' ):       os.remove( f'{self.name}.amdl.zip' )
            if os.path.exists( f'{self.name}-adipls.log.zip' ): os.remove( f'{self.name}-adipls.log.zip' )
            if os.path.exists( f'{self.name}.adipls.zip' ):     os.remove( f'{self.name}.adipls.zip' )
            if os.path.exists( f'{self.name}.rkr.zip' ):        os.remove( f'{self.name}.rkr.zip' )

        if self.params.out_h5:
            if write_amdl:
                self.osc_amdl(osc=self.osch5)

            steps=[]
            if self.params.nom_output_.endswith('_adiah5'):
                suf = '-ad'
                ind = 9
            elif self.params.nom_output_.endswith('_nadiah5'):
                suf = '-nad'
                ind = 10
            elif self.params.nom_output_.endswith('_platoh5'):
                suf = '-plato'
                ind = 12
            else:
                suf = '-inv'
                ind = 10

            if mods is not None:
                nmax - len(mods) if nmax is None else nmax
                steps = [f'{i:05d}-{self.name}{suf}.osch5' for i in mods[:nmax]]
            else:
                with h5py.File(self.osch5, 'r') as f:
                    if not steps:
                        steps = list(f.keys())[:nmax]

            agsms  = []
            iagsms = []

            with zipfile.ZipFile( f'{self.name}.agsm.zip', 'a' )       if self.fparams.to_zip else nullcontext() as agsmz, \
                 zipfile.ZipFile( f'{self.name}.ssm.zip', 'a' )        if self.fparams.to_zip else nullcontext() as ssmz, \
                 zipfile.ZipFile( f'{self.name}.amdl.zip', 'a' )       if self.fparams.to_zip else nullcontext() as amdlz, \
                 zipfile.ZipFile( f'{self.name}-adipls.log.zip', 'a' ) if self.fparams.to_zip else nullcontext() as alogz, \
                 zipfile.ZipFile( f'{self.name}.adipls.zip', 'a' )     if self.fparams.to_zip else nullcontext() as adipz, \
                 zipfile.ZipFile( f'{self.name}.rkr.zip', 'a' )        if self.fparams.to_zip else nullcontext() as rkrz:
                for i, s in enumerate( steps ):
                    imod = mods[i] if mods is not None else i
                    iagsms.append( imod )
                    step = s[:-ind]
                    if str2bool( self.fparams.remesh_ ): self.redistrb(step + '.amdl')

                    if i>0 and self.fparams.reduce_ == 'gauss_norad':
                        if os.path.exists( f"{steps[i-1]}_rad.adipls" ):
                            os.remove( f"{steps[i-1]}_rad.adipls" )
                        if os.path.exists( f"{steps[i-1]}_nrad.adipls" ):
                            os.remove( f"{steps[i-1]}_nrad.adipls" )

                    agsms = self._CModel__reduce_adipls( step, agsms, i, verbose=verbose )
                    if self.fparams.to_zip:
                        if self.fparams.reduce_ == 'gauss_norad':
                            ff = [agsmz, agsmz, ssmz, ssmz, adipz, adipz, amdlz, alogz]
                            nn = ['_rad.agsm', '_nrad.agsm', '_rad.ssm', '_nrad.ssm', '_rad.adipls',
                                '_nrad.adipls', '.amdl', '-adipls.log']
                        else:
                            ff = [agsmz, ssmz, amdlz, alogz, adipz]
                            nn = ['.agsm', '.ssm', '.amdl', '-adipls.log', '.adipls']

                        if self.fparams.rotkr:
                            ff.append( rkrz )
                            nn.append( '.rkr' )

                        for (f,n) in zip(ff, nn):
                            f.write( f'{step}{n}' )
                            os.remove( f'{step}{n}' )

        else:
            if oscs is None:
                # amdl files should be written ?
                if self.fparams.adipls_new_amdl:
                    amdl_given = True
                    oscs = [self.fparams.adipls_amdl_name]
                elif write_amdl:
                    if oscs is None:
                        if self.all_osc:
                            oscs = np.array(self.osc)
                            if mods:
                                oscs = oscs[mods]
                            elif self.fparams.step > 1:
                                oscs = oscs[::self.fparams.step]
                        else:
                            oscs = [self.osc]
                # the list of amdl files is given ?
                elif amdls is not None:
                    amdl_given = True
                    oscs = amdls
                # of the above
                else:
                    oscs = [self.osc]
            # if oscs and amdls are given: abort.
            else:
                if amdls is not None:
                    raise CESAMError("oscs and amdls should not be given at the same time.")

            agsms  = []
            iagsms = []

            with zipfile.ZipFile( f'{self.name}.agsm.zip', 'a' )       if self.fparams.to_zip else nullcontext() as agsmz, \
                 zipfile.ZipFile( f'{self.name}.ssm.zip', 'a' )        if self.fparams.to_zip else nullcontext() as ssmz, \
                 zipfile.ZipFile( f'{self.name}.amdl.zip', 'a' )       if self.fparams.to_zip else nullcontext() as amdlz, \
                 zipfile.ZipFile( f'{self.name}-adipls.log.zip', 'a' ) if self.fparams.to_zip else nullcontext() as alogz, \
                 zipfile.ZipFile( f'{self.name}.adipls.zip', 'a' )     if self.fparams.to_zip else nullcontext() as adipz, \
                 zipfile.ZipFile( f'{self.name}.rkr.zip', 'a' )        if self.fparams.to_zip else nullcontext() as rkrz:
                for i, osc in enumerate( oscs ):
                    imod = int( osc[:5] ) if self.all_osc else -1
                    iagsms.append( imod )
                    if write_amdl and not self.fparams.adipls_new_amdl:
                        # Make amdl file
                        self.osc_amdl(osc=osc)
                        name = osc[:6] + self.name if self.all_osc else self.name
                    elif amdl_given:
                        name = osc[:-5]
                    else:
                        name = self.name

                    g_grav = self.ctes.ggrav if amdl_given else self.read_osc_glob( osc )[13]


                    if str2bool( self.fparams.remesh_ ) and not amdl_given: self.redistrb(name + '.amdl')


                    agsms = self._CModel__reduce_adipls( name, agsms, i, ggrav=g_grav, verbose=verbose )
                    if self.fparams.to_zip:
                        if self.fparams.reduce_ == 'gauss_norad':
                            ff = [agsmz, agsmz, ssmz, ssmz, adipz, adipz, amdlz, alogz]
                            nn = ['_rad.agsm', '_nrad.agsm', '_rad.ssm', '_nrad.ssm', '_rad.adipls',
                                '_nrad.adipls', '.amdl', '-adipls.log']
                        else:
                            ff = [agsmz, ssmz, amdlz, alogz, adipz]
                            nn = ['.agsm', '.ssm', '.amdl', '-adipls.log', '.adipls']

                        if self.fparams.rotkr:
                            ff.append( rkrz )
                            nn.append( '.rkr' )

                        for (f,n) in zip(ff, nn):
                            f.write( f'{name}{n}' )
                            os.remove( f'{name}{n}' )

        if read:
            if self.fparams.reduce_ != 'none':
                self.read_agsm( imods=iagsms, reduce=True )
            elif from_gui:
                self.read_agsm( )
            else:
                self.read_agsm( mods=agsms )

    def calc_freqs_acor(self, mods=None, quiet=False):
        """
        <CModel.calc_freqs_acor( mods=None, quiet=False )>

        Computes frequencies using the ACOR oscillation code (R.-M. Ouazzani et al., 2015).

        :kparam mods: Indices of the time steps for which frequencies must be computed.
        :ktype mods: list of integers

        :kparam read: If True, frequencies are read after computation.
        :ktype read: boolean
        """

        if self.fparams.acor_params_select_ == 1:
            if self.fparams.acor_m > self.fparams.lmax:
                raise AttributeError("For computing frequencies, m must be lower than self.fparams.lmax.")

        oscs   = []
        osc2ds = []

        if self.all_osc:
            if self.fparams.input_file_[:2]=='d1':
                oscs = np.array(self.osc)
                if mods:
                    oscs = oscs[mods]
                elif self.fparams.only_last:
                    oscs = oscs[-1]
                elif self.fparams.step1d > 1:
                    oscs = oscs[::self.fparams.step1d]
            if self.fparams.input_file_[-1]=='2':
                osc2ds = np.array(self.osc2d)
                if mods:
                    osc2ds = osc2ds[mods]
                elif self.fparams.only_last:
                    osc2ds = osc2ds[-1]
                elif self.fparams.step2d > 1:
                    osc2ds = osc2ds[::self.fparams.step2d]
        else:
            if self.fparams.input_file_[:2]=='d1':
                oscs = [self.osc]
            if self.fparams.input_file_[-1]=='2':
                osc2ds = [self.osc2d]


        for osc in oscs:
            if not self.fparams.mod_suff:
                self.fparams.mod_suff = "1d"
            if self.fparams.mod_flag_ != 'custom':
                if self.all_osc:
                    self.fparams.mod_flag_custom = osc[:6] + self.name
                else:
                    self.fparams.mod_flag_custom = self.name

            self.calc_freq_lmax_acor('cestam', osc, False, quiet=quiet)

        for osc2d in osc2ds:
            if not self.fparams.mod_suff:
                self.fparams.mod_suff = "2d"
            if self.fparams.mod_flag_ != 'custom':
                if self.all_osc:
                    self.fparams.mod_flag_custom = osc2d[:6] + self.name
                else:
                    self.fparams.mod_flag_custom = self.name

            self.calc_freq_lmax_acor('cesam2d', osc2d, True, quiet=quiet)


    def calc_freq_lmax_acor(self, evol_code, nmod, d2, quiet=True):
        """
        <CModel.calc_freq_lmax_acor( evol_code, nmod, d2, quiet=True )>

        Computes frequencies using the ACOR oscillation code (R.-M. Ouazzani et al., 2015),
        up to a certain degree lmax. ACOR can only find modes with a given parity between
        l and m. Therfore, in general, we need to run ACOR twice for even and odd (l,m).
        The difficulty is to find the number of needed spherical harmonics giving l and m.
        There may be a more elegant way than what I did here.

        :param evol_code: Name of the stellar evolution code that was used to compute the
            structure model.
        :type evol_code: string

        :param nmod: Name of the model that contains the structure.
        :type nmod: string

        :param d2: If True, the structure model is a 2D model.
        :type d2: boolean

        :kparam quiet: If True, standard output will not be printed when frequencies are computed.
        :ktype quiet: boolean
        """
        procs = []
        if self.fparams.acor_params_select_ == 1:
            lfreq_name = "liste_frequences_"+self.fparams.mod_flag_custom
            if self.fparams.input_file_[:2]=='d1':
                if self.fparams.acor_rot:
                    lfreq_name += "_"+self.fparams.mod_suff
                else:
                    lfreq_name += "_norot"
            # If self.fparams.lmax and m are specified
            self.afreq = CFacor(f'{lfreq_name}_M{self.fparams.acor_nsh:d}_m{self.fparams.acor_m:d}' \
                                        f'_par{not self.fparams.acor_parity:d}')
            if (1-self.fparams.acor_m % 2) and (1-self.fparams.lmax % 2):
                if self.fparams.lmax == 0:
                    # parameters for (l,m) of same parity
                    self.fparams.acor_parity = True
                    self.fparams.acor_nsh    = 1
                    procs = self._CModel__run_acor( evol_code=evol_code, nmod=nmod, d2=d2, quiet=quiet, procs=procs )
                else:
                    # parameters for (l,m) of same parity
                    self.fparams.acor_parity = True
                    self.fparams.acor_nsh    = self.fparams.lmax - self.fparams.acor_m//2
                    procs = self._CModel__run_acor( evol_code=evol_code, nmod=nmod, d2=d2, quiet=quiet, procs=procs )
                    # parameters for (l,m) of opposite parity
                    self.fparams.acor_parity = False
                    self.fparams.acor_nsh    = self.fparams.lmax - (self.fparams.acor_m//2 + 1)
                    procs = self._CModel__run_acor( evol_code=evol_code, nmod=nmod, d2=d2, quiet=quiet, procs=procs )

            elif (1-self.fparams.acor_m % 2) and self.fparams.lmax % 2:
                self.fparams.acor_parity = True
                self.fparams.acor_nsh    = self.fparams.lmax - (self.fparams.acor_m//2 + 1)
                procs = self._CModel__run_acor( evol_code=evol_code, nmod=nmod, d2=d2, quiet=quiet, procs=procs )

                self.fparams.acor_parity = False
                self.fparams.acor_nsh    = self.fparams.lmax - (self.fparams.acor_m//2 + 1)
                procs = self._CModel__run_acor( evol_code=evol_code, nmod=nmod, d2=d2, quiet=quiet, procs=procs )

            elif self.fparams.acor_m % 2 and (1-self.fparams.lmax % 2):
                self.fparams.acor_parity = True
                self.fparams.acor_nsh    = self.fparams.lmax - (self.fparams.acor_m//2 + 1)
                procs = self._CModel__run_acor( evol_code=evol_code, nmod=nmod, d2=d2, quiet=quiet, procs=procs )

                self.fparams.acor_parity = False
                self.fparams.acor_nsh    = self.fparams.lmax - self.fparams.acor_m//2
                procs = self._CModel__run_acor( evol_code=evol_code, nmod=nmod, d2=d2, quiet=quiet, procs=procs )

            else:
                self.fparams.acor_parity = True
                self.fparams.acor_nsh    = self.fparams.lmax - (self.fparams.acor_m//2 + 1)
                procs = self._CModel__run_acor( evol_code=evol_code, nmod=nmod, d2=d2, quiet=quiet, procs=procs )

                self.fparams.acor_parity = False
                self.fparams.acor_nsh    = self.fparams.lmax - ((self.fparams.acor_m+1)//2 + 1)
                procs = self._CModel__run_acor( evol_code=evol_code, nmod=nmod, d2=d2, quiet=quiet, procs=procs )

            for proc in procs:
                proc.join()

            self.afreq = procs[0].cfacor
            self.afreq.merge( procs[1].cfacor )


        elif self.fparams.acor_params_select_ == 2:
            # If number of spherical harmonics, m, and parity are specified.
            print(f"Calculating freqs for model {self.fparams.mod_flag_custom}...")
            run_acor = self._CModel__run_acor(evol_code=evol_code, nmod=nmod, quiet=quiet)
            self.afreq = run_acor.cfacor

        self.clean_freq_acor(self.afreq)
        self.write_acor_freq(self.fparams.mod_flag_custom+"_"+self.fparams.mod_suff+".acor", self.afreq)
        os.system("rm -f liste_frequences_* acor_freq_*")


    def osc_amdl(self, osc=None):
        """
        <CModel.osc_amdl( osc=None )>

        Wrapper for the Fortran executable that converts an osc file to an amdl file.

        :kparam osc: Name of the osc file to be converted.
        :ktype osc: string
        """
        cmd = f'{self.name}\n'
        if osc is not None:
            cmd += f'{osc}\n'
        else:
            cmd += f'{self.osc}\n'

        s = subp.Popen('osc-amdl.x', stdin=subp.PIPE, stdout=subp.PIPE)

        sortie, err = s.communicate(input=cmd.encode())
        sortie = sortie.decode().split('\n')
        for line in sortie:
            if line.startswith(' Wrote'):
                print( GREEN + line + NO_COLOR )

    def osc_losc(self, osc=None):
        """
        <CModel.osc_losc( osc=None )>

        Wrapper for the Fortran executable that converts an osc file to an losc file.

        :kparam osc: Name of the osc file to be converted.
        :ktype osc: string
        """
        s = os.popen('osc-losc.x','w')
        s.write(f'{self.name}\n')
        if osc is not None:
            s.write(f'{osc}\n')
        else:
            s.write(f'{self.osc}\n')

        s.close()


    ## TODO: Do the same for ACOR. The routine is already written into ACOR


    def redistrb(self, amdl, out_name=None):
        """
        <CModel.redistrb( amdl, out_name=None )>

        Runs the ADIPLS redistribution algorithm.

        :param amdl: Name of .amdl file.
        :type amdl: string

        :kparam out_name: If given, the amdl file is not overwritten. Default: None.
        :ktype out_name: string
        """

        # file redistrib.in for p or g modes
        if self.fparams.remesh_ in ('p', 'g'):
            fname = 'redistrib_' + self.name + '.in'
            with open(fname, 'w') as f:
                f.write(f'2 \'{amdl}\'    @\n')
                f.write(f'3 \'{amdl}.1\'    @\n')

                f.write('-1 \'\'        @\n')
                f.write('nn,icnmsh\n')
                f.write(f'{self.fparams.npoints},,,  @\n')
                f.write('icase,icvzbn\n')
                if self.fparams.remesh_ == 'p':
                    f.write('11,,,  @\n')
                else:
                    f.write('12,,,  @\n')

                f.write('cg,cx,ca,cdgr,cddgr,alphsf\n')
                f.write(',,,,,,,,,,,,,,  @\n')

                f.write('nout,cn,irsu,unew  \n')
                f.write('30,,,,,,,,,  @\n')

                f.write('nmodel,kmodel,itsaml,ioldex\n')
                f.write(',,,,,,,,,,,  @\n')

            res = subp.Popen(['redistrb.d', fname],stdout=subp.PIPE).communicate()[0]
            new_amdl = out_name if out_name is not None else amdl
            os.rename( amdl + '.1', new_amdl )

        elif self.fparams.remesh_ == 'r':
            fname = 'redistrib_' + self.name + '.c.in'
            with open(fname, 'w') as f:
                f.write(f'2 \'{amdl}\'    @\n')
                f.write(f'3 \'{amdl}.1\'    @\n')
                f.write('-1 \'\'        @\n')

                f.write('nn,icnmsh\n')
                f.write(f'{self.fparams.npoints},,,  @\n')

                f.write('icase,icvzbn,nsmth,ndisc,dlxdsc,dlgrmx,cacvzb\n')
                f.write('201,,,0.013,,5.,,,,,,,,,,,, @\n')

                f.write('cg,cx,ca,sig1,sig2,lmax,alphsf,adda4,accrm\n')
                f.write('1.,0.05,0.05,0.,0.,2,,0.02,0.01,,,,,,,,,,,  @\n')

                f.write('nout,cn,irsu,unew\n')
                f.write('60,,,,,,,,,  @\n')

                f.write('nmodel,kmodel,itsaml,ioldex\n')
                f.write(',,,,,,,,,,,  @\n')

            res = subp.Popen(['redistrb.c.d', fname], stdout=subp.PIPE).communicate()[0]
            new_amdl = out_name if out_name is not None else amdl
            os.rename( amdl + '.1', new_amdl )


    def remesh(self, rep=None, don=None):
        """
        <CModel.remesh( rep=None, don=None )>

        Remesh a model with points at a suitable location for accurate computation of
        p-modes (`'p'`), g-modes (`'g'`), modes in RGB stars (`'r'`), customized settings (`'o'`).

        :kparam rep: Name of the rep file of the model that needs to be remeshed.
        :ktype rep: string

        :kparam don: Name of the don file that was used to compute the model.
        :ktype don: string
        """

        if not rep:
            rep = self.file_dat
        if not don:
            don = self.name

        if self.fparams.remesh_ == 'n':
            cmd = f'{rep}\n{don}\nn\ny\n'
        elif self.fparams.remesh_ in ('p', 'g', 'r', 'o'):
            cmd = f'{rep}\n{don}\ny\no\n{self.fparams.cts[0]}, {self.fparams.cts[1]}, '\
                    f'{self.fparams.cts[2]}\n{self.fparams.npoints}\ny\n'
        else:
            cmd = f'{rep}\n{don}\ny\na\n{self.fparams.npoints}\ny\n'

        print("Re-meshing model and converting to amdl...", end='')
        s = subp.Popen('rep-osc_v3.x', stdin=subp.PIPE,
                             stdout=subp.PIPE, stderr=subp.PIPE)
        res, err = s.communicate(input=cmd.encode())
        print(f"{GREEN}[Done]{NO_COLOR}\n")
        messg = res.decode().split('\n')[-2]
        print(f"{GREEN_U}{messg}{NO_COLOR}")

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def __setexec(self, nmod=None):

        fname = 'setexec_' + self.name + '.in'
        if not nmod:
            nmod = self.name

        if os.path.exists(fname):
            os.remove(fname)

        agsm = f'{nmod}.agsm.1'
        miss = f'{nmod}.ssm.miss'

        with open(fname, 'w') as f:
            f.write(f'2 {agsm}   @\n')
            f.write(f'3 {miss}   @\n')
            f.write('4 \'ttt.setex.exec\'   @\n')
            f.write('-1 ''   @\n')
            f.write('1,,,,0,1000,,,,,,,,,,,,,  @\n')
            f.write('0,0  @\n')
            f.write('2,,,,,,,,,,,,,,,,,,  @\n')

        res = subp.Popen(['setexec.d', fname], stdout=subp.PIPE).communicate()[0]


    def __selsum(self, filein, fileout_gsm, fileout_ssm, n1, n2):
        """
        <CModel.__selsum( filein, fileout_gsm, fileout_ssm, n1, n2 )>
        Private wrapper that calls ADIPLS's selsum program.
        It select desired modes in a summary file.

        :param filein: Input file that contains modes characteristics.
        :type filein: str

        :param fileout_gsm: File name to output grand summary.
        :type fileout_gsm: str

        :param fileout_ssm: File name to output short summary.
        :type fileout_ssm: str

        :param n1: Minimum order of in the selected modes set.
        :type n1: int

        :param n1: Maximum order of in the selected modes set.
        :type n1: int
        """

        fname = f'selsum_{self.name}.in'
        if os.path.exists(fname): os.remove(fname)

        with open(fname, 'w') as f:
            f.write(f'2 {filein}   @\n')
            f.write(f'11 {fileout_gsm}   @\n')
            f.write(f'12 {fileout_ssm}   @\n')
            f.write('13 \'0\'   @\n')
            f.write('21 \'0\'   @\n')
            f.write('22 \'0\'   @\n')
            f.write('23 \'0\'   @\n')
            f.write('-1 \'\'   @\n')
            f.write(',,,   @\n')
            f.write('1,1   @\n')
            f.write(f'{n1},{n2},,,,,10010,,,,,,,,,,,,,,,,,,,,    @\n')
            f.write(f'0,,,{self.fparams.lmax},,,,,,,,,,,,,,,,,,    @\n' )
            f.write(f'0,{-self.fparams.lmax},,,,,,,,,,,,,,,,,,,,,,,,,,,,,    ,,,,,,,,,,,,,,,,,,,,  @\n' )
            f.write('1    @\n')

        res = subp.Popen(['selsum.d', fname], stdout=subp.PIPE).communicate()[0]


    def run_adipls( self, nn, nmod=None, ggrav=None, verbose=True ):
        """
        <CModel.run_adipls( nn, nmod=None )>

        Wrapper of `CModel.__run_adipls` for public use.

        :param nn:

        :kparam nmod: Name of the stellar track.
        :ktype nmod: string

        :kparam ggrav: Gravitational constant to be used. Default: uses the globally defined ggrav.
        :ktype ggrav: float

        :kparam verbose: If True, information is written for the user. Default: True.
        :ktype verbose: boolean
        """
        ggrav = self.ctes.ggrav if ggrav is None else ggrav
        self.__run_adipls( nn, nmod=nmod, ggrav=ggrav, verbose=verbose )

    def __select_modes( self, gsm, ssm ):
        # select modes with n=0-50000 (p modes) and -500000 - -1 (g modes)
        if self.fparams.modes_ != 'a':
            if self.fparams.modes_ == 'p':
                n1 = 0
                n2 = 5000000
            elif self.fparams.modes_ == 'g':
                n1 = -50000000
                n2 = -1
            os.rename( gsm, gsm+'.1' )
            self._CModel__selsum( gsm+'.1', gsm, ssm, n1, n2 )
            os.remove( gsm+'.1' )

    def __reduce_adipls( self, step, agsms, i, ggrav=None, verbose=True ):
        """
        :kparam verbose: If True, information is written for the user. Default: True.
        :ktype verbose: boolean
        """
        ggrav = self.ctes.ggrav if ggrav is None else ggrav
        # calculate first try
        if self.fparams.reduce_ == 'none':
            gsm  = step + '.agsm'
            ssm  = step + '.ssm'
            agsms.append( gsm )
            print(f"Calculating freqs for model {step}...")
            self._CModel__run_adipls( 1, nmod=step, ggrav=ggrav, verbose=verbose )
            self.__select_modes( gsm, ssm )

        elif self.fparams.reduce_ == 'gauss_norad':
            lmin = self.fparams.lmin
            nsel = self.fparams.nsel
            if lmin == 0:
                # radial modes
                self.fparams.nsel = 1
                gsm  = step + '.agsm'
                gsmr = step + '_rad.agsm'
                ssm  = step + '.ssm'
                ssmr = step + '_rad.ssm'
                infile = step + '_rad.adipls'
                agsms.append( gsmr )
                print(f"Calculating freqs for model {step}...")
                self._CModel__run_adipls( 1, nmod=step, ggrav=ggrav, infile=infile, verbose=verbose )
                os.rename( gsm, gsmr )
                os.rename( ssm, ssmr )


                self.__select_modes( gsmr, ssmr )

            if nsel > 1 or lmin > 0:
                # radial modes
                self.fparams.nsel = nsel-1 if lmin == 0 else nsel
                self.fparams.lmin = max( 1, lmin )
                sig1 = self.fparams.sig1
                sig2 = self.fparams.sig2
                # sigma from mosser 2012
                numax = 3104.0 * (self.mstar[i]) / (self.rstar[i])**2 * np.sqrt(5777.0/10.0**self.log_teff[i])
                delta_nuenv = 0.66 * numax * 0.88 / (2*np.sqrt(2 * np.log(2)))
                nu_min = max( 0.0, numax - self.fparams.nred * delta_nuenv )
                nu_max = numax + self.fparams.nred * delta_nuenv

                self.fparams.sig1 = max( 3.0, sigma2( nu_min*1e-3, self.mstar[i]*self.ctes.msun, self.rstar[i]*self.ctes.rsun, ggrav ) )
                self.fparams.sig2 = sigma2( nu_max*1e-6, self.mstar[i]*self.ctes.msun, self.rstar[i]*self.ctes.rsun, ggrav )
                gsm   = step + '.agsm'
                gsmnr = step + '_nrad.agsm'
                ssm   = step + '.ssm'
                ssmnr = step + '_nrad.ssm'
                infile = step + '_nrad.adipls'
                agsms.append( gsmr )
                print(f"Calculating freqs for model {step}...")
                self._CModel__run_adipls( 1, nmod=step, ggrav=ggrav, infile=infile, verbose=verbose )
                os.rename( gsm, gsmnr )
                os.rename( ssm, ssmnr )


                self.__select_modes( gsmnr, ssmnr )

            self.fparams.sig1 = sig1
            self.fparams.sig2 = sig2
            self.fparams.lmin = lmin
            self.fparams.nsel = nsel



        return agsms

    def __run_adipls( self, nn, infile=None, nmod=None, ggrav=None, verbose=True ):
        """
        <CModel.__run_adipls( nn, nmod=None )>

        Private method called internally to run ADIPLS oscillation computations.

        :param nn:

        :kparam nmod: Name of the stellar track.
        :ktype nmod: string

        :kparam ggrav: Gravitational constant to be used. Default: uses the globally defined ggrav.
        :ktype ggrav: float

        :kparam verbose: If True, information is written for the user. Default: True.
        :ktype verbose: boolean
        """
        ggrav = self.ctes.ggrav if ggrav is None else ggrav
        if not nmod:
            nmod = self.name

        nfmode = int(self.fparams.amde)
        irotkr = int(self.fparams.rotkr)
        igm1kr = int(self.fparams.gm1kr)

        infile = infile if infile is not None else f'{self.name}.adipls'

        if os.path.exists(infile): os.remove(infile)

        if nn == 1:
            suf = ''
        else:
            suf = '.2'

        amde = f'{nmod}.amde{suf}'
        agsm = f'{nmod}.agsm{suf}'
        ssm  = f'{nmod}.ssm{suf}'
        log  = f'{nmod}-adipls.log{suf}'
        miss = f'{nmod}.ssm.miss'
        amdl = f'{nmod}.amdl'
        rkr  = f'{nmod}.rkr{suf}'
        gkr  = f'{nmod}.gkr{suf}'

        with open(infile, 'w') as f:

            f.write(f'2  {amdl}   @\n')
            if self.fparams.amde: f.write(f'4  {amde}    @\n')
            f.write(f'9  {log}   @\n')

            if nn == 1:
                f.write('10 \'0\'   @\n')
            else:
                f.write(f'10 {miss}   @\n')

            f.write(f'11 {agsm}   @\n')
            if self.fparams.rotkr: f.write(f'12 {rkr}   @\n')
            if self.fparams.gm1kr: f.write(f'13 {gkr}   @\n')
            f.write(f'15 {ssm}   @\n')

            f.write('-1 \'\'   @\n')
            f.write('\ncntrd,\n')
            f.write('mod.osc.cst.int.out.dgn     @\n')

            f.write('\nmod:\n')
            f.write('ifind,xmod,imlds,in,irname,nprmod,\n')
            f.write(',,,,,,,   @\n')
            f.write('ntrnct,ntrnsf,imdmod,\n')
            f.write(',,,,,,,,,,,,,,,,,,,,,, @\n')

            f.write('\nosc:\n')
            f.write('el,nsel,els1,dels,dfsig1,dfsig2,nsig1,nsig2,\n')
            if nn == 1:
                f.write(f'0,{self.fparams.nsel},{self.fparams.lmin},{self.fparams.dels},,,,,,,,,,,,,,,,,,,,,,     @\n')
            else:
                f.write(f',0,0,{self.fparams.dels},,,,,,,,,,,,,,,,     @\n')

            f.write('itrsig,sig1,istsig,inomde,itrds,\n')
            if nn == 1:
                f.write(f'1,{self.fparams.sig1},,,,,,,,,,,,,,,   @\n')
                f.write('dfsig,nsig,iscan,sig2,\n')
                f.write(f'0,{self.fparams.nsig_},{self.fparams.iscan},{self.fparams.sig2},,,,,,,,,,,,,,,,,,,     @\n')
            else:
                f.write('4,,6000,1,10,,,,,,,,   @\n')
                f.write('dfsig,nsig,iscan,sig2,\n')
                f.write('0,,,,,,,,,,,,,,,,,,,,,,,,     @\n')

            f.write('eltrw1,eltrw2,sgtrw1,sgtrw2,\n')
            f.write(',,,,,,,,,,,,,,,,    @\n')

            f.write('\ncst:\n')

            f.write('cgrav,\n')
            f.write(f'{ggrav}               @\n')

            f.write('\nint:\n')

            f.write('iplneq,iturpr,icow,alb,\n')
            f.write(',,,,,             @\n')
            f.write('istsbc,fctsbc,ibotbc,fcttbc,\n')
            f.write(f'{self.fparams.istsbc_},{self.fparams.fsbc},,,,,,,,,,,,,,  @\n')
            f.write('mdintg,iriche,xfit,fcnorm,eps,epssol,itmax,\n')
            if self.fparams.remesh_ == 'r':
                f.write(f'{self.fparams.mdintg_},{int(self.fparams.iriche)},0.1,,1.0d-9,,15,,,,,,,,,,,  @\n')
            else:
                f.write(f'{self.fparams.mdintg_},{int(self.fparams.iriche)},0.9,,,,,,,,,,,,,,,  @\n')
            f.write('fsig,dsigmx,irsevn,\n')
            f.write(',,,,,,,,,,,,,,,,  @\n')

            f.write('\nout:\n')

            f.write('istdpr,nout,nprcen,irsord,iekinr,\n')
            f.write(f'9,,,{self.fparams.irsord},{self.fparams.iekinr_},,,,,,,,,,,     @\n')
            f.write('iper,ivarf,kvarf,npvarf,nfmode,\n')
            f.write(f'1,{self.fparams.ivarf},2,0,{nfmode} @\n')
            f.write('irotkr,nprtkr,igm1kr,npgmkr,ispcpr,,\n')
            f.write(f'{irotkr},1,{igm1kr},1,,,,,,,,,,     @\n')
            icaswn = 10 + 10000*self.fparams.istsbc_
            f.write('icaswn,sigwn1,sigwn2,frqwn2,frqwn2,iorwn1,iorwn2,frlwn1,frlwn2,ewnmax\n')
            f.write(f'{icaswn},,,,,,,,,,,,,,,,,,,,,,,,,,,,,   @\n')

            f.write('\ndgn:\n')

            f.write(',,,,,,,,,     @\n')
            f.write(',,,,,,,,,,,,    @\n')

        stdout = subp.PIPE if verbose else subp.DEVNULL
        res = subp.Popen( ['adipls.c.d', infile], stdin=subp.PIPE, stdout=stdout, stderr=subp.PIPE ).communicate()[0]



    def __run_acor( self, evol_code="2k202k20", nmod=None, d2=False, aparam_name='ACOR_parameters',
        quiet=False, procs=None ):
        """
        <CModel.__run_acor( evol_code="CESAM2k20", nmod=None, d2=False, aparam_name='ACOR_parameters',
            quiet=False, procs=None )>

        Private method called internally to run ACOR oscillation computations.

        :kparam evol_code: Name of the evolution code that computed the structure model.
        :ktype  evol_code: string

        :kparam nmod: Name of the model that contains the structure.
        :ktype  nmod: string

        :kparam d2: If True, the structure model is a 2D model.
        :ktype  d2: boolean

        :kparam aparam_name: Name of the ACOR parameter file. Default: `'ACOR_parameters'`.
        :ktype  aparam_name: string

        :kparam quiet: If True, standard output will not be printed when frequencies are computed.
        :ktype  quiet: boolean

        :kparam procs: List containing the ACOR instances that need to be ran in parallel.
        :ktype  procs: list of CRunAcor instances
        """

        pard = {True:'even', False:'odd'}
        par  = not self.fparams.acor_parity
        if self.fparams.acor_nsh > 0:
            print(f"Calculating freqs for model {self.fparams.mod_flag_custom} for {pard[par]} parities...")
            if procs is None:
                procs = CRunAcor( self.fparams, evol_code=evol_code, nmod=nmod, d2=d2, aparam_name=aparam_name )
                procs.start( )
                procs.join( )
            else:
                procs.append( CRunAcor( self.fparams, evol_code=evol_code, nmod=nmod, d2=d2, aparam_name=aparam_name+str(len(procs)), quiet=quiet ) )
                procs[-1].start()

        return procs


    def read_agsm(self, hdf5=False, filename=None, imods=None, mods=None, dnu=True, reduce=False ):
        """
        <CModel.read_agsm( hdf5=False, filename=None, mods=None, dnu=True )>

        Wrapper that reads agsm files, either in classical binary format, or in an HDF5 file.

        :kparam hdf5: If True, agsm file is an HDF5 file. Default False
        :ktype hdf5: boolean

        :kparam filename: Name of the file we want to read.
        :ktype filename: string

        :kparam mods: Indices of the agsm files to be read.
        :ktype mods: list of integers

        :kparam imods: Indices of the agsm files to be read.
        :ktype  imods: list of integers

        :kparam dnu: If True, frequencies are post-processed, i.e. we compute frequency ratios,
            large and small seperation, etc. Default: True.
        :ktype dnu: boolean

        :kparam reduce: If True, will read agsm in two seperate files for radial and non-radial modes. Default: False.
        :ktype reduce: boolean
        """
        if hdf5:
            self.read_agsmh5( filename=filename, mods=mods, dnu=dnu )
        elif self.fparams.to_zip:
            self.read_agsm_zip( filename=filename, mods=mods, imods=imods, dnu=dnu, reduce=reduce )
        else:
            if reduce:
                self.read_agsm_rad_nrad( imods=imods, dnu=dnu )
            else:
                self.read_agsm_original( filename=filename, mods=mods, dnu=dnu )



    def read_agsm_original( self, filename=None, mods=None, dnu=True ):
        """
        <CModel.read_agsm_original( filename=None, mods=None, dnu=True )>

        Read frequencies in an original binary agsm file.

        :kparam filename: Name of the file we want to read.
        :ktype  filename: string

        :kparam mods: Indices of the agsm files to be read.
        :ktype  mods: list of integers or list of strings

        :kparam dnu: If True, frequencies are post-processed, i.e. we compute frequency ratios,
            large and small seperation, etc. Default: True.
        :ktype  dnu: boolean
        """

        if filename is not None:
            if not os.path.exists( filename ):
                raise CESAMError('The filename you entered does not exist.\n \
                    Are you in the right directory?\n \
                    Did you forget the .agsm extension?')
            self.agsm = [filename]
        elif mods is not None:
            if isinstance( mods[0], (int, np.integer) ):
                self.agsm = ['%05d-%s.agsm' % (i, self.name) for i in mods]
            elif type( mods[0] ) == str:
                self.agsm = mods
        else:
            self.agsm = [i for i in os.listdir('.') if i[6:-5] == self.name and
                         i[:5].isdigit() and i[-5:] == '.agsm']

            self.agsm.sort()

            if len(self.agsm) == 0:
                self.agsm = self.name + '.agsm'
                if os.path.exists(self.name + '.agsm'):
                    self.agsm = [self.name + '.agsm']
                else:
                    print('Freqs. not yet calculated...')
                    return 0

        self.freq, self.per, self.beta, self.freq_var, self.freq_ri, self.sigma2, self.l, self.n, \
            self.freq_ord, self.per_ord, self.ekin_ord, self.ekin, self.ls, self.nuls, \
            self.dd01, self.nudd01, self.d02, self.rd02, self.nud02, self.d13, self.rd13, \
            self.nud13, self.dd01, self.rdd01, self.nudd01, self.dd10, self.rdd10, \
            self.nudd10 = [[] for _ in range(28)]

        nfiles = len(self.agsm)
        for i in tqdm( range( nfiles ), desc='Reading .agsm files' ):
            with FortranBinaryFile( self.agsm[i], mode='r' ) as f:
                self.freq.append(np.array([]))
                self.sigma2.append(np.array([]))
                self.per.append(np.array([]))
                self.beta.append(np.array([]))
                self.freq_var.append(np.array([]))
                self.freq_ri.append(np.array([]))
                self.ekin.append(np.array([]))
                self.l.append([])
                self.n.append([])

                g = self.ctes.ggrav

                while True:
                    if not self.__read_agsm_record( f ): break

        if nfiles: print( '' )

        if dnu:
            self.process_freqs( )

    def read_agsm_rad_nrad( self, imods, dnu=True ):
        """
        <CModel.read_agsm_rad_nrad( filename=None, imods=None, dnu=True )>

        Read frequencies in an original binary agsm file.

        :param imods: Indices of the agsm files to be read.
        :type  imods: list of integers

        :kparam dnu: If True, frequencies are post-processed, i.e. we compute frequency ratios,
            large and small seperation, etc. Default: True.
        :ktype  dnu: boolean
        """
        if isinstance( imods[0], (int, np.integer) ):
            if imods[0] == -1:
                self.agsm_rad  = [f'{self.name}_rad.agsm']
                self.agsm_nrad = [f'{self.name}_nrad.agsm']
            else:
                self.agsm_rad  = ['%05d-%s_rad.agsm' % (i, self.name) for i in imods]
                self.agsm_nrad = ['%05d-%s_nrad.agsm' % (i, self.name) for i in imods]
            self.agsm = self.agsm_rad + self.agsm_nrad
            self.agsm.sort()

        else:
            raise CESAMError('imods should be a list of integers.')

        self.freq, self.per, self.beta, self.freq_var, self.freq_ri, self.sigma2, self.l, self.n, \
            self.freq_ord, self.per_ord, self.ekin_ord, self.ekin, self.ls, self.nuls, \
            self.dd01, self.nudd01, self.d02, self.rd02, self.nud02, self.d13, self.rd13, \
            self.nud13, self.dd01, self.rdd01, self.nudd01, self.dd10, self.rdd10, \
            self.nudd10 = [[] for _ in range(28)]

        nfiles = len(self.agsm_rad)
        for i in tqdm( range( nfiles ), desc='Reading .agsm files' ):
            with FortranBinaryFile( self.agsm_rad[i], mode='r' ) as f_r, FortranBinaryFile( self.agsm_nrad[i], mode='r' ) as f_nr:
                self.freq.append(np.array([]))
                self.sigma2.append(np.array([]))
                self.per.append(np.array([]))
                self.beta.append(np.array([]))
                self.freq_var.append(np.array([]))
                self.freq_ri.append(np.array([]))
                self.ekin.append(np.array([]))
                self.l.append([])
                self.n.append([])

                g    = self.ctes.ggrav

                while True:
                    if not self.__read_agsm_record( f_r ): break

                while True:
                    if not self.__read_agsm_record( f_nr ): break

        if nfiles: print( '' )

        if dnu:
            self.process_freqs( )

    def read_agsm_zip( self, filename=None, mods=None, imods=None, reduce=False, dnu=True ):
        """
        <CModel.read_agsm_zip( filename=None, mods=None, imods=None, reduce=False, dnu=True )>

        Read frequencies in an original binary agsm file.

        :kparam filename: Name of the file we want to read.
        :ktype  filename: string

        :kparam mods: Names of the agsm files to be read.
        :ktype  mods: string or list of strings

        :param imods: Indices of the agsm files to be read.
        :type  imods: list of integers

        :kparam dnu: If True, frequencies are post-processed, i.e. we compute frequency ratios,
            large and small seperation, etc. Default: True.
        :ktype  dnu: boolean
        """
        if filename is not None:
            self.agsm = [filename]
        elif mods is not None and isinstance( mods, (str) ):
            self.agsm = [mods]
        elif mods is not None and isinstance( mods, (np.ndarray, list) ):
            self.agsm = mods
        elif imods is not None and reduce and isinstance( imods[0], (int, np.integer) ):
            n = self.name
            self.agsm_rad  = [f'{n}_rad.agsm']  if imods[0] == -1 else ['%05d-%s_rad.agsm'  % (i, n) for i in imods]
            self.agsm_nrad = [f'{n}_nrad.agsm'] if imods[0] == -1 else ['%05d-%s_nrad.agsm' % (i, n) for i in imods]
            self.agsm      = self.agsm_rad + self.agsm_nrad
            self.agsm.sort()
        else:
            with zipfile.ZipFile( f"{self.name}.agsm.zip", 'r') as zipf: self.agsm = zipf.namelist()

        nfiles = len( self.agsm )

        self.freq, self.per, self.beta, self.freq_var, self.freq_ri, self.sigma2, self.l, self.n, \
            self.freq_ord, self.per_ord, self.ekin_ord, self.ekin, self.ls, self.nuls, \
            self.dd01, self.nudd01, self.d02, self.rd02, self.nud02, self.d13, self.rd13, \
            self.nud13, self.dd01, self.rdd01, self.nudd01, self.dd10, self.rdd10, \
            self.nudd10 = [[] for _ in range(28)]

        modulo = 2 if reduce else 1

        with zipfile.ZipFile( f"{self.name}.agsm.zip", 'r') as zipf:
            for i in tqdm( range( nfiles ), desc='Reading .agsm files' ):
                if self.agsm[i] in zipf.namelist():
                    with FortranBinaryFile( self.agsm[i], zipf=zipf, mode='r' ) as file:
                        if not (i % modulo):
                            self.freq.append(np.array([]))
                            self.sigma2.append(np.array([]))
                            self.per.append(np.array([]))
                            self.beta.append(np.array([]))
                            self.freq_var.append(np.array([]))
                            self.freq_ri.append(np.array([]))
                            self.ekin.append(np.array([]))
                            self.l.append([])
                            self.n.append([])

                        while True:
                            if not self.__read_agsm_record( file ): break

        if nfiles: print( '' )

        if dnu:
            self.process_freqs( )

    def read_agsmh5( self, filename=None, mods=None, dnu=True ):
        """
        <CModel.read_agsmh5( filename=None, mods=None, dnu=True )>

        Read frequencies in an HDF5 file.

        :kparam filename: Name of the file we want to read.
        :ktype  filename: string

        :kparam mods: Indices of the agsm files to be read.
        :ktype  mods: list of integers

        :kparam dnu: If True, frequencies are post-processed, i.e. we compute frequency ratios,
            large and small seperation, etc. Default: True.
        :ktype  dnu: boolean
        """
        steps    = []

        if filename is not None:
            setps = ['/'+filename]
            file = filename[6:]

        if mods is not None:
            steps = ['/%05d-%s.agsmh5' % (i, self.name) for i in mods]
            file = self.name + ".agsmh5"
        else:
            file = self.name + ".agsmh5"
            if not os.path.exists( file ):
                print(f'Freqs. not yet calculated or {file} missing...')
                return 0

        self.freq, self.per, self.beta, self.freq_var, self.sigma2, self.freq_ri, self.l, self.n, \
            self.freq_ord, self.per_ord, self.ekin_ord, self.ekin, self.ls, self.nuls, \
            self.dd01, self.nudd01, self.d02, self.rd02, self.nud02, self.d13, self.rd13, \
            self.nud13, self.dd01, self.rdd01, self.nudd01, self.dd10, self.rdd10, \
            self.nudd10   = [[] for _ in range(28)]

        print(f'Reading {file}...', end='')
        with h5py.File(file, 'r') as f:
            if not steps:
                steps = list(f.keys())
            for step in steps:
                g = self.ctes.ggrav

                self.l.append(        f[step+'/l'][:].astype(int) )
                self.n.append(        f[step+'/n'][:].astype(int) )
                self.ekin.append(     f[step+'/E'][:] )
                self.freq_var.append( f[step+'/nuVar'][:]*1.0e3 )
                self.beta.append(     f[step+'/betanl'][:] )
                self.freq_ri.append(  f[step+'/nuRi'][:]*1.0e3 )
                self.sigma2.append(   f[step+'/sig2'][:] )
                mstar = f[step+'/Ma'][:]
                rstar = f[step+'/Ra'][:]
                omega2 = self.sigma2[-1]*g*mstar/(rstar**3)
                nu = np.sqrt(omega2)/(2.*np.pi)*1000.0 # in mHz

                self.freq.append(     nu*1000.0 ) # in muHz
                self.per.append(      1.0/nu*1000.0) # in muHz

        print(f"{GREEN}[Done]{NO_COLOR}")
        if dnu:
            self.process_freqs( )

    def __read_agsm_record( self, file ):
        try:
            cs  = file.readRecordNative('d')
            self.l[-1].append(int(cs[17]))
            self.n[-1].append(int(cs[18]))
            self.ekin[-1]     = np.append( self.ekin[-1], cs[23] )
            self.freq_var[-1] = np.append( self.freq_var[-1], cs[26]*1.0e3 )
            self.beta[-1]     = np.append( self.beta[-1], cs[35] )
            self.freq_ri[-1]  = np.append( self.freq_ri[-1], cs[36]*1.0e3 )
            self.sigma2[-1]   = np.append( self.sigma2[-1], cs[19]  )
            mstar             = cs[1]
            rstar             = cs[2]
            omega2            = cs[19]*self.ctes.ggrav*mstar/(rstar**3)
            nu                = np.sqrt(omega2)/(2.*np.pi)*1000.0 # in mHz

            self.freq[-1]     = np.append(self.freq[-1], nu*1000.0) # in muHz
            self.per[-1]      = np.append(self.per[-1], 1.0/nu*1000.0) # in muHz

            return True

        except IndexError:
            return False


    def process_freqs( self ):
        """
        <CModel.process_freqs( )>

        Compute data products from frequencies such as frequency ratios, large and small
        seperation, etc.
        """

        comb_cond = lambda k, n, xx : bool( np.prod([(k[i] in xx) for i in range(n)]))
        comb2    = lambda k, xx : xx[k[0]] - xx[k[1]]
        comb5    = lambda k, xx : xx[k[0]] - 4*xx[k[1]] + 6*xx[k[2]] - 4*xx[k[3]] + xx[k[4]]

        nfiles = len( self.l )
        for kk in range( nfiles ):
            dim = len( self.l[kk] )
            if dim > 0:
                self.fparams.lmax = np.max( self.l[kk] )
                nmax = np.max( np.abs( np.array( self.n[kk] ) ) )
            else:
                self.fparams.lmax = 0
                nmax = 0

            self.freq_ord.append( {} )
            self.ekin_ord.append( {} )
            self.per_ord.append( {} )
            self.ls.append( [] )
            self.nuls.append( [] )

            for l, n, nu, ekin, per in zip( self.l[kk], self.n[kk], self.freq[kk], self.ekin[kk], self.per[kk] ):
                self.freq_ord[kk][(l, n)] = nu
                self.ekin_ord[kk][(l, n)] = ekin
                self.per_ord[kk][(l, n)]  = per

            if nmax > 0:
                # large separation
                self.ls[kk]   = np.zeros( (self.fparams.lmax+1, nmax+1) )
                self.nuls[kk] = np.zeros( (self.fparams.lmax+1, nmax+1) )
                for i in range( self.fparams.lmax+1 ):
                    self.ls[kk][i]   = np.zeros(nmax+1)
                    self.nuls[kk][i] = np.zeros(nmax+1)
                    for j in range( 1, nmax+1 ):
                        if ((i,j) in self.freq_ord[kk].keys()) and ((i,j-1) in self.freq_ord[kk].keys()):
                            self.ls[kk][i][j]   = self.freq_ord[kk][(i,j)] - self.freq_ord[kk][(i,j-1)]
                            self.nuls[kk][i][j] = self.freq_ord[kk][(i,j)]
                        else:
                            self.ls[kk][i][j]   = np.nan
                            self.nuls[kk][i][j] = np.nan

                self.d02.append(    np.zeros( nmax+1 ) )
                self.rd02.append(   np.zeros( nmax+1 ) )
                self.nud02.append(  np.zeros( nmax+1 ) )
                self.d13.append(    np.zeros( nmax+1 ) )
                self.rd13.append(   np.zeros( nmax+1 ) )
                self.nud13.append(  np.zeros( nmax+1 ) )
                self.dd01.append(   np.zeros( nmax+1 ) )
                self.rdd01.append(  np.zeros( nmax+1 ) )
                self.nudd01.append( np.zeros( nmax+1 ) )
                self.dd10.append(   np.zeros( nmax+1 ) )
                self.rdd10.append(  np.zeros( nmax+1 ) )
                self.nudd10.append( np.zeros( nmax+1 ) )

                # small separation d02
                for j in range(1,nmax):
                    k = [(0,j), (2,j-1)]

                    if comb_cond( k, 2, self.freq_ord[kk].keys() ):
                        self.d02[kk][j]   = comb2( k, self.freq_ord[kk] )
                        self.nud02[kk][j] = self.freq_ord[kk][k[0]]
                        if np.isnan( self.ls[kk][1][j] ):
                            self.rd02[kk][j] = np.nan
                        else:
                            self.rd02[kk][j] = self.d02[kk][j]/self.ls[kk][1][j]
                    else:
                        self.d02[kk][j]   = np.nan
                        self.nud02[kk][j] = np.nan
                        self.rd02[kk][j]  = np.nan

                # small separation d13
                for j in range(1,nmax):
                    k1 = [(1,j), (3,j-1)]
                    if comb_cond( k, 2, self.freq_ord[kk].keys() ):
                        self.d13[kk][j]   = comb2( k, self.freq_ord[kk] )
                        self.nud13[kk][j] = self.freq_ord[kk][k[0]]

                        if (j == nmax) or (np.isnan( self.ls[kk][0][j+1] )):
                            self.rd13[kk][j] = np.nan
                        else:
                            self.rd13[kk][j] = self.d13[kk][j]/self.ls[kk][0][j+1]
                    else:
                        self.d13[kk][j]   = np.nan
                        self.nud13[kk][j] = np.nan
                        self.rd13[kk][j]  = np.nan

                if nmax > 0:
                    for j in range(1,nmax):
                        # 5 point separtion dd01
                        k = [(0,j-1), (1,j-1), (0,j), (1,j), (0,j+1)]
                        if comb_cond( k, 5, self.freq_ord[kk].keys() ):
                            self.dd01[kk][j]   = comb5( k, self.freq_ord[kk] )
                            self.dd01[kk][j]   = self.dd01[kk][j]/8.0
                            self.nudd01[kk][j] = self.freq_ord[kk][k[2]]
                            if np.isnan(self.ls[kk][1][j]):
                                self.rdd01[kk][j] = np.nan
                            else:
                                self.rdd01[kk][j] = self.dd01[kk][j]/self.ls[kk][1][j]

                        else:
                            self.dd01[kk][j]   = np.nan
                            self.nudd01[kk][j] = np.nan
                            self.rdd01[kk][j]  = np.nan

                        # 5 point separtion dd10
                        k = [(1,j-1), (0,j), (1,j), (0,j+1), (1,j+1)]
                        if comb_cond( k, 5, self.freq_ord[kk].keys() ):
                            self.dd10[kk][j]   = comb5( k, self.freq_ord[kk] )
                            self.dd10[kk][j]   =-self.dd10[kk][j]/8.0
                            self.nudd10[kk][j] = self.freq_ord[kk][k[2]]

                            if np.isnan(self.ls[kk][0][j+1]):
                                self.rdd10[kk][j] = np.nan
                            else:
                                self.rdd10[kk][j] = self.dd01[kk][j]/self.ls[kk][0][j+1]
                        else:
                            self.dd10[kk][j]   = np.nan
                            self.nudd10[kk][j] = np.nan
                            self.rdd10[kk][j]  = np.nan

                    self.nud02[kk][self.nud02[kk]==0.0]   = np.nan
                    self.nud13[kk][self.nud13[kk]==0.0]   = np.nan
                    self.nudd01[kk][self.nudd01[kk]==0.0] = np.nan
                    self.nudd10[kk][self.nudd10[kk]==0.0] = np.nan


        if len(self.agsm)== 1 :
            self.l        = self.l[0]
            self.n        = self.n[0]
            self.freq     = self.freq[0]
            self.sigma2   = self.sigma2[0]
            self.freq_var = self.freq_var[0]
            self.freq_ri  = self.freq_ri[0]
            self.per      = self.per[0]
            self.beta     = self.beta[0]
            self.ekin     = self.ekin[0]
            self.freq_ord = self.freq_ord[0]
            self.ekin_ord = self.ekin_ord[0]
            if nmax > 0:
                self.ls     = self.ls[0]
                self.nuls   = self.nuls[0]
                self.d02    = self.d02[0]
                self.rd02   = self.rd02[0]
                self.nud02  = self.nud02[0]
                self.d13    = self.d13[0]
                self.rd13   = self.rd13[0]
                self.nud13  = self.nud13[0]
                self.dd01   = self.dd01[0]
                self.rdd01  = self.rdd01[0]
                self.nudd01 = self.nudd01[0]
                self.dd10   = self.dd10[0]
                self.rdd10  = self.rdd10[0]
                self.nudd10 = self.nudd10[0]

        self.l        = np.array(self.l,        dtype=object)
        self.n        = np.array(self.n,        dtype=object)
        self.freq     = np.array(self.freq,     dtype=object)
        self.freq_var = np.array(self.freq_var, dtype=object)
        self.freq_ri  = np.array(self.freq_ri,  dtype=object)
        self.per      = np.array(self.per,      dtype=object)


    def read_acor( self, filename, suf='1d' ):
        """
        <CModel.read_acor( filename, suf='1d' )>
        Reads a reconstructed ACOR frequency file.

        :param filename: Name of the model.
        :type  filename: string

        :kparam suf: Suffix appended to the filename. Default: '1d'.
        :ktype  suf: string
        """
        lfreq_name = filename + '_' + suf + '.acor'

        print(f'Reading {lfreq_name}...', end='')
        data = np.loadtxt( lfreq_name, skiprows=1)
        if data.shape[1] == 21:
            self.number, self.m, self.l, self.nt, self.npt, self.ngt, self.n, self.np, \
                self.ng = np.transpose( data[:,0:9].astype(int) )

            self.Sigma_in, self.Sigma_out, self.Sigma_out2, self.freq, self.ekin, self.I, \
                self.Knl, self.beta, self.I_core, self.Knl_core, self.beta_core, \
                self.Omk = np.transpose( data[:,9:21] )
        else:
            self.number, self.m, self.l, self.n, self.np, self.ng = np.transpose( data[:,0:6].astype(int) )
            self.nt, self.npt, self.ngt = [ -42 * np.ones( self.number.size, dtype=int ) for _ in range(3)]

            self.Sigma_in, self.Sigma_out, self.Sigma_out2, self.freq, self.ekin, self.I, \
                self.Knl, self.beta, self.I_core, self.Knl_core, self.beta_core, \
                self.Omk = np.transpose( data[:,6:18] )
        print(f"{GREEN}[Done]{NO_COLOR}")



    def clean_freq_acor(self, cfacor):
        """
        <CModel.clean_freq_acor( cfacor )>
        Removes the modes that have been found several times.

        :param cfacor: Instance of class CFacor that represents the frequencies.
        :type  cfacor: Instance of class CFacor
        """
        nmax = np.max(cfacor.n)
        for i in np.arange(len(cfacor.n)):
            cfacor.index = np.append( cfacor.index, nmax*(cfacor.l[i]+1)+(cfacor.n[i] + cfacor.nt[i]) )
        indices = np.argsort( cfacor.index )
        cfacor.number     = cfacor.number[indices]
        cfacor.m          = cfacor.m[indices]
        cfacor.l          = cfacor.l[indices]
        if hasattr( cfacor, 'nt'):
            cfacor.nt         = cfacor.nt[indices]
            cfacor.npt        = cfacor.npt[indices]
            cfacor.ngt        = cfacor.ngt[indices]
        cfacor.n          = cfacor.n[indices]
        cfacor.np         = cfacor.np[indices]
        cfacor.ng         = cfacor.ng[indices]
        cfacor.Sigma_in   = cfacor.Sigma_in[indices]
        cfacor.Sigma_out  = cfacor.Sigma_out[indices]
        cfacor.Sigma_out2 = cfacor.Sigma_out2[indices]
        cfacor.freq       = cfacor.freq[indices]
        cfacor.ekin       = cfacor.ekin[indices]
        cfacor.I          = cfacor.I[indices]
        cfacor.Knl        = cfacor.Knl[indices]
        cfacor.Beta       = cfacor.Beta[indices]
        cfacor.I_core     = cfacor.I_core[indices]
        cfacor.Knl_core   = cfacor.Knl_core[indices]
        cfacor.Beta_core  = cfacor.Beta_core[indices]
        cfacor.Omk        = cfacor.Omk[indices]

        dup = self.__duplicate_mode(np.array([cfacor.n, cfacor.nt, cfacor.l]))
        dups = np.array([], dtype='int')
        for i in dup:
            dups = np.append(dups, i[1:])

        cfacor.number     = np.delete(cfacor.number, dups)
        cfacor.m          = np.delete(cfacor.m, dups)
        cfacor.l          = np.delete(cfacor.l, dups)
        if hasattr( cfacor, 'nt'):
            cfacor.nt         = np.delete(cfacor.nt, dups)
            cfacor.npt        = np.delete(cfacor.npt, dups)
            cfacor.ngt        = np.delete(cfacor.ngt, dups)
        cfacor.n          = np.delete(cfacor.n, dups)
        cfacor.np         = np.delete(cfacor.np, dups)
        cfacor.ng         = np.delete(cfacor.ng, dups)
        cfacor.Sigma_in   = np.delete(cfacor.Sigma_in, dups)
        cfacor.Sigma_out  = np.delete(cfacor.Sigma_out, dups)
        cfacor.Sigma_out2 = np.delete(cfacor.Sigma_out2, dups)
        cfacor.freq       = np.delete(cfacor.freq, dups)
        cfacor.ekin       = np.delete(cfacor.ekin, dups)
        cfacor.I          = np.delete(cfacor.I, dups)
        cfacor.Knl        = np.delete(cfacor.Knl, dups)
        cfacor.Beta       = np.delete(cfacor.Beta, dups)
        cfacor.I_core     = np.delete(cfacor.I_core, dups)
        cfacor.Knl_core   = np.delete(cfacor.Knl_core, dups)
        cfacor.Beta_core  = np.delete(cfacor.Beta_core, dups)
        cfacor.Omk        = np.delete(cfacor.Omk, dups)
        cfacor.index      = np.arange(len(cfacor.n))


    def __duplicate_mode(self, data, minoccur=2):
        """
        <CModel.__duplicate_mode( data, minoccur=2 )>

        Find duplicate modes in acor frequency files.

        :param data: Array of integers n and l.
        :type data: 2d array

        :kparam minoccur: Threshold at which an tuple (n,l) is considered to be a duplicate.
        :ktype minoccur: integer
        """
        ind    = np.lexsort(data)
        diff   = np.any(data.T[ind[1:]] != data.T[ind[:-1]], axis=1)
        edges  = np.where(diff)[0] + 1
        result = np.split(ind, edges)
        return [group for group in result if len(group) >= minoccur]

    def write_acor_freq(self, filename, cfacor):
        """
        <CModel.write_acor_freq( filename, cfacor )>

        Write acor frequencies in a file .acor.

        :param filename: Name of thde frequency file.
        :type filename: string

        :param cfacor: CFacor instance representing acor frequencies.
        :type cfacor: CFacor instance
        """
        fformat = "{:6d}{:4d}{:4d}{:4d}{:4d}{:4d}{:4d}{:4d}{:4d}{:22.12e}{:22.12e}{:22.12e}" \
            "{:22.12e}{:22.12e}{:22.12e}{:22.12e}{:22.12e}{:22.12e}{:22.12e}" \
            "{:22.12e}{:22.12e}\n"
        with open(filename, 'w') as f:
            f.write("numero   m   l  nt npt ngt   n  np  ng    Sigma_out/D_nu    "
                "    Sigma_out             Sigma_out^2           Sigma_out (muHz)"
                "      Ec(mode)              I(mode)               Knl(mode)     "
                "        Beta(mode)            I_core(mode)          Knl_core(mod"
                "e)        Beta_core(mode)       Omk * 1d-6 / 2 Pi\n")
            for i in cfacor.index:
                f.write(fformat.format(i, cfacor.m[i], cfacor.l[i], cfacor.nt[i],
                    cfacor.npt[i], cfacor.ngt[i], cfacor.n[i], cfacor.np[i],
                    cfacor.ng[i], cfacor.Sigma_in[i], cfacor.Sigma_out[i],
                    cfacor.Sigma_out2[i], cfacor.freq[i], cfacor.ekin[i],
                    cfacor.I[i], cfacor.Knl[i], cfacor.Beta[i], cfacor.I_core[i],
                    cfacor.Knl_core[i], cfacor.Beta_core[i], cfacor.Omk[i]))

    def pfreqs(self, p_or_g, nmod=-1):

        try:
            self.freq_ord
        except (AttributeError, NameError):
            self.read_agsm()

        if isinstance(self.freq_ord, dict):
            self.fparams.lmax = max(self.l)
            nmax = max(self.n)
            for l in range(self.fparams.lmax):
                print('-------------------')
                print(' l    n   nu (muHz)   ')
                print('-------------------')
                for n in range(nmax):
                    if (l, n) in self.freq_ord.keys():
                        print(f'{l:2d}  {n:3d} {self.freq_ord[(l, n)]:10.2f}')
        else:
            self.fparams.lmax = max(self.l[nmod])
            nmax = max(self.n[nmod])
            for l in range(self.fparams.lmax):
                print('-------------------')
                print(' l    n   nu (muHz)   ')
                print('-------------------')
                for n in range(nmax):
                    if (l, n) in self.freq_ord[nmod].keys():
                        print(f'{l:2d}  {n:3d} {self.freq_ord[nmod][(l, n)]:10.2f}')

        print('-------------------')

    def read_amdl(self, filename=None):
        """
        <CModel.read_amdl( filename=None )>

        Reads an amdl file. Stores the result in `self.file_data[8]` and `self.aa[6,npoints] .

        :kparam filename: Name of th file to be read. If no filename is specified, reads `self.amdl`.
        :ktype filename: str
        """
        if filename is None:
            filename = self.amdl

        if not os.path.exists( filename ):
            print(f'File {filename} does not exist.')
            return

        cmd = f'{filename}\n'

        s = subp.Popen('convert-amdl.x', stdin=subp.PIPE,
                             stdout=subp.PIPE)
        res, err = s.communicate(input=cmd.encode())

        res = res.decode().split('\n')

        if res:
            self.n_amdl = int(res[0])
            self.aa     = np.zeros((6, self.n_amdl))
            self.file_data   = list(map(float, res[1].split()))
            for i in range(self.n_amdl):
                self.aa[:, i] = list(map(float, res[i+2].split()))
        else:
            print(f'Error reading file {filename}')


    def read_amde(self, l, n, freq=None, filename=None, mod=None):
        r"""
        <CModel.read_amde( l, n, freq=None, filename=None, mod=None )>

        Reads an amde file. Reads eigenfunctions for given l and n from file filename

        :param l: Degree of the mode.
        :type l:  int

        :param n: Radial order of the mode.
        :type n:  int

        :kparam freq: If given, look for mode with closest frequency to freq. Default: None
        :ktype freq: float

        :kparam filename: Filename in which data should be read.
            If None, look for `mod-name.amde` if mod is given and otherwise, for `self.amde`.
            Default None.
        :ktype filename: str

        :kparam mod: Number of the time step. Default: None.
        :ktype mod: int


        :member xeig: Dimension: npoints
        :mtype xeig: np.ndarray

        :member yeig: Dimension: 6, npoints.
        :mtype yeig: np.ndarray

        :member yeig[0]: $\xi_r/R$ for radial modes or $\xi_r/R$ for non-radial modes.
        :member yeig[1]: $p'/(w^2 R^2 \rho)$ for radial modes or $l(l+1) \xi_h/R$ for non-radial modes.
        :member yeig[2]: $(4 \pi r^3 \rho/M)^(1/2)*self.yeig[0]$ for radial modes or $-x \phi'/(g r)$ for non-radial modes.
        :member yeig[3]: $x^2 d/dx(self.yeig[2]/x)$ for non-radial modes.
        :member yeig[4]: $(4 \pi r^3 \rho/M)^(1/2)*self.yeig[0]$ for non-radial modes.
        :member yeig[5]: $(4 \pi r^3 \rho/M)^(1/2)/\sqrt{l(l+1)}*self.yeig[1]$ for non-radial modes.
        """

        if mod is not None and filename is not None:
            raise ValueError( "`mod` and `filename` cannot both be given a non None value." )
        elif mod is not None:
            if type( mod ) != int:
                raise TypeError( "`mod` should be an integer." )
            filename = '%05d-%s.amde' % (mod, self.name)
        elif filename is not None:
            if type(filename) not in [str, np.str_]:
                raise TypeError( "`filename` should be an str." )
        else:
            filename = self.amde

        if not os.path.exists( filename ):
            raise CESAMError( f"The file {self.amde} does not exist here. "
                "Please specify a correct filename through `filename` "
                "or a time step index through `mod`." )

        if freq is not None:
            cmd = '%s\n%d\n%s\n%f\n' % (filename, l, 'f', freq)
        else:
            cmd = '%s\n%d\n%s\n%d\n' % (filename, l, 'n', n)

        s = subp.Popen('convert-amde.x', stdin=subp.PIPE, stdout=subp.PIPE)
        res, err = s.communicate(input=cmd.encode())
        res = res.decode().split('\n')

        if res:
            if res[0]:
                nn = int(res[0])
                print(f'#layers = {nn}')
                self.xeig = np.zeros(nn)
                self.yeig = np.zeros((6, nn))
                for i in range(nn):
                    dat = list(map(float, res[i+1].split()))
                    self.xeig[i]   = dat[0]
                    self.yeig[:,i] = dat[1:]
                return 1

        try:
            del self.xeig
            del self.yeig
        except (NameError, AttributeError):
            pass
        print('Mode not calculated...')
        return 0


    def read_rotkr(self, l, n, freq=None, filename=None, mod=None):
        """
        <CModel.read_rotkr( l, n, freq=None, filename=None, mod=None )>
        Reads an rotational kernel file.
        """

        if mod is not None:
            filename = f'{mod:05d}-{self.name}.rkr'
        elif filename is not None:
            filename = f'{self.name}.rkr'

        if freq is not None:
            cmd = f"{filename}\n{l:d}\nf\n{freq}\n"
        else:
            cmd = f"{filename}\n{l:d}\nn\n{freq}\n"


        s = subp.Popen('read_rotk.x', stdin=subp.PIPE,
                             stdout=subp.PIPE)
        res, err = s.communicate(input=cmd.encode())
        res = res.decode().split('\n')

        if res:
            nn = int(res[0])
            self.xrot = np.zeros(nn)
            self.krot = np.zeros(nn)
            for i in range(nn):
                dat = list(map(float, res[i+1].split()))
                self.xrot[i] = dat[0]
                self.krot[i] = dat[1]
        else:
            print('Mode not calculated...')


    p1       = lambda self, a, z, i       :  a[2*i+1] * z + a[2*i+0]
    p2       = lambda self, a, z, i       : ((a[3*i+2] * z) + a[3*i+1]) * z + a[3*i+0]
    ludwig_s = lambda self, a, x, y: a[0] + a[5] * x + a[6] * y + a[1] * np.exp( a[3] * x + a[4] * y )
    tanner_s = lambda self, fs, x, y, z: fs[2](z) + fs[4](z) * np.exp( (fs[0](z)*np.log10(x) + fs[1](z)*y-fs[3](z)) / fs[5](z) )
    magic_s  = lambda self, a, x, y, z: self.p2(a,z,0) + x*self.p2(a,z,1) + y*self.p2(a,z,2) + self.p2(a,z,3) * np.exp( x*self.p2(a,z,4) + y*self.p2(a,z,5))
    magic2   = lambda self, a, x, y, z: self.p2(a,z,0) + x*self.p2(a,z,1) + y*self.p2(a,z,2) + x**2*self.p2(a,z,3) + y**2*self.p2(a,z,4) + x*y*self.p2(a,z,5) + self.p2(a,z,6) * np.exp( x*self.p2(a,z,7) + y*self.p2(a,z,8))
    X2       = lambda self, p, x, y, z :  self.p2(p,z,0)*np.log10(x) + self.p2(p,z,1)*y
    cus2_fun = lambda self, p, x, y, z : (self.p2(p,z,2) * np.exp(self.X2(p,x,y,z) - self.p2(p,z,3)) ) \
        * (self.p2(p,z,4) + self.p2(p,z,5) * np.cos( self.p2(p,z,6)*self.X2(p,x,y,z) - self.p2(p,z,7) ))

    cus2_fun = lambda self, p, x, y, z : (self.p2(p,z,2) + self.p2(p,z,3) * np.exp(self.X2(p,x,y,z) - self.p2(p,z,4)) ) \
        * (self.p2(p,z,5) + self.p2(p,z,6) * np.cos( self.p2(p,z,7)*self.X2(p,x,y,z) - self.p2(p,z,8) ))

    def init_s( self, grid='original' ):
        """
        <CModel.init_s( grid='original' )>

        Stores prescription parameters.

        :kparam grid: Grid on which the coefficients are calibrated, Either `'original'` (default),
            `'cifist'`, `'cif_mdw'` or `'cif_red'` (recommended).
        :ktype grid: string

        :raises ValueError: If law does not matched one of the available option.
        """
        if grid == 'original':
            self.lud_par = np.array([1.6488, 0.0740, 0.0, 1.7860, -1.6762, 0.1274, -0.1412])
            self.mu_lud = 1.0/(2.0*0.704 + 0.75*0.28 + 0.5*0.016)

            feh = [-4.0, -2.0, -1.0, 0.0, 0.5]
            #                            FeH     -4.0     -2.0     -1.0      0.0      0.5
            A      = interpolate.interp1d( feh, [ 0.9985,  0.9981,  0.9974,  0.9967,  0.9961], kind=3 )
            B      = interpolate.interp1d( feh, [-0.0553, -0.0623, -0.0720, -0.0811, -0.0884], kind=3 )
            s0     = interpolate.interp1d( feh, [ 1.104,   1.254,   1.304,   1.336,   1.396] , kind=3 )
            x0     = interpolate.interp1d( feh, [ 3.606,   3.603,   3.540,   3.485,   3.435] , kind=3 )
            beta   = interpolate.interp1d( feh, [ 1.216,   1.439,   1.127,   1.051,   0.929] , kind=3 )
            tau    = interpolate.interp1d( feh, [ 0.0985,  0.0899,  0.0973,  0.1056,  0.1009], kind=3 )

            self.mag_par    = np.array([ 1.5789,  0.0455,  0.0111,   # a
                                    0.0784, -0.0183,  0.0071,   # b
                                   -0.1076, -0.0028, -0.0042,   # c
                                    0.1602,  0.0618,  0.0062,   # d
                                    1.2867, -0.0824,  0.0970,   # e
                                   -1.2136, -0.0338, -0.0764])  # f

            self.g_mu = interpolate.splrep( [-4.0, -3.0, -2.0, -1.0, -0.5, 0.0, 0.5],
                [0.59374534, 0.59375636, 0.59385896, 0.59488762, 0.59554686, 0.5994291, 0.61156366])

        elif grid == 'cifist':
            self.lud_par = np.array([1.67199769, 0.10100425, 0.0, 1.53800252, -1.41784705, 0.10501243, -0.14985532])
            self.mu_lud = 1.0 / (2.0*0.7373416 + 0.75*0.2492016 + 0.5*1.3456777e-2)

            feh = [-3.0, -2.0, -1.0, 0.0]

            A      = interpolate.interp1d( feh, [ 0.98961169,  0.99493933,  1.00499735,  1.00081353], kind=3 )
            B      = interpolate.interp1d( feh, [-0.06093353, -0.06698903, -0.07727121, -0.08185904], kind=3 )
            s0     = interpolate.interp1d( feh, [ 1.34685727,  1.26991333,  1.39872341,  1.45743295], kind=3 )
            x0     = interpolate.interp1d( feh, [ 3.56012518,  3.55049407,  3.51834665,  3.51622312], kind=3 )
            beta   = interpolate.interp1d( feh, [ 1.1962127,   1.02777147,  0.73165717,  1.23378744], kind=3 )
            tau    = interpolate.interp1d( feh, [ 0.06962119,  0.10314898,  0.09256627,  0.08633943], kind=3 )

            self.mag_par = np.array([1.67901310,     5.16687533e-2,  8.04353235e-03,
                                     1.06696845e-1, -1.40526587e-2,  4.94382783e-03,
                                    -1.55069343e-1, -1.85089003e-2, -2.83874574e-03,
                                     9.49214263e-2,  3.95751680e-2,  6.50800196e-03,
                                     1.59280342,    -1.95275001e-1, -1.53999199e-03,
                                    -1.46081040,     9.44902574e-2,  2.42135097e-03])

            self.g_mu = interpolate.splrep(
                [-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5],
                [0.5937374992328304, 0.5937400021879967, 0.5937478979739941,
                 0.5937727554097134, 0.5938511865043716, 0.5940994036676523,
                 0.5948836051069895, 0.5962405885928419, 0.5994079808613778,
                 0.6115162974442507])

            self.cus_par = np.array([-2.18378202e+00, -7.64622218e-01, -7.74188328e-02,
                                      -1.55087596e+00, -2.81033993e-01, -5.85868061e-02,
                                       1.24094347e+00,  2.05095969e-01,  2.33473027e-02,
                                      -2.66512069e-01, -4.11068835e-02, -2.40004599e-02,
                                      -1.57434322e-01,  1.27097690e-02,  3.67150778e-03,
                                       3.97172348e-01,  8.23168318e-03,  9.09349035e-03,
                                       3.95326488e+00,  8.58301306e-01,  9.09009468e-02,
                                       4.56841394e-01, -1.18689768e-02,  1.86234974e-02,
                                      -3.90940603e-01,  1.43312769e-02, -4.21773510e-03])
        elif grid == 'cif_mdw':
            self.lud_par = np.array([1.78893741, 0.01393728, 0.0, 2.83739722, -2.73529477, 0.20939864, -0.22146688])
            self.mu_lud = 1.0 / (2.0*0.7373416 + 0.75*0.2492016 + 0.5*1.3456777e-2)


            feh = [-3.0, -2.0, -1.5, -1.0, 0.0]

            A      = interpolate.interp1d( feh, [ 0.99226251, 0.95128643, 1.01024002, 1.00978999, 0.98941775], kind=3 )
            B      = interpolate.interp1d( feh, [-0.06109710,-0.06426326,-0.05042124,-0.08347966,-0.08521810], kind=3 )
            s0     = interpolate.interp1d( feh, [ 1.34685687,-0.32288953, 0.53781023,-0.07734037, 0.96002090], kind=3 )
            x0     = interpolate.interp1d( feh, [ 3.57150949, 3.57574060, 3.50458747, 3.49088598, 3.55465243], kind=3 )
            beta   = interpolate.interp1d( feh, [ 1.22832343, 3.36887271, 1.45804976, 2.04932174, 2.33885905], kind=3 )
            tau    = interpolate.interp1d( feh, [ 0.06980791, 0.54629542, 0.15673423, 0.59327275, 0.21258317], kind=3 )

            self.mag_par = np.array([ 1.80148322000, 0.0577016827, 8.38371649e-06,
                                      0.21215910200,-0.0388993900,-0.00623365372,
                                     -0.21819514000,-0.0242559349,-0.00487801174,
                                      0.00953529951, 0.0109658947, 0.00349894985,
                                      3.15817802000,-1.4452995500,-0.484573558,
                                     -3.06805112000, 1.4728516500, 0.561482463])

            self.g_mu = interpolate.splrep(
                [-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5],
                [0.5937374992328304, 0.5937400021879967, 0.5937478979739941,
                 0.5937727554097134, 0.5938511865043716, 0.5940994036676523,
                 0.5948836051069895, 0.5962405885928419, 0.5994079808613778,
                 0.6115162974442507])
            self.cus_par = np.array([ 0.80688874, -0.45354039, -0.07114942,
                                      -0.60076283, -0.35785253, -0.04721316,
                                       0.35466526,  0.24423163,  0.0305956,
                                      -0.19129024, -0.06727497, -0.01128207,
                                      -0.04142519, -0.01938396, -0.00278367,
                                       0.19342332,  0.08328439,  0.0107597,
                                       0.96125589,  0.53837118,  0.08279529,
                                       0.84485214, -0.01106269,  0.11195951,
                                      -0.65246183,  0.02974554, -0.05169905])

            # Mine old:
            #self.cus2_par = np.array([ 1.37837160e+01, -3.93297922e-02,  3.21206362e-02,
            #                          -1.12193744e+00, -8.16563381e-02,  2.26348485e-03,
            #                           6.26008357e-01,  2.18466872e-01,  2.67144435e-02,
            #                           4.94719879e+01, -1.21569234e+00, -9.92111313e-02,
            #                           5.21897704e+00,  1.19007110e+00, -1.19143651e-02,
            #                           6.19464798e-01, -1.72172888e-01, -1.10513661e-02,
            #                           3.72358413e+01, -8.41108302e+00, -3.20390452e-01,
            #                           1.51400829e+00, -9.35476383e-01,  4.39244101e-02,
            #                          -7.17828480e-01, -8.76002336e-02, -7.17316496e-02])
            # Justin's:
            #self.cus2_par = np.array( [ 6.33472012e+00, -1.17538820e-01,  4.40484200e-02,
            #                           -5.07914480e-01, -1.54680000e-02,  3.79104000e-03,
            #                           -3.02789046e+01, -1.73112714e+01, -3.42456659e+00,
            #                            2.61083536e+01,  5.19526673e-01,  1.02822600e-02,
            #                           -2.31738524e+01, -1.43207230e+01, -2.25982117e+00,
            #                            1.89186773e+01,  1.18191306e+01,  1.88566104e+00,
            #                            6.10432910e-01,  8.61525000e-03,  4.45608000e-03,
            #                           -2.41684163e+01, -1.98890730e-01,  1.96936030e-01] )

            self.cus2_par = np.array([ 1.12687311e+01 ,-1.21246578e-01  ,3.07647049e-02
                            ,-9.11646146e-01,-5.02661739e-02  ,4.66674888e-03
                            ,6.02383952e-01 ,-2.37084521e-01,1.28870808e-01
                            ,1.66358977e-06 ,-1.15444771e-06 ,-3.91989841e-07
                            ,2.57867014e+01 ,-9.70787983e-01 ,-8.62200232e-02
                            ,1.67935956e+00,7.22110516e-01  ,9.97674348e-02
                            ,5.02679263e-01  ,2.73597548e-01,4.25156273e-02
                            ,9.08857894e-01 ,-1.05265248e-01  ,1.51149962e-02
                            ,1.82447092e+00 ,-4.59881972e+00  ,9.04642780e-01])

        elif grid == 'cif_mdw2':
            feh = [-3.0, -2.0, -1.5, -1.0, 0.0]

            A      = interpolate.interp1d( feh, [ 0.99226251, 0.95128643, 1.01024002, 1.00978999, 0.98941775], kind=3 )
            B      = interpolate.interp1d( feh, [-0.06109710,-0.06426326,-0.05042124,-0.08347966,-0.08521810], kind=3 )
            s0     = interpolate.interp1d( feh, [ 1.34685687,-0.32288953, 0.53781023,-0.07734037, 0.96002090], kind=3 )
            x0     = interpolate.interp1d( feh, [ 3.57150949, 3.57574060, 3.50458747, 3.49088598, 3.55465243], kind=3 )
            beta   = interpolate.interp1d( feh, [ 1.22832343, 3.36887271, 1.45804976, 2.04932174, 2.33885905], kind=3 )
            tau    = interpolate.interp1d( feh, [ 0.06980791, 0.54629542, 0.15673423, 0.59327275, 0.21258317], kind=3 )

            self.g_mu = interpolate.splrep(
                [-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5],
                [0.5937374992328304, 0.5937400021879967, 0.5937478979739941,
                 0.5937727554097134, 0.5938511865043716, 0.5940994036676523,
                 0.5948836051069895, 0.5962405885928419, 0.5994079808613778,
                 0.6115162974442507])

            self.cus2_par = np.array([
                             1.31148433e+01, -9.31674092e-01,  1.93532251e-01,
                            -1.06215735e+00,  1.17293354e-02,  3.51311087e-03,
                             5.99546973e-01, -2.46621933e-01,  9.40483666e-02,
                             1.49511848e-07, -9.95452052e-08, -8.24121705e-10,
                             2.99049287e+01, -4.15038623e+00,  5.50774294e-01,
                             1.69345428e+00,  6.18088827e-01,  6.93286608e-02,
                             6.22503487e-01,  3.31156596e-01,  5.71817227e-02,
                             7.07317418e-01, -1.76873537e-02,  1.06887715e-03,
                            -1.18750703e+00, -3.46996675e+00,  6.61750211e-01  ])


        elif grid == 'cif_red':
            self.lud_par = np.array([1.78893741, 0.01393728, 0.0, 2.83739722, -2.73529477, 0.20939864, -0.22146688])
            self.mu_lud = 1.0 / (2.0*0.7373416 + 0.75*0.2492016 + 0.5*1.3456777e-2)


            feh = [-3.0, -2.0, -1.5, -1.0, 0.0]

            A      = interpolate.interp1d( feh, [ 0.99226251, 0.95128643, 1.01024002, 1.00978999, 0.98941775], kind=3 )
            B      = interpolate.interp1d( feh, [-0.06109710,-0.06426326,-0.05042124,-0.08347966,-0.08521810], kind=3 )
            s0     = interpolate.interp1d( feh, [ 1.34685687,-0.32288953, 0.53781023,-0.07734037, 0.96002090], kind=3 )
            x0     = interpolate.interp1d( feh, [ 3.57150949, 3.57574060, 3.50458747, 3.49088598, 3.55465243], kind=3 )
            beta   = interpolate.interp1d( feh, [ 1.22832343, 3.36887271, 1.45804976, 2.04932174, 2.33885905], kind=3 )
            tau    = interpolate.interp1d( feh, [ 0.06980791, 0.54629542, 0.15673423, 0.59327275, 0.21258317], kind=3 )

            self.mag_par = np.array([ 1.67263298,  0.0352205,   0.00386878,
                                      0.10501256, -0.06697099, -0.05256333,
                                     -0.14985512,  0.00498434,  0.01468436,
                                      0.10200094,  0.02628645, -0.01584492,
                                      1.53800265,  1.06588406,  1.3428408,
                                     -1.41784815, -0.95096173, -1.09586513])

            self.g_mu = interpolate.splrep(
                [-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5],
                [0.5937374992328304, 0.5937400021879967, 0.5937478979739941,
                 0.5937727554097134, 0.5938511865043716, 0.5940994036676523,
                 0.5948836051069895, 0.5962405885928419, 0.5994079808613778,
                 0.6115162974442507])
            self.cus_par = np.array([ 0.80688874, -0.45354039, -0.07114942,
                                      -0.60076283, -0.35785253, -0.04721316,
                                       0.35466526,  0.24423163,  0.0305956,
                                      -0.19129024, -0.06727497, -0.01128207,
                                      -0.04142519, -0.01938396, -0.00278367,
                                       0.19342332,  0.08328439,  0.0107597,
                                       0.96125589,  0.53837118,  0.08279529,
                                       0.84485214, -0.01106269,  0.11195951,
                                      -0.65246183,  0.02974554, -0.05169905])

        elif grid == 'extended':
            self.mag_par = np.array([ 1.56354660e+00, -3.00831464e-02, -7.45627059e-04,
                                      4.20187069e-02, -5.19911057e-02,  4.89058780e-03,
                                     -1.08923156e-01,  2.37246571e-02,  1.71333180e-03,
                                      2.00465658e-01,  1.26026681e-01,  2.07122883e-02,
                                      1.17133266e+00, -2.27989737e-01,  1.79083620e-01,
                                     -1.05694510e+00,  1.44257388e-01, -9.57347694e-02] )
            self.g_mu = interpolate.splrep(
                [-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5],
                [0.5937374992328304, 0.5937400021879967, 0.5937478979739941,
                 0.5937727554097134, 0.5938511865043716, 0.5940994036676523,
                 0.5948836051069895, 0.5962405885928419, 0.5994079808613778,
                 0.6115162974442507])

            # change this !
            A      = np.array( [1.0, 1.0, 1.0, 1.0, 1.0] )
            B      = np.array( [1.0, 1.0, 1.0, 1.0, 1.0] )
            s0     = np.array( [1.0, 1.0, 1.0, 1.0, 1.0] )
            x0     = np.array( [1.0, 1.0, 1.0, 1.0, 1.0] )
            beta   = np.array( [1.0, 1.0, 1.0, 1.0, 1.0] )
            tau    = np.array( [1.0, 1.0, 1.0, 1.0, 1.0] )

        else:
            raise ValueError( "the grid you ask does not exist. Available grids: \
                'original', 'cifist', 'cif_mdw', 'cif_red'.")

        self.tan_par = [A, B, s0, x0, beta, tau]

    def fmu( self, mu1d, murhd, srhd=None ):
        if srhd is None:
            return murhd / mu1d
        else:
            return murhd / mu1d - (mu1d - murhd) / (srhd)


    def get_ds( self, law, grid='cifist', i=-1, fmu_simple=True ):
        """
        <CModel.get_ds( law, grid='cifist', i=-1 )>

        Compute the entropy offset between a given prescription on a given grid, and
        the current model.

        :param law: Name of the prescription. Either `'ludwig'`, `'tanner'` or `'magic'` (recommended).
        :type law: string

        :kparam grid: Grid on which the coefficients are calibrated, Either `'original'`,
            `'cifist'` (default), `'cif_mdw'` or `'cif_red'`.
        :ktype grid: string

        :kparam i: index of the model taken as reference (default: last, i=-1).
        :ktype i: integer

        :return: Entropy offset.
        :rtype: float

        :raises CESAMError: If you asked to use a model that was not computed or not used.

        :raises ValueError: If law does not matched one of the available option.
        """
        self.init_s( grid=grid )

        if self.params.nom_output_.startswith('all_'):
            try:
                mu  = self.mu[i]
                s1d = self.senv[i]
            except AttributeError:
                # mu extracted from my optimal standard model:
                mu  = self.glob[i][18]
                s1d = self.glob[i][17]
        else:
            try:
                mu  = self.mu[i]
                s1d = self.senv[i]
            except AttributeError:
                if i != -1:
                    raise CESAMError( "Only last model available" )
                else:
                    # mu extracted from my optimal standard model:
                    mu  = self.glob[18]
                    s1d = self.glob[17]

        # Teff and log g extracted from my optimal standard model:
        teff = 10**self.log_teff[i]
        logg = self.log_g[i]
        if law == 'ludwig':
            # Evaluation of the functional
            x      = (teff - 5770.0) * 0.001
            y      = logg - 4.439333
            s_rhd  = self.ludwig_s( self.lud_par, x, y )

            # computation of fmu:
            if fmu_simple:
                fmu = self.fmu( mu, self.mu_lud )
            else:
                fmu = self.fmu( mu, self.mu_lud, srhd=s_rhd )

        elif law == 'tanner':
            z      = np.log10( self.ab_s['Fe56'][i] / np.abs( self.ab_s['H1'][i] ) ) - np.log10(56.0) + 4.5

            s_rhd  = self.tanner_s( self.tan_par, teff, logg, z )
            if fmu_simple:
                fmu = self.fmu( mu, interpolate.splev(z, self.g_mu) )
            else:
                fmu = self.fmu( mu, interpolate.splev(z, self.g_mu), s_rhd )

        elif law == 'magic':
            x      = (teff - 5777)/1000
            y      = logg - 4.44
            z      = np.log10( self.ab_s['Fe56'][i] / np.abs( self.ab_s['H1'][i] ) ) - np.log10(56.0) + 4.5

            s_rhd  = self.magic_s( self.mag_par, x, y, z )
            if fmu_simple:
                fmu = self.fmu( mu, interpolate.splev(z, self.g_mu) )
            else:
                fmu = self.fmu( mu, interpolate.splev(z, self.g_mu), s_rhd )
        elif law == 'custom':
            x      = (teff - 5777)/1000
            y      = logg - 4.44
            z      = np.log10( self.ab_s['Fe56'][i] / np.abs( self.ab_s['H1'][i] ) ) - np.log10(56.0) + 4.5

            s_rhd  = self.magic2( self.cus_par, x, y, z )
            if fmu_simple:
                fmu = self.fmu( mu, interpolate.splev(z, self.g_mu) )
            else:
                fmu = self.fmu( mu, interpolate.splev(z, self.g_mu), s_rhd )
        elif law == 'custom2':
            z      = np.log10( self.ab_s['Fe56'][i] / np.abs( self.ab_s['H1'][i] ) ) - np.log10(56.0) + 4.5

            s_rhd  = self.cus2_fun( self.cus2_par, teff, logg, z )
            if fmu_simple:
                fmu = self.fmu( mu, interpolate.splev(z, self.g_mu) )
            else:
                fmu = self.fmu( mu, interpolate.splev(z, self.g_mu), s_rhd )
        else:
            raise ValueError("The law '" + law + "' is not defined. Choose either 'ludwig', 'tanner', or 'magic'")

        return s1d - s_rhd*fmu

    def get_N2( self ):
        """
        <CModel.get_N2( )>

        Computes the value of the Brunt-Vaisala frequencies for all time step associated to the model.
        """
        in2 = 14 if not self.params.nom_output_.endswith('_platoh5') else 5
        if self.all_osc:
            N2 = []
            for i in range(self.nmod):
                ra = self.var[i][0]
                ma = self.var[i][1]
                grav = self.ctes.ggrav*ma*self.glob[i][0]/ra**3
                N2.append(self.var[i][in2] * grav)
        else:
            ra = self.var[0]
            ma = self.var[1]
            grav = self.ctes.ggrav*ma*self.glob[0]/ra**3
            N2 = self.var[in2] * grav

        return N2

    def get_I( self, imod=None, imin=0, imax=-1):
        """
        <CModel.get_I( )>

        Computes moment of inertia for all time step associated to the model.
        :kparams imin: Minimum layer index from which moment of inertia must be computed
        :ktype imin: int

        :kparams imax: Maximum layer index from which moment of inertia must be computed
        :ktype imax: int
        """
        if self.all_osc:
            I = []
            rrange = range(self.nmod) if imod is None else [imod]
            for i in rrange:
                ra = self.var[i][0]
                ma = self.var[i][1]
                ra2 = ra**2
                M  = self.mstar[i]*self.ctes.msun
                I.append( M * np.sum(0.5*(ra2[imin:imax-1] + ra2[imin+1:imax]) * (ma[imin:imax-1] - ma[imin+1:imax])) )
        else:
            ra = self.var[0]
            ma = self.var[1]
            ra2 = ra**2
            M  = self.mstar[-1]*self.ctes.msun
            I  = M * np.sum(0.5*(ra2[imin:imax-1] + ra2[imin+1:imax]) * (ma[imin:imax-1] - ma[imin+1:imax]))

        return I

    def get_J( self, imod=None, imin=0, imax=-1):
        """
        <CModel.get_J( )>

        Computes angular momentum for all time step associated to the model.
        :kparams imin: Minimum layer index from which angular momentum must be computed
        :ktype imin: int

        :kparams imax: Maximum layer index from which angular momentum must be computed
        :ktype imax: int
        """
        if self.all_osc:
            J = []
            rrange = range(self.nmod // self.params.osc_step) if imod is None else [imod]
            for i in rrange:
                ra = self.var[i][0]
                ma = self.var[i][1]
                om = self.var[i][15]
                omra2 = om*ra**2
                M  = self.mstar[i]*self.ctes.msun
                J.append( M * np.sum(0.5*(omra2[imin:imax-1] + omra2[imin+1:imax]) * (ma[imin:imax-1] - ma[imin+1:imax])) )
        else:
            ra = self.var[0]
            ma = self.var[1]
            om = self.var[15]
            omra2 = om*ra**2
            M  = self.mstar[-1]*self.ctes.msun
            J  = M * np.sum(0.5*(omra2[imin:imax-1] + omra2[imin+1:imax]) * (ma[imin:imax-1] - ma[imin+1:imax]))

        return J
