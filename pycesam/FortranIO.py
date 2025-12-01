#!/usr/bin/env python3
# coding: utf-8

# Code d'Evolution Stellaire, Adaptatif et Modulaire for the 2020 decade.
# Copyright (c) 1997-2023 The Cesam2k20 authors
# SPDX-License-Identifier : GPL-3.0-or-later


#============================================================
#
#           This file is part of FEval, a module for the
#           evaluation of Finite Element results
#
# Licencse: FEval is provided under the GNU General Public License (GPL)
#
# Authors:  Martin Lthi, tinu@tnoo.net
#
# Homepage: http://feval.sourceforge.net
#
# History:  2001.09.21 (ml):  Code cleaned up for intial release
#           2002.08.15 (ml):  Speed improvement of 20% when byteswapping is
#                             hardcoded (but less elegant)
#
#============================================================

# Fortran binary IO (at present only input of array)
# The handling of zipped files is based on Konrad Hinsen's TextFile

import os, string, sys
import struct
import numpy as np

MODE = 'rb'

# this is taken from the Numeric.py module, where it is defined globally
# LittleEndian = Numeric.fromstring("\001"+"\000"*7, 'i')[0] == 1

def myFromstring(endian):
    """
    Handle the problems caused by endian types.
    o Big endian:    Sun Sparc, HP etc.
    o Little endian: Intel machines (PC)
    Give here the endian type of the file, the endian type of the machine
    is automatically determined by the constant LittleEndian above.
    The conventions are as those given in 'struct' module:
    '@','=' : File is native (was produced on the same machine)
    '<'     : File is little endian
    '>','!' : File is big endian (network)

    :param endian: Endianness
    :type endian: str
    """
    littleendian = False
    try:
        if np.LittleEndian:  littleendian = True
    except:
        if np.little_endian: littleendian = True

    if littleendian:
        if endian in ['>','!']:
            return 'swap', lambda s, t: np.frombuffer(s,t).byteswap()
        else:
            return 'noswap', np.frombuffer
    elif endian in ['<']:
        return 'swap', lambda s, t: np.frombuffer(s,t).byteswap()
    else:
        return 'noswap', np.frombuffer


class FortranBinaryFile:

    def __init__(self, filename, mode='r', zipf=None, endian='@', verbose=False):
        """
        Fortran binary file support
        It is assumed that the file consists of Fortran records only.
        Constructor: FortranBinaryFile(filename, mode)

        :param filename: Name of the file to be read.
        :type filename: str

        :param mode: Openning mode ('r': read; 'w': write).
        :type mode: str
        """
        self.filename              = os.path.expanduser(filename)
        self.endian                = endian
        self.mode                  = mode
        self.swap, self.fromstring = myFromstring(endian)
        self.verbose               = verbose
        self.zipf                  = zipf

        if self.verbose:
            print('the byte-swapping is', self.swap)
        # this accelerates the file reading by about 20%
        if self.swap=='swap':
            setattr( FortranBinaryFile, 'readRecord', FortranBinaryFile.__dict__['readRecordByteswapped'] )
        else:
            setattr( FortranBinaryFile, 'readRecord', FortranBinaryFile.__dict__['readRecordNative'] )
        self.readRecord = self.readRecordByteswapped

    def __del__(self):
        self.file.close()

    def __enter__(self):
        if self.zipf is not None:
            self.file = self.zipf.open( self.filename, mode=self.mode )
        else:
            self.file = open( self.filename, self.mode+'b' )
        return self

    def __exit__(self, exception_type, exception_value, exception_traceback):
        self.file.close()

    def open( self, mode='r' ):
        """
        :param mode: Openning mode ('r': read; 'w': write).
        :type mode: str
        """
        if mode == 'r':
            if not os.path.exists(self.filename):
                raise IOError(2, 'No such file or directory: '+ self.filename)
            if self.filename[-2:] == '.Z':
                self.file = os.popen("uncompress -c " + self.filename, 'rb')
            elif self.filename[-3:] == '.gz':
                self.file = os.popen("gunzip -c " + self.filename, 'rb')
            else:
                try:
                    print( 'here1')
                    if zipf is not None:
                        self.file = self.zipf.open( self.filename, 'rb' )
                    else:
                        self.file = open(self.filename, 'rb')
                except(IOError, details):
                    if type(details) == type(()):
                        details = details + (self.filename,)
                    raise(IOError, details)

            self.filesize = os.path.getsize(self.filename)

        elif mode == 'w':
            if self.filename[-2:] == '.Z':
                self.file = os.popen("compress > " + self.filename, mode)
            elif self.filename[-3:] == '.gz':
                self.file = os.popen("gzip > " + self.filename, mode)
            else:
                try:
                    print( 'here2')
                    if zipf is not None:
                        self.file = self.zipf.open( self.filename, 'wb' )
                    else:
                        self.file = open(self.filename, 'wb')
                except(IOError, details):
                    if type(details) == type(()):
                        details = details + (self.filename,)
                    raise(IOError, details)
        else:
            raise IOError(0, 'Illegal mode: ' + repr(mode))

    def close(self):
        if not self.file.closed:
            self.file.close()

    def flush(self):
        self.file.flush()

    def readRecordNative(self, dtype=None):
        a = self.file.read(4)   # record size in bytes
        recordsize = np.frombuffer(a,'i')
        record = self.file.read(recordsize[0])
        self.file.read(4)   # record size in bytes

        if dtype in ('f', 'i', 'I', 'b', 'B', 'h', 'H',  'l', 'L', 'd'):
            return np.frombuffer(record,dtype)
        elif dtype in ('c', 'x'):
            return struct.unpack(self.endian+'1'+dtype, record)
        else:
            return (None, record)

    def readRecordByteswapped(self, dtype=None):
        a = self.file.read(4)   # record size in bytes
        recordsize = np.frombuffer(a,'i').byteswap()
        record = self.file.read(recordsize[0])
        self.file.read(4)   # record size in bytes

        if dtype in ('f', 'i', 'I', 'b', 'B', 'h', 'H',  'l', 'L', 'd'):
            return np.frombuffer(record,dtype).byteswap()
        elif dtype in ('c', 'x'):
            return struct.unpack(self.endian+'1'+dtype, record)
        else:
            return (None, record)

    def readBytes(self, recordsize, dtype, offset = None):
        if offset:
            self.file.seek(offset)
            record = self.file.read(recordsize*struct.calcsize(dtype))
        if dtype in ('b', 'B', 'h', 'H', 'i', 'I', 'l', 'L', 'f', 'd'):
            return self.fromstring(record, dtype)
        elif dtype in ('c', 'x'):
            return struct.unpack(self.endian+'1'+ttype, record)
        else:
            return (None, record)

### Test

if __name__ == '__main__':
    pass


