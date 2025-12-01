#!/usr/bin/env python3
# coding: utf-8

# Code d'Evolution Stellaire, Adaptatif et Modulaire for the 2020 decade.
# Copyright (c) 1997-2023 The Cesam2k20 authors
# SPDX-License-Identifier : GPL-3.0-or-later


import signal
import numpy as np
import os, sys
import subprocess as subp
from traits.trait_list_object import TraitListObject as ListType

def get_hash( ):
    cwd = os.getcwd( )
    os.chdir( os.environ['CESDIR'] )
    githash = subp.check_output(['git', 'rev-parse', 'HEAD']).decode('ascii').strip()
    os.chdir( cwd )
    return githash

class Data:

    def __init__(self, header):

        self.header = header

class CESAMError( Exception ):

    __module__ = 'pycesam'

    def __init__( self, message ):
        self.message = message
        super().__init__( self.message )

class CESAMGUIError( Exception ):

    __module__ = 'pycesam.gui'

    def __init__( self, message ):
        self.message = message
        super().__init__( self.message )

def savitzky_golay(data, kernel=11, order=4):
    """
    <savitzky_golay(data, kernel=11, order=4)>
    Applies a Savitzky-Golay filter

    :param data: Input data
    :type data: np.array

    :kparam kernel: A positive integer > 2*order giving the kernel size
    :ktype kernel: int

    :kparam order: Order of the polynomal
    :ktype order: int

    :return: Smoothed data.
    :rtype: np.array

    :example:
    >>> smoothed = savitzky_golay(<rough>, [kernel = value], [order = value]

    :raise TypeError: `kernel` and `order` must be of type int.
    :raise ValueError: kernel must be positive, odd and less than order + 2.
    """
    try:
        kernel = abs(int(kernel))
        order  = abs(int(order))
    except(ValueError, msg):
        raise TypeError("kernel and order have to be of type int (floats will be converted).")
    if kernel % 2 != 1 or kernel < 1:
        raise ValueError(f"kernel size must be a positive odd number, was: {kernel}")
    if kernel < order + 2:
        raise ValueError("kernel is to small for the polynomals\nshould be > order + 2")

    # a second order polynomal has 3 coefficients
    order_range = range(order+1)
    half_window = (kernel -1) // 2
    b           = np.asmatrix([[k**i for i in order_range]
        for k in range(-half_window, half_window+1)])
    # since we don't want the derivative, else choose [1] or [2], respectively
    m           = np.linalg.pinv(b).A[0]
    window_size = len(m)
    half_window = (window_size-1) // 2

    # precompute the offset values for better performance
    offsets     = range(-half_window, half_window+1)
    offset_data = zip(offsets, m)

    smooth_data = list()

    # temporary data, with padded zeros (np.since we want the same length after
    # smoothing)
    data = np.concatenate((np.zeros(half_window), data, np.zeros(half_window)))
    for i in range(half_window, len(data) - half_window):
        value = 0.0
        for offset, weight in offset_data:
            value += weight * data[i + offset]

        smooth_data.append(value)

    return np.array(smooth_data)


def isbool( value ):
    """
    Function that determines if a string represents a boolean.
    Warning: this function does not evaluate the boolean value of the string if the string
    is indeed a boolean. See `tools.str2bool`.

    :param value: The input value
    :type value: str

    :returns: True if the string is a boolean.
    :rtype: bool
    """
    if not isinstance( value, str ):
        return False
    return value.lower() in ['true', '.true.', 't', 'y', 'yes',
        'false', '.false.', 'f', 'n', 'no']

def str2bool( value ):
    """
    Converts a string to boolean.

    :param value: The input value.
    :type value: str

    :returns: The boolean value of the string.
    :rtype: bool
    """
    if   value.lower() in ['true',  '.true.',  't', 'y', 'yes', '1']:
        return True
    elif value.lower() in ['false', '.false.', 'f', 'n', 'no', '0']:
        return False

def str2list( value ):
    """
    Converts the string representation of a int/float list to an int/float list.

    :param value: The input value.
    :type value: str

    :returns: A list corresponding to the string representation.
    :rtype: list
    """
    value = value.split( '[' )[1].split( ']' )[0]
    value = value.split(',')
    isfloat = False
    for v in value:
        if '.' in v:
            isfloat = True

    if isfloat:
        return [float( v ) for v in value]
    else:
        return [int( v ) for v in value]

def str2type( attribute, value ):
    """
    Converts a string representation to a variable of the correct type.

    :param attribute: The attribute or the variable to which the variable must be
        stored. This attributes defines the output type
    :type attribute: int, str, float, bool or list-like

    :param value: The string to be converted.
    :type value: str

    :returns: A value corresponding to the string representation.
    :ktype: int, str, float, bool or list-like
    """
    atype = type( attribute )
    if atype == str:
        return value
    elif atype == int:
        return int( value )
    elif atype == float:
        return float( value )
    elif atype == bool:
        return str2bool( value )
    elif atype in [ListType, list]:
        return str2list( value )

def swap( a, b ):
    """
    Swaps the value of two variables.

    :param a,b: The value that should be exchanged.
    :type a, b: any type

    :returns: A tuple of the two values in the reverse order.
    :rtype: any type
    """
    return b, a

getattr_l  = lambda ll, a: np.array([getattr( l, a ) for l in ll])
inline_str = lambda ll, fmt='{:s}': " ".join( fmt.format(l) for l in ll )
inline     = lambda ll, fmts='{:>10}', fmtf='{:>8g}': " ".join( fmts.format( fmtf.format(l) ) for l in ll)

# P2: 2nd order legendre polynomial
lp2cos  = lambda x: 0.5 * (3.0 * np.cos(x)**2 - 1.0)
# dP2 / dcos(theta)
dlp2cos = lambda x:-3.0 * np.cos(x) * np.sin(x)
# compute sigma**2 from frequency
sigma2 = lambda nu, mass, rad, grav: (nu*2*np.pi)**2 * rad**3 / grav / mass



class AlarmException(Exception):
    pass

def alarmHandler(signum, frame):
    raise AlarmException

def nonBlockingInput(prompt, timeout, default=''):
    signal.signal(signal.SIGALRM, alarmHandler)
    signal.alarm(timeout)
    try:
        text = input(prompt)
        signal.alarm(0)
        return text
    except AlarmException:
        pass
    signal.signal(signal.SIGALRM, signal.SIG_IGN)
    return default

def email_attachement( file ):
    from email.mime.base      import MIMEBase
    from email.encoders       import encode_base64
    part = MIMEBase('text', "plain")
    part.set_payload( open( file, "rb" ).read() )
    encode_base64(part)
    part.add_header('Content-Disposition', f'attachment; filename="{os.path.basename( file )}"')

    return part

def email( mdl, email_adress, text ):
    from email.mime.text      import MIMEText
    from email.mime.multipart import MIMEMultipart
    import smtplib
    sys.stdout.flush()
    # Create the container (outer) email message.
    s = 'Cesam2k20 bug report'
    msg            = MIMEMultipart()
    msg['Subject'] = 'Cesam2k20 bug report'
    msg['From']    = email_adress
    msg['To']      = 'troubleshootingcesam2k20@disroot.org'
    msg.preamble   = text

    msg.attach( email_attachement( mdl.don ) )
    msg.attach( email_attachement( mdl.err_file ) )

    # Send the email via our own SMTP server.
    s = smtplib.SMTP('localhost')
    try:
        s.sendmail( email_adress, 'troubleshootingcesam2k20@disroot.org', msg.as_string())
        s.quit()
    except:
        pass