#!/usr/bin/env python3
# coding: utf-8

# Code d'Evolution Stellaire, Adaptatif et Modulaire for the 2020 decade.
# Copyright (c) 1997-2023 The Cesam2k20 authors
# SPDX-License-Identifier : GPL-3.0-or-later


from pycesam import *
from pycesam.tools import *
from threading import Thread
import matplotlib.pyplot as plt

class RunParallel(Thread):
    def __init__(self, model):
        self.model = model
        self.model.reinit = True
        Thread.__init__(self)

    def run(self):
        self.model(mkdon=False)
        if not self.model.finished:
            self.model('rep', dt0=1.0)
            if not self.model.finished:
                print('Error, model not calculated!')

        self.finished = self.model.finished
                

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def find_models(*args, **kargs):
    """find_models(*args, **kargs)

Finds all models with a calculated HR file; reads HR file if required.
Returns a list of CModel instances with models found.

Parameters
----------
args : strings, optional
    Selects models; Example: find_models('1msun') only finds models
    with '1msun' in name.
kargs : optional keyword arguments
    read=True/False: read/doesn't read HR file
    obs=True/False: compute observational CMD
        
Examples
--------
Finding all models in current folder, m will contain models:

>>> m = find_models()

Finding models whose name contains '1.0Msun' and '_ov0.2':

>>> m = find_models('1.0Msun', '_ov0.2')

Finding models whose name contains '1.0Msun' and computing CMD:

>>> m = find_models('1.0Msun', obs=True)

    """
    
    names = []
    if len(args) > 0:
        for j in args:
            names.extend([i[:-4] for i in os.listdir('.')
                          if i.endswith('.don') and j in i])
    else:
        names.extend([i[:-4] for i in os.listdir('.')
                      if i.endswith('.don')])

    names.sort()
    
    m = [CModel(i, read=False) for i in names]

    if len(kargs) > 0:
        if 'read' in kargs:
            if kargs['read']:
                for i in m:
                    if 'obs' in kargs:
                        i.read_hr(obs=kargs['obs'])

        else:
            i.read_hr()

    return m



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def run_parallel(model_names):
    procs = []
    status = []
    for model_name in model_names:
        model = CModel(model_name, read=False)
        procs.append(RunParallel(model))
        procs[-1].start()
        
    for proc in procs:
        proc.join()

    for proc in procs: status.append(proc.finished)

    if not all(procs):
        print('There was an error... Not all processes completed sucessfully.')

    return status



