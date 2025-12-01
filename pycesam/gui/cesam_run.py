#!/usr/bin/env python3
# coding: utf-8

# Code d'Evolution Stellaire, Adaptatif et Modulaire for the 2020 decade.
# Copyright (c) 1997-2023 The Cesam2k20 authors
# SPDX-License-Identifier : GPL-3.0-or-later


from pycesam.cesam_run import *

from traits.api import HasTraits, Trait, Float, File, Bool
from traitsui.api import View, Group, Item, FileEditor
from traitsui.menu import OKCancelButtons

class CRunGUI(CRun):

     view = View(Item('job'),
                 Item('nohup'),
                 Group(Item('type_file', style='custom',
                            show_label=False,
                            enabled_when="job != 'From previous model'"),
                       Item('mod_init',
                            editor=FileEditor(filter=['*.pms']),
                            visible_when="job == 'From PMS'"),
                       Item('mod_init',
                            editor=FileEditor(filter=['*.zams']),
                            visible_when="job == 'From ZAMS'"),
                       Item('mod_init',
                            editor=FileEditor(filter=['*.rep', '*.dat']),
                            visible_when="job == 'From previous model'"),
                       visible_when="job != 'Frequencies'"),
                 Group( Item('dt0',
                      visible_when="job in ['From ZAMS', 'From previous model']") ),
                 Group( Item('li_zams',
                      visible_when="job == 'From ZAMS'") ),
                 Group( Item('c_iben',
                      visible_when="job == 'From PMS'") ),
                 Group( Item('pre_pms',
                      visible_when="job == 'From PMS'") ),
                 buttons = OKCancelButtons, kind='live',
                 title='Cesam2k20 run Options',width=700)

