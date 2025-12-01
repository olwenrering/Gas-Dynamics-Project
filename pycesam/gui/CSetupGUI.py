#!/usr/bin/env python3
# coding: utf-8

# Code d'Evolution Stellaire, Adaptatif et Modulaire for the 2020 decade.
# Copyright (c) 1997-2023 The Cesam2k20 authors
# SPDX-License-Identifier : GPL-3.0-or-later

import os

from traits.api    import File, HasTraits
from traitsui.api   import FileEditor, View, Item
from traitsui.menu  import OKCancelButtons

class CSetupGUI(HasTraits):
    don = File(label='Data file')

    view = View( Item( 'don', editor=FileEditor( filter=['*.don'] ) ),
                     buttons=OKCancelButtons )

    def __init__(self):
        self.don = os.getcwd() + '/model.don'