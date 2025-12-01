#!/usr/bin/env python3
# coding: utf-8

# Code d'Evolution Stellaire, Adaptatif et Modulaire for the 2020 decade.
# Copyright (c) 1997-2023 The Cesam2k20 authors
# SPDX-License-Identifier : GPL-3.0-or-later


from pyface.qt import QtGui, QtCore

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT

from traits.api import Any, Instance
try:
    from traitsui.qt.editor import Editor
except ModuleNotFoundError:
    from traitsui.qt4.editor import Editor
from traitsui.basic_editor_factory import BasicEditorFactory
from traitsui.api import Handler

class _MPLFigureEditor(Editor):

   scrollable  = True

   def init(self, parent):
        matplotlib.use('Qt5Agg')  # wxPython does not work...

        self.control = self._create_canvas(parent)
        self.set_tooltip()

   def update_editor(self):
        pass

   def _create_canvas(self, parent):
        """ Create the MPL canvas. """
        # matplotlib commands to create a canvas
        frame = QtGui.QWidget()
        mpl_canvas = FigureCanvas(self.value)
        mpl_canvas.setParent(frame)
        mpl_toolbar = NavigationToolbar2QT(mpl_canvas,frame)

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(mpl_canvas)
        vbox.addWidget(mpl_toolbar)
        frame.setLayout(vbox)

        return frame

class MPLFigureEditor(BasicEditorFactory):

    klass = _MPLFigureEditor

class MPLInitHandler(Handler):
    """
    Handler calls mpl_setup() to initialize mpl events
    """

    def init(self, info):
        """
        This method gets called after the controls have all been
        created but before they are displayed.
        """
        info.object.mpl_setup()
        return True


class MPLHandler(Handler):
    """
    Handler for closing matplotlib window when embedded in a TraitsUI object.
    """

    def close(self, info, is_ok):
        plt.close(info.object.figure)

        return True
