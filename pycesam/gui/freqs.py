#!/usr/bin/env python3
# coding: utf-8

# Code d'Evolution Stellaire, Adaptatif et Modulaire for the 2020 decade.
# Copyright (c) 1997-2023 The Cesam2k20 authors
# SPDX-License-Identifier : GPL-3.0-or-later


import os
from pycesam.freqs import *
from traits.api import HasTraits, Trait, Float, Int, Bool, File, Range
from traitsui.api import View, Group, Item, HGroup
from traitsui.menu import OKCancelButtons

class FreqsGUI(Freqs):

    traits_view = \
        View(Item('osc_code_name'),
            Group(
                HGroup(
                    Item('to_zip'),
                    Item( 'adipls_new_amdl' ),
                    Item( 'adipls_amdl_name', visible_when="adipls_new_amdl" )),
                Item('modes', style='custom'),
                HGroup( Item('reduce' ), Item('nred', visible_when="reduc_!='none" ) ),
                Item('step', visible_when='all_osc'),
                HGroup(
                    Item('lmin'),
                    Item('lmax'),
                    Item('dels'),
                    Item('nsel', style='readonly'),
                    show_border=True, label='Options for l'),
                HGroup(
                    Item('sig1'),
                    Item('sig2'),
                    Item('iscan'),
                    Item('nsig'), show_border=True,
                    label='Frequency interval'),
                Item('remesh'),
                Item('cts', visible_when="remesh=='Customize'"),
                Item('npoints', visible_when="remesh!='Do not re-mesh'"),
                Item('irsord'),
                HGroup(Item('mdintg'), Item('iriche')),
                Item('iekinr'),
                Item('rotkr'),
                Item('gm1kr'),
                Item('amde'),
                Item('istsbc', tooltip="See Eqs. from 2.16-2.20 of Notes on adiabatic oscillation programme, 8th ed., JCD."),
                Item('fctsbc', visible_when="istsbc_ == 0 "),
                Item('fsbc', visible_when="fctsbc_ == 2"),label='ADIPLS run Options', visible_when="osc_code_name=='ADIPLS'"),
            Group(
                Group(
                    HGroup(
                        Item('input_file'),
                        Item('step1d', visible_when="input_file_ in ['d1', 'd12'] and not only_last"),
                        Item('step2d', visible_when="input_file_ in ['d2', 'd12'] and not only_last"),
                        Item('only_last')),
                    Item('acor_star_type'),
                    HGroup(
                        Item('mod_flag'),
                        Item('mod_flag_custom', visible_when="mod_flag=='Custom'"),
                        Item('mod_suff')),
                    Item('acor_exe_path'),
                    Item('out_freq_path'),
                    HGroup(
                        Item('acor_eigen'),
                        Item('out_eigen_path', visible_when="acor_eigen")),
                    Item('acor_rm_input'),
                    show_border=True, label='Paths and stellar model'),
                Group(
                    HGroup(
                        Item('acor_rot'),
                        Item('acor_rot_grad', visible_when="acor_rot"),
                        Item('acor_rot_toro', visible_when="acor_rot"),
                        Item('acor_rot_prof', visible_when="acor_rot")),
                    Item('acor_profile_ext'),
                    Group(
                        Item('profile_filename'),
                        Item('unif_om', visible_when="acor_rot_prof=='Uniform'"),
                        Item('diff_om', visible_when="acor_rot_prof=='Differential'"),
                        HGroup(
                            Item('biz_rad'),
                            Item('biz_oms'),
                            Item('biz_omc'),
                            visible_when="acor_rot_prof in ['Bizonal', 'Bizonal with smooth gradient']"),
                        HGroup(
                            Item('lat_omeq'),
                            Item('lat_dom'), visible_when="acor_rot_prof=='Latitudinal diff: Om = Om_eq - DOm np.cos^2 theta'"),
                        visible_when='acor_profile_ext'),
                    show_border=True, label='Rotation caracteristics'),
                Group(
                    Item('acor_freq_target'),
                    HGroup(
                        Item('acor_freq_units'),
                        Item('acor_freq_start'),
                        Item('acor_freq_step'),
                        Item('acor_nsteps'),
                        visible_when="acor_freq_target=='Scanning'"),
                    HGroup(
                        Item('acor_freq_file_type'),
                        Item('acor_freq_filename'),
                        visible_when="acor_freq_target=='File'"),
                    show_border=True, label='Frequency target'),
                Group(
                    Item('acor_params_select'),
                    HGroup(
                        Item('acor_nsh'),
                        Item('acor_m'),
                        Item('acor_parity'), visible_when="acor_params_select_==2"),
                    HGroup(
                        Item('lmax'),
                        Item('acor_m'), visible_when="acor_params_select_==1"),
                    show_border=True, label='Mode caracteristics'),
                Group(
                    HGroup(
                        Item('acor_rad_res'),
                        Item('acor_newt_it'),
                        Item('acor_follow')),
                    show_border=True, label='Computations parameters'),
                label='ACOR run Options', visible_when="osc_code_name=='ACOR'"),
            buttons = OKCancelButtons, kind='live', width=900, scrollable=True )

    remesh_view = \
        View(
            Item('remesh'),
            Item('cts', visible_when="remesh=='Customize'"),
            Item('npoints', visible_when="remesh!='Do not re-mesh'"),
            buttons=OKCancelButtons, kind='live', title='Remesh options',
            width=700, resizable=True)
