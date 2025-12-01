#!/usr/bin/env python3
# coding: utf-8

# Code d'Evolution Stellaire, Adaptatif et Modulaire for the 2020 decade.
# Copyright (c) 1997-2023 The Cesam2k20 authors
# SPDX-License-Identifier : GPL-3.0-or-later


from traitsui.api import View, Item, CodeEditor, \
    Tabbed, ListEditor, Group, DirectoryEditor, Include, HGroup, VGroup
from traitsui.menu import Action, LiveButtons, OKButton, OKCancelButtons

from pycesam.cparams import *

class CParametersGUI(CParameters):

    literal = Trait('Description of options', 'Function names',
                    label='Show options')
    basic   = Trait('Advanced', 'Basic', label='Show options')

    # NL Cesam ----------------------------------------------------
    nl_cesam = Group(Item('nom_chemin', editor=DirectoryEditor(),
                          visible_when="basic=='Advanced'"),
                     Item('nom_ctes',visible_when="basic=='Advanced'"),
                     HGroup(Item('trunc_nom_output'),
                           Group(
                           HGroup(Item('all_rep'),
                                  Item('acor')),
                           HGroup(Item('out_2d'),
                                  Item('out_h5')))),
                     Item('nom_osc',visible_when="basic=='Advanced'"),
                     Item('n_max',visible_when="basic=='Advanced'"),
                     Item('precision'),
                     label='Main')

    # NL Evol ----------------------------------------------------
    nl_evol = Group(
                VGroup(
                    Item('arret'),
                    Item('disclaim_evol', show_label=False, style='readonly')),
                VGroup(
                    Item('agemax'),
                    Item('log_teff'),
                    Item('he_core'),
                    Item('t_stop'),
                    Item('x_stop'),
                    Item('xcx0'),
                    Item('rhoc'),
                    Item('l_stop', tooltip="Given in solar unit (not log).\n \
                                            If l_stop < 0: stop when luminosity decreases below l_stop\n  \
                                            If l_stop > 0: stop when luminosity increases above l_stop"),
                    Item('r_stop', tooltip="Given in solar unit (not log).\n \
                                            If r_stop < 0: stop when luminosity decreases below r_stop\n  \
                                            If r_stop > 0: stop when luminosity increases above r_stop"),
                    Item('nb_max_modeles'),
                    Item('dtlist') ),label='Evolution')

    # NL Mass ----------------------------------------------------
    nl_mass = Group(
                    VGroup(
                        Item('mtot'),
                        Item('disclaim_mass', show_label=False, style='readonly')),
                    Item('nom_pertm_solar'),
                    Item('nom_pertm_rgb'),
                    Item('nom_pertm_agb'),
                    Item('nom_pertm_hot'),
                    Item('mdot'), label='Mass')

    # NL Chim ----------------------------------------------------
    nl_chim = Group(Item('garde_xish',visible_when="basic=='Advanced'"),
                    Item('grille_fixe',visible_when="basic=='Advanced'"),
                    Item('modif_chim',visible_when="basic=='Advanced'"),
                    Item('nom_abon'),
                    HGroup(Item('x_or_y', show_label=False, width=250),
                           Item('x0', show_label=False, visible_when="x_or_y=='Initial X'",
                                width=250),
                           Item('y0', show_label=False, visible_when="x_or_y=='Initial Y'",
                                width=250)),
                    HGroup(Item('z_or_zsx', show_label=False, width=250),
                           Item('z0', show_label=False, visible_when="z_or_zsx=='Initial Z'",
                                width=250),
                           Item('zsx0', show_label=False,
                                visible_when="z_or_zsx=='Initial Z/X'", width=250)),
                           Item('x0', style='readonly'),
                           Item('y0', style='readonly'),
                           Item('z0', style='readonly'),
                           Item('zsx0', style='readonly'),label='Chemical composition')

    # NL Conv ----------------------------------------------------
    nl_conv = Group(Item('nom_conv'),
                    HGroup(VGroup(Item('nom_alpha'),
                                  Item('nom_ovshts')),
                           VGroup(Item('alpha', visible_when="nom_alpha_ in ['alpha_f', 'entropy_ludwig99', 'entropy_magic13', 'entropy_tanner16', 'entropy_custom']"),
                                  Item('ovshts', visible_when="nom_ovshts_ in ['ovshts_f','ovshts_f_diff']"),
                                  Item('ovshti'))),
                    Item('jpz',visible_when="basic=='Advanced'"),
                    Item('cpturb',visible_when="basic=='Advanced'"),
                    Group( HGroup(
                        Item('ledoux'),
                        Item('20'),
                        Item('solberg_hoiland'),
                        Item('20'),
                        Item('gough_tayler')), label='Convection criterion (default: Schwarzschild):'),label='Convection')
    # NL Diff ----------------------------------------------------
    nl_nodiff = Group(Item('diffusion'),label='Diffusion')
    nl_diff = Group(Item('diffusion'),
                    Group( Item('nom_diffm'),
                           Item('nom_difft'),
                           Item('nom_difft_fing'),
                           Item('nom_difft_smc'),
                           Item('d_turb'),
                           Item('re_nu'),
                           Item('nom_frad'), visible_when='diffusion',
                           ), label='Diffusion')

    # NL Rot ----------------------------------------------------
    nl_rot = Group(VGroup(
                   Item('w_rot'),
                   Item('unit'),
                   Item('nom_thw')),
                   Group(HGroup(VGroup(Item('rot_zc'),
                                       Item('nom_diffmg_ts')),
                                VGroup(Item('nom_diffw'),
                                       Item('nom_gsf'))),
                         HGroup(Item('init_rot'),
                         Item('g_modes'),
                         Item('mixed_modes')),
                         visible_when="nom_thw not in ['None', 'Solid body', 'Local conservation of angular momentum']"),
                   Group(Item('nom_pertw'),
                         Item('p_pertw'),
                         Item('w_sat'),
                         label='Loss of angular momentum',show_border=True,
                         visible_when="nom_thw != 'None'"),
                   Group(Item('tau_disk'),
                         Item('p_disk'),
                         label='Initial conditions',show_border=True,
                         visible_when="nom_thw != 'None'")
                   , label='Rotation', scrollable=True)
    # NL Etat -------------------------------------------------
    nl_etat = Group(Item('nom_etat'),
                    Item('f_eos'), label='EoS')

    # NL Opa -------------------------------------------------
    nl_opa = Group(Item('nom_opa'),
                   Item('nom_opa_cond'),
                   Item('f_opa'),label='Opacity')

    # NL Nuc -------------------------------------------------
    nl_nuc = Group(Item('nom_nuc'),
                   Item('nom_nuc_cpl'),
                   Item('nom_neutrinos'),
                   Item('mitler'),label='Nuclear')

    # NL Atm -------------------------------------------------
    nl_atm = Group(Item('nom_atm'),
                   Group(Item('nom_tdetau'),
                         Item('tau_max', visible_when="basic=='Advanced'"),
                         Item('lim_ro', visible_when="basic=='Advanced'"),
                         visible_when="nom_atm == u'Restitution of the atmosphere from a T(\u03C4) law'"),
                   label='Atmosphere')
    # NL Mag -------------------------------------------------
    nl_mag = Group(Item('nom_mag'),
                   Group(Item('mag_inhib', visible_when="nom_mag_ in ['mag_inhib', 'mag_inhib_rmax', 'mag_inhib_brmax']"),
                         VGroup(Item('br'),
                                Item('bphi'),
                                Item('btheta'), visible_when="nom_mag_ in ['mag_f', 'mag_inhib_rmax']"),
                         VGroup(Item('xmax'),
                                Item('xmax_unit'), visible_when="nom_mag_=='mag_inhib_rmax'"),
                        VGroup(Item('brmax'), visible_when="nom_mag_=='mag_inhib_brmax'"),
                         visible_when="nom_mag_ != 'mag_0'"),
                   label='Magnetic field')

    # NL regl -------------------------------------------------
    nl_regl = Group(
        HGroup( Item('disclaim_num', show_label=False, style='readonly') ),
        HGroup( Item('precision') ),
        Group(
            HGroup(
            VGroup(
                Item( "ini0"       ), Item( "m_qs"       ), Item( "m_rot"      ), Item( "n_atm" ),
                Item( "ordre"      ), Item( "yld_rep"    ), Item( "m_qs2d"     )),
            VGroup(
                Item( 'l0'         ), Item( "m_ch"       ), Item( "m_tds"      ), Item( "30"    ),
                Item( '30'         ), Item( "osc2d_step" ), Item( "m_rot2d"    )),
            VGroup(
                Item( '30'         ), Item( "m_vth"      ), Item( "m_ptm"      ), Item( "n_min" ), Item( "n_min_post" ),
                Item( '30'         ), Item( "osc_step"   ), Item( "m_legendre" ))), label='Integers' ),

        Group(
            HGroup(
            VGroup(
                Item( "ctel"       ), Item( "cter"       ), Item( "d_grav"     ), Item( "coeff_d_grav_post"     ),
                Item( "dlntc"     ), Item( "dlogg"      ), Item( "dw"         ), Item( "evolved"    ),
                Item( "dt0"       ), Item( "dt_m"       ), Item( "fmin_abon"  ), Item( "loc_zc"     ),
                Item( "psi0"      ), Item( "yld_dx"     ), Item( "x_tams"     ) ),
            VGroup(
                Item( "ctep"       ), Item( "ctet"       ), Item( "d_grav_int" ), Item( "drhoc"     ),
                Item( "dalpha"     ), Item( "drhoca"     ), Item( "dn_fixe"    ), Item( "dtmax"     ),
                Item( '30'         ), Item( '30'         ), Item( "precix"     ), Item( "ro_test"   ),
                Item( "yld_rac"    ), Item( "dx_tams"    ) ),
            VGroup(
                Item( "ctem"       ), Item( '30'         ), Item( '30'         ), Item( "dteff"     ),
                Item( "dlum"       ), Item( "dsenv"      ), Item( "dpsi"       ), Item( "dtmin"     ),
                Item( '30'         ), Item( '30'         ), Item( "precit"     ), Item( "q0"        ),
                Item( "tau_min"    ), Item( "y_agb"      ))), label="Floats" ),
        Group(
            HGroup(
            VGroup( Item( "exact_stop" ), Item( "new_bv"    ), Item( "simple"         ), Item( "out_pms"      )),
            VGroup( Item( "en_masse"   ), Item( "yld_hr"    ), Item( "extended_osc2d" ), Item( "g_modes_ab"   )),
            VGroup( Item( "kipp"       ), Item( "yld_l23"   ), Item( "vdiff_out"      ), Item( "clean"   )),
            VGroup( Item( "lisse"      ), Item( "yld_fgong" ), Item( "grad_out"       )),
            VGroup( Item( "mu_saha"    ), Item( "d2"        ), Item( "th_out"         )),
            VGroup( Item( "mvt_dis"    ), Item( "general"   ), Item( "eos_inplace"    )) ), label="Boolean"),
        Group( Item( "ltau_rep" ), Item( "ec_s" ), Item( "lim_zc" ), label="Strings" ),
        label='Numerical', orientation='vertical', scrollable=True )

    # NL extra  -------------------------------------------------
    nl_extra = Group(
        HGroup( Item('disclaim_extra', show_label=False, style='readonly') ),
        Group( HGroup(
            # VGroup( Item( 'source_cts' ) ),
            VGroup( Item( '30'      ) ),
            VGroup( Item( '30'      ) ),
            VGroup( Item( '30'      ) )),
            label='Physical constants'),
        Group(
            HGroup( Item('disclaim_osm', show_label=False, style='readonly') ),
            HGroup(
            VGroup( Item( 'osm_p_pertw' ) ),
            VGroup( Item( '30'      ) ),
            VGroup( Item( '30'      ) ),
            VGroup( Item( '30'      ) )),
            label='OSM renormalization factor'),
        Group( HGroup(
            VGroup( Item( 'debug'          ), Item( 'ts_out'   ), Item( 'hse_out' ) ),
            VGroup( Item( 'convection_out' ), Item( 'rota_out' ) ),
            VGroup( Item( 'lim_zc_out'     ), Item( '30'       ) ),
            VGroup( Item( 'plato_track'    ), Item( '30'       ) )),
            label='Extra outputs'),
        Group( HGroup( Item( 'full_hr' ), Item( 'writes_tau_conv' ), Item( 'writes_conv' ), Item( 'writes_nuc_neutr' ),
                Item( 'writes_burn_zones' ), Item( 'writes_2d' ), Item( 'writes_l_nuc' ), Item( 'writes_mdot' ) ),
            Item( 'n_contours_burn' ), Item( 'burn_contours'   ),
            label='HRnew output'),
        Group( HGroup(
            VGroup( Item( 'rho_ext'    ) ),
            VGroup( Item( 'easy_krisw' ) ),
            VGroup( Item( '30'         ) ),
            VGroup( Item( '30'         ) )),
            label='Atmosphere'),
        Group(HGroup(
            VGroup( Item( 'eta_vdif'   ) ),
            VGroup( Item( 'no_diff_post'   ) ),
            VGroup( Item( 'b_strength' ) ),
            VGroup( Item( '30'         ) ),
            VGroup( Item( '30'         ) )),
            label='Atomic diffusion'),
        Group(HGroup(
            VGroup( Item( 'eta_grad' ) ),
            VGroup( Item( 'goop' ) ),
            VGroup( Item( '30' ) ),
            VGroup( Item( '30' ) )),
            label='Radiative diffusion'),
        Group(HGroup(
            VGroup( Item( 't_mix' ),    Item( 'd_expo_2' ), Item( 'eta_vdif' ), Item( 'goop' ),     Item( 'teff_zams_sun' ) ),
            VGroup( Item( 'm_mix' ),    Item( 'l_expo' ),   Item( '30' ), Item( 'errabop' ),  Item( 'teff_zams_2' )   ),
            VGroup( Item( 'expo_mix' ), Item( 'l_expo_2' ), Item( '30' ),       Item( 't_op_sup' ), Item( 't0_sun' )        ),
            VGroup( Item( 'cst_mix' ),  Item( '30' ),       Item( '30' ),       Item( 't_op_inf' ), Item( 't0_2' )          )),
            label='Turbulent diffusion'),
        Group( HGroup(
            VGroup( Item( 'n_fake_ce'   ) ),
            VGroup( Item( 'd_conv'      ) ),
            VGroup( Item( 'mass_tol_cz' ) ),
            VGroup( Item( 'conv_block'  ) ),
            VGroup( Item( 'no_tot_rad_post'  ) )),
            label='Convection'),
        Group( HGroup(
            VGroup( Item( 'ds' ),              Item( 'ec_extrap' ), Item( 'alpha_m10' ), Item( 'a_m10_z005' ) ),
            VGroup( Item( 'ec_grid' ),         Item( '30' ),        Item( 'alpha_m15' ), Item( 'a_m15_z005' ) ),
            VGroup( Item( 'ec_corr' ),         Item( '30' ),        Item( '30' ),        Item( 'a_m10_z02' ) ),
            VGroup( Item( 'alpha_autoguess' ), Item( '30' ),        Item( '30' ),        Item( 'a_m15_z02' ) )),
            label='Entropy calibration'),
        Group( HGroup(
            VGroup( Item( 't_ov_inf' ),          Item( 'alpha_jcd11' ), Item( 'claret18_shift' ), Item( 'ovs_type' ) ),
            VGroup( Item( 'start_ovshti_zams' ), Item( 'beta_jcd11' ),  Item( '30' ),             Item( 'ovi_type' ) ),
            VGroup( Item( '30' ),                Item( '30' ),          Item( '30' ),             Item( '30' ) ),
            VGroup( Item( '30' ),                Item( '30' ),          Item( '30' ),             Item( '30' ) )),
            label='Overshoot and overshooting regions'),
        Group( HGroup(
            VGroup( Item( 'mass_loss_eta' ),   Item( 'v_inf_on_v_esc_hot' ) ),
            VGroup( Item( 'mass_loss_zeta' ),  Item( 'v_inf_on_v_esc_cold' ) ),
            VGroup( Item( 'log_lum_min_jvh' ), Item( 'eta_r_rgb_1' ) ),
            VGroup( Item( 'log_lum_max_jvh' ), Item( 'eta_r_rgb_2' ) )),
            label='Mass loss'),
        Group( HGroup(
            VGroup( Item( 'errabop' ) ),
            VGroup( Item( 't_op_sup' ) ),
            VGroup( Item( 't_op_inf' ) ),
            VGroup( Item( '30' ) )),
            label='Opacities OP'),
        Group( HGroup(
            VGroup( Item( 't_oplib_sup' ) ),
            VGroup( Item( 't_oplib_inf' ) ),
            VGroup( Item( 'use_findne' ) )),
            label='Opacities OPLIB'),
        Group( HGroup(
            VGroup( Item( 'lt_compton_l' ) ),
            VGroup( Item( 'lt_compton_h' ) ),
            VGroup( Item( 't_compton' ) ),
            VGroup( Item( '30' ) )),
            label='Compton opacities'),
        Group( HGroup(
            VGroup( Item( 'npt_central_cz' ), Item( 'nu_v_fact' ), Item( 'nu_h_fact' ), Item( 'omega_sun' ), Item( 're_c' ) ),
            VGroup( Item( 'nmore_shells' ),   Item( 'nu_v_add' ),  Item( 'nu_h_add' ),  Item( 'fast_rot' ),  Item( 'ri_c' )  ),
            VGroup( Item( '30' ),             Item( 'dv_fact' ),   Item( 'deff_fact' ), Item( '30' ),        Item( '30' ) ),
            VGroup( Item( '30' ),             Item( 'dv_add' ),    Item( 'deff_add' ),  Item( '30' ),        Item( '30' )  )),
            label='Rotation'),
        Group( HGroup(
            VGroup( Item( 'pertw_profile'  ) ),
            VGroup( Item( 'p_power_matt15' ) ),
            VGroup( Item( 'tau_conv_sun'   ) ),
            VGroup( Item( '30' ) )),
            label='Angular momentum loss'),
        Group( HGroup(
            VGroup( Item( 'ts_dt_max', width=0.2 ),         Item( '30'           ),           Item( 'ts_weak_b' ),    Item( 'ts_extrapolate' ) ),
            VGroup( Item( 'lsmooth_fact', width=0.2 ),      Item( '30'           ),           Item( 'ts_hysteresis' ) ),
            VGroup( Item( 'ts_average', width=0.2 ),        Item( 'lsmooth_min', width=0.2 ), Item( 'ts_patch_regime' ) ),
            VGroup( Item( 'ts_average_length', width=0.2 ), Item( 'alpha_fuller' ),           Item( 'ts_qlim' ) )),
            label='Tayler-Spruit instability'),
        Group( HGroup(
            VGroup( Item( 'ell_max' ),      Item( 't0_corsaro' ), Item( 'ratio_work_max' ) ),
            VGroup( Item( 'ell_min' ),      Item( 'eta_sun' ),    Item( '30' ) ),
            VGroup( Item( 'ndelta' ),       Item( 'eta_0' ),      Item( '30' ) ),
            VGroup( Item( 'mode_damping' ), Item( '30' ),         Item( '30' ) )),
            label='Mixed modes'),
        Group(HGroup(
            VGroup( Item( 'r_transition' ) ),
            VGroup( Item( 'width_transition' ) ),
            VGroup( Item( 'omega_c' ) ),
            VGroup( Item( 'omega_s' ) )),
            label='Setup for ad hoc rotation profile'),
        Group(HGroup(
            VGroup( Item( 'struct2d' ) ),
            VGroup( Item( 'struct2d_file' ) )),
            label='Structure model to be deformed'),
        Group( HGroup( 
             VGroup(Item( 'screen' ),  Item( 'eta_nuc' ))),
            label='Nuclear reactions'),
        Group( HGroup(
            VGroup( Item( 'r_h2h' ) ),
            VGroup( Item( 'r_14n15o' ) ),
            VGroup( Item( '30' ) ),
            VGroup( Item( '30' ) )),
            label='Modification of individual nuclear reaction rates'),
        Group( HGroup(
            VGroup( Item( 'inj_a' ) ),
            VGroup( Item( 'inj_mu' ) ),
            VGroup( Item( '30' ) ),
            VGroup( Item( '30' ) )),
            label='Custom energy injection below CZ'),
        label='Fine tuning', orientation='vertical', scrollable=True )

    traits_view = View(Item('basic', style='custom', show_label=False),
                       Tabbed(Include('nl_cesam'),
                              Include('nl_mass'),
                              Include('nl_evol'),
                              Include('nl_chim'),
                              Include('nl_conv'),
                              Include('nl_diff'),
                              Include('nl_rot'),
                              Include('nl_etat'),
                              Include('nl_opa'),
                              Include('nl_nuc'),
                              Include('nl_atm'),
                              Include('nl_mag'),
                              Include('nl_regl'),
                              Include('nl_extra'),
                              show_border=True,scrollable=True),
                       buttons = LiveButtons, kind='live',
                       title='Cesam2k20 Options')

    def __init__( self, name ):
        """
        Class that represents the input parameters of a Cesam2k20 model.

        :param name: Name of the model.
        :type name: string
        """
        self.initialized = False
        super().__init__( name )


        if not self.nom_chemin.endswith('/'):
            self.nom_chemin = self.nom_chemin + '/'

        self.f_eos[0] = self.nom_chemin + '/EOS/eos_opal2005/' + self.f_eos[0]
        self.f_opa[0] = self.nom_chemin + '/OPA/opa_yveline/' + self.f_opa[0]
        self.f_opa[1] = self.nom_chemin + self.f_opa[1]
        self.f_opa[2] = self.nom_chemin + self.f_opa[2]
        self.f_opa[3] = self.nom_chemin + self.f_opa[3]

        self.initialized = True

