#!/usr/bin/env python3
# coding: utf-8

# Code d'Evolution Stellaire, Adaptatif et Modulaire for the 2020 decade.
# Copyright (c) 1997-2023 The Cesam2k20 authors
# SPDX-License-Identifier : GPL-3.0-or-later


import os
import re
import subprocess as subp
from copy import deepcopy
from traits.api import HasTraits, Int, Trait, Float, Instance, Code, \
    String, Bool, Range, List, TraitError, File, Directory

from pycesam.constants import *
from pycesam.defaults import *
from pycesam.tools import *

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

dict_bool = {True:"T", False:"F"}
pat_char = r'=\s*(\b[TFtf]\b)'
pat_bool = r'=\s*(\.(true|false)\.)'
pat_str = r'\'([^,\']+)\''
pat_flt = r'=\s*([+-]?\d+((\.)?\d*([eEdD][+-]?\d+)?)?)'

FS = ['f_eos', 'f_opa']
FI = ['iter_qs']
FF = ['burn_contours']


class CParameters(HasTraits):

    # Contents of namelists ------------------------------------------
    nlist = {'nl_cesam' : ['nom_chemin', 'nom_ctes',
                           'nom_output', 'nom_osc', 'acor', 'n_max', 'precision'],
             'nl_mass'  : ['mtot', 'nom_pertm', 'nom_pertm_solar', 'nom_pertm_rgb',
                           'nom_pertm_agb', 'nom_pertm_hot', 'mdot'],
             'nl_evol'  : ['agemax', 'arret', 'dtlist', 'log_teff',
                           'nb_max_modeles', 'he_core', 't_stop',
                           'x_stop', 'xcx0', 'rhoc', 'l_stop', 'r_stop'],
             'nl_chim'  : ['grille_fixe', 'nom_abon', 'modif_chim',
                           'garde_xish', 'x0', 'y0', 'zsx0'],
             'nl_conv'  : ['nom_conv', 'nom_alpha', 'alpha', 'nom_ovshts', 'ovshts', 'ovshti',
                           'jpz', 'cpturb', 'ledoux', 'solberg_hoiland', 'gough_tayler'],
             'nl_diff'  : ['diffusion', 'nom_diffm', 'nom_difft', 'nom_difft_fing', 'nom_difft_smc',
                           'd_turb', 're_nu', 'nom_frad'],
             'nl_rot'   : ['w_rot', 'unit', 'nom_diffw', 'nom_diffmg_ts', 'nom_gsf',
                           'nom_thw', 'rot_zc', 'nom_pertw', 'p_pertw',
                           'lim_jpz', 'g_modes', 'mixed_modes', 'init_rot',
                           'tau_disk', 'p_disk', 'w_sat'],
             'nl_etat'  : ['nom_etat', 'f_eos'],
             'nl_opa'   : ['nom_opa', 'nom_opa_cond', 'f_opa'],
             'nl_nuc'   : ['nom_nuc', 'nom_nuc_cpl', 'nom_neutrinos', 'mitler'],
             'nl_atm'   : ['nom_atm', 'nom_tdetau', 'tau_max',
                           'lim_ro'],
             'nl_mag'   : ['nom_mag', 'mag_inhib', 'br', 'bphi', 'btheta', 'xmax', 'xmax_unit', 'brmax']}

    nlsublist = {
    'num': ['ini0', 'l0', 'm_qs', 'm_ch', 'm_vth', 'm_rot', 'm_tds', 'm_ptm', 'n_atm', 'n_min', 'n_min_post',
        'ordre', 'yld_rep', 'osc2d_step', 'osc_step', 'm_qs2d', 'm_rot2d', 'm_legendre',
        'n_theta', 'iter_qs', 'ctel', 'ctep', 'ctem', 'cter', 'ctet', 'd_grav', 'coeff_d_grav_post', 'd_grav_int',
        'dlntc', 'drhoc', 'dteff', 'dlogg', 'dalpha', 'dlum', 'dw', 'drhoca', 'dsenv',
        'evolved', 'dn_fixe', 'dpsi', 'dt0', 'dtmax', 'dtmin', 'fmin_abon', 'loc_zc', 'precix',
        'precit', 'psi0', 'ro_test', 'q0', 'yld_dx', 'yld_rac', 'tau_min', 'dt_m', 'x_tams',
        'y_agb', 'dx_tams', 'mcz_ext_min', 'exact_stop', 'en_masse', 'kipp', 'lisse',
        'mu_saha', 'mvt_dis', 'new_bv', 'yld_hr', 'yld_l23', 'yld_fgong', 'd2', 'general',
        'simple', 'extended_osc2d', 'vdiff_out', 'grad_out', 'th_out', 'eos_inplace', 'out_pms',
        'clean', 'ltau_rep', 'ec_s', 'lim_zc', 'g_modes_ab'],
    'extra' : ['a_m10_z005', 'a_m10_z02', 'a_m15_z005', 'a_m15_z02', 'alpha_autoguess', 'alpha_fuller',
        'alpha_jcd11', 'alpha_m10', 'alpha_m15', 'b_strength', 'beta_jcd11', 'burn_contours',
        'claret18_shift', 'conv_block', 'convection_out', 'hse_out', 'cst_mix', 'd_conv', 'd_expo_2', 'debug',
        'deff_add', 'deff_fact', 'ds', 'dv_add', 'dv_fact', 'easy_krisw', 'ec_corr', 'ec_extrap',
        'ec_grid', 'ell_max', 'ell_min', 'errabop', 'eta_0', 'eta_grad', 'eta_r_rgb_1', 'eta_r_rgb_2',
        'eta_sun', 'eta_vdif', 'expo_mix', 'fast_rot', 'full_hr', 'goop', 'inj_a', 'inj_mu', 'l_expo',
        'l_expo_2', 'lim_zc_out', 'log_lum_max_jvh', 'log_lum_min_jvh', 'lsmooth_fact', 'lsmooth_min',
        'lt_compton_h', 'lt_compton_l', 'm_mix', 'mass_loss_eta', 'mass_loss_zeta', 'mass_tol_cz', 'mode_damping',
        'n_contours_burn', 'n_fake_ce', 'ndelta', 'nmore_shells', 'no_diff_post', 'no_tot_rad_post' ,'npt_central_cz', 'nu_h_add', 'nu_h_fact',
        'nu_v_add', 'nu_v_fact', 'omega_c', 'omega_s', 'omega_sun', 'opt_tau_turb', 'osm_p_pertw',
        'ovi_type', 'ovs_type', 'p_power_matt15', 'pertw_profile', 'plato_track', 'r_14n15o', 'r_h2h',
        'r_transition', 'ratio_work_max', 're_c', 'rho_ext', 'ri_c', 'rota_out', 'source_cts', 'start_ovshti_zams',
        'struct2d', 'struct2d_file', 't0_2', 't0_corsaro', 't0_sun', 't_compton', 't_mix', 't_op_inf',
        't_op_sup', 't_oplib_inf', 't_oplib_sup', 'use_findne', 't_ov_inf', 'tau_conv_sun', 'teff_zams_2', 'teff_zams_sun', 'ts_average',
        'ts_average_length', 'ts_dt_max', 'ts_extrapolate', 'ts_hysteresis', 'ts_out', 'ts_patch_regime',
        'ts_qlim', 'ts_weak_b', 'v_inf_on_v_esc_cold', 'v_inf_on_v_esc_hot', 'width_transition',
        'writes_burn_zones', 'writes_conv', 'writes_2d', 'writes_l_nuc', 'writes_mdot', 'writes_nuc_neutr',
        'writes_tau_conv']}


    # Dictionaries for mapped lists --------------------------------------
    _dict_val = PARAMS_DICT

    initialized = Bool(False)

    # NL Cesam ----------------------------------------------------
    nom_chemin = String(os.environ['CESDIR'] + '/SUN_STAR_DATA/',
                        label='Path to data')

    nom_ctes = Trait('ctes_iau2015_codata2018', _dict_val['nom_ctes'], label='Physical constants')


    nom_output = Trait('Final model (non-adiabatic)', _dict_val['nom_output'],
                       label='Type of output file')
    nom_osc = Trait('None', _dict_val['nom_osc'],
                       label='Code to compute frequencies on-the-fly')
    trunc_nom_output = Trait('Final model (non-adiabatic)', _dict_val['trunc_nom_output'],
                       label='Type of output file')

    acor   = Bool(False, label='floats written to ACOR precision ?')
    out_2d = Bool(False, label='Save 2D outputs ?')
    out_h5 = Bool(False, label='Save HDF5 outputs ?')

    all_rep = Bool(False, label='Save all rep files')

    n_max = Int(2500, label='Max number of shells', tooltip="Number of layers in the internal structure")

    precision = Trait('CoRoT', _dict_val['precision'])

    # NL Evol ----------------------------------------------------
    agemax         = Float(4570.0, label='Maximum age (Myrs)')
    disclaim_evol  = String( "At present Cesam2k20 has been extensively tested and validated for low and intermediate mass stars up to the RGB phase.\n"
        "For the later phases and more massive stars, a lot of testing and development is still needed.",
                        label='DISCLAIMER')

    arret          = Trait('Other', _dict_val['arret'], label='Stop condition')

    dtlist         = Float(1.0e12, label='List time interval (Myrs)',
        tooltip='time interval between models written in .lis')
    log_teff       = Float(10.0, label='log Teff', tooltip="Stop condition: logTeff crosses log_teff. log_teff in solar unit.\n"
        "\tif log_teff < 0: stop when logTeff decreases below log_teff\n" \
        "\tif log_teff > 0: stop when logTeff increases above log_teff" )
    nb_max_modeles = Int(10000, label='Max number of mods')
    he_core        = Float(-1.0, label='Mass of He core')
    t_stop         = Float(5.0e10, label='Central temperature')
    x_stop         = Float(-1.0, label='Central X')
    xcx0           = Float(-1.0, label='Fraction of initial X in core.')
    rhoc           = Float(-1.0, label='Central density')
    l_stop         = Float(1.0e10, label='Luminosity crosses', tooltip="Stop condition: luminosity crosses l_stop. l_stop in solar unit (not log).\n"
        "\tif l_stop < 0: stop when luminosity decreases below l_stop\n" \
        "\tif l_stop > 0: stop when luminosity increases above l_stop")
    r_stop         = Float(1.0e10, label='Radius crosses', tooltip="Stop condition: total radius crosses r_stop. r_stop in solar unit.\n"
        "\tif r_stop < 0: stop when temperature decreases below r_stop\n" \
        "\tif r_stop > 0: stop when temperature increases above r_stop")


    # NL Mass ----------------------------------------------------
    disclaim_mass  = String( "Disclaimer: At present Cesam2k20 has been extensively tested and validated for low and intermediate mass stars.\n"
        "For more massive stars, a lot of testing and development is still needed.",
        label='DISCLAIMER')

    mtot = Float(1.0, label=u'Initial mass (M\u2299)')
    nom_pertm = Trait('Constant mass loss rate',
                        _dict_val['nom_pertm'], label='Mass loss prescription')
    nom_pertm_solar = Trait('Constant mass loss rate',
                        _dict_val['nom_pertm_solar'], label='Mass loss prescription for solar-type')
    nom_pertm_rgb = Trait('Constant mass loss rate',
                        _dict_val['nom_pertm_rgb'], label='Mass loss prescription for RGB')
    nom_pertm_agb = Trait('Constant mass loss rate',
                        _dict_val['nom_pertm_agb'], label='Mass loss prescription for AGB')
    nom_pertm_hot = Trait('Constant mass loss rate',
                        _dict_val['nom_pertm_hot'], label='Mass loss prescription for Teff > 10000 K')
    mdot = Float(0.0, label=u'Mass loss rate (M\u2299/yr)')


    # NL Chim ----------------------------------------------------
    grille_fixe = Bool(False, label='Fixed grid for diffusion')
    nom_abon = Trait('Solar (AGS 2009) + S10 recom.',
                     _dict_val['nom_abon'], label='Initial abundances')

    modif_chim = Bool(False, label='Takes into account custom files')
    garde_xish = Bool(False, label='Z/X from mixture')
    x0 = Float(0.73115705882274107, label='Initial X')
    y0 = Float(0.25544294117725902, label='Initial Y')
    z0 = Float(0.0134, label='Initial Z')
    zsx0 = Float(0.018327115683702433, label='Initial Z/X')
    z_or_zsx = Trait('Initial Z', 'Initial Z/X')
    x_or_y = Trait('Initial Y', 'Initial X')


    # NL Conv ----------------------------------------------------
    nom_conv        = Trait(u'MLT with l=\u03B1 Hp', _dict_val['nom_conv'], label='Convection prescription')
    nom_alpha       = Trait(u'Fixed alpha', _dict_val['nom_alpha'], label='Alpha prescription')
    alpha           = Float(1.6420906995607241, label=u'\u03B1 (convection)')
    nom_ovshts      = Trait(u'Fixed core overshoot', _dict_val['nom_ovshts'], label='Core overshoot prescription')
    ovshts          = Float(0.0, label='Ovsht. from core', tooltip='If step overshoot, use value in Hp. If diffusive overshoot, use f0 value (~Hp/11.4)')
    ovshti          = Float(0.0, label='Ovsht. from envelope', tooltip='If step overshoot, use value in Hp. If diffusive overshoot, use f0 value (~Hp/11.4)')
    jpz             = Bool(False, label='JPZ prescription', tooltip='Overshoot param. according to Zahn 1991. Should only be used with step overshoot.')
    cpturb          = Float(0.0, label='Constant for Pturb')
    ledoux          = Bool(False, label='Ledoux criterion')
    solberg_hoiland = Bool(False, label="Solberg-Hoiland criterion")
    gough_tayler    = Bool(False, label="Gough-Tayler criterion")


    # NL Diff ----------------------------------------------------

    diffusion = Bool(False)
    nom_diffm = Trait('Michaud and Proffit formalism' ,
                      _dict_val['nom_diffm'], label='Microscopic diffusion')
    nom_difft = Trait('Constant turbulent coeff.',
                      _dict_val['nom_difft'], label='Turbulent diffusion')

    nom_difft_fing = Trait('None',
                      _dict_val['nom_difft_fing'], label='Thermohaline convection')

    nom_difft_smc = Trait('None',
                      _dict_val['nom_difft_smc'], label='Semi-convection')

    d_turb = Float(0.0, label='Isotropic turbulent diffusion coefficient')
    re_nu = Float(0.0, label='Radiative diffusivity coefficient')
    nom_frad = Trait('None', _dict_val['nom_frad'],
                     label='Radiative accelerations')


    # NL Rot ----------------------------------------------------
    w_rot = Float(0.0, label=u'Initial \u03a9')
    unit = Trait('kms/s', _dict_val['unit'])
    nom_diffw = Trait('None', _dict_val['nom_diffw'],
                      label='Turbulent viscosity coefficients')
    nom_diffmg_ts = Trait('None', _dict_val['nom_diffmg_ts'],
                       label='Magnetic fields')
    nom_thw = Trait('None', _dict_val['nom_thw'],
                    label='Transport of angular momentum')
    rot_zc = Trait('Solid body', _dict_val['rot_zc'],
                   label='Rotation in CZ')
    nom_pertw = Trait('None', _dict_val['nom_pertw'], label='Prescription')
    p_pertw = Float(6.5e47, label='Constant K')
    lim_jpz = Bool(True)
    init_rot = Bool(False,label='Use new rotation profile')
    mixed_modes = Bool(False)
    g_modes = Bool(False)
    tau_disk = Float(0.0, label='Disk lifetime (Myr)')
    p_disk = Float(0.0, label='Disk period (days)')
    w_sat = Float(8.0, label=u'\u03a9 sat ( \u03a9 \u2299)')
    nom_gsf = Trait('None', _dict_val['nom_gsf'],
            label="Prescription for GSF instability")

   # NL Etat -------------------------------------------------

    nom_etat = Trait(u'OPAL 2005, Z\u22600', _dict_val['nom_etat'],
                     label='Equation of state')

    f_eos = List(File, ['eos5_sunags09_nodiff.bin', '', '', '', '', '', '', ''],
                 label='EoS tables', comparison_mode=2)


    # NL Opa -------------------------------------------------
    nom_opa = Trait('OPAL and AF, interp. of Y. Lebreton',
                    _dict_val['nom_opa'], label='Opacity')
    nom_opa_cond = Trait("Pothekin 1999 (table of 2021)",
                    _dict_val['nom_opa_cond'], label='Conductive opacity')

    f_opa = List(File, ['opa_yveline_AGSS09_met.bin', 'TOPS_MW',
                     'TOPS_LMC', 'TOPS_SMC', '', '', '', ''],
                 label='Opacity tables', comparison_mode=2)


    # NL Nuc -------------------------------------------------
    nom_nuc = Trait('PP+CNO, 13 elems + Ne, Na, Mg, Al, Si, P, S, Ca, Fe, Ni in eq., 1<T6<80',
                    _dict_val['nom_nuc'], label='Nuclear reaction rates')

    nom_nuc_cpl = Trait('NACRE + LUNA', _dict_val['nom_nuc_cpl'],
                        label='Compilation of reaction rates')
    nom_neutrinos = Trait('Haft et al. (1994) + Weigert (1966)' ,
                        _dict_val['nom_neutrinos'], label='Neutrino losses')

    mitler = Bool(False, label='Screening according to Mitler (1997)')


    # NL Atm -------------------------------------------------
    nom_atm = Trait(u'Restitution of the atmosphere from a T(\u03C4) law',
                    _dict_val['nom_atm'], label='External limit')
    nom_tdetau = Trait(u'Fully radiative T(\u03C4), Hopf',
                       _dict_val['nom_tdetau'],
                       label=u'T(\u03C4) prescription')
    tau_max = Range(0.0, 100.0, 20.0, label=u'Maximum \u03C4')
    lim_ro = Bool(True, label='external condition on density')

    # NL Mag -------------------------------------------------
    nom_mag   = Trait('No magnetic field', _dict_val['nom_mag'],
        label='Prescription for the magnetic field')
    mag_inhib = Float(0.0, label=u'Magnetic inhibition parameter \u03b4')
    br        = Float(0.0, label=u'Radial component Br')
    bphi      = Float(0.0, label=u'Azimuthal component Bphi')
    btheta    = Float(0.0, label=u'Latitudinal component Btheta')
    xmax    = Float(0.0, label=u'Depth-like quantity above which \u03b4 is constant')
    xmax_unit   = Trait('Fraction of the star radius', _dict_val['xmax_unit'], label='Unit of Rmax')
    brmax    = Float(0.0, label=u'Ceiling value of br')

    # NL regl -------------------------------------------------
    _dict_val_num = NUM_DICT

    disclaim_num  = String( "DISCLAIMER: The documentation of each option is accessible by simply hovering the name of each variable.\n"
        "Changing the first widget 'precision' will revert the settings to preset values. You can then change the parameters individually.")

    ini0           = Int( DEFAULT_NUM['ini0'],    label='ini0',       tooltip="Number of Newton-Raphson iterations with re-estimation of chemical composition, and limits RZ/CZ" )
    l0             = Int( DEFAULT_NUM['l0'],    label='l0',         tooltip="[Deprecated] In output ASCII files, number of points added in the neighbourhood of limits ZR/CZ" )
    m_qs           = Int( DEFAULT_NUM['m_qs'],    label='m_qs',       tooltip="order of the B-Splines for interpolation of quasi-static variables" )
    m_ch           = Int( DEFAULT_NUM['m_ch'],    label='m_ch',       tooltip="order of the B-Splines for interpolation of the chemical composition" )
    m_vth          = Int( DEFAULT_NUM['m_vth'],    label='m_vth',      tooltip="order of the B-Splines for interpolation of thermodynamic variblaes (only used in `rep.osc_v3.f90`)" )
    m_rot          = Int( DEFAULT_NUM['m_rot'],    label='m_rot',      tooltip="order of the B-Splines for interpolation of rotation quantities" )
    m_tds          = Int( DEFAULT_NUM['m_tds'],    label='m_tds',      tooltip="order of the B-Splines for interpolation of T dS" )
    m_ptm          = Int( DEFAULT_NUM['m_ptm'],    label='m_ptm',      tooltip="ordre of the B-Splines for interpolation of mass (used when there is mass loss." )
    n_atm          = Int( DEFAULT_NUM['n_atm'],   label='n_atm',      tooltip="Number of layers in the atmosphere" )
    n_min          = Int( DEFAULT_NUM['n_min'],  label='n_min',      tooltip="Minimum number of layers" )
    n_min_post     = Int( DEFAULT_NUM['n_min_post'],  label='n_min_post',      tooltip="Minimum number of layers after MS" )
    ordre          = Int( DEFAULT_NUM['ordre'],    label='ordre',      tooltip="Order of the scheme used for the interpolation of nuclear reaction rates with `rk_imps`" )
    yld_rep        = Int( DEFAULT_NUM['yld_rep'],    label='yld_rep',    tooltip="write one `.rep` file every `yld_rep` (needs `nbmax_modeles<0`, `all_rep=true`)" )
    osc2d_step     = Int( DEFAULT_NUM['osc2d_step'],    label='osc2d_step', tooltip="write one `.osc2d` file every `osc2d_step`" )
    osc_step       = Int( DEFAULT_NUM['osc_step'],    label='osc_step',   tooltip="write one `.osc` file every `osc_step` (redundant with `yld_osc`)" )
    # For 2D computations:
    m_qs2d         = Int( DEFAULT_NUM['m_qs2d'], label='m_qs2d', tooltip="order of the spline for 2D interpolation of 1D quantities on 2D grid" )
    m_rot2d        = Int( DEFAULT_NUM['m_rot2d'], label='m_rot2d', tooltip="order of the spline for 2D interpolation of rotation" )
    m_legendre     = Int( DEFAULT_NUM['m_legendre'], label='m_legendre', tooltip="Order of the decomposition on Legendre polynomials" )
    n_theta        = Int( DEFAULT_NUM['n_theta'], label='n_theta', tooltip="Number of points in the angular mesh" )

    iter_qs        = List( Int, DEFAULT_NUM['iter_qs'], label='iter_qs', tooltip="Controls quasi-static variables" )

    ctel           = Float( DEFAULT_NUM['ctel'], label='ctel', tooltip="Repartition constant for the luminosity" )
    ctep           = Float( DEFAULT_NUM['ctep'], label='ctep', tooltip="Repartition constant for the pressure" )
    ctem           = Float( DEFAULT_NUM['ctem'], label='ctem', tooltip="Repartition constant for the mass" )
    cter           = Float( DEFAULT_NUM['cter'], label='cter', tooltip="Repartition constant for the radius" )
    ctet           = Float( DEFAULT_NUM['ctet'], label='ctet', tooltip="Repartition constant for the temperature" )
    coeff_d_grav_post = Float( DEFAULT_NUM['coeff_d_grav_post'], label='d_grav_post', tooltip="Coefficient of maximum variation of TdS after MS" )
    d_grav         = Float( DEFAULT_NUM['d_grav'], label='d_grav', tooltip="Maximum variation of TdS" )
    d_grav_int     = Float( DEFAULT_NUM['d_grav_int'], label='d_grav_int', tooltip="Maximum variation of TdS inside `static`" )
    dlntc          = Float( DEFAULT_NUM['dlntc'], label='dlntc', tooltip="Maximum variation of core temperature between two consecutive time steps." )
    drhoc          = Float( DEFAULT_NUM['drhoc'], label='drhoc', tooltip="Maximum relative variation of cenral density between two consecutive time steps" )
    dteff          = Float( DEFAULT_NUM['dteff'], label='dteff', tooltip="Maximum relative variation of Teff between two consecutive time steps" )
    dlogg          = Float( DEFAULT_NUM['dlogg'], label='dlogg', tooltip="Maximum dex variation of log g between two consecutive time steps" )
    dalpha         = Float( DEFAULT_NUM['dalpha'], label='dalpha', tooltip=u"Maximum relative variation of \u03b1 between two consecutive time steps" )
    dlum           = Float( DEFAULT_NUM['dlum'], label='dlum', tooltip="Maximum dex variation of luminosity between two consecutive time steps" )
    dw             = Float( DEFAULT_NUM['dw'], label='dw', tooltip=u"Maximum relative variation of central \u03a9 between two consecutive time steps" )
    drhoca         = Float( DEFAULT_NUM['drhoca'], label='drhoca', tooltip=u"Maximum variation of central \u03c1 to average \u03c1 between two consecutive time steps" )
    dsenv          = Float( DEFAULT_NUM['dsenv'], label='dsenv', tooltip="Maximum variation of adiabatic entropy between two consecutive time steps" )
    evolved        = Float( DEFAULT_NUM['evolved'], label='evolved', tooltip=u"Fraction of the initial X content below which variations of Teff, luminosity, log g and \u03c1_c are limited" )
    dn_fixe        = Float( DEFAULT_NUM['dn_fixe'], label='dn_fixe', tooltip="Modification rate for a fixed grid of chemical composition" )
    dpsi           = Float( DEFAULT_NUM['dpsi'], label='dpsi', tooltip=u"Maximum variation of \u03c8 between two consecutive time steps before `n_qs` is modify" )
    dt0            = Float( DEFAULT_NUM['dt0'], label='dt0', tooltip="Initial time step when starting from ZAMS" )
    dtmax          = Float( DEFAULT_NUM['dtmax'], label='dtmax', tooltip="Maximum time step" )
    dtmin          = Float( DEFAULT_NUM['dtmin'], label='dtmin', tooltip="Minimum time step" )
    fmin_abon      = Float( DEFAULT_NUM['fmin_abon'], label='fmin_abon', tooltip="Factor linking `ab_min` et `abon_ini`" )
    loc_zc         = Float( DEFAULT_NUM['loc_zc'], label='loc_zc', tooltip="Precision of locating RC/CZ boundaries" )
    precix         = Float( DEFAULT_NUM['precix'], label='precix', tooltip="Precision on Newton-Raphson iterations for spatial integrations." )
    precit         = Float( DEFAULT_NUM['precit'], label='precit', tooltip="Maximum relative variation for temporal integration of chemical composition." )
    psi0           = Float( DEFAULT_NUM['psi0'], label='psi0', tooltip="Distribution constant to ensure." )
    ro_test        = Float( DEFAULT_NUM['ro_test'], label='ro_test', tooltip="Test for variation of gravitational energy if `rho > ro_test`." )
    q0             = Float( DEFAULT_NUM['q0'], label='q0', tooltip="In output ASCII files, a point is placed at q0>0 times the spacing between the first two points." )
    yld_dx         = Float( DEFAULT_NUM['yld_dx'], label='yld_dx', tooltip="Maximum relative variation of X in the core between two consecutive time steps" )
    yld_rac        = Float( DEFAULT_NUM['yld_rac'], label='yld_rac', tooltip="Optical depth of the matching between interior and atmosphere" )
    tau_min        = Float( DEFAULT_NUM['tau_min'], label='tau_min', tooltip="Minimum optical depth" )
    dt_m           = Float( DEFAULT_NUM['dt_m'], label='dt_m', tooltip="Parameter needed to modulate the value of `dtmax` with the value of the mass: `dtmax = dtmax / mtot**(dt_m)`" )
    x_tams         = Float( DEFAULT_NUM['x_tams'], label='x_tams', tooltip="Fraction of hydrogen at center that defines TAMS" )
    y_agb          = Float( DEFAULT_NUM['y_agb'], label='y_agb', tooltip="Fraction of helium at center that defines AGB" )
    dx_tams        = Float( DEFAULT_NUM['dx_tams'], label='dx_tams', tooltip="precision on the fraction of hydrogen at center that defines TAMS" )
    mcz_ext_min    = Float( DEFAULT_NUM['mcz_ext_min'], label='mcz_ext_min', tooltip="external CZ must have at least this fraction of mstar" )

    exact_stop     = Bool( DEFAULT_NUM['exact_stop'], label='exact_stop', tooltip="If true, a last model is computed to match exactly the stop criteria (e.g. to reach given core temperature)" )
    en_masse       = Bool( DEFAULT_NUM['en_masse'], label='en_masse', tooltip="If true, use lagrangian coordinates" )
    kipp           = Bool( DEFAULT_NUM['kipp'], label='kipp', tooltip="If true, use Kippenhahn's approximation TdS = dE + PdV" )
    lisse          = Bool( DEFAULT_NUM['lisse'], label='lisse', tooltip="If true, chemical composition and rotation quantities are smoothed" )
    mu_saha        = Bool( DEFAULT_NUM['mu_saha'], label='mu_saha', tooltip="If true, we use Saha equation to compute mean molecular weight \u03bc" )
    mvt_dis        = Bool( DEFAULT_NUM['mvt_dis'], label='mvt_dis', tooltip="If true, adjustment of chemical composition due to discontinuity displacements." )
    new_bv         = Bool( DEFAULT_NUM['new_bv'], label='new_bv', tooltip=u"If true [not recommended], calculate the Brunt-Vaissala freq using \u03c6=d ln \u03c1 / ln \u03bc" )
    yld_hr         = Bool( DEFAULT_NUM['yld_hr'], label='yld_hr', tooltip="If true, writes a special file `.HRYL`" )
    yld_l23        = Bool( DEFAULT_NUM['yld_l23'], label='yld_l23', tooltip="If true, the unknown is L^{2/3} instead of L" )
    yld_fgong      = Bool( DEFAULT_NUM['yld_fgong'], label='yld_fgong', tooltip="If true, fgong files are generated" )
    # For 2D computations:
    d2             = Bool( DEFAULT_NUM['d2'], label='d2', tooltip="If true, we compute the 2D structure of the star according to Roxburgh 2006 method" )
    general        = Bool( DEFAULT_NUM['general'], label='general', tooltip="If True, general structure equation are used" )
    simple         = Bool( DEFAULT_NUM['simple'], label='simple', tooltip="If true, `Omega_2` is not accounted for" )
    extended_osc2d = Bool( DEFAULT_NUM['extended_osc2d'], label='extended_osc2d', tooltip="If true, more quantities are added to osc2d files" )
    # For diffusion:
    vdiff_out      = Bool( DEFAULT_NUM['vdiff_out'], label='vdiff_out', tooltip="If true, write diffusion velocities to an ASCII file" )
    grad_out       = Bool( DEFAULT_NUM['grad_out'], label='grad_out', tooltip="If true, write radiative accelerations to an ASCII file" )
    th_out         = Bool( DEFAULT_NUM['th_out'], label='th_out', tooltip="If true, write information about thermohaline convection to an ASCII file" )
    eos_inplace    = Bool( DEFAULT_NUM['eos_inplace'], label='eos_inplace', tooltip="If true, EoS table must be stored in the directory when the model is computed" )
    # Other
    out_pms        = Bool( DEFAULT_NUM['out_pms'], label='out_pms', tooltip="Output .osc on PMS" )
    clean          = Bool( DEFAULT_NUM['clean'], label='clean', tooltip="If true, remove dublicate or non stricly increasing layer before" )
    g_modes_ab     = Bool( DEFAULT_NUM['g_modes_ab'], label='g_modes_ab', tooltip="If True, uses Adams-Bashforth scheme to couple last 3 times step when interpolating g-mode flux" )

    ltau_rep       = Trait( DEFAULT_NUM['ltau_rep'], _dict_val_num['ltau_rep'], label="ltau_rep", tooltip=u"log \u03c4 follows either a linear trend or a piecewise cubic trend" )
    ec_s           = Trait( DEFAULT_NUM['ec_s'], _dict_val_num['ec_s'], label="ec_s", tooltip="Method of determination of adiabatic entropy" )
    lim_zc         = Trait( DEFAULT_NUM['lim_zc'], _dict_val_num['lim_zc'], label="lim_zc", tooltip="Switch between original `lim_zc`, and `lim_zc` that follows Gabriel (2014)'s recommendations" )


    # NL extra -------------------------------------------------
    # extra = Code()
    _dict_val_extra = EXTRA_DICT

    disclaim_extra  = String( "DISCLAIMER: The documentation of each option is accessible by simply hovering the name of each variable." )

    # Atmosphere
    rho_ext             = Float( DEFAULT_EXTRA['rho_ext'], label="rho_ext", tooltip="External density" )
    easy_krisw          = Bool( DEFAULT_EXTRA['easy_krisw'], label='easy_krisw', tooltip="if true, and nom_atm='krisw', first model is computed with eddington" )

    # OSM renormalization factors
    disclaim_osm        = String( "DISCLAIMER: Avoid ill-conditionned hessian matrix by renormalizaing a given parameter" )
    osm_p_pertw         = Float( DEFAULT_EXTRA['osm_p_pertw'], label="osm_p_pertw", tooltip="Renormalize the constant p_pertw (should be around 1e31 for Matt15)" )


    # Rotation
    npt_central_cz      = Int( DEFAULT_EXTRA['npt_central_cz'], label='npt_central_cz', tooltip='Number of points in the core CZ (used when rotation is taken into account).'
        '`npt_central_cz = 0`: radiative central boundary condition.' )
    nmore_shells        = Int( DEFAULT_EXTRA['nmore_shells'], label='nmore_shells', tooltip='Factor that increase the number of layers in the radiative zone (used when rotation is taken into account).'
        '`nmore_shells = 1` no extra points' )
    nu_v_fact           = Float( DEFAULT_EXTRA['nu_v_fact'], label='nu_v_fact', tooltip="Multiplicative factor for $\nu_v$" )
    nu_v_add            = Float( DEFAULT_EXTRA['nu_v_add'], label='nu_v_add', tooltip="Offset for $\nu_v$" )
    nu_h_fact           = Float( DEFAULT_EXTRA['nu_h_fact'], label='nu_h_fact', tooltip="Multiplicative factor for $\nu_h$" )
    nu_h_add            = Float( DEFAULT_EXTRA['nu_h_add'], label='nu_h_add', tooltip="Offset for $\nu_h$" )
    deff_fact           = Float( DEFAULT_EXTRA['deff_fact'], label='deff_fact', tooltip="Multiplicative factor for deff" )
    deff_add            = Float( DEFAULT_EXTRA['deff_add'], label='deff_add', tooltip="Offset for deff" )
    dv_fact             = Float( DEFAULT_EXTRA['dv_fact'], label='dv_fact', tooltip="Multiplicative factor pour dv" )
    dv_add              = Float( DEFAULT_EXTRA['dv_add'], label='dv_add', tooltip="Offset pour dv" )
    re_c                = Float( DEFAULT_EXTRA['re_c'], label='re_c', tooltip="Critical Ryenolds number" )
    ri_c                = Float( DEFAULT_EXTRA['ri_c'], label='ri_c', tooltip="Critical Richardson number" )
    omega_sun           = Float( DEFAULT_EXTRA['omega_sun'], label='omega_sun', tooltip="Surface angular velocity of the sun" )
    fast_rot            = Float( DEFAULT_EXTRA['fast_rot'], label='fast_rot', tooltip="if `fast_rot` > 0.0, the first `fast_rot` Myrs are "
        "computed with a slow rotation rate. The rotation is increased linearly with age until it reaches the "
        "desired initial rotation rate at an age of `fast_rot`. This avoids the case of breakup velocity reached "
        "at the beginning." )

    opt_tau_turb        = Int( DEFAULT_EXTRA['opt_tau_turb'], label='opt_tau_turb', tooltip='Caracteristic time of turbulence.' )

    # Mixed modes
    ell_max             = Int( DEFAULT_EXTRA['ell_max'], label='ell_max', tooltip='degree $\ell_{max}$ of mixed modes' )
    ell_min             = Int( DEFAULT_EXTRA['ell_min'], label='ell_min', tooltip='degree $\ell_{min}$ of mixed modes' )
    ndelta              = Int( DEFAULT_EXTRA['ndelta'], label='ndelta', tooltip='Number of large speration to the left or to the right of $\nu_{max}$.' )
    eta_sun             = Float( DEFAULT_EXTRA['eta_sun'], label='eta_sun', tooltip="Damping rate for the Sun at nu_max" )
    t0_corsaro          = Float( DEFAULT_EXTRA['t0_corsaro'], label='t0_corsaro', tooltip="factor T0 Corsaro et al. 2012, eq (13)" )
    eta_0               = Float( DEFAULT_EXTRA['eta_0'], label='eta_0', tooltip="eta_0 for the mixed modes" )
    ratio_work_max      = Float( DEFAULT_EXTRA['ratio_work_max'], label='ratio_work_max', tooltip="To limit ell_max" )
    mode_damping        = Bool( True, label='mode_damping', tooltip="includes damping of mixed modes" )

    # Convection
    n_fake_ce           = Int( DEFAULT_EXTRA['n_fake_ce'], label='n_fake_ce', tooltip='Minimum number of mesh point to be convective in the envelope (important for numerical reasons)' )
    d_conv              = Float( DEFAULT_EXTRA['d_conv'], label='d_conv', tooltip="Turbulent diffusion coefficient in convective zone" )
    mass_tol_cz         = Float( DEFAULT_EXTRA['mass_tol_cz'], label='mass_tol_cz', tooltip="determination of the limits of convection zones (10percent by default)" )
    conv_block          = Bool( False, label='conv_block', tooltip="if true, convection is blocked because of the surface magnetic field (for Ap stars)" )
    no_tot_rad_post     = Bool( False, label='no_tot_rad_post', tooltip="if true, after tams tot_rad model get a fake convective envelope with n_fake_ce points (important for numerical reasons)" )

    # Setup for ad hoc rotation profile:
    r_transition        = Float( DEFAULT_EXTRA['r_transition'], label='r_transition', tooltip="transition radius (r/Rstar) between core and envelope" )
    width_transition    = Float( DEFAULT_EXTRA['width_transition'], label='width_transition', tooltip="width of transition region between core and envelope" )
    omega_c             = Float( DEFAULT_EXTRA['omega_c'], label='omega_c', tooltip="(1 muHz) central angular velocity" )
    omega_s             = Float( DEFAULT_EXTRA['omega_s'], label='omega_s', tooltip="(55 nHz) surface angular velocity" )

    # Structure model to be deformed
    struct2d            = Trait( DEFAULT_EXTRA['struct2d'], _dict_val_extra['struct2d'], label='struct2d', tooltip="variable used to deform a model from an old structure. Values can be the following: '', 'cles', 'cesam2k20'" )
    struct2d_file       = File( DEFAULT_EXTRA['struct2d_file'], label='struct2d_file', tooltip="Name of the file where structure is written" )

    # For turbulent diffusion:
    t_mix               = Float( DEFAULT_EXTRA['t_mix'], label='t_mix', tooltip="Reference temperature for the mixing in the Montreal code (difft_montreal.f90)" )
    m_mix               = Float( DEFAULT_EXTRA['m_mix'], label='m_mix', tooltip="Reference mass for the mixing in the Montreal code (difft_montreal.f90)" )
    expo_mix            = Float( DEFAULT_EXTRA['expo_mix'], label='expo_mix', tooltip="Density exponent for the mixing in the Montreal code (difft_montreal.f90)" )
    cst_mix             = Float( DEFAULT_EXTRA['cst_mix'], label='cst_mix', tooltip="Numerical stability constant for the mixing in the Montreal code (difft_xxx.f90)" )

    d_expo_2            = Float( DEFAULT_EXTRA['d_expo_2'], label='d_expo_2', tooltip="Amplitude of the second exponential (difft_expo.f90; typical value: 200.0)" )
    l_expo              = Float( DEFAULT_EXTRA['l_expo'], label='l_expo', tooltip="Width of the first exponential, in percent of the radius (difft_expo.f90; typical value: 0.01)" )
    l_expo_2            = Float( DEFAULT_EXTRA['l_expo_2'], label='l_expo_2', tooltip="Width of the second exponential, in percent of the radius (difft_expo.f90; typical value: 0.23)" )

    # For radiative diffusion:
    eta_grad            = Float( DEFAULT_EXTRA['eta_grad'], label='eta_grad', tooltip="Allows to in/decrease the efficiency of g_rad, in order to account for uncertainties (~30%, 0.70<eta_gra<1.3 )" )
    goop                = Float( DEFAULT_EXTRA['goop'], label='goop', tooltip="Age from which g_rad are taken into account" )

    # For atomic diffusion:
    eta_vdif            = Float( DEFAULT_EXTRA['eta_vdif'], label='eta_vdif', tooltip="Allows to in/decrease the efficiency of atomic diffusion, in order to account for uncertainties (~20%, 0.80<eta_vdiff<1.2 )" )
    no_diff_post        = Bool( DEFAULT_EXTRA['no_diff_post'], label='no_diff_post', tooltip='Allows to cut atomic diffusion after MS' )
    b_strength          = Float( DEFAULT_EXTRA['b_strength'], label='b_strength', tooltip="Magnetic field" )

    errabop             = Float( DEFAULT_EXTRA['errabop'], label='errabop', tooltip="Maximum abundance difference to recompute opacities" )
    t_op_sup            = Float( DEFAULT_EXTRA['t_op_sup'], label='t_op_sup', tooltip="Temperature range in which OP monochromatic opacities are use to compute the Rosseland mean" )
    t_op_inf            = Float( DEFAULT_EXTRA['t_op_inf'], label='t_op_inf', tooltip="Temperature range in which OP monochromatic opacities are use to compute the Rosseland mean" )

    t_oplib_sup         = Float( DEFAULT_EXTRA['t_oplib_sup'], label='t_oplib_sup', tooltip="Temperature range in which OPLIB monochromatic opacities are use to compute the Rosseland mean" )
    t_oplib_inf         = Float( DEFAULT_EXTRA['t_oplib_inf'], label='t_oplib_inf', tooltip="Temperature range in which OPLIB monochromatic opacities are use to compute the Rosseland mean" )
    use_findne          = Bool( DEFAULT_EXTRA['use_findne'], label='use_findne', tooltip="Use findne subroutine in which OPLIB monochromatic opacities are use to compute the Rosseland mean")

    teff_zams_sun       = Float( DEFAULT_EXTRA['teff_zams_sun'], label='teff_zams_sun', tooltip="difft_grille.f90" )
    teff_zams_2         = Float( DEFAULT_EXTRA['teff_zams_2'], label='teff_zams_2', tooltip="difft_grille.f90" )
    t0_sun              = Float( DEFAULT_EXTRA['t0_sun'], label='t0_sun', tooltip="difft_grille.f90" )
    t0_2                = Float( DEFAULT_EXTRA['t0_2'], label='t0_2', tooltip="difft_grille.f90" )

    # Overshoot and overshooting regions
    t_ov_inf            = Float( DEFAULT_EXTRA['t_ov_inf'], label='t_ov_inf', tooltip="overshoot inférieur si T > t_ov_inf à la limite de la zone convective" )
    alpha_jcd11         = Float( DEFAULT_EXTRA['alpha_jcd11'], label='alpha_jcd11', tooltip="Controls jcd11 temperature gradient" )
    beta_jcd11          = Float( DEFAULT_EXTRA['beta_jcd11'], label='beta_jcd11', tooltip="Controls jcd11 temperature gradient" )
    claret18_shift      = Float( DEFAULT_EXTRA['claret18_shift'], label='claret18_shift', tooltip="Shift Eq. 2 of Claret et al. 2018 of a constante value" )
    ovs_type            = Trait( DEFAULT_EXTRA['ovs_type'], _dict_val_extra['ovs_type'], label='ovs_type', tooltip="gradient de température dans les zones d'overshoot (adia, radia, jcd11= expression hybride - JCD, Monteiro et al. 2011, MNRAS.414.1158C)" )
    ovi_type            = Trait( DEFAULT_EXTRA['ovi_type'], _dict_val_extra['ovi_type'], label='ovi_type', tooltip="gradient de température dans les zones d'overshoot (adia, radia, jcd11= expression hybride - JCD, Monteiro et al. 2011, MNRAS.414.1158C)" )
    start_ovshti_zams   = Bool( DEFAULT_EXTRA['start_ovshti_zams'], label='start_ovshti_zams', tooltip='include overshoot below the enveloppe only after the ZAMS' )

    # Compton opacities
    lt_compton_l        = Float( DEFAULT_EXTRA['lt_compton_l'], label='lt_compton_l', tooltip="Compton opacity is calculated when logT >= lt_compton_l" )
    lt_compton_h        = Float( DEFAULT_EXTRA['lt_compton_h'], label='lt_compton_h', tooltip="Above lt_compton_h opacity is equal to Compton opacity" )
    t_compton           = Float( DEFAULT_EXTRA['t_compton'], label='t_compton', tooltip="If T > t_compton chemical elements are fully ionized and we use Compton opacity" )

    # Entropy calibration
    ds                  = Float( DEFAULT_EXTRA['ds'], label='ds', tooltip="entropy offset needed for solar calibration" )
    ec_grid             = Trait( DEFAULT_EXTRA['ec_grid'], _dict_val_extra['ec_grid'], label='ec_grid', tooltip="Coefficient of the entropy prescriptions" )
    ec_extrap           = Trait( DEFAULT_EXTRA['ec_extrap'], _dict_val_extra['ec_extrap'], label='ec_extrap', tooltip="What to do when we are outside of the grid?" )
    alpha_autoguess     = Bool( DEFAULT_EXTRA['alpha_autoguess'], label='alpha_autoguess', tooltip="If True, a guess is found automatically for alpha" )
    ec_corr             = Bool( DEFAULT_EXTRA['ec_corr'], label='ec_corr', tooltip="If False, entropies are not corrected for offset and mean molecular weight variations" )
    alpha_m10           = Float( DEFAULT_EXTRA['alpha_m10'], label='alpha_m10', tooltip="(use for 1st guess) alpha of 1.0Msun model with your entropy prescription" )
    alpha_m15           = Float( DEFAULT_EXTRA['alpha_m15'], label='alpha_m15', tooltip="(use for 1st guess) alpha of 1.5Msun model with your entropy prescription" )
    a_m10_z005          = Float( DEFAULT_EXTRA['a_m10_z005'], label='a_m10_z005', tooltip="(use for 1st guess) alpha of 1.0Msun and mu = 0.59 model with your entropy prescription" )
    a_m15_z005          = Float( DEFAULT_EXTRA['a_m15_z005'], label='a_m15_z005', tooltip="(use for 1st guess) alpha of 1.5Msun and mu = 0.59 model with your entropy prescription" )
    a_m10_z02           = Float( DEFAULT_EXTRA['a_m10_z02'], label='a_m10_z02', tooltip="(use for 1st guess) alpha of 1.0Msun and mu = 0.61 model with your entropy prescription" )
    a_m15_z02           = Float( DEFAULT_EXTRA['a_m15_z02'], label='a_m15_z02', tooltip="(use for 1st guess) alpha of 1.5Msun and mu = 0.61 model with your entropy prescription" )

    # Mass loss
    mass_loss_eta       = Float( DEFAULT_EXTRA['mass_loss_eta'], label='mass_loss_eta', tooltip="eta from Reimers law. Between 0.5 and 1.0" )
    mass_loss_zeta      = Float( DEFAULT_EXTRA['mass_loss_zeta'], label='mass_loss_zeta', tooltip="zeta from Reimers law: metalicity scaling, in [0.5 , 0.8]" )
    log_lum_min_jvh     = Float( DEFAULT_EXTRA['log_lum_min_jvh'], label='log_lum_min_jvh', tooltip="minimal luminosity consistent with the model (in log)" )
    log_lum_max_jvh     = Float( DEFAULT_EXTRA['log_lum_max_jvh'], label='log_lum_max_jvh', tooltip="maximal luminosity consistent with the model (in log)" )
    v_inf_on_v_esc_hot  = Float( DEFAULT_EXTRA['v_inf_on_v_esc_hot'], label='v_inf_on_v_esc_hot', tooltip="Value of V_inf/V_esc on the hot side of bi-stability jump" )
    v_inf_on_v_esc_cold = Float( DEFAULT_EXTRA['v_inf_on_v_esc_cold'], label='v_inf_on_v_esc_cold', tooltip="Value of V_inf/V_esc on the cool side of bi-stability jump" )
    eta_r_rgb_1         = Float( DEFAULT_EXTRA['eta_r_rgb_1'], label='eta_r_rgb_1', tooltip="Value of eta_r for mass value of 1M_sol" )
    eta_r_rgb_2         = Float( DEFAULT_EXTRA['eta_r_rgb_2'], label='eta_r_rgb_2', tooltip="Value of eta_r for mass value of [3,4,5,7]M_sol" )

    # Angular momentum loss
    p_power_matt15      = Float( DEFAULT_EXTRA['p_power_matt15'], label='p_power_matt15', tooltip="power in Matt 15" )
    tau_conv_sun        = Float( DEFAULT_EXTRA['tau_conv_sun'], label='tau_conv_sun', tooltip="Solar convective turnover time (in days)" )
    pertw_profile       = Trait( DEFAULT_EXTRA['pertw_profile'], _dict_val_extra['pertw_profile'], label='pertw_profile', tooltip="How we deal with angular momentum loss at low tau_conv" )

    # Nuclear
    screen              = Bool( DEFAULT_EXTRA['screen'], label='screen', tooltip="If False, totaly removes the screening in nuclear reactions" )
    eta_nuc             = Float( DEFAULT_EXTRA['eta_nuc'], label='eta_nuc', tooltip="Efficiency of pp and cno nuclear reactions" )

    # Modification of individual nuclear reaction rates
    r_h2h               = Float( DEFAULT_EXTRA['r_h2h'], label='r_h2h', tooltip="change in the rate of p + p -> 2H" )
    r_14n15o            = Float( DEFAULT_EXTRA['r_14n15o'], label='r_14n15o', tooltip="change in the rate of N14 + p -> O15" )

    # Custom energy injection below CZ:
    # injection distribution is a truncated sinus
    inj_a               = Float( DEFAULT_EXTRA['inj_a'], label="inj_a", tooltip="Total quantity of energy injected [erg/s]" )
    inj_mu              = Float( DEFAULT_EXTRA['inj_mu'], label="inj_mu", tooltip="location extra of the sin: mu * m_zc (bottom env) [%]" )

    # Tayler Spruit instability
    ts_dt_max           = Float( DEFAULT_EXTRA['ts_dt_max'], label="ts_dt_max", tooltip="maximum time step used when Tayler-Spruit is included" )
    lsmooth_fact        = Float( DEFAULT_EXTRA['lsmooth_fact'], label="lsmooth_fact", tooltip="Factor to in/decrease caracteristic smoothing scale" )
    lsmooth_min        = Float( DEFAULT_EXTRA['lsmooth_min'], label="lsmooth_min", tooltip="Minimum caracteristic smoothing scale" )
    alpha_fuller        = Float( DEFAULT_EXTRA['alpha_fuller'], label="alpha_fuller", tooltip="alpha coefficient used in Fuller+2019's formulation" )
    ts_average          = Trait( DEFAULT_EXTRA['ts_average'], _dict_val_extra['ts_average'], label='ts_average', tooltip="Kernel used to average the TS instability diffusion coefficients" )
    ts_average_length   = Trait( DEFAULT_EXTRA['ts_average_length'], _dict_val_extra['ts_average_length'], label='ts_average_length', tooltip="Characteristic width for sliding average of Tayler-Spruit coefficient" )
    ts_weak_b           = Bool( DEFAULT_EXTRA['ts_weak_b'], label='ts_weak_b', tooltip="If True, assume small B field and TS dynamo sets on only when shear induced turbulence is high enough" )
    ts_hysteresis       = Bool( DEFAULT_EXTRA['ts_hysteresis'], label='ts_hysteresis', tooltip="If True, when shear decreases and Ri goes from < Ric to > Ric, if" )
    ts_patch_regime     = Bool( DEFAULT_EXTRA['ts_patch_regime'], label='ts_patch_regime', tooltip="If true, uses 'patched formula' from Spruit 2002 that combined the regimes" )
    ts_qlim             = Bool( DEFAULT_EXTRA['ts_qlim'], label='ts_qlim', tooltip="If true, account for the stability criterion (almost always verified)" )
    ts_extrapolate      = Bool( DEFAULT_EXTRA['ts_extrapolate'], label='ts_extrapolate', tooltip="Extrapolate TS instability diffusion coeff using Adams-Bashforth scheme" )

    # Constants
    source_cts          = Trait( DEFAULT_EXTRA['source_cts'], _dict_val_extra['source_cts'], label='source_cts', tooltip="Source for the constants (this may have no effect on some choice of constants)" )

    # extra outputs
    debug               = Bool( DEFAULT_EXTRA['debug'], label='debug', tooltip="write several files with intermediate variables" )
    convection_out      = Bool( DEFAULT_EXTRA['convection_out'], label='convection_out', tooltip="Writes additional information in a file about convection" )
    hse_out      = Bool( DEFAULT_EXTRA['hse_out'], label='hse_out', tooltip="Writes additional information in a file about hydrostatic equilibrium" )
    rota_out            = Trait( DEFAULT_EXTRA['rota_out'], _dict_val_extra['rota_out'], label='rota_out', tooltip="Writes rota.dat files" )
    ts_out              = Bool( DEFAULT_EXTRA['ts_out'], label='ts_out', tooltip="If true, writes additional information in a file about Tayler-Spruit instability" )
    lim_zc_out          = Bool( DEFAULT_EXTRA['lim_zc_out'], label='lim_zc_out', tooltip="Write debug info from lim_zc" )

    # HRnew output
    full_hr             = Bool( DEFAULT_EXTRA['full_hr'], label='full_hr', tooltip="If True, writes full output to the HR file" )
    writes_tau_conv     = Bool( DEFAULT_EXTRA['writes_tau_conv'], label='writes_tau_conv', tooltip="If True, calculate and writes tau_conv" )
    writes_nuc_neutr    = Bool( DEFAULT_EXTRA['writes_nuc_neutr'], label='writes_nuc_neutr', tooltip="If True, writes nuclear neutrinos" )
    writes_l_nuc        = Bool( DEFAULT_EXTRA['writes_l_nuc'], label='writes_l_nuc', tooltip="If True, write nuclear and grav. energies" )
    writes_mdot         = Bool( DEFAULT_EXTRA['writes_mdot'], label='writes_mdot', tooltip="If True, writes mass loss rate" )
    writes_conv         = Bool( DEFAULT_EXTRA['writes_conv'], label='writes_conv', tooltip="If true, writes alpha and quantities related to the adiabat" )
    writes_2d           = Bool( DEFAULT_EXTRA['writes_2d'], label='writes_2d', tooltip="If true, writes quantities related to the 2D strcuture" )
    writes_burn_zones   = Bool( DEFAULT_EXTRA['writes_burn_zones'], label='writes_burn_zones', tooltip="If True, writes where epsilon_nuc greater than burn_contours" )
    n_contours_burn     = Int( DEFAULT_EXTRA['n_contours_burn'], label='n_contours_burn', tooltip='Number of contours for burning zones' )
    burn_contours       = List( Float, DEFAULT_EXTRA['burn_contours'], label='burn_contours', tooltip="contours for burning zones only burn_contours(1:n_contours_burn) taken into account" )

    plato_track         = Bool( DEFAULT_EXTRA['plato_track'], label='plato_track', tooltip="If True, the HR file is replaced by a HDF5 file containing all information needed by the PLATO consortium. It takes precedence on all other kind of output: neither HR file nor osc files will be written." )

    dict_eos_dir  = {'etat_opal':   'eos_opalplus/', 'etat_opalX':  'eos_opal2001/',
                     'etat_opalZ':  'eos_opal2001/', 'etat_opal5Z': 'eos_opal2005/',
                     'etat_opal5Z_spline': 'eos_opal2005/'}
    dict_make_eos = {'etat_opalX':  'make_eos.x',    'etat_opalZ':  'make_eos.x',
                     'etat_opal5Z': 'make_eos5.x', 'etat_opal5Z_spline': 'make_eos5.x'}

    eos_dir = 'eos_opal2005'


    def __init__( self, name, verbose=True ):
        """
        <CParameters.__init__( name )>

        Class that represents the input parameters of a Cesam2k20 model, to be manipulated
        through GUI.

        :param name: Name of the model.
        :type name: string

        :kparam verbose: If True, prints additional information.
        :ktype verbose: dict
        """
        self.initialized  = False
        self.name         = name
        self.don          = name + '.don'
        self.verbose      = verbose
        self.write_rgl    = False
        self.write_extra  = False
        self._subnlist_changed = {'num':[], 'extra':[]}
        self._num_changed = []
        self._extra_changed = []
        if os.path.exists(self.don): # read data file
            self.read_params( )
        self._precision_changed( )

        if len( self.nom_chemin ) > 120:
            raise CESAMError( 'The string variable `nom_chemin` is too long.'
                'Only 120 characters allowed.' )


        if not self.nom_chemin.endswith('/'):
            self.nom_chemin = self.nom_chemin + '/'

        self.f_eos[0] = self.nom_chemin + self.eos_dir + '/' + self.f_eos[0]
        self.f_opa[0] = self.nom_chemin + 'OPA/opa_yveline/' + self.f_opa[0]
        self.f_opa[1] = self.nom_chemin + self.f_opa[1]
        self.f_opa[2] = self.nom_chemin + self.f_opa[2]
        self.f_opa[3] = self.nom_chemin + self.f_opa[3]

        self.initialized = True


    def __dict_keys(self, dct, value, ref_dict=_dict_val):
        """
        <CParameters.__dict_keys( dct, value )>

        Fetch the key associated to a value in a dictionnary of dictionnaries.

        :param dct: Dictionary of dictionnaries form which the key should be fetched.
        :type dct: dict of dict

        :param value: Value that correspond to the kwy that should be fetched.
        :type value: any object
        """
        if dct not in ref_dict:
            return value
        else:
            for key in ref_dict[dct]:
                if ref_dict[dct][key] == value:
                    return key

        raise TraitError(f"Value {value} was not found in {dct}.")


    def __make_eos(self, curr_dir, eosfile, exec_path):
        """
        Compute new EoS table.

        :param curr_dir: Directory where the computation of the model should take place.
        :type curr_dir: string

        :param eosfile: Name of the EoS file.
        :type eosfile: string

        :param exec_path: Path to `make_eos5.x` executable.
        :type exec_path: string
        """
        cmd = f"{self.z0}\n{eosfile.replace('//', '/')}\n{exec_path}\n"

        proc = subp.Popen(exec_path + self.dict_make_eos[self.nom_etat_], stdin=subp.PIPE, stdout=subp.PIPE)

        try:
            messg, err = proc.communicate(input=cmd.encode())
        except TypeError:
            print(f'Could not make file {self.f_eos[0]}\n')
        finally:
            if curr_dir:
                os.path.join(curr_dir)


    def __process_value(self, name):

        if name in self._dict_val.keys() or \
           name in self._dict_val_num.keys() or \
           name in self._dict_val_extra.keys():
            value = self.__getattribute__(name+'_')
        else:
            value = self.__getattribute__(name)

        if name == 'nom_ctes': value = self.__change_ctes( value )

        if isinstance(value, str): value = str(value)

        if name in FS:
            for i, v in enumerate(value):
                value[i] = v.replace( ' ', '' )
            value = map(str, value)
            value = map(repr, value)
            value = ', '.join(value)
        elif name in FF or name in FI:
            value = map(repr, value)
            value = ', '.join(value)
        else:
            if isinstance(value, bool):
                value = dict_bool[value]
            else:
                value = repr(value)

        return value

    def __change_ctes( self, value ):
        if value in ['ctes_94', 'ctes_85', 'ctes_94m', 'ctes_94_asplund', 'ctes_gs98',
            'ctes_aag21', 'ctes_aag21phot', 'ctes_mg22phot']:
            print(f"{RED} [[DEPRECATED]] {NO_COLOR}", end='')
            if value in ['ctes_85', 'ctes_94', 'ctes_94_asplund']:
                print(f"{RED} You choose nom_ctes={value}. It now corresponds to ctes_lide1994_nacre. {NO_COLOR}")
                value = 'ctes_lide1994_nacre'
            elif value in ['ctes_94m']:
                print(f"{RED} You choose nom_ctes={value}. It now corresponds to ctes_lide1994_nacre_int. {NO_COLOR}")
                value = 'ctes_lide1994_nacre_int'
            elif value in ['ctes_gs98', 'ctes_aag21', 'ctes_aag21phot', 'ctes_mg22phot']:
                if self.__getattribute__('source_cts') == 'IAU2015_CODATA2018':
                    print(f"{RED} You choose nom_ctes={value}. It now corresponds to ctes_iau2015_codata2018. {NO_COLOR}")
                    value = 'ctes_iau2015_codata2018'
                elif self.__getattribute__('source_cts') == 'Lide1994':
                    print(f"{RED} You choose nom_ctes={value}. It now corresponds to ctes_lide1994_nacre. {NO_COLOR}")
                    value = 'ctes_lide1994_nacre'
                print(f"{RED} extra%source_cts is also DEPRECATED. {NO_COLOR}", end='')

            print(f"{RED} The change was made automatically {NO_COLOR}")

        return value


    def _nom_chemin_changed(self):
        """
        <CParameters._nom_chemin_changed( )>
        Automatically every time `nom_chemin` is changed.

        :raises CESAMError: Length should not exceed 120 characters.
        """
        if len( self.nom_chemin ) > 120:
            raise CESAMError( 'The string variable `nom_chemin` is too long.'
                'Only 120 characters allowed.' )

    def _nom_pertm_changed(self):

        err_mess = """No Value was specified for {:s} and we tried to use the value of NOM_PERTM = {:s},
            which is invalid. In your .don input file, please specify a value for NOM_PERTM_SOLAR
            among: {:s}."""


        if self.nom_pertm_ in list( self._dict_val['nom_pertm_solar'].values() ):
            self.nom_pertm_solar = self.nom_pertm
        else:
            raise CESAMError( err_mess.format('nom_pertm_solar', self.nom_pertm,
                inline_str( list( self._dict_val['nom_pertm_solar'].values() ) ) ) )

        if self.nom_pertm_ in list( self._dict_val['nom_pertm_rgb'].values() ):
            self.nom_pertm_rgb = self.nom_pertm
        else:
            raise CESAMError( err_mess.format('nom_pertm_rgb', self.nom_pertm,
                inline_str( list( self._dict_val['nom_pertm_rgb'].values() ) ) ) )

        if self.nom_pertm_ in list( self._dict_val['nom_pertm_agb'].values() ):
            self.nom_pertm_agb = self.nom_pertm
        else:
            raise CESAMError( err_mess.format('nom_pertm_agb', self.nom_pertm,
                inline_str( list( self._dict_val['nom_pertm_agb'].values() ) ) ) )

        if self.nom_pertm_ in list( self._dict_val['nom_pertm_hot'].values() ):
            self.nom_pertm_hot = self.nom_pertm
        else:
            raise CESAMError( err_mess.format('nom_pertm_hot', self.nom_pertm,
                inline_str( list( self._dict_val['nom_pertm_hot'].values() ) ) ) )

    def _precision_changed(self):
        """
        <CParameters._precision_changed( )>
        Automatically every time `precision` is changed.

        :raises CESAMError: Length should not exceed 120 characters.
        """
        for k in PRECISIONS[self.precision_]:
            if k not in self._subnlist_changed['num']:
                self.__setattr__( k, PRECISIONS[self.precision_][k] )

    def _check_changes( self ):
        for n in self.nlsublist['num']:
            if not n in self._subnlist_changed['num']:
                n_ = self.__getattribute__( n )
                if n_ != DEFAULT_NUM[n]:
                    self._subnlist_changed['num'].append( n )

        for e in self.nlsublist['extra']:
            if not e in self._subnlist_changed['extra']:
                e_ = self.__getattribute__( e )
                if e_ != DEFAULT_EXTRA[e]:
                    self._subnlist_changed['extra'].append( e )


    def _x0_changed(self):
        if self.initialized:
            if self.x_or_y == 'Initial X':
                if self.z_or_zsx == 'Initial Z':  # Keeps initial Z
                    self.zsx0 = self.z0/self.x0
                else:
                    self.z0 = self.zsx0*self.x0

                self.y0 = 1.0 - self.x0 - self.z0


    def _y0_changed(self):
        if self.initialized:
            if self.x_or_y == 'Initial Y':
                if self.z_or_zsx == 'Initial Z':  # Keeps initial Z
                    self.x0 = 1.0 - self.y0 - self.z0
                    self.zsx0 = self.z0/self.x0
                else:
                    self.x0 = (1.0 - self.y0)/(1.0 + self.zsx0)
                    self.z0 = self.zsx0*self.x0


    def _z0_changed(self):
        if self.initialized:
            if self.z_or_zsx == 'Initial Z':  # Keeps initial Z
                if self.x_or_y == 'Initial X':
                    self.y0 = 1.0 - self.x0 - self.z0
                    self.zsx0 = self.z0/self.x0
                else:
                    self.x0 = 1.0 - self.y0 - self.z0
                    self.zsx0 = self.z0/self.x0


    def _zsx0_changed(self):
        if self.initialized:
            if self.z_or_zsx == 'Initial Z/X':  # Keeps initial Z/X
                if self.x_or_y == 'Initial X':
                    self.z0 = self.zsx0*self.x0
                    self.y0 = 1.0 - self.x0 - self.z0
                else:
                    self.x0 = (1.0 - self.y0)/(1.0 + self.zsx0)
                    self.z0 = self.zsx0*self.x0

    def change_nom_output(self):
        if ((self.trunc_nom_output_ == 'no_output') or self.trunc_nom_output_.endswith('ascii')) :
            self.nom_output_ = self.trunc_nom_output_
        else:
            if self.out_2d:
                self.nom_output_ = self.trunc_nom_output_.replace('2d', '') + '2d'
            else:
                self.nom_output_ = self.trunc_nom_output_.replace('2d', '')
            if self.out_h5:
                self.nom_output_ = self.nom_output_.replace('h5', '') + 'h5'
            else:
                self.nom_output_ = self.nom_output_.replace('h5', '')

    def _trunc_nom_output_changed(self):
        self.change_nom_output()

    def _out_2d_changed(self):
        self.change_nom_output()
        if self.out_2d:
            self.acor = True
        else:
            self.acor = False

    def _out_h5_changed(self):
        self.change_nom_output()


    def init_params(self):
        if os.path.exists(self.don): # read data file
            self.read_params()

    def copy(self, other):
        name = self.name
        self.params = self.copy_traits(other)
        self.name = name
        self.don = name + '.don'

    def set_param( self, line, namelist, ref_dict=_dict_val ):
        if namelist in FS:
            files = re.findall(pat_str, line)
            n = len(files)
            if n < 8: files.extend(' '*(8 - n))
            self.__setattr__(namelist, files)
        elif namelist in FF:
            line = line.split('=')[1]
            floats = [eval(w.lower().replace( 'd', 'e' )) for w in line.split(',')[:-1]]
            self.__setattr__(namelist, floats)
        elif namelist in FI:
            line = line.split('=')[1]
            ints = [int(w.lower().replace( 'd', 'e' )) for w in line.split(',')[:-1]]
            self.__setattr__(namelist, ints)
        else:
            if re.search(pat_char, line):
                value = re.findall(pat_char, line)[0].upper()
                if value == 'T':
                    self.__setattr__(namelist, True)
                else:
                    self.__setattr__(namelist, False)
            elif re.search(pat_bool, line, re.I):
                value = re.findall(pat_bool, line, re.I)[0][0].upper()
                if value == '.TRUE.':
                    self.__setattr__(namelist, True)
                else:
                    self.__setattr__(namelist, False)
            elif re.search(pat_str, line):
                value = re.findall(pat_str, line)[0]
                value = value.replace( ' ', '' )
                if value == 'difft_nu': value = 'difft_cte'
                if namelist == 'nom_output':
                    self.__setattr__('out_2d', '2d' in value)
                    self.__setattr__('out_h5', 'h5' in value)
                    tmp_value = value
                key = self.__dict_keys(namelist, value, ref_dict=ref_dict )
                self.__setattr__(namelist, key)
                if namelist == 'nom_output':
                    tmp_value = self.__dict_keys(namelist, tmp_value.replace('2d', '').replace('h5', ''), ref_dict=ref_dict )
                    self.__setattr__('trunc_nom_output', tmp_value)
            elif re.search(pat_flt, line):
                value = re.findall(pat_flt, line)[0][0]
                value = re.sub('[DdE]', 'e', value)
                self.__setattr__(namelist, eval(value))

    def read_params(self):
        """
        <CParameters.read_params( )>

        Reads a .don file so to obtain model parameters and physical inputs.
        """
        if self.verbose: print(f"Reading {self.don}...", end='')
        with open(self.don, 'r') as f:
            inl = f.readlines()

        self.num          = []
        self.extra        = []

        n_lines           = len(inl)

        for i in range(n_lines):
            if 'NL_RLG' in inl[i].upper():
                j = i + 1
                while '/' not in inl[j]:
                    self.num.append( inl[j] )
                    j += 1
            if 'NL_EXTRA' in inl[i].upper():
                j = i + 1
                while '/' not in inl[j]:
                    self.extra.append( inl[j] )
                    j += 1


        for line in inl:
            if re.search(r'^\s*\b\w+\b\s*=', line):
                namelist = re.findall(r'^\s*\b(\w+)\b\s*=', line)[0].lower()
                self.set_param( line, namelist )

        for num in self.num:
            var = num.split( '%' )[1].split( '=' )[0].strip()
            self.set_param( num, var.lower(), ref_dict=self._dict_val_num )
            self._subnlist_changed['num'].append( var )



        for extra in self.extra:
            var = extra.split( '%' )[1].split( '=' )[0].strip()
            self.set_param( extra, var.lower(), ref_dict=self._dict_val_extra )
            self._subnlist_changed['extra'].append( var )


        self.z0 = self.zsx0*self.x0
        self.all_rep = self.nb_max_modeles < 0.0
        if self.verbose: print(f"{GREEN}[Done]{NO_COLOR}")
        self.nom_ctes = self.__change_ctes( self.nom_ctes )



    def mkdon( self, overwrite=True, replace_eos=True, verbose=True, eos_inplace=False ):
        """
        <CParameters.mkdon( overwrite=True, replace_eos=True, verbose=True, eos_inplace=False )>

        Creates a .don file with the model parameters. Creates an EoS OPAL file
        for the right metalicity.

        :param overwrite: If True, overwrites an existing file (default: True, optional)
        :type overwrite: boolean

        :param replace_eos: If True, overwrites EoS OPAL file (default: True, optional).
        :type replace_eos: boolean

        :param verbose: If True, print some details (default: True, optional).
        :type verbose: boolean

        :param eos_inplace: If True, the eos file is stored in the folder where the computation tkes place.
            (default: True, optional)
        :type eos_inplace: boolean
        """


        if self.all_rep:
            self.nb_max_modeles =-abs(self.nb_max_modeles)
        else:
            self.nb_max_modeles = abs(self.nb_max_modeles)

        filenames = deepcopy(self.f_eos)
        filenames[0] = os.path.basename(filenames[0])
        self.trait_set(f_eos=filenames, trait_change_notify=False)

        filenames = deepcopy(self.f_opa)
        filenames[0] = os.path.basename(filenames[0])
        filenames[1] = os.path.basename(filenames[1])
        filenames[2] = os.path.basename(filenames[2])
        filenames[3] = os.path.basename(filenames[3])
        self.trait_set(f_opa=filenames, trait_change_notify=False)

        eos_inplace = self.eos_inplace

        self._check_changes( )

        model_dir = os.getcwd()
        if 'opal' in self.nom_etat_:
            eosdir = f"{os.environ.get('CESDIR')}/SUN_STAR_DATA/EOS/{self.dict_eos_dir[self.nom_etat_]}".replace( '//', '/' )
            exec_path = f"{eosdir}/".replace( '//', '/' )
            if eos_inplace:
                eosfile   = f"{model_dir}/{self.f_eos[0]}".replace( '//', '/' )
            else:
                eosfile   = f"{self.nom_chemin}/EOS/{self.dict_eos_dir[self.nom_etat_]}/{self.f_eos[0]}".replace( '//', '/' )

            if replace_eos:
                if not eos_inplace:
                    os.path.join( eosdir )

                if os.path.exists(eosfile):
                    if self.nom_etat_ not in ['etat_opal', 'etat_opalX']:
                        cmd = f"{eosfile}\n"
                        proc = subp.Popen(f'{exec_path}read_Z.x', stdin=subp.PIPE,
                            stdout=subp.PIPE)

                        try:
                            messg, err = proc.communicate(input=cmd.encode())
                            messg = messg.decode().split('\n')
                            try:
                                z_file = max( abs( float(messg[0]) ), 1e-10 )
                            except ValueError:
                                print(f'Error in file {self.f_eos[0]}, rebuilding.')
                                self.__make_eos(model_dir, eosfile, exec_path)
                                z_file = self.z0

                            if abs((z_file - self.z0)/z_file) > 1.0e-8:
                                if verbose:
                                    print(f'Z = {z_file} in file {self.f_eos[0]}')
                                    print(f'Z = {self.z0} in model; rebuilding file {self.f_eos[0]}')
                                self.__make_eos(model_dir, eosfile, exec_path)
                            else:
                                if verbose:
                                    print(f'Z = {z_file} in file {self.f_eos[0]}')
                                    print(f'Z = {self.z0} in model; agrees!')

                        except TypeError:
                            print(f'Error in file {self.f_eos[0]}, rebuilding.')
                            self.__make_eos(model_dir, eosfile, exec_path)
                else:
                    self.__make_eos(model_dir, eosfile, exec_path)

        os.path.join( model_dir )


        nlists = ['nl_cesam', 'nl_mass', 'nl_evol', 'nl_chim', 'nl_conv', 'nl_diff',
                  'nl_rot', 'nl_etat', 'nl_opa', 'nl_nuc', 'nl_atm', 'nl_mag']
        nsublist_dict = {'num':'RLG', 'extra':'EXTRA'}

        if overwrite:
            cont = []
            for nlist in nlists:
                cont.append(f' &{nlist.upper()}\n')
                for opt in self.nlist[nlist][:-1]:
                    value = self.__process_value(opt)
                    cont.append(f' {opt.upper()} = {str(value)},\n')

                # end of namelist
                opt = self.nlist[nlist][-1]
                value = self.__process_value(opt)
                cont.append(f' {opt.upper()} = {str(value)}\n')
                cont.append(' /\n')

            for key in nsublist_dict.keys():
                cont.append(f' &NL_{nsublist_dict[key]}\n')
                if len( self._subnlist_changed[key] ) > 1:
                    lnl  = self._subnlist_changed[key][:-1]
                    lnl1 = self._subnlist_changed[key][-1]
                elif len( self._subnlist_changed[key] ) > 0:
                    lnl = []
                    lnl1 = self._subnlist_changed[key][0]
                else:
                    lnl = []
                    lnl1 = None


                for opt in lnl:
                    value = self.__process_value(opt)
                    cont.append(f' {key}%{opt.lower()} = {str(value)},\n')
                # end of namelist
                if lnl1 is not None:
                    opt   = lnl1
                    value = self.__process_value(opt)
                    cont.append(f' {key}%{opt.lower()} = {str(value)}\n')
                cont.append(' /\n')

            with open(self.don, 'w') as f:
                f.writelines(cont)
