#!/usr/bin/env python3
# coding: utf-8

# Code d'Evolution Stellaire, Adaptatif et Modulaire for the 2020 decade.
# Copyright (c) 1997-2025 The Cesam2k20 authors
# SPDX-License-Identifier : GPL-3.0-or-later


# lsun     = 3.846e33              # solar luminosity
# msun     = 1.98919e33            # solar mass
mearth   = 5.977e27              # earth mass
# rsun     = 6.9599e10             # solar radius
# gmsun    = 1.32712438e26
# ggrav    = gmsun/msun            # grav. constant G
mbol_sun = 4.75                  # solar bolometric magnitude
w_sun    = 2.86e-6               # solar angular velocity
bv_sun   = 0.64                  # solar B-V (Neckel (1986): +0.650
amu      = 1.6605402e-24         # atomic mass unit
avogadro = 1.0/amu               # avogadro number
aradia   = 7.565912199839849e-15 # radiation constant a
clight   = 2.99792458e10         # light speed
echarg   = 4.8032068e-10         # electron charge
eve      = 1.60217733e-12        # electron volt
hpl      = 6.6260755e-27         # planck constant
kbol     = 1.380658e-16          # boltzman constant
rgas     = kbol*avogadro         # perfect gas constant
mel      = 9.1093897e-28         # electron mass
sigma    = aradia*clight/4.0     # Stefan's constant
ua       = 149.59787066e6*1.0e5  # astronomical unit

# solar metalicity
zsx_sun = {'solaire_aag21' : 0.0187,
           'solaire_ags09_0' : 0.0181, #M. Deal: opacity tables not yet available in SUN_STAR_DATA
           'solaire_ags09' : 0.01781,  #M. Deal
           'solaire_ags09_alpha' : 0.01781, #M. Deal
           'solaire_gn' : 0.0244,
           'solaire_gs98': 0.0231,
           'solaire_ags05' : 0.0165}

Nfeh_sun = {'solaire_aag21' : -4.54,
           'solaire_ags09_0' : -4.5,
           'solaire_ags09' : -4.55,
           'solaire_ags09_alpha' : -4.55,
           'solaire_gn' : -4.5,
           'solaire_gs98': -4.5,
           'solaire_ags05' : -4.55}

ggrav_dict = {
    'ctes_lide1994_nacre' : 6.67168e-08,
    'ctes_lide1994_nacre_int' : 6.67168e-08,
    'ctes_iau2015_codata2018' : 6.67430e-8}

msun_dict = {
    'ctes_lide1994_nacre' : 1.98919e33,
    'ctes_lide1994_nacre_int' : 1.98919e33,
    'ctes_iau2015_codata2018' : 1.98841e+33}

rsun_dict = {
    'ctes_lide1994_nacre' : 6.9599e10,
    'ctes_lide1994_nacre_int' : 6.9599e10,
    'ctes_iau2015_codata2018' : 6.957e+10}

lsun_dict = {
    'ctes_lide1994_nacre' : 3.846e33,
    'ctes_lide1994_nacre_int' : 3.846e33,
    'ctes_iau2015_codata2018' : 3.828e33}


nucleo = {'H1'  : 1.007825,
          'H2'  : 2.0141018,
          'He3' : 3.0160293,
          'He4' : 4.0026033,
          'Li6' : 6.015121,
          'Li7' : 7.0160040,
          'Be7' : 7.0169292,
          'Be9' : 9.0121821,
          'B11' : 11.0093055,
          'C12' : 12.0,
          'C13' : 13.0033548,
          'N13' : 13.0057386,
          'N14' : 14.003074,
          'N15' : 15.001089,
          'O16' : 15.9949146,
          'O17' : 16.9991315,
          'O18' : 17.9991604,
          'Fe56': 55.847,
          'F18' : 18.0009377,
          'F19' : 18.9984032,
          'Ne20': 19.9924402,
          'Ne21': 20.9938467,
          'Ne22': 21.9913855,
          'Na23': 22.9897697,
          'Mg23': 22.9941249,
          'Mg24': 23.9850419,
          'Mg25': 24.985837,
          'Mg26': 25.982593,
          'Si28': 28.0855,
          'Al27': 26.9854,
          'S32' : 31.972070,
          'P31' : 30.973762,
          'n'   : 1.008665}

RUN_DICT = {'job' : {
        'From PMS' : 'pms',
        'From ZAMS' : 'zams',
        'From previous model' : 'rep',
        'Frequencies' : 'freqs'}}

PARAMS_DICT = {   'nom_ctes' : {
                               'ctes_lide1994_nacre' : 'ctes_lide1994_nacre',
                               'ctes_lide1994_nacre_int' : 'ctes_lide1994_nacre_int',
                               'ctes_iau2015_codata2018' : 'ctes_iau2015_codata2018',
                               'ctes_94' : 'ctes_94',
                               'ctes_85' : 'ctes_85',
                               'ctes_94m' : 'ctes_94m',
                               'ctes_94_asplund' : 'ctes_94_asplund',
                               'ctes_gs98' : 'ctes_gs98',
                               'ctes_aag21' : 'ctes_aag21',
                               'ctes_aag21phot' : 'ctes_aag21phot'},
                  'nom_des' : {'In r' : 'des_r',
                               'In m' : 'des_m',
                               'Zoom' : 'zoom',
                               'None' : 'no_des'},
                  'nom_output' : {'Final model (adiabatic)'                 : 'osc_adia',
                                  'Final model (adiabatic; 2D)'             : 'osc_adia2d',
                                  'Final model (adiabatic; HDF5)'           : 'osc_adiah5',
                                  'Final model (adiabatic; 2D + HDF5)'      : 'osc_adia2dh5',
                                  'Final model (inversions)'                : 'osc_invers',
                                  'Final model (inversions; 2D)'            : 'osc_invers2d',
                                  'Final model (inversions; HDF5)'          : 'osc_inversh5',
                                  'Final model (inversions; 2D + HDF5)'     : 'osc_invers2dh5',
                                  'Final model (non-adiabatic)'             : 'osc_nadia',
                                  'Final model (non-adiabatic; 2D)'         : 'osc_nadia2d',
                                  'Final model (non-adiabatic; HDF5)'       : 'osc_nadiah5',
                                  'Final model (non-adiabatic; 2D + HDF5)'  : 'osc_nadia2dh5',
                                  'Final model (PLATO format)'              : 'osc_plato',
                                  'Final model (PLATO; 2D)'                 : 'osc_plato2d',
                                  'Final model (PLATO; HDF5)'               : 'osc_platoh5',
                                  'Final model (PLATO; 2D + HDF5)'          : 'osc_plato2dh5',
                                  'All models (adiabatic)'                  : 'all_adia',
                                  'All models (adiabatic; 2D)'              : 'all_adia2d',
                                  'All models (adiabatic; HDF5)'            : 'all_adiah5',
                                  'All models (adiabatic; 2D + HDF5)'       : 'all_adia2dh5',
                                  'All models (inversions)'                 : 'all_invers',
                                  'All models (inversions; 2D)'             : 'all_invers2d',
                                  'All models (inversions; HDF5)'           : 'all_inversh5',
                                  'All models (inversions; 2D; HDF5)'       : 'all_invers2dh5',
                                  'All models (non-adiabatic)'              : 'all_nadia',
                                  'All models (non-adiabatic; 2D)'          : 'all_nadia2d',
                                  'All models (non-adiabatic; HDF5)'        : 'all_nadiah5',
                                  'All models (non-adiabatic; 2D; HDF5)'    : 'all_nadia2dh5',
                                  'All models (PLATO format)'               : 'all_plato',
                                  'All models (PLATO format; 2D)'           : 'all_plato2d',
                                  'All models (PLATO format; HDF5)'         : 'all_platoh5',
                                  'All models (PLATO format; 2D + HDF5)'    : 'all_plato2dh5',
                                  'None'                                    : 'no_output'},
                  'nom_osc' : {'ADIPLS' : 'adipls',
                               'None' : 'none'},
                  'trunc_nom_output' : {'Final model (adiabatic)'           : 'osc_adia',
                                        'Final model (inversions)'          : 'osc_invers',
                                        'Final model (non-adiabatic)'       : 'osc_nadia',
                                        'Final model (PLATO format)'        : 'osc_plato',
                                        'All models (adiabatic)'            : 'all_adia',
                                        'All models (inversions)'           : 'all_invers',
                                        'All models (non-adiabatic)'        : 'all_nadia',
                                        'All models (PLATO format)'         : 'all_plato',
                                        'None'                              : 'no_output',
                                        'Final model (ASCII file)'          : 'ascii',
                                        'All models (ASCII file)'           : 'all_ascii'},
                  'precision' : {'Realistic precision' : 'pr',
                                 'Realistic Eulerian precision' : 'er',
                                 'Realistic precision with adjustments' : 'aj',
                                 'High precision (Lagrangean)' : 'sp',
                                 'High precision (Eulerian)' : 'sr',
                                 'Long run' : 'av',
                                 'CoRoT' : 'co',
                                 'PLATO' : 'pl',
                                 'Mariejo' : 'mj',
                                 'Evolved low mass' : 'lm',
                                 'Solar accuracy' : 'sa',
                                 'Ultra high precision' : 'hp',
                                 'Normal precision (Lagrangean)' : 'np',
                                 'Normal precision (Eulerian)' : 'nr',
                                 'Maximum number of shells' : 'mx',
                                 'External file rg' : 'rg'},
                  'arret' : {'Other' : 'else',
                             'ZAMS' : 'zams',
                             'TAMS' : 'post',
                             'He burning' : 'cohe',
                             'C burning' : 'coca',
                             'O burning' : 'coox'},
                  'nom_pertm' : {'No mass loss'                            : 'pertm_0',
                                 'Constant mass loss rate'                 : 'pertm_ext',
                                 'Constant mass loss rate, stops at 1Msun' : 'pertm_msol',
                                 'Mass loss mc^2 (nuclear burning)'        : 'pertm_tot',
                                 'Reimers (1975)'                          : 'pertm_reimers',
                                 'Waldron (1990)'                          : 'pertm_waldron',
                                 'Blöcker (1994)'                          : 'pertm_B_1994',
                                 'Cranmer and Saar (2011)'                 : 'pertm_CS_2011',
                                 'van Loon et al. (2005)'                  : 'pertm_vanLoon',
                                 'de Jager et al. (1987)'                  : 'pertm_JNH_1987',
                                 'Nieuwenhuijzen and de Jager (1990)'      : 'pertm_NJ_1990',
                                 'Vink et al. (2001)'                      : 'pertm_vink'},
                  'nom_pertm_solar' : {'No mass loss'                      : 'pertm_0',
                                 'Constant mass loss rate'                 : 'pertm_ext',
                                 'Constant mass loss rate, stops at 1Msun' : 'pertm_msol',
                                 'Reimers (1975)'                          : 'pertm_reimers',
                                 'Mass loss mc^2 (nuclear burning)'        : 'pertm_tot'},
                  'nom_pertm_rgb' : {'No mass loss'                        : 'pertm_0',
                                 'Constant mass loss rate'                 : 'pertm_ext',
                                 'Reimers (1975)'                          : 'pertm_reimers',
                                 'van Loon et al. (2005)'                  : 'pertm_vanLoon'},
                  'nom_pertm_agb' : {'No mass loss'                        : 'pertm_0',
                                 'Constant mass loss rate'                 : 'pertm_ext',
                                 'Reimers (1975)'                          : 'pertm_reimers',
                                 'van Loon et al. (2005)'                  : 'pertm_vanLoon'},
                  'nom_pertm_hot' : {'No mass loss'                        : 'pertm_0',
                                 'Constant mass loss rate'                 : 'pertm_ext',
                                 'Vink et al. (2001)'                      : 'pertm_vink'},
                  'nom_abon' : {'Solar (AAG 2021)' : 'solaire_aag21',
                                'Solar (AGS 2009)' : 'solaire_ags09_0',
                                'Solar (AGS 2009) + S10 recom.' : 'solaire_ags09',
                                'Solar (AGS 2009) + alpha enhancement' : 'solaire_ags09_alpha',
                                'Solar (GN 1993)' : 'solaire_gn',
                                'Solar (GS 1998; with OPAL/Wichita op. tables)' : 'solaire_gs98',
                                'Solar (AGS 2005)' : 'solaire_ags05',
                                'Solar (AGS 2005) with Li' : 'solaire_ags05i',
                                'Meteorites (AG 1989)' : 'meteorites_ag',
                                'Meteorites (GS 1998)' : 'meteorites_gs',
                                'Custom' : 'mixture'},
                  'nom_conv' : {u'MLT with l\u27F60 at CZ/RZ' : 'conv_a0',
                                u'CM 1991 with l=\u03B1 Hp' : 'conv_cm',
                                u'CM 1991 with l=\u03B1 Hp with ' : 'conv_cm_reza',
                                u'CGM 1996 with l=\u03B1 Hp' : 'conv_cgm_reza',
                                u'MLT with l=\u03B1 Hp' : 'conv_jmj',
                                u'MLT (Ireland2018)' : 'conv_ireland2018',
                                u'MLT with l=\u03B1 Hp, thick' : 'conv_jmj_thick'},
                  'nom_alpha': {u'Fixed alpha' : 'alpha_f',
                                u'Prescription of Sonoi+19' : 'alpha_sonoi19',
                                u'Entropy-calibrated model with Magic+13'        : 'entropy_magic13',
                                u'Entropy-calibrated model with Ludwig+99'       : 'entropy_ludwig99',
                                u'Entropy-calibrated model with Tanner+16'       : 'entropy_tanner16',
                                u'Entropy-calibrated model with custom function' : 'entropy_custom'},
                  'nom_ovshts': {u'Fixed core overshoot' : 'ovshts_f',
                                 u'Fixed core overshoot with diffusive transport' : 'ovshts_f_diff',
                                u'[step overshoot only] Prescription of Claret \& Torres (2016)' : 'ovshts_claret16',
                                u'Prescription of Claret \& Torres (2018)' : 'ovshts_claret18',
                                u'Prescription of Claret \& Torres (2018) with diffusive transport' : 'ovshts_claret18_diff'},
                  'nom_diffm' : {'Burgers formalism' : 'diffm_br',
                                 'Michaud and Proffit formalism' : 'diffm_mp',
                                 'None' : 'diffm_0'},
                  'nom_difft' : {'None' : 'difft_c0',
                                 'Constant turbulent coeff.' : 'difft_cte',
                                 'Exponential  mixing' : 'difft_expo',
                                 'Montreal mixing' : 'difft_montreal',
                                 'Montreal mixing at cz' : 'difft_montreal_cz',
                                 'Mixing grille: 7Li solaire -> calibration K. Verma+19' : 'difft_grille',
                                 'Precription of M. Gabriel' : 'difft_gab',
                                 'Coefficients bellow solar CZ according to Gabriel' : 'difft_sun'},
                  'nom_difft_fing' :{'None' : 'no_fing',
                                     'Prescription of Brown et al. 2013' : 'brown2013'},
                  'nom_difft_smc' :{'None' : 'no_smc',
                                     'Prescription of Langer et al. 1983' : 'langer1983'},
                  'nom_frad' : {'None' : 'no_frad',
                                'Formalism of G. Alecian 2004' : 'alecian2004',
                                'Formalism of G. Alecian 2020' : 'alecian2020'},
                  'unit' : {'kms/s' : 'kms/s',
                            'rad/s' : 'rad/s'},
                  'nom_thw' : {'None' : 'rot_0',
                               'Solid body' : 'cons_glob_mnt_cin',
                               'Local conservation of angular momentum' : 'cons_loc_mnt_cin',
                               'Turbulent diffusion and m. circulation (Talon & Zahn 1997)' : 'diff_tz97',
                               'Turbulent diffusion and m. circulation (MZ 2004)' : 'diff_mz04',
                               'Turbulent diffusion and m. circulation (MZ 2004) in 2D' : 'diff_mz04_2d',
                               'Rotation profile read from file (used when using libcestam.so)': 'file'
                               },
                  'rot_zc' : {'Solid body' : 'solide',
                              'Local conservation of angular momentum' : 'local_j'},
                  'nom_pertw' : {'None' : 'pertw_0',
                                 'Kawaler (1988)' : 'pertw_kaw',
                                 'Matt et al. (2015)' : 'pertw_matt15',
                                 'Reiners and Mohanty (2012)' : 'pertw_rm',
                                 'Schumanich(1972)' : 'pertw_sch',
                                 'Proportional to local rotation kinatic energy' : 'pertw_loc',
                                 'Resulting from mass loss' : 'pertw_ptm'},
                  'nom_diffw' : {'None' : 'diffw_0',
                                 'Mathis et al. (2018)' : 'diffw_mathis2018',
                                 'Mathis, Palacios and Zahn (2004)' : 'diffw_mpz',
                                 'Palacios (2003)' : 'diffw_p03',
                                 'Constant (Deff=300, Dh=1e6, Dv=250)' : 'diffw_cte'},
                  'nom_gsf' : {'None' : 'no_gsf',
                               'Hirschi+2010' : 'hirschi10',
                               'Barker+2019, equator' : 'barker19_eq'},
                  'nom_diffmg_ts' : {'None' : 'diffmg_ts_0',
                                  'Taylor-Spruit dynamo (Spruit+2002)' : 'spruit2002',
                                  'Taylor-Spruit dynamo (Maeder & Meynet 2004)' : 'maeder2004',
                                  'Taylor-Spruit dynamo (Fuller+2019)' : 'fuller2019',
                                  'Taylor-Spruit dynamo (Daniel+2023)' : 'daniel2023'},
                  'nom_etat'  : {u'OPAL 2005, Z\u22600' : 'etat_opal5Z',
                                 u'OPAL 2005 Spline, Z\u22600' : 'etat_opal5Z_spline',
                                 'CEFF' : 'etat_ceff',
                                 'EFF' : 'etat_eff',
                                 'Simplified EoS, H and He fully ionised (GONG1)' : 'etat_gong1',
                                 'Simplified EoS, H and He4 only (GONG2)' : 'etat_gong2',
                                 'MHD' : 'etat_mhd',
                                 'EoS SAHA-S' : 'etat_saha',
                                 'CMS19' : 'etat_cms19',
                                 'Default CMS 2019' : 'etat_cms19_tab',
                                 'OPAL/CMS Merger' : 'merged_etat_cms19',
                                 'IVL CMS19' : 'mix_etat_cms19',
                                 'Ideal Mixed CMS 2019' : 'mix_ideal_etat_cms19',
                                 'IVL OPAL/CMS Merger' : 'merged_mix_etat_cms19',
                                 'Merged and Ideal Mixed CMS 2019' : 'merged_mix_ideal_etat_cms19',
                                 'CD21' : 'etat_cd21',
                                 'Default CD 2021' : 'etat_cd21_tab'},
                  'nom_opa' : {'OPAL and AF, interp. of Y. Lebreton' : 'opa_yveline',
                               'Improved Kramers (GONG)' : 'opa_gong',
                               'OPAL for GN 1993, Kurucz' : 'opa_int_zsx',
                               'OPAL with correction for Z>0.1 for C, O' : 'opa_opalCO',
                               'OPAL2 with correction for Z>0.1 for C, N, O' : 'opa_opal2_cno',
                               'OPAL2 with correction for Z>0.1 for C, O' : 'opa_opal2_co',
                               'OPAL and AF, interp. of Y. Lebreton + TOPS (D. Cordier)' : 'opa_yveline_daniel',
                               'OPAL and AF, interp. of Y. Lebreton with smoothing' : 'opa_yveline_lisse',
                               'OPAL, interp. of CLES' : 'opa_cles',
                               'OP (monochromatic)' : 'opa_mono_OP',
                               'OPLIB (monochromatic)' : 'opa_mono_OPLIB',
                               'OPAL GN93 (smooth)' : 'opa_smooth',
                               'OPAL and AF, interp. of G. Houdek' : 'opa_houdek9',
                               'Simple Opacity Model' : 'opa_simple',
                               'Freedman et al.' : 'opa_freedmanFIT_ro'},
                  'nom_opa_cond' : {"Pothekin 1999 (table of 2006)" : 'cond_yveline',
                                    "Pothekin 1999 (table of 2021)" : 'cond_yveline21',
                                    'Iben (1975)' : 'cond_iben75'},
                  'nom_nuc' : {'PP+CNO, 9 elems, H2, Li7 and Be7 in eq., 1<T6<40' : 'ppcno9',
                               'PMS, epsilon proportional to T' : 'iben',
                               'Simplified PP' : 'pp1',
                               'Deuterium only' : 'deuterium',
                               'PP with H, He3, He4 with H2, Li7 and Be7 in eq., 1<T6<80' : 'pp3',
                               'PP+CNO, 9 elems + Fe, H2, Li7 and Be7 in eq., 1<T6<80' : 'ppcno9Fe',
                               'PP+CNO, 9 elems + Ne, Na, Mg, Al, Si, P, S, Ca, Fe in eq., 1<T6<80' : 'ppcno9grad',
                               'PP+CNO, 10 elems, H2 and Be7 in eq., 0.5<T6<80' : 'ppcno10',
                               'PP+CNO, 10 elems + Fe, H2 and Be7 in eq., 0.5<T6<80' : 'ppcno10Fe',
                               'PP+CNO, 10 elems + K, H2 and Be7 in eq., 0.5<T6<80' : 'ppcno10K',
                               'PP+CNO, 10 elems + Li6, Be9, B11 and Fe, H2 and Be7 in eq., 0.5<T6<80' : 'ppcno10BeBFe',
                               'PP+CNO, 11 elems, Be7 in eq., 0.5<T6<80' : 'ppcno11',
                               'PP+CNO, 12 elems, 1<T6<80' : 'ppcno12',
                               'PP+CNO, 12 elems + Be9, 0.5<T6<80' : 'ppcno12Be',
                               'PP+CNO, 12 elems + Li6, Be9, B11 and Fe, 0.5<T6<80' : 'ppcno12BeBFe',
                               'PP+CNO, 12 elems + Li6, 0.5<T6<80' : 'ppcno12Li',
                               'PP+CNO, 13 elems + Ne, Na, Mg, Al, Si, P, S, Ca, Fe, Ni in eq., 1<T6<80' : 'ppcno13gradNi',
                               u'PP+CNO+3\u03B1, 12 elems (Ne22), H2, Li7 and Be7 in eq., 1<T6<800' : 'ppcno3a12Ne',
                               u'PP+CNO+3\u03B1, 9 elems, H2, Li7 and Be7 in eq., 1<T6<800' : 'ppcno3a9',
                               u'PP+CNO+3\u03B1+C+O, 17 elems, H2, Li7 and Be7 in eq., 1<T6<3000' : 'ppcno3aco',
                               u'PP+CNO+3\u03B1+C+O, 23 elems, H2, Li7 and Be7 in eq., 0.5<T6<3000' : 'nrj31_YL'},
                  'nom_nuc_cpl' : {'NACRE + LUNA' : 'NACRE_LUNA',
                                   'Caughlan & Fowler (1988)' : 'Cau-Fow',
                                   'Adelberger et al. (1998)' : 'Adelb',
                                   'Adelberger et al. (1998) + LUNA' : 'Adelb_LUNA',
                                   'NACRE' : 'NACRE',
                                   'NACRE II' : 'NACRE2',
                                   'Reaclib' : 'Reaclib'},
                  'nom_neutrinos' : {'Haft et al. (1994) + Weigert (1966)' : 'neutrinos_h94w66',
                                     'Itoh (1996)' : 'neutrinos_itoh96',
                                     'None' : 'neutrinos0'},
                  'nom_atm' : {u'Restitution of the atmosphere from a T(\u03C4) law' : 'lim_atm',
                               'Simplified restitution of the atmosphere (GONG1)' : 'lim_gong1',
                               'Simplified restitution of the atmosphere (one shell)' : 'lim_tau1'},
                  'nom_tdetau' : {u'Fully radiative T(\u03C4), Hopf' : 'hopf',
                                  u'Fully radiative T(\u03C4), Eddington' : 'edding',
                                  u'Fully radiative T(\u03C4), Eddington w. irradiation' : 'eddingG10',
                                  'hsra' : 'hsra',
                                  'Kurucz solar atmosphere for Teff=5750 K' : 'k5750',
                                  'Kurucz solar atmosphere for Teff=5777 K' : 'k5777',
                                  'MARCS' : 'MARCS',
                                  'Kurucz atmosphere, Roger' : 'roger',
                                  'Kurucz' : 'kuruc',
                                  'Roger + YL' : 'rogYL',
                                  'Krishna Swamy' : 'krisw',
                                  'Vernazza (Model C)' : 'verna',
                                  'Trampedach et al. (2014)' : 'tramp',
                                  'Ball+21 fit of Trampedach+14 Sun' : 'ball_s',
                                  'Ball+21 fit of Trampedach+14 Grid' : 'ball_g'},
                  'nom_mag' : {'No magnetic field' : 'mag_0',
                               'Influence through fixed magnetic convective-inhibition parameter' : 'mag_inhib',
                               'Constant magnetic field' : 'mag_f',
                               'Magnetic convective-inhibition parameter saturating for r > Rmax' : 'mag_inhib_rmax',
                               'Magnetic convective-inhibition parameter saturating for Br > Brmax' : 'mag_inhib_brmax'},
                  'xmax_unit':{
                        'Fraction of the solar radius': 'rsun',
                        'Fraction of the star radius': 'rstar',
                        'log10( rho )': 'lrho'}
                  }

FPARAMS_DICT = {
    'osc_code_name' : {
        'ADIPLS' : 'adipls',
        'ACOR'   : 'acor'},
    'modes' : {
        'p-modes' : 'p',
        'g-modes' : 'g',
        'all' : 'a'},
    'reduce' : {
        'None': 'none',
        'Gaussian envelope for non-radial modes': 'gauss_norad'},
    'mdintg' : {
        'Shooting method' : 1,
        'Relaxation method' : 2,
        '4th-order shooting method' : 5},
    'iekinr' : {
        'With surface vertical displacement' : 0,
        'With total photospheric displacement' : 1},
    'nsig' : {
        u'Frequency\u00B2' : 1,
        'Frequency' : 2,
        '1/Frequency' : 3,
        '1/Frequency at low freq, frequency at high freq' : 4,
        'Red giants' : 10},
    'remesh' : {
        'Do not re-mesh' : 'n',
        'Adapt for solar p-modes' : 'p',
        'Adapt for solar g-modes' : 'g',
        'Just increase mesh density' : 'a',
        'Adapt for red giants' : 'r',
        'Customize' : 'o'},
    'istsbc' : {
        'Simple surface condition' : 0,
        'Exponentially decaying solution' : 1,
        u'Use \u03B4r = 0 at surface' : 9},
    'fctsbc' : {
        u'\u03B4p = 0' : 0,
        "p' = 0" : 1,
        'Custom' : 2},
    'input_file' : {
        'Compute frequencies from 1D output files' : 'd1',
        'Compute frequencies from 2D output files' : 'd2',
        'Compute frequencies from 1D and 2D output files' : 'd12'},
    'mod_flag' : {
        'Model name' : 'modnam',
        'Custom' : 'custom'},
    'acor_star_type' : {
        'Low mass star' : 1,
        'Massive star' : 2,
        '\u03b3 Dor' : 3},
    'acor_rot_prof' : {
        'Uniform' : 1,
        'Differential' : 2,
        'Bizonal' : 3,
        'Bizonal with smooth gradient' : 4,
        'Latitudinal diff: Om = Om_eq - DOm np.cos^2 theta' : 5},
    'acor_freq_target' : {
        'Scanning' : 1,
        'File' : 2,
        '\u0394\u03a0' : 3,
        '\u0394\u03bd' : 4},
    'acor_freq_units' : {
        'microHz' : 1,
        'Omk' : 2},
    'acor_freq_file_type' : {
        'ADIPLS' : 1,
        'ACOR' : 2,
        'Marie-Jo' : 3,
        'SebS' : 4},
    'acor_params_select' : {
        'lmax and m' : 1,
        'Number of sph. harm., m and parity' : 2}
    }

NUM_DICT = {
   'ltau_rep' : {
        'Linear': 'lin',
        'Cubic': 'cub'
        },
   'ec_s' : {
        'Bottom' : 'bottom',
        'Average': 'average'
        },
   'lim_zc' : {
        'Gabriel (2014)':'gabriel14',
        'Legacy'        :'original',
        'Planets'       :'planet'
        },
}

EXTRA_DICT = {
    'ovs_type': {
        'Adiabatic gradient': 'adia',
        'Radiative gradient': 'radia'
        },
    'ovi_type': {
        'Adiabatic gradient': 'adia',
        'Radiative gradient': 'radia',
        "Follows JCD+2011's prescription": 'jcd11',
        "Follows Deal+2013's prescription": 'md23'
        },
    'struct2d': {
        "None" : '',
        "CLES": 'cles',
        "Cesam2k20": 'cesam2k20'
        },
    'source_cts': {
        "IAU2015_CODATA2018": 'IAU2015_CODATA2018',
        'Lide1994': 'Lide1994'
        },
    'ec_grid': {
        'Reduced CIFIST grid [discouraged]' : 'cif_red',
        'CIFIST + Mdwarfs grid' : 'cif_mdw',
        'CIFIST + updated Mdwarfs grid': 'cif_mdw2',
        'CIFIST grid' : 'cifist',
        'CIFIST grid + RGB + [Fe/H] = +0.5' : 'extended',
        'Extended CIFIST grid + hot stars' : 'cifist_hot',
        'Coefficient from original papers' : 'original'
        },
    'pertw_profile': {
        'Angular momentum loss as computed' : 'constant',
        'Angular momentum loss smoothly set to 0 as convective turnover time goes to 0': 'logistic'
        },
    'ts_average': {
        'Dirac function': 'dirac',
        'Positive piecewise parabola': 'parabola'
        },
    'ts_average_length': {
        'None': 'none',
        'Tayler scale': 'lts',
        'Pressure scale height': 'hp',
        'Omega scale height': 'hw'
        },
    'rota_out' : {
        'None' : 'none',
        'End of evolution' : 'end_evol',
        'All time steps' : 'all_mod',
        'Each iteration' : 'all_iter'},
    'ec_extrap' : {
        "Keep previous value": 'keep',
        "Keep previous value, but smoothly swith to EC when reentering domain.": 'keep_smooth',
        "Eval fitting function with current values of (teff,logg,z)": 'eval',
        "Project (teff,logg,z) onto convex hull, and evaluate fitting function there": 'proj_coord',
        "Project (teff,logg,z) onto convex hull, and extrapolate entropy there)": 'proj_s',
        "Extrapolate at fixed Teff": 'tfixed',
        "Fall back on Ludwig99": 'ludwig',
        "Fall back on Sonoi's prescription" : 'sonoi'}
}

