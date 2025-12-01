DEFAULT_NUM = {'ini0' :4, 'l0' :0, 'm_qs' :2, 'm_ch' :2, 'm_vth' :3, 'm_rot' :3, 'm_tds' :2,
    'm_ptm' :2, 'n_atm' :75, 'n_min' :150, 'n_min_post' :150, 'ordre' :2, 'yld_rep' :1, 'osc2d_step' :5,
    'osc_step' :1, 'm_qs2d' :3, 'm_rot2d' :3, 'm_legendre' :2, 'n_theta' :100,
    'iter_qs' : [0,0,0,0,0,0,0], 'ctel' : 0.0, 'ctep' : -1.0, 'ctem' : 15.0, 'cter' : 0.0,
    'ctet' : -1.0, 'coeff_d_grav_post' : 1.0, 'd_grav' : 5.0e3, 'd_grav_int' : 5.0e3, 'dlntc' : 0.07, 'drhoc' : -1.0,
    'dteff' : -1.0, 'dlogg' : -1.0, 'dalpha' : 0.03, 'dlum' : -1.0, 'dw' : -1.0,
    'drhoca' : -1.0, 'dsenv' : 5.0e-4, 'evolved' : 0.1, 'dn_fixe' : 0.05, 'dpsi' : 0.05,
    'dt0' : 1.0, 'dtmax' : 200.0, 'dtmin' : 1.0e-15, 'fmin_abon' : 0.05, 'loc_zc' : 1.0e-3,
    'precix' : 1.0e-3, 'precit' : 0.15, 'psi0' : 0.08, 'ro_test' : 0.1, 'q0' : 0.05,
    'yld_dx' : -1.0, 'yld_rac' : 1.0, 'tau_min' : 1.0e-4, 'dt_m' : -2.5, 'x_tams' : 0.01,
    'y_agb' : 0.001, 'dx_tams' : 1.0e-4, 'mcz_ext_min' : 1.0e-10, 'exact_stop' : True,
    'en_masse' : True, 'kipp' : False, 'lisse' : False, 'mu_saha' : True, 'mvt_dis' : True,
    'new_bv' : False, 'yld_hr' : False, 'yld_l23' : False, 'yld_fgong' : False, 'd2' : True,
    'general' : False, 'simple' : True, 'extended_osc2d' : False, 'vdiff_out' : False,
    'grad_out' : False, 'th_out' : False, 'eos_inplace' : False, 'g_modes_ab' : False, 'out_pms' : True,
    'clean' : True, 'ltau_rep' : 'Cubic', 'ec_s' : 'Average', 'lim_zc' : 'Gabriel (2014)'}


DEFAULT_EXTRA = {'rho_ext' : 3.55e-9, 'easy_krisw' : False, 'osm_p_pertw' : 1.0,
    'npt_central_cz' : 0, 'nmore_shells' : 1, 'nu_v_fact' : 1.0, 'nu_v_add' : 0.0,
    'nu_h_fact' : 1.0, 'nu_h_add' : 0.0, 'deff_fact' : 1.0, 'deff_add' : 0.0, 'dv_fact' : 1.0,
    'dv_add' : 0.0, 're_c' : 10.0, 'ri_c' : 1.0/6.0, 'omega_sun' : 2.86e-6, 'fast_rot' : -1.0,
    'opt_tau_turb' : 0, 'ell_max' : 10, 'ell_min' : 1, 'ndelta' : 5, 'eta_sun' : 4.36e-6,
    't0_corsaro' : 601.0, 'eta_0' : -1.0, 'ratio_work_max' : 1.0e5, 'mode_damping' : True,
    'n_fake_ce' : 1, 'd_conv' : 1.0e13, 'mass_tol_cz' : 0.10, 'conv_block' : False, 'no_tot_rad_post' : False,
    'r_transition' : 1.009e-2, 'width_transition' : 1.211e-3, 'omega_c' : 6.28e-6,
    'omega_s' : 3.45e-7, 'struct2d' : 'None', 'struct2d_file' : '', 't_mix' : 0.0,
    'm_mix' : 0.0, 'expo_mix' : 3.0, 'cst_mix' : 0.0, 'd_expo_2' : 0.0, 'l_expo' : 0.0,
    'l_expo_2' : 0.0, 'eta_grad' : 1.0, 'goop' : 0.0, 'eta_vdif' : 1.0, 'no_diff_post' : False, 'b_strength' : 0.0,
    'errabop' : 0.1, 't_op_sup' : 6.5e6, 't_op_inf' : 1.0e4, 't_oplib_sup' : 1.16e9, 't_oplib_inf' : 1.0e4,
    'use_findne' : True, 'teff_zams_sun' : 5700,
    'teff_zams_2' : 7000, 't0_sun' : 10**6.43, 't0_2' : 10**5.90, 't_ov_inf' : 5.0e5,
    'alpha_jcd11' : 0.3, 'beta_jcd11' : 0.2, 'claret18_shift' : 0.0,
    'ovs_type' : 'Radiative gradient', 'ovi_type' : 'Adiabatic gradient',
    'start_ovshti_zams' : True, 'lt_compton_l' : 8.2, 'lt_compton_h' : 8.7,
    't_compton' : 7e7, 'ds' : 0.0, 'ec_grid' : 'CIFIST grid', 'ec_extrap' : 'Keep previous value',
    'alpha_autoguess' : False, 'ec_corr' : True, 'alpha_m10' : -1.0, 'alpha_m15' : -1.0,
    'a_m10_z005' : -1.0, 'a_m15_z005' : -1.0, 'a_m10_z02' : -1.0, 'a_m15_z02' : -1.0,
    'mass_loss_eta' : 0.5, 'mass_loss_zeta' : 0.5, 'log_lum_min_jvh' : 2.5,
    'log_lum_max_jvh' : 6.7, 'v_inf_on_v_esc_hot' : 2.6, 'v_inf_on_v_esc_cold' : 1.3,
    'eta_r_rgb_1' : 0.5, 'eta_r_rgb_2' : 1.0, 'p_power_matt15' : 2.0, 'tau_conv_sun' : 35.0,
    'pertw_profile' : 'Angular momentum loss smoothly set to 0 as convective turnover time goes to 0',
    'screen' : True, 'eta_nuc' : 1.0, 'r_h2h' : 1.0, 'r_14n15o' : 1.0, 'inj_a' : -1.0, 'inj_mu' : 0.0001,
    'ts_dt_max' : 0.05, 'lsmooth_fact' : 1.0, 'lsmooth_min' : 0.0, 'alpha_fuller' : 1.0,
    'ts_average' : 'Positive piecewise parabola', 'ts_average_length' : 'None',
    'ts_weak_b' : False, 'ts_hysteresis' : False, 'ts_patch_regime' : True,
    'ts_qlim' : False, 'ts_extrapolate' : False, 'source_cts' : 'IAU2015_CODATA2018',
    'debug' : False, 'convection_out' : False, 'hse_out' : False, 'rota_out' : 'None', 'ts_out' : False,
    'lim_zc_out' : False, 'full_hr' : True, 'writes_tau_conv' : True, 'writes_nuc_neutr' : True,
    'writes_l_nuc' : True, 'writes_mdot' : True, 'writes_conv' : True, 'writes_2d' : True, 'writes_burn_zones' : True,
    'n_contours_burn' : 4, 'burn_contours' : [1.0, 1.0e1, 1.0e3, 1.0e5,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0],
    'plato_track' : False}


PRECISIONS = {
    'er':{
        'm_qs' : 1, 'm_rot' : 1, 'en_masse' : False, 'ctem' : 0.0, 'cter' : 15.0},
    'sp':{
        'precix' : 1.0e-4, 'precit' : 0.05, 'psi0' : 0.06, 'd_grav' : 0.5, 'd_grav_int' : 5.0,
        'loc_zc' : 1.0e-4, 'dtmax' : 50.0, 'ini0' : 5, 'n_atm' : 100, 'fmin_abon' : 0.01,
        'dlntc' : 0.05},
    'sr':{
        'm_qs' : 1, 'm_rot' : 1, 'en_masse' : False, 'ctem' : 0.0, 'cter' : 15.0,
        'precix' : 1.0e-4, 'precit' : 0.05, 'psi0' : 0.06, 'd_grav' : 0.5, 'd_grav_int' : 5.0,
        'loc_zc' : 1.0e-4, 'dtmax' : 50.0, 'ini0' : 5, 'n_atm' : 100, 'fmin_abon' : 0.01,
        'dlntc' : 0.05},
    'av':{
        'm_qs' : 1, 'kipp' : True, 'mu_saha' : False, 'exact_stop' : False, 'lisse' : True,
        'q0' : 0.0
        },
    'co':{
        'ini0' : 6, 'precix' : 1.0e-4, 'dtmax' : 50.0, 'm_ch' : 3, 'precit' : 0.05,
        'psi0' : 0.06, 'd_grav' : 0.5, 'd_grav_int' : 2.5, 'loc_zc' : 1.0e-4, 'n_atm' : 100,
        'q0' : 0.0, 'fmin_abon' : 0.01, 'dlntc' : 0.05
        },
    'pl':{
        'ini0' : 6, 'precix' : 1.0e-5, 'precit' : 0.01, 'dtmax' : 10.0, 'ordre' : 2,
        'm_ch' : 2, 'psi0' : 0.02, 'd_grav' : 0.2, 'd_grav_int' : 1.0, 'loc_zc' : 1.0e-5,
        'n_atm' : 200, 'n_min' : 3000, 'n_max' : 4000, 'q0' : 0.0, 'l0':0, 'fmin_abon' : 0.01,
        'dlntc' : 0.05, 'drhoc' : 0.05, 'dteff' : 0.005, 'dlogg' : 0.1, 'dlum' : 0.03,
        'dalpha' : 0.01, 'yld_dx' : 0.01, 'dt_m' : -3.0, 'mvt_dis' : True
        },
    'mj':{
        'mu_saha' : True, 'ini0' : 6, 'precix' : 1.0e-4, 'dtmax' : 50.0, 'm_ch' : 3,
        'precit' : 0.05, 'psi0' : 0.03, 'd_grav' : 0.5, 'loc_zc' : 1.0e-4, 'n_atm' : 100,
        'q0' : 0.01, 'l0' : 5, 'lisse' : False, 'fmin_abon' : 0.01, 'dlntc' : 0.05
        },
    'lm':{
        'precit' : 0.2, 'dtmax' : 300.0, 'ini0' : 5, 'kipp' : True, 'dn_fixe' : 0.1,
        'mu_saha' : False, 'q0' : 0.0, 'l0' : 0
        },
    'sa':{
        'precix' : 1.0e-5, 'precit' : 0.02, 'psi0' : 0.06, 'd_grav' : 0.5, 'loc_zc' : 1.0e-5,
        'dtmax' : 50.0, 'ini0' : 5, 'n_atm' : 100, 'lisse' : False, 'q0' : 0.01, 'l0' : 5,
        'fmin_abon' : 0.01, 'new_bv' : False, 'dlntc' : 0.05
        },
    'hp':{
        'm_qs' : 3, 'm_ch' : 4, 'ordre' : 4, 'precix' : 1.0e-6, 'precit' : 0.01,
        'ro_test' : 0.05, 'psi0' : 0.05, 'd_grav' : 0.05, 'loc_zc' : 1.0e-6, 'dtmax' : 50.0,
        'ini0' : 6, 'n_atm' : 150, 'n_max' : 2500, 'q0' : 0.01, 'l0' : 7, 'fmin_abon' : 0.01,
        'dlntc' : 0.03
        },
    'np':{
        'm_qs' : 1, 'm_ch' : 2, 'ordre' : 1, 'precix' : 5.0e-3, 'precit' : 0.3, 'psi0' : 0.1,
        'loc_zc' : 5.0e-3, 'dtmax' : 300.0, 'ini0' : 2, 'd_grav' : 0.5, 'd_grav_int' : 2.5,
        'n_atm' : 50, 'kipp' : True, 'mvt_dis' : False, 'mu_saha' : False, 'exact_stop' : False,
        'q0' : 0.0, 'l0' : 0, 'dlntc' : 0.1
        },
    'nr':{
        'm_qs' : 1, 'ordre' : 1, 'precix' : 5.0e-3, 'precit' : 0.3, 'psi0' : 0.1,
        'en_masse' : False, 'ctem' : 0.0, 'cter' : 15.0, 'loc_zc' : 5.0e-3,
        'dtmax' : 300.0, 'ini0' : 3, 'n_atm' : 50, 'mvt_dis' : False, 'mu_saha' : False,
        'exact_stop' : False, 'kipp' : True, 'q0' : 0.0, 'l0' : 0, 'dlntc' : 0.1
        },
    'mx':{
        'precix' : 1.0e-5, 'precit' : 0.02, 'psi0' : 1.0e-8, 'd_grav' : 0.5, 'loc_zc' : 1.0e-5,
        'dtmax' : 50.0, 'ini0' : 5, 'n_atm' : 100, 'lisse' : False, 'q0' : 0.01, 'l0' : 4,
        'fmin_abon' : 0.01, 'dlntc' : 0.05
        }
}
