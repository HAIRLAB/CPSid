function xdot=AP_simulation(time,x)
% End Matlab code

% Start Octave code
%function xdot=f(x,t)
% End Octave code

% Compartment: id = Compartment, name = Compartment, constant
	compartment_Compartment=1.0;
% Parameter:   id =  V_membrane, name = V
% Parameter:   id =  R, name = R
	global_par_R=8.3143;
% Parameter:   id =  T, name = T
	global_par_T=310.0;
% Parameter:   id =  F, name = F
	global_par_F=96.4867;
% Parameter:   id =  Cm, name = Cm
	global_par_Cm=100.0;
% Parameter:   id =  i_st, name = i_st
% Parameter:   id =  stim_start, name = stim_start
	global_par_stim_start=100.0;
% Parameter:   id =  stim_end, name = stim_end
	global_par_stim_end=50000.0;
% Parameter:   id =  stim_period, name = stim_period
	global_par_stim_period=1000.0;
% Parameter:   id =  stim_duration, name = stim_duration
	global_par_stim_duration=2.0;
% Parameter:   id =  stim_amplitude, name = stim_amplitude
	global_par_stim_amplitude=-2000.0;
% Parameter:   id =  i_Na, name = i_Na
% Parameter:   id =  E_Na, name = E_Na
% Parameter:   id =  g_Na, name = g_Na
	global_par_g_Na=7.8;
% Parameter:   id =  m, name = m
% Parameter:   id =  alpha_m, name = alpha_m
% Parameter:   id =  beta_m, name = beta_m
% Parameter:   id =  m_inf, name = m_inf
% Parameter:   id =  tau_m, name = tau_m
% Parameter:   id =  h, name = h
% Parameter:   id =  alpha_h, name = alpha_h
% Parameter:   id =  beta_h, name = beta_h
% Parameter:   id =  h_inf, name = h_inf
% Parameter:   id =  tau_h, name = tau_h
% Parameter:   id =  j, name = j
% Parameter:   id =  alpha_j, name = alpha_j
% Parameter:   id =  beta_j, name = beta_j
% Parameter:   id =  j_inf, name = j_inf
% Parameter:   id =  tau_j, name = tau_j
% Parameter:   id =  i_K1, name = i_K1
% Parameter:   id =  E_K, name = E_K
% Parameter:   id =  g_K1, name = g_K1
	global_par_g_K1=0.09;
% Parameter:   id =  i_to, name = i_to
% Parameter:   id =  K_Q10, name = K_Q10
	global_par_K_Q10=3.0;
% Parameter:   id =  g_to, name = g_to
	global_par_g_to=0.1652;
% Parameter:   id =  oa, name = oa
% Parameter:   id =  alpha_oa, name = alpha_oa
% Parameter:   id =  beta_oa, name = beta_oa
% Parameter:   id =  tau_oa, name = tau_oa
% Parameter:   id =  oa_infinity, name = oa_infinity
% Parameter:   id =  oi, name = oi
% Parameter:   id =  alpha_oi, name = alpha_oi
% Parameter:   id =  beta_oi, name = beta_oi
% Parameter:   id =  tau_oi, name = tau_oi
% Parameter:   id =  oi_infinity, name = oi_infinity
% Parameter:   id =  i_Kur, name = i_Kur
% Parameter:   id =  g_Kur, name = g_Kur
% Parameter:   id =  ua, name = ua
% Parameter:   id =  alpha_ua, name = alpha_ua
% Parameter:   id =  beta_ua, name = beta_ua
% Parameter:   id =  tau_ua, name = tau_ua
% Parameter:   id =  ua_infinity, name = ua_infinity
% Parameter:   id =  ui, name = ui
% Parameter:   id =  alpha_ui, name = alpha_ui
% Parameter:   id =  beta_ui, name = beta_ui
% Parameter:   id =  tau_ui, name = tau_ui
% Parameter:   id =  ui_infinity, name = ui_infinity
% Parameter:   id =  K_Q10_ultrarapid_delayed_rectifier_K_current_ui_gate, name = K_Q10
	global_par_K_Q10_ultrarapid_delayed_rectifier_K_current_ui_gate=NaN;
% Parameter:   id =  i_Kr, name = i_Kr
% Parameter:   id =  g_Kr, name = g_Kr
	global_par_g_Kr=0.029411765;
% Parameter:   id =  xr, name = xr
% Parameter:   id =  alpha_xr, name = alpha_xr
% Parameter:   id =  beta_xr, name = beta_xr
% Parameter:   id =  tau_xr, name = tau_xr
% Parameter:   id =  xr_infinity, name = xr_infinity
% Parameter:   id =  i_Ks, name = i_Ks
% Parameter:   id =  g_Ks, name = g_Ks
	global_par_g_Ks=0.12941176;
% Parameter:   id =  xs, name = xs
% Parameter:   id =  alpha_xs, name = alpha_xs
% Parameter:   id =  beta_xs, name = beta_xs
% Parameter:   id =  tau_xs, name = tau_xs
% Parameter:   id =  xs_infinity, name = xs_infinity
% Parameter:   id =  i_Ca_L, name = i_Ca_L
% Parameter:   id =  g_Ca_L, name = g_Ca_L
	global_par_g_Ca_L=0.12375;
% Parameter:   id =  d, name = d
% Parameter:   id =  d_infinity, name = d_infinity
% Parameter:   id =  tau_d, name = tau_d
% Parameter:   id =  f, name = f
% Parameter:   id =  f_infinity, name = f_infinity
% Parameter:   id =  tau_f, name = tau_f
% Parameter:   id =  f_Ca, name = f_Ca
% Parameter:   id =  f_Ca_infinity, name = f_Ca_infinity
% Parameter:   id =  tau_f_Ca, name = tau_f_Ca
% Parameter:   id =  i_NaK, name = i_NaK
% Parameter:   id =  Km_Na_i, name = Km_Na_i
	global_par_Km_Na_i=10.0;
% Parameter:   id =  Km_K_o, name = Km_K_o
	global_par_Km_K_o=1.5;
% Parameter:   id =  i_NaK_max, name = i_NaK_max
	global_par_i_NaK_max=0.59933874;
% Parameter:   id =  f_NaK, name = f_NaK
% Parameter:   id =  sigma, name = sigma
% Parameter:   id =  i_B_Na, name = i_B_Na
% Parameter:   id =  i_B_Ca, name = i_B_Ca
% Parameter:   id =  i_B_K, name = i_B_K
% Parameter:   id =  g_B_Na, name = g_B_Na
	global_par_g_B_Na=6.744375E-4;
% Parameter:   id =  g_B_Ca, name = g_B_Ca
	global_par_g_B_Ca=0.001131;
% Parameter:   id =  g_B_K, name = g_B_K
	global_par_g_B_K=0.0;
% Parameter:   id =  E_Ca, name = E_Ca
% Parameter:   id =  i_NaCa, name = i_NaCa
% Parameter:   id =  I_NaCa_max, name = I_NaCa_max
	global_par_I_NaCa_max=1600.0;
% Parameter:   id =  K_mNa, name = K_mNa
	global_par_K_mNa=87.5;
% Parameter:   id =  K_mCa, name = K_mCa
	global_par_K_mCa=1.38;
% Parameter:   id =  K_sat, name = K_sat
	global_par_K_sat=0.1;
% Parameter:   id =  gamma, name = gamma
	global_par_gamma=0.35;
% Parameter:   id =  i_CaP, name = i_CaP
% Parameter:   id =  i_CaP_max, name = i_CaP_max
	global_par_i_CaP_max=0.275;
% Parameter:   id =  i_rel, name = i_rel
% Parameter:   id =  Fn, name = Fn
% Parameter:   id =  K_rel, name = K_rel
	global_par_K_rel=30.0;
% Parameter:   id =  u, name = u
% Parameter:   id =  tau_u, name = tau_u
% Parameter:   id =  u_infinity, name = u_infinity
% Parameter:   id =  v, name = v
% Parameter:   id =  tau_v, name = tau_v
% Parameter:   id =  v_infinity, name = v_infinity
% Parameter:   id =  w, name = w
% Parameter:   id =  tau_w, name = tau_w
% Parameter:   id =  w_infinity, name = w_infinity
% Parameter:   id =  i_tr, name = i_tr
% Parameter:   id =  tau_tr, name = tau_tr
	global_par_tau_tr=180.0;
% Parameter:   id =  I_up_max, name = I_up_max
	global_par_I_up_max=0.005;
% Parameter:   id =  i_up, name = i_up
% Parameter:   id =  K_up, name = K_up
	global_par_K_up=9.2E-4;
% Parameter:   id =  i_up_leak, name = i_up_leak
% Parameter:   id =  Ca_up_max, name = Ca_up_max
	global_par_Ca_up_max=15.0;
% Parameter:   id =  CMDN_max, name = CMDN_max
	global_par_CMDN_max=0.05;
% Parameter:   id =  TRPN_max, name = TRPN_max
	global_par_TRPN_max=0.07;
% Parameter:   id =  CSQN_max, name = CSQN_max
	global_par_CSQN_max=10.0;
% Parameter:   id =  Km_CMDN, name = Km_CMDN
	global_par_Km_CMDN=0.00238;
% Parameter:   id =  Km_TRPN, name = Km_TRPN
	global_par_Km_TRPN=5.0E-4;
% Parameter:   id =  Km_CSQN, name = Km_CSQN
	global_par_Km_CSQN=0.8;
% Parameter:   id =  Ca_CMDN, name = Ca_CMDN
% Parameter:   id =  Ca_TRPN, name = Ca_TRPN
% Parameter:   id =  Ca_CSQN, name = Ca_CSQN
% Parameter:   id =  Na_i, name = Na_i
% Parameter:   id =  Ca_i, name = Ca_i
% Parameter:   id =  K_i, name = K_i
% Parameter:   id =  Ca_rel, name = Ca_rel
% Parameter:   id =  Ca_up, name = Ca_up
% Parameter:   id =  V_cell, name = V_cell
	global_par_V_cell=20100.0;
% Parameter:   id =  V_i, name = V_i
% Parameter:   id =  V_rel, name = V_rel
% Parameter:   id =  V_up, name = V_up
% Parameter:   id =  B1, name = B1
% Parameter:   id =  B2, name = B2
% Parameter:   id =  Na_o, name = Na_o
	global_par_Na_o=140.0;
% Parameter:   id =  Ca_o, name = Ca_o
	global_par_Ca_o=1.8;
% Parameter:   id =  K_o, name = K_o
	global_par_K_o=5.4;
% rateRule: variable = V_membrane
global_par_V_membrane = x(1);
% rateRule: variable = m
global_par_m = x(2);
% rateRule: variable = h
global_par_h = x(3);
% rateRule: variable = j
global_par_j = x(4);
% rateRule: variable = oa
global_par_oa = x(5);
% rateRule: variable = oi
global_par_oi = x(6);
% rateRule: variable = ua
global_par_ua = x(7);
% rateRule: variable = ui
global_par_ui = x(8);
% rateRule: variable = xr
global_par_xr = x(9);
% rateRule: variable = xs
global_par_xs = x(10);
% rateRule: variable = d
global_par_d = x(11);
% rateRule: variable = f
global_par_f = x(12);
% rateRule: variable = f_Ca
global_par_f_Ca = x(13);
% rateRule: variable = u
global_par_u = x(14);
% rateRule: variable = v
global_par_v = x(15);
% rateRule: variable = w
global_par_w = x(16);
% rateRule: variable = Na_i
global_par_Na_i = x(17);
% rateRule: variable = K_i
global_par_K_i = x(18);
% rateRule: variable = Ca_i
global_par_Ca_i = x(19);
% rateRule: variable = Ca_up
global_par_Ca_up = x(20);
% rateRule: variable = Ca_rel
global_par_Ca_rel = x(21);
% assignmentRule: variable = i_st
% time = 200;
	global_par_i_st=piecewise(global_par_stim_amplitude, (time >= global_par_stim_start) && (time <= global_par_stim_end) && ((time-global_par_stim_start-floor((time-global_par_stim_start)/global_par_stim_period)*global_par_stim_period) <= global_par_stim_duration), 0);
% assignmentRule: variable = i_Na
	global_par_E_Na=global_par_R*global_par_T/global_par_F*log(global_par_Na_o/global_par_Na_i);

	global_par_i_Na=global_par_Cm*global_par_g_Na*global_par_m^3*global_par_h*global_par_j*(global_par_V_membrane-global_par_E_Na);
% assignmentRule: variable = E_Na
% assignmentRule: variable = alpha_m
	global_par_alpha_m=piecewise(3.2, (47.13+global_par_V_membrane) == 0, 0.32*(global_par_V_membrane+47.13)/(1-exp((-0.1)*(global_par_V_membrane+47.13))));
% assignmentRule: variable = beta_m
	global_par_beta_m=0.08*exp((-global_par_V_membrane)/11);
% assignmentRule: variable = m_inf
	global_par_m_inf=global_par_alpha_m/(global_par_alpha_m+global_par_beta_m);
% assignmentRule: variable = tau_m
	global_par_tau_m=1/(global_par_alpha_m+global_par_beta_m);
% assignmentRule: variable = alpha_h
	global_par_alpha_h=piecewise(0.135*exp((global_par_V_membrane+80)/(-6.8)), (40+global_par_V_membrane) < 0, 0);
% assignmentRule: variable = beta_h
	global_par_beta_h=piecewise(3.56*exp(0.079*global_par_V_membrane)+310000*exp(0.35*global_par_V_membrane), (40+global_par_V_membrane) < 0, 1/(0.13*(1+exp((global_par_V_membrane+10.66)/(-11.1)))));
% assignmentRule: variable = h_inf
	global_par_h_inf=global_par_alpha_h/(global_par_alpha_h+global_par_beta_h);
% assignmentRule: variable = tau_h
	global_par_tau_h=1/(global_par_alpha_h+global_par_beta_h);
% assignmentRule: variable = alpha_j
	global_par_alpha_j=piecewise(((-127140)*exp(0.2444*global_par_V_membrane)-3.474E-5*exp((-0.04391)*global_par_V_membrane))*(global_par_V_membrane+37.78)/(1+exp(0.311*(global_par_V_membrane+79.23))), (40+global_par_V_membrane) < 0, 0);
% assignmentRule: variable = beta_j
	global_par_beta_j=piecewise(0.1212*exp((-0.01052)*global_par_V_membrane)/(1+exp((-0.1378)*(global_par_V_membrane+40.14))), (40+global_par_V_membrane) < 0, 0.3*exp((-2.535E-7)*global_par_V_membrane)/(1+exp((-0.1)*(global_par_V_membrane+32))));
% assignmentRule: variable = j_inf
	global_par_j_inf=global_par_alpha_j/(global_par_alpha_j+global_par_beta_j);
% assignmentRule: variable = tau_j
	global_par_tau_j=1/(global_par_alpha_j+global_par_beta_j);
% assignmentRule: variable = E_K
	global_par_E_K=global_par_R*global_par_T/global_par_F*log(global_par_K_o/global_par_K_i);
% assignmentRule: variable = i_K1
	global_par_i_K1=global_par_Cm*global_par_g_K1*(global_par_V_membrane-global_par_E_K)/(1+exp(0.07*(global_par_V_membrane+80)));
% assignmentRule: variable = i_to
	global_par_i_to=global_par_Cm*global_par_g_to*global_par_oa^3*global_par_oi*(global_par_V_membrane-global_par_E_K);
% assignmentRule: variable = alpha_oa
	global_par_alpha_oa=0.65*(exp((global_par_V_membrane+10)/(-8.5))+exp((global_par_V_membrane-10-40)/(-59)))^(-1);
% assignmentRule: variable = beta_oa
	global_par_beta_oa=0.65*(2.5+exp((global_par_V_membrane+10+72)/17))^(-1);
% assignmentRule: variable = tau_oa
	global_par_tau_oa=1/((global_par_alpha_oa+global_par_beta_oa)*global_par_K_Q10);
% assignmentRule: variable = oa_infinity
	global_par_oa_infinity=(1+exp((global_par_V_membrane+10+10.47)/(-17.54)))^(-1);
% assignmentRule: variable = alpha_oi
	global_par_alpha_oi=(18.53+1*exp((global_par_V_membrane+10+103.7)/10.95))^(-1);
% assignmentRule: variable = beta_oi
	global_par_beta_oi=(35.56+1*exp((global_par_V_membrane+10-8.74)/(-7.44)))^(-1);
% assignmentRule: variable = tau_oi
	global_par_tau_oi=1/((global_par_alpha_oi+global_par_beta_oi)*global_par_K_Q10);
% assignmentRule: variable = oi_infinity
	global_par_oi_infinity=(1+exp((global_par_V_membrane+10+33.1)/5.3))^(-1);
% assignmentRule: variable = g_Kur
	global_par_g_Kur=0.005+0.05/(1+exp((global_par_V_membrane-15)/(-13)));
% assignmentRule: variable = i_Kur
	global_par_i_Kur=global_par_Cm*global_par_g_Kur*global_par_ua^3*global_par_ui*(global_par_V_membrane-global_par_E_K);
% assignmentRule: variable = alpha_ua
	global_par_alpha_ua=0.65*(exp((global_par_V_membrane-10)/(-8.5))+exp((global_par_V_membrane+10-40)/(-59)))^(-1);
% assignmentRule: variable = beta_ua
	global_par_beta_ua=0.65*(2.5+exp((global_par_V_membrane+10+72)/17))^(-1);
% assignmentRule: variable = tau_ua
	global_par_tau_ua=1/((global_par_alpha_ua+global_par_beta_ua)*global_par_K_Q10);
% assignmentRule: variable = ua_infinity
	global_par_ua_infinity=(1+exp((global_par_V_membrane+10+20.3)/(-9.6)))^(-1);
% assignmentRule: variable = alpha_ui
	global_par_alpha_ui=(21+1*exp((global_par_V_membrane+10-195)/(-28)))^(-1);
% assignmentRule: variable = beta_ui
	global_par_beta_ui=1/exp((global_par_V_membrane+10-168)/(-16));
% assignmentRule: variable = tau_ui
	global_par_tau_ui=1/((global_par_alpha_ui+global_par_beta_ui)*global_par_K_Q10);
% assignmentRule: variable = ui_infinity
	global_par_ui_infinity=(1+exp((global_par_V_membrane+10-109.45)/27.48))^(-1);
% assignmentRule: variable = i_Kr
	global_par_i_Kr=global_par_Cm*global_par_g_Kr*global_par_xr*(global_par_V_membrane-global_par_E_K)/(1+exp((global_par_V_membrane+15)/22.4));
% assignmentRule: variable = alpha_xr
	global_par_alpha_xr=piecewise(0.0015, abs(global_par_V_membrane+14.1) < 1E-10, 0.0003*(global_par_V_membrane+14.1)/(1-exp((global_par_V_membrane+14.1)/(-5))));
% assignmentRule: variable = beta_xr
	global_par_beta_xr=piecewise(0.00037836118, abs(global_par_V_membrane-3.3328) < 1E-10, 7.3898E-5*(global_par_V_membrane-3.3328)/(exp((global_par_V_membrane-3.3328)/5.1237)-1));
% assignmentRule: variable = tau_xr
	global_par_tau_xr=1/(global_par_alpha_xr+global_par_beta_xr);
% assignmentRule: variable = xr_infinity
	global_par_xr_infinity=(1+exp((global_par_V_membrane+14.1)/(-6.5)))^(-1);
% assignmentRule: variable = i_Ks
	global_par_i_Ks=global_par_Cm*global_par_g_Ks*global_par_xs^2*(global_par_V_membrane-global_par_E_K);
% assignmentRule: variable = alpha_xs
	global_par_alpha_xs=piecewise(0.00068, abs(global_par_V_membrane-19.9) < 1E-10, 4E-5*(global_par_V_membrane-19.9)/(1-exp((global_par_V_membrane-19.9)/(-17))));
% assignmentRule: variable = beta_xs
	global_par_beta_xs=piecewise(0.000315, abs(global_par_V_membrane-19.9) < 1E-10, 3.5E-5*(global_par_V_membrane-19.9)/(exp((global_par_V_membrane-19.9)/9)-1));
% assignmentRule: variable = tau_xs
	global_par_tau_xs=0.5*(global_par_alpha_xs+global_par_beta_xs)^(-1);
% assignmentRule: variable = xs_infinity
	global_par_xs_infinity=(1+exp((global_par_V_membrane-19.9)/(-12.7)))^(-0.5);
% assignmentRule: variable = i_Ca_L
	global_par_i_Ca_L=global_par_Cm*global_par_g_Ca_L*global_par_d*global_par_f*global_par_f_Ca*(global_par_V_membrane-65);
% assignmentRule: variable = d_infinity
	global_par_d_infinity=(1+exp((global_par_V_membrane+10)/(-8)))^(-1);
% assignmentRule: variable = tau_d
	global_par_tau_d=piecewise(4.579/(1+exp((global_par_V_membrane+10)/(-6.24))), abs(global_par_V_membrane+10) < 1E-10, (1-exp((global_par_V_membrane+10)/(-6.24)))/(0.035*(global_par_V_membrane+10)*(1+exp((global_par_V_membrane+10)/(-6.24)))));
% assignmentRule: variable = f_infinity
	global_par_f_infinity=exp((-(global_par_V_membrane+28))/6.9)/(1+exp((-(global_par_V_membrane+28))/6.9));
% assignmentRule: variable = tau_f
	global_par_tau_f=9*(0.0197*exp((-0.0337^2)*(global_par_V_membrane+10)^2)+0.02)^(-1);
% assignmentRule: variable = f_Ca_infinity
	global_par_f_Ca_infinity=(1+global_par_Ca_i/0.00035)^(-1);
% assignmentRule: variable = tau_f_Ca
	global_par_tau_f_Ca=2;
% assignmentRule: variable = sigma
	global_par_sigma=1/7*(exp(global_par_Na_o/67.3)-1);
% assignmentRule: variable = f_NaK
	global_par_f_NaK=(1+0.1245*exp((-0.1)*global_par_F*global_par_V_membrane/(global_par_R*global_par_T))+0.0365*global_par_sigma*exp((-global_par_F)*global_par_V_membrane/(global_par_R*global_par_T)))^(-1);
% assignmentRule: variable = i_NaK
	global_par_i_NaK=global_par_Cm*global_par_i_NaK_max*global_par_f_NaK*1/(1+(global_par_Km_Na_i/global_par_Na_i)^1.5)*global_par_K_o/(global_par_K_o+global_par_Km_K_o);
% assignmentRule: variable = E_Ca
	global_par_E_Ca=global_par_R*global_par_T/(2*global_par_F)*log(global_par_Ca_o/global_par_Ca_i);
% assignmentRule: variable = i_B_Na
	global_par_i_B_Na=global_par_Cm*global_par_g_B_Na*(global_par_V_membrane-global_par_E_Na);
% assignmentRule: variable = i_B_Ca
	global_par_i_B_Ca=global_par_Cm*global_par_g_B_Ca*(global_par_V_membrane-global_par_E_Ca);
% assignmentRule: variable = i_B_K
	global_par_i_B_K=global_par_Cm*global_par_g_B_K*(global_par_V_membrane-global_par_E_K);
% assignmentRule: variable = i_NaCa
	global_par_i_NaCa=global_par_Cm*global_par_I_NaCa_max*(exp(global_par_gamma*global_par_F*global_par_V_membrane/(global_par_R*global_par_T))*global_par_Na_i^3*global_par_Ca_o-exp((global_par_gamma-1)*global_par_F*global_par_V_membrane/(global_par_R*global_par_T))*global_par_Na_o^3*global_par_Ca_i)/((global_par_K_mNa^3+global_par_Na_o^3)*(global_par_K_mCa+global_par_Ca_o)*(1+global_par_K_sat*exp((global_par_gamma-1)*global_par_V_membrane*global_par_F/(global_par_R*global_par_T))));
% assignmentRule: variable = i_CaP
	global_par_i_CaP=global_par_Cm*global_par_i_CaP_max*global_par_Ca_i/(0.0005+global_par_Ca_i);
% assignmentRule: variable = Fn
global_par_V_rel = 96.48;
global_par_i_rel=global_par_K_rel*global_par_u^2*global_par_v*global_par_w*(global_par_Ca_rel-global_par_Ca_i);

if time >= 100 && time <=102
    global_par_i_rel =6;
end

global_par_Fn=1000*(1E-15*global_par_V_rel*global_par_i_rel-1E-15/(2*global_par_F)*(0.5*global_par_i_Ca_L-0.2*global_par_i_NaCa));

% if time>101
%     global_par_i_rel
% end
% assignmentRule: variable = i_rel
% assignmentRule: variable = tau_u
	global_par_tau_u=8;
% assignmentRule: variable = u_infinity
	global_par_u_infinity=(1+exp((-(global_par_Fn-3.4175E-13))/1.367E-15))^(-1);
% assignmentRule: variable = tau_v
	global_par_tau_v=1.91+2.09*(1+exp((-(global_par_Fn-3.4175E-13))/1.367E-15))^(-1);
% assignmentRule: variable = v_infinity
	global_par_v_infinity=1-(1+exp((-(global_par_Fn-6.835E-14))/1.367E-15))^(-1);
% assignmentRule: variable = tau_w
	global_par_tau_w=piecewise(6*0.2/1.3, abs(global_par_V_membrane-7.9) < 1E-10, 6*(1-exp((-(global_par_V_membrane-7.9))/5))/((1+0.3*exp((-(global_par_V_membrane-7.9))/5))*1*(global_par_V_membrane-7.9)));
% assignmentRule: variable = w_infinity
	global_par_w_infinity=1-(1+exp((-(global_par_V_membrane-40))/17))^(-1);
% assignmentRule: variable = i_tr
	global_par_i_tr=(global_par_Ca_up-global_par_Ca_rel)/global_par_tau_tr;
% assignmentRule: variable = i_up
	global_par_i_up=global_par_I_up_max/(1+global_par_K_up/global_par_Ca_i);
% assignmentRule: variable = i_up_leak
	global_par_i_up_leak=global_par_I_up_max*global_par_Ca_up/global_par_Ca_up_max;
% assignmentRule: variable = Ca_CMDN
	global_par_Ca_CMDN=global_par_CMDN_max*global_par_Ca_i/(global_par_Ca_i+global_par_Km_CMDN);
% assignmentRule: variable = Ca_TRPN
	global_par_Ca_TRPN=global_par_TRPN_max*global_par_Ca_i/(global_par_Ca_i+global_par_Km_TRPN);
% assignmentRule: variable = Ca_CSQN
	global_par_Ca_CSQN=global_par_CSQN_max*global_par_Ca_rel/(global_par_Ca_rel+global_par_Km_CSQN);
% assignmentRule: variable = V_i
	global_par_V_i=global_par_V_cell*0.68;
% assignmentRule: variable = V_rel
	global_par_V_rel=0.0048*global_par_V_cell;
% assignmentRule: variable = V_up
	global_par_V_up=0.0552*global_par_V_cell;
% assignmentRule: variable = B1
	global_par_B1=(2*global_par_i_NaCa-(global_par_i_CaP+global_par_i_Ca_L+global_par_i_B_Ca))/(2*global_par_V_i*global_par_F)+(global_par_V_up*(global_par_i_up_leak-global_par_i_up)+global_par_i_rel*global_par_V_rel)/global_par_V_i;
% assignmentRule: variable = B2
	global_par_B2=1+global_par_TRPN_max*global_par_Km_TRPN/(global_par_Ca_i+global_par_Km_TRPN)^2+global_par_CMDN_max*global_par_Km_CMDN/(global_par_Ca_i+global_par_Km_CMDN)^2;

	xdot=zeros(21,1);
	% rateRule: variable = V_membrane
	xdot(1) = (-(global_par_i_Na+global_par_i_K1+global_par_i_to+global_par_i_Kur+global_par_i_Kr+global_par_i_Ks+global_par_i_B_Na+global_par_i_B_Ca+global_par_i_NaK+global_par_i_CaP+global_par_i_NaCa+global_par_i_Ca_L+global_par_i_st))/global_par_Cm;
	% rateRule: variable = m
	xdot(2) = (global_par_m_inf-global_par_m)/global_par_tau_m;
	% rateRule: variable = h
	xdot(3) = (global_par_h_inf-global_par_h)/global_par_tau_h;
	% rateRule: variable = j
	xdot(4) = (global_par_j_inf-global_par_j)/global_par_tau_j;
	% rateRule: variable = oa
	xdot(5) = (global_par_oa_infinity-global_par_oa)/global_par_tau_oa;
	% rateRule: variable = oi
	xdot(6) = (global_par_oi_infinity-global_par_oi)/global_par_tau_oi;
	% rateRule: variable = ua
	xdot(7) = (global_par_ua_infinity-global_par_ua)/global_par_tau_ua;
	% rateRule: variable = ui
	xdot(8) = (global_par_ui_infinity-global_par_ui)/global_par_tau_ui;
	% rateRule: variable = xr
	xdot(9) = (global_par_xr_infinity-global_par_xr)/global_par_tau_xr;
	% rateRule: variable = xs
	xdot(10) = (global_par_xs_infinity-global_par_xs)/global_par_tau_xs;
	% rateRule: variable = d
	xdot(11) = (global_par_d_infinity-global_par_d)/global_par_tau_d;
	% rateRule: variable = f
	xdot(12) = (global_par_f_infinity-global_par_f)/global_par_tau_f;
	% rateRule: variable = f_Ca
	xdot(13) = (global_par_f_Ca_infinity-global_par_f_Ca)/global_par_tau_f_Ca;
	% rateRule: variable = u
	xdot(14) = (global_par_u_infinity-global_par_u)/global_par_tau_u;
	% rateRule: variable = v
	xdot(15) = (global_par_v_infinity-global_par_v)/global_par_tau_v;
	% rateRule: variable = w
	xdot(16) = (global_par_w_infinity-global_par_w)/global_par_tau_w;
	% rateRule: variable = Na_i
	xdot(17) = ((-3)*global_par_i_NaK-(3*global_par_i_NaCa+global_par_i_B_Na+global_par_i_Na))/(global_par_V_i*global_par_F);
	% rateRule: variable = K_i
	xdot(18) = (2*global_par_i_NaK-(global_par_i_K1+global_par_i_to+global_par_i_Kur+global_par_i_Kr+global_par_i_Ks+global_par_i_B_K))/(global_par_V_i*global_par_F);
	% rateRule: variable = Ca_i
	xdot(19) = global_par_B1/global_par_B2;
	% rateRule: variable = Ca_up
	xdot(20) = global_par_i_up-(global_par_i_up_leak+global_par_i_tr*global_par_V_rel/global_par_V_up);
	% rateRule: variable = Ca_rel
	xdot(21) = (global_par_i_tr-global_par_i_rel)*(1+global_par_CSQN_max*global_par_Km_CSQN/(global_par_Ca_rel+global_par_Km_CSQN)^2)^(-1);
end
