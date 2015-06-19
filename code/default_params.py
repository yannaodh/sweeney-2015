from brian import *
import numpy

default_params = {

### Cell parameters ###
   'tau_m'   : 20. * ms, # Membrane time constant
    'c_m'     : 0.2 * nF, # Capacitance
    'tau_exc' :  3. * ms, # Synaptic time constant (excitatory)
    'tau_inh' :  7. * ms, # Synaptic time constant (inhibitory)
    'tau_ref' :  5. * ms, # Refractory period
    'El'      : -80 * mV, # Leak potential
    'Ee'      :  0. * mV, # Reversal potential (excitation)
    'Ei'      : -70.* mV, # Reversal potential (inhibition)
    'Vt_init' : -50 * mV, # Spike Threhold
    'Vr'      : -60 * mV, # Spike Reset
    'Ca_tau'  : 10 * msecond,            #tau for Ca decay
    'nNOS_tau': 100 * msecond,            #tau for nNOS activation from Ca conc.
    'ca_spike': 1.0,                     #Ca influx after spike
    ### Hill function (nNOS =  hillfunction(Ca)) coefficients ###
    'hill_K'  : 1.0,
    'hill_n'  : 3.0,
    ## range for modulation of threshold
    'Vt_min'  : - 52.0 *mV,
    'Vt_max'  : - 42.0 *mV,
    'bounded_Vt': True,
    'Vt_init_nomod' : -47.8 * mV,
    'mod_tau' : 100, #50 # should depend on dt_modulation [mod_tau( in ms) = dt_modulation*mod_tau]

    'weight_dist' : 'uniform',  ### 'uniform, lognormal, seperate_groups_uniform'

    'modulate_only_input_groups' : False,
    'single_group' : False,
    'main_group_input_type' : 'uniform', ###!!normal doesn't work yet!!!  #types: uniform, normal [add lognormal option?]
    'NO_monitor' : True,
    'diffusing': True   ,
    'modulating' : True ,
    'NO_mod'  : True,
    'no_groupsize' : True,
    'centre_group_and_ring' : True,
    'mod_targetNOconc' : False  ,
    'init_homogenous_inputs' : False,
    'restore_inputs': False,
    'reconfigure_inputs': False,
    'n_cells'       : 300, #12500                   # Total number of cells
    #if remote:
    #    n_cells : 12500,
    'n_exc_ratio'   : 0.8 ,      # 4:1 ratio for exc/inh
    'size'         : 1.0,                      # Size of the network - in mm (roughly... it is a torus)
    #ds            : min(0.0025*size,0.01)       # spatial resolution of network (needed for saving NO concentrations at each point)
    'ds'            : 0.0025,
    'sampletime'    : 1500 * ms,               # Simulation time
    'NOsettletime'  : 500 * ms,
    'settletime'    : 500 * ms,
    'autothresh_time': 2000 *ms,
    'sim_step'      : 500 * ms ,                 # Display snapshots every sim_step ms
    'epsilon'       : 0.02     ,               # Probability density
    'C'             : 1000 ,                    # mean number of synaptic inputs per neuron (makes epsilon redundant) [possibly way too higher for small networks (N<=1000)]
    'use_C_over_epsilon': True,
    's_lat'         : 1.0       ,              # Spread of the lateral connections
    'g_exc'         : 4.  * nS   ,#4.             # Excitatory conductance
    'g_inh'         : 64.0 * nS  , #64.             # Inhibitory conductance
    'g_cov'         : 0.5 ,
    'g_ext'         : 80 * nS,   #200.             # External drive
    'velocity'      : 0.3 * mm/ms,             # velocity
    'ext_rate'      : 0.0 * Hz, #200 * Hz               # Rate of the external source
    'ext_mean'      : 3.5 ,
    'ext_std'      : 2.5 ,
    'dt_diffusion'  : 1*ms , #was 50 ms
    'NO_diff'       : 0.1 ,                    # diffusion constant for NO
    'NOdecay_rate'  : 10.0 ,                    # decay rate for NO (due to e.g. oxidisation etc.) -> should result in half-life of ~ 5s for free NO according to philippedes
    'local_NOdecay_factor': 5.0,
    'dt_modulation' : 5*ms,

    'target_NO' : 6e-9,

    'maxthresh_NO_auto': True,
    'NOthresh_scaling': 1.5,
    'scale_to_mean': True,
    'scale_to_rate': True,
    'Nsample_rate_scale': 20,
    'maxthresh_NO' : 3e-7,         #value for NO concentration at which Vt -> Vt_max
    'nNOS_inhibitory_only': False,
    'nNOS_ampa_activated': False,
    'global_NOreadout': False,
    'neighbour_readout': False,
    'k_neighbours': 1,

    'num_input_groups' : 2,
    'input_group_size' : 0.2,    # 2.0
    'input_group_pop' : 100, #50
    'input_rates' : [10*Hz,0.0*Hz],

    'ind_ext_input': True,
    'tau_OU': 1*ms,
    'sigma_OU': 1*mV,
    'nu_ext': 2*mV,

    'no_input_change_autothresh': False,
    'connect_groups_after_settle': False,
    'connect_group_settletime': 10000*ms,
    'connection_group_strength_factor': 2.0,
    'dynamic_connectivity': False,

    'present_new_input_stats': False,
    'get_network_dynamic_range': False,
    'do_orientation_tasks': False,

    #### STDP params (originally taken from STDP2.py brian example, Adapted from Song, Miller and Abbott (2000), Song and Abbott (2001)
    #and van Rossum et al (2000).
    'tau_pre': 20 * ms,
    'tau_post': 20 *ms,
    'gmax': 20.0 *nS,
    'dA_pre': .01,
    'dA_post': -.01 * 20*ms / 20*ms * 2.5,

#    'ring_r2' : input_group_size/2*1.0,
#    'ring_r1' : input_group_size/2*0.8,

    'size_display_group' : 25,
    'sample_neurons_N': 5,

    'record_pop_rate': False,

    'log': True,
    'remote': False,
    'interactive_plotting': True,
    'savefig_id': False,
    'dump_each_simstep': False,
    'flush_after_dump': True,
    'save_raster_spikes': False,

    'reload_sim': False,
    'reload_sim_from_dir': False,

    'sim_id': '',
    'logdir': 'output/'
    }


runsim_test = {
    'n_cells': 1000,
    'input_group_pop': 100,
    'sim_step': 1000 *ms,
    'sampletime': 25000 *ms,
    'NOsettletime': 50000 *ms,
    'NOthresh_scaling': 1.0,
    'autothresh_time': 60000 *ms,
    'settletime': 500 *ms,
    'maxthresh_NO_auto': True,
    'scale_to_mean': True,
    'scale_to_rate': True,
    'Nsample_rate_scale': 20,
    'sample_neurons_N': 0,
    'record_pop_rate': True,
    'mod_targetNOconc': True,
    'init_homogenous_inputs': False,
    'single_group': True,
    'main_group_input_type': 'normal',
    'num_input_groups': 0,
    'interactive_plotting': False,
    'log': False,
    'input_rates' : [10*Hz,5*Hz],
    'ext_rate': 2.5 * Hz,
    'ext_mean': 10.0,
    'ext_std': 5.0,
    'g_ext': 80.0 * nS, ## originally 80.0 * nS
    'g_exc': 4.0 * nS,
    'g_inh': 90.0 * nS,
    'C': 50,
#    'epsilon': 0.02,
    'NO_diff': 0.1,
    'NOdecay_rate': 10.0,
    'local_NOdecay_factor': 5.0,
    'ds': 0.0020,
    'maxthresh_NO': 4.0e-6,
    'Vt_init_nomod' : -50.0 * mV,
    'Vt_min': -58.0 *mV,
    'Vt_max': -40.0 *mV,
    'bounded_Vt': True,
    's_lat': 0.1,
    'weight_dist': 'uniform',
    'centre_group_and_ring': False,
    'savefig_id': True,
    'diffusing': True,
    'mod_tau': 250,
    'modulating': True,
    'nNOS_inhibitory_only': False,
    'nNOS_ampa_activated': False,
    'global_NOreadout': False,
    'restore_inputs': False,
    'reconfigure_inputs': False,
    'ind_ext_input': True,
    'tau_OU': 1*ms,
    'sigma_OU': 1.0*mV,
    'nu_ext': 1.0*mV,
    'size': 1.0,
    'use_C_over_epsilon': True,
    'variable_targets': False,
    'target_variance': 3.3,
    'decay_weights': False,
    'do_scaling': False,
    'only_stdp_first': False,
    'only_stdp_second': True,
    'stdp_only_time': 30000*ms,
    'homeostasis_only_first': False,
    'homeostasis_only_second': False,
    'homeostasis_only_time':60000*ms,
    'weight_sample_time': 10000*ms,
    'get_transfer_funs': False,
    'N_transf_presentations': 100,
    't_transf_presentation': 5000*ms,
    'do_discrim_task': False,
    'get_pattern_responses': False,
    't_discrim_presentation': 1000*ms,
    'N_discrim_presentations': 200,
    'discrim_stim_difference': 0.1,
    'N_discrim_stim': 50,
    'discrim_pattern_size': 25,
    'subgroup_presentations_during_homeostasis': False
}

runsim_test.update({'diffusing':False,'main_group_input_type':'normal', 'target_NO': 0.00040, 'mod_tau':5000})
runsim_test.update({'global_NOreadout': False})

runsim_test.update({'gmax':200.0*nS, 'autothresh_time':80000*ms,'mod_tau':1000, 'ext_rate':10.0*Hz, 'C':20})

runsim_test.update({'NO_diff':0.025, 'size':0.5,'ds':0.00045})


sequence_params = {
    'N_seq': 2,
    'L_seq': 5,
    'seq_ISI': 10 *ms,
    'g_seq_ext': 20*nS,
    'g_seq_exc': 10*nS,
    'N_training_per': 250,
    'N_test_per': 15,
    'C': 20,
    'dA_pre': 0.001,
    'dA_post': -.0025,
    'gmax': 20.0*nS,
    'n_cells': 100,
    'mod_tau': 100,
    'interactive_plotting': True,
    'weight_sample_time': 1000*ms,
    'perturb_sample_weights': False,
    'manual_STDP_weight_dist': False,
    'combo_readout': True,
    'alpha': 0.0
}

pattern_params = {
    ### should have target rate pretty high in order to get good information from spike counts over exposure time
    'homeostasis_only_time':60000*ms,
    'homeostasis_only_first': True,
    'modulating': True,
    'plasticity': True,
    'freeze_recall_homeostasis': False,
    'ext_mean': 25.0,
    'C': 50,
    'n_cells': 200,
    'dA_post': -.012,
    'N_patterns': 4,
    'pattern_group_size': 32,
    'N_random_test_patterns': 3,
    'N_groups_per_pattern': 2,
    'presentations_per_epoch': 10,
    'N_epoch_no_stdp': 50,
    'one_shot_learning': True,
    'N_epoch': 100,
    'pattern_exposure_time': 2500*ms,
    'pattern_response_time': 250*ms,
    'sequential_pattern_learning': False
}

test_params = {
    'autothresh_time':1000*ms,
    'NOsettletime':1000*ms,
    'N_epoch':1,
    'homeostasis_only_time':1000*ms,
    'sampletime':1000*ms,
    'homeostasis_only_first': False,
    'ext_mean': 10.0,
    'nu_ext': 1.0*mV ,
    'sigma_OU': 2.0*mV ,
    'dt_diffusion': 10*ms,
    'N_discrim_presentations': 20,
    'discrim_stim_difference': 0.1,
    'N_discrim_stim': 2
}

discrim_params = {
    'autothresh_time':30000*ms,
    'NOsettletime': 100000*ms,
    'sampletime':25000*ms,
    'ds': 0.0015,
    'size': 1.0,
    'n_cells': 2500,
    'g_exc': 5.5*nS,
    'bounded_Vt': True,
    'epsilon': 0.02,
    'ext_std': 10,
    'ext_mean':10,
    'ext_rate': 5.0*Hz,
    'n_cells': 2500,
    'NO_diff': 0.1,
    'use_C_over_epsilon': False,
    'reconfigure_inputs': False,
}


rate_distribution_params = {
    'autothresh_time':30000*ms,
    'NOsettletime': 100000*ms,
    'sampletime':25000*ms,
    'ds': 0.002,
    'size': 1.0,
    'main_group_input_type':'normal',
    'n_cells': 5000,
    'bounded_Vt': True,
    'epsilon': 0.02,
    'ext_std': 10,
    'ext_mean':10,
    'ext_rate': 5.0*Hz,
    'use_C_over_epsilon': False,
    'reconfigure_inputs': True,
    'g_exc': 5.5*nS,
}

reconf_params = {
    'autothresh_time':30000*ms,
    'NOsettletime': 100000*ms,
    'sampletime':25000*ms,
    'ds': 0.002,
    'size': 1.0,
    'main_group_input_type':'normal',
    'n_cells': 2500,
    'g_exc': 5.5*nS,
    'bounded_Vt': True,
    'epsilon': 0.02,
    'ext_std': 10,
    'ext_mean':10,
    'ext_rate': 5.0*Hz,
    'use_C_over_epsilon': False,
    'reconfigure_inputs': True

}

new_input_stats_params = {
    'present_new_input_stats': True,
    'scale_current_input_stats': True, ### conserves rank of distribution
    'new_main_group_input_type': 'normal',
    'new_ext_std': 10.0,
    'new_ext_mean': 30.0,
    'NOsettletime': 1000*ms
}

get_network_dynamic_range_params = {
    'get_network_dynamic_range': True,
    'dynamic_range_sampletime': 2000*ms,
    'sim_step': 1000*ms,
    'network_input_range': range(0,100,1),
    'NOsettletime': 1000*ms,
    'sampletime': 1000*ms,
    'just_save_dynamic_range': True
}

comb_pars_reconf_exp = []
#comb_pars_reconf_exp.append({'diffusing':True})
#comb_pars_reconf_exp.append({'diffusing':True, 'ds':0.002})

comb_pars_reconf_exp.append({'diffusing':True, 'ds':0.0015})
comb_pars_reconf_exp.append({'diffusing':False})
comb_pars_reconf_exp.append({'diffusing':False,'variable_targets':True,'var_targets_premod_std':5.0,'var_targets_premod_mean':2.0})

reconf_exp_plasticity = {
    'reconfigure_inputs': True,
    'g_exc': 5.5*nS,
    'ext_std': 10.0,
    'homeostasis_only_first': False,
    'only_stdp_first': False,
    'only_stdp_second': False,
    'homeostasis_only_second': False,
    'mod_tau': 2500
}

input_group_exp = {
    'reconfigure_inputs': False,
    'n_cells': 2500,
    'input_group_pop': 250,
    'single_group': False,
    'main_group_input_type': 'uniform',
    'num_input_groups': 2,
    'input_rates': [10*Hz,5*Hz],
    'ext_rate': 2.5*Hz,
    'NOsettletime': 400000*ms,
    'autothresh_time': 100000*ms,
    'sampletime': 50000*ms,
    'ind_ext_input': False,
    'use_C_over_epsilon': True,
    'C': 100,
    'sim_step': 1000*ms,
    'g_exc': 5.5*nS,
}

toy_network = {
    'n_cells': 50,
    'epsilon': 0.4,
    'weight_sample_time': 100*ms,
    'connect_groups': False,
    'N_connection_groups': 1,
    'connection_group_size': 20,
    'connection_group_strength': 15.0*nS,
    'NOsettletime': 10000000*ms,
    'autothresh_time': 20000*ms,
    'sampletime': 20000*ms,
    'sample_neurons_N': 40,
    'dynamic_connectivity': True,
    'reconfigure_inputs': False,
    'ext_rate':5.0*Hz,
    'nu_ext':1.0*mV,
    'sigma_OU':1.0*mV,
    'tau_OU':1.0*ms,
    'sim_step':10000*ms,
    'ind_ext_input': False,
    'maxthresh_NO_auto': True,
    'target_NO':0.000175,
    'dump_each_simstep': True
}

orientation_task_pars = {
    'do_orientation_tasks': True,
    'just_save_orientation_responses': False,
    'sim_step': 1000*ms,
    'NOsettletime': 300000*ms,
    'sampletime': 1000*ms,
    'orientation_width': 0.025,
    'N_discrim_stim': 100,
    'subgroup_pattern_mag': 2.5,
    'no_orient_pref_bias': 00.0,
    't_discrim_presentation': 1000*ms
}


STDP_net_2500 = {
    'n_cells': 2500,
    'epsilon': 0.1,
    'bounded_Vt': False,
    'mod_tau': 15000*ms
}


