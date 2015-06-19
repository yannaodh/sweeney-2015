### WARNING: this script launches a process for each combination
### of parameters specified, so either run on a multicore machine
### or run a single case at a time

from brian import *
import network_simulation
import os
import datetime
import common
import numpy
import collections
import multiprocessing

script_pars = {
    'log': True,
    'remote': os.environ['USER']=='s1144270',
    'get_NOtarget_ext_rate': False
    }

now = datetime.datetime.now()

if script_pars['log']:
    script_pars['logdir'] = 'simulation-results/fig-4/'
    os.system('mkdir simulation-results')
    os.system('mkdir simulation-results/fig-4')

## default simulation parameters. Will be updated below depending on experiment
simulation_pars = {
    'n_cells': 5000,
    'input_group_pop': 100,
    'sim_step': 1000 *ms,
    'sampletime': 50000 *ms,
    'NOsettletime': 300000 *ms,
    'NOthresh_scaling': 1.0,
    'autothresh_time': 200000 *ms,
    'settletime': 500 *ms,
    'maxthresh_NO_auto': True,
    'scale_to_mean': True,
    'scale_to_rate': False,
    'Nsample_rate_scale': 20,
    'sample_neurons_N': 0,
    'record_pop_rate': False,
    'mod_targetNOconc': True,
    'init_homogenous_inputs': False,
    'single_group': True ,
    'main_group_input_type': 'uniform',
    'num_input_groups': 0,
    'interactive_plotting': False,
    'log': False,
    'input_rates' : [10*Hz,5*Hz],
    'ext_rate': 20.0 * Hz,
    'ext_mean': 25.0,
    'ext_std': 5.0,
    'g_ext': 80.0 * nS, ## originally 80.0 * nS
    'g_exc': 4.0 * nS,
    'g_inh': 64.0 * nS,
    'C': 100,
    'NO_diff': 0.1,
    'NOdecay_rate': 10.0,
    'local_NOdecay_factor': 5.0,
    'ds': 0.0020,
    'maxthresh_NO': 4.0e-6,
    'Vt_min': -58.0 *mV,
    'Vt_max': -40.0 *mV,
    'bounded_Vt': True,
    's_lat': 0.1,
    'weight_dist': 'uniform',
    'centre_group_and_ring': False,
    'savefig_id': False,
    'logdir': script_pars['logdir'],
    'diffusing': True,
    'mod_tau': 2500,
    'modulating': True,
    'nNOS_inhibitory_only': False,
    'nNOS_ampa_activated': False,
    'global_NOreadout': False,
    'restore_inputs': False,
    'reconfigure_inputs': False,
    'ind_ext_input': False,
    'tau_OU': 1*ms,
    'sigma_OU': 1*mV,
    'nu_ext': 2*mV,
    'size': 1.0,
    'use_C_over_epsilon': True,
    'variable_targets': False,
    'target_variance': 0.0,
    'var_targets_premod_std': 5.0
    }

comb_pars = []

from default_params import comb_pars_reconf_exp

comb_pars = comb_pars_reconf_exp

comb_pars_reconf_exp = []
comb_pars_reconf_exp.append({'diffusing':True, 'ds':0.0015})

from default_params import input_group_exp
simulation_pars.update(input_group_exp)

simulation_pars['dump_each_simstep'] = False

simulation_pars['reload_sim'] = False
simulation_pars['reload_sim_from_dir'] = False

sweep_pars = {
      'trial_ID': [0]
}

if script_pars['log']:
    common.save_pickle_safe(script_pars['logdir']+'simulations_pars.pickle',
            {'simulation_pars':simulation_pars,'comb_pars':comb_pars,'sweep_pars':sweep_pars})

results = []

def run_iter(iter_pars):
    print iter_pars
    logfile = iter_pars['logdir']+"{0}{1}.pickle".format('results',iter_pars['sim_id'])
    iter_pars['logfile']=logfile
    results.append(network_simulation.main(iter_pars))
    if script_pars['log']:
        common.save_pickle_safe(iter_pars['logdir']+"{0}{1}.pickle".format('results',iter_pars['sim_id']),results[-1])

iter_pass = []

simulation_pars_orig = simulation_pars.copy()
for combo in comb_pars:
    simulation_pars = simulation_pars_orig.copy()
    simulation_pars.update(combo)
    for parkey in sweep_pars.keys():
        if isinstance(sweep_pars[parkey],collections.Iterable):
            for parval in sweep_pars[parkey]:
                simulation_pars.update({parkey:parval})
                simulation_pars['sim_id'] = common.create_info_string(combo) + common.create_info_string(dict(filter(lambda i:i[0] in sweep_pars, simulation_pars.iteritems())))
                iter_pass.append(simulation_pars)

                multiprocessing.Process(target=run_iter,args=[simulation_pars]).start()

                ## comment out above line, and use below line for single core process
                #run_iter(simulation_pars)


