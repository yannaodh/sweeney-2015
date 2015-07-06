from brian import *
import numpy, os
from itertools import product
from random import shuffle
import datetime
import pickle
import sys
import common
from glob import glob

import scipy.stats
import scipy.sparse

from default_params import default_params, runsim_test

sim_properties = default_params

sim_properties.update(runsim_test)

#### returns random unique x,y positions for each cell
def random_pos(N,L):
    allpairs = list(product(*[xrange(int(L/ds))]*2))
    shuffle(allpairs)
    return numpy.array(allpairs[:N])
#### WARNING - unique positions not guaranteed here: check to see if overlaps are common??

def generate_pos():
    positions = numpy.zeros(shape=(n_cells,2))
    if not sim_properties['no_groupsize']:
        for i in xrange(num_input_groups):
            input_loc = random_pos(1,size)
            scatter = random_pos(input_group_pop,input_group_size)
            scatter = subtract(scatter,input_group_size/(2*ds))
            scattered = add(scatter,input_loc)
            positions[i*input_group_pop:((i+1)*input_group_pop)] = scattered
            positions[num_input_groups*input_group_pop:] = random_pos(n_cells-(num_input_groups*input_group_pop),size)
    else:
        positions = random_pos(n_cells,size)

    positions[positions>=size/ds]=positions[positions>=size/ds]-size/ds
    positions[positions<0]=positions[positions<0]+size/ds

    return positions

def threshold_fun(states):
    return states >= all_cells.Vt


### check if x11 session is available (for plotting)
def X_is_running():
    from subprocess import Popen, PIPE
    p = Popen(["xset", "-q"], stdout=PIPE, stderr=PIPE)
    p.communicate()
    return p.returncode == 0

### Function to get the distance between one position and an array of positions
### This is needed to used the vectorized form of the connections in the brian.Connection objects
### Note that the plane is wrapped, to avoid any border effects.
def get_distance(x, y):
    #print 'get_distance called'
    d1       = abs(x - y)
    min_d    = numpy.minimum(d1, size - d1)
    #print 'min_d', min_d
    #print 'size', len(min_d)
    return numpy.sqrt(numpy.sum(min_d**2))

### Function returning the probabilities of connections as a functions of distances
def probas(i, j, x, y,epsilon):
    #print 'probas called: i = ', i
    #print 'j = ', j
    distance = get_distance(x[i], y[j])
    #print 'distance call finished'
    s_lat = sim_properties['s_lat']
    return epsilon * numpy.exp(-distance**2/(2*s_lat**2))

### Function returning linear delays as function of distances
def delays(i, j, x, y,velocity):
    distance = get_distance(x[i], y[j])
    return 0.1*ms + (distance * mm )/ velocity

## this is called every timestep with a list of all neurons that spiked
def ca_influx(spikes):
    if sim_properties['nNOS_inhibitory_only']:   ## THIS IS DONE INEFFICIENTLY: since nNOS and Ca values are still tracked for
    ## excitatroy cells: best would be to have different eqs for exc cells
        inh_cells.Ca[spikes] = inh_cells.Ca[spikes] + sim_properties['ca_spike']
    else:
        all_cells.Ca[spikes] = all_cells.Ca[spikes] + sim_properties['ca_spike']

def feed_nNOSstates_to_grid():
    NO_rates[(all_cells.position[:]/(size*ds))[:].astype(int)[:,0],(all_cells.position[:]/(size*ds))[:].astype(int)[:,1]] = all_cells.nNOS

def feed_NOgrid_to_states():
    all_cells.NO = NO_concs[(all_cells.position[:]/(size*ds))[:].astype(int)[:,0],(all_cells.position[:]/(size*ds))[:].astype(int)[:,1]]
    #print numpy.mean(all_cells.NO)

def diffuse_NOgrid(): # taken from Philippedes et al. (2000, J Neurosci). Also code snippet from http://goo.gl/18lc8
    #	This function uses a numpy expression to evaluate the derivatives in the Laplacian, and calculates u[i,j] based on ui[i,j].
    global NO_concs

    NO_diff = sim_properties['NO_diff']
    NOdecay_rate = sim_properties['NOdecay_rate']

    ui = NO_concs.copy()
    nx = int(size/ds)
    ny = int(size/ds)
    leftX = arange(nx)-1
    rightX = (arange(nx)+1)%nx
    leftY = arange(ny)-1
    rightY = (arange(ny)+1)%ny
    # assuming no sinks so far
    NO_concs = ui + dt_diff*(NO_diff*((ui[rightX,:]-2*ui + ui[leftX,:])/dx2 + (ui[:,rightY] - 2*ui
        + ui[:,leftY])/dy2)+ NO_rates - NOdecay_rate*ui)

    # make sure all values for NO_conc >= 0
    NO_concs.clip(min=0,out=NO_concs)


def nodiffuse_NOgrid():
    #	This function uses a numpy expression to evaluate the derivatives in the Laplacian, and calculates u[i,j] based on ui[i,j].
    global NO_concs#, maxthresh_NO

    ui = NO_concs.copy()
    NO_concs = ui + dt_diff*(NO_rates - sim_properties['local_NOdecay_factor']*sim_properties['NOdecay_rate']*ui)

    # make sure all values for NO_conc >= 0
    NO_concs.clip(min=0,out=NO_concs)

def show_NOgrid(show):
    NO_max_t.append(numpy.max(NO_concs))
    NO_min_t.append(numpy.min(NO_concs))
    NO_mean_t.append(numpy.mean(NO_concs))
    NO_std_t.append(numpy.std(NO_concs))
    NO_rate_max_t.append(numpy.max(NO_rates))
    NO_rate_min_t.append(numpy.min(NO_rates))

    if show:
        fig3.cla()
        im = fig3.pcolor(NO_concs.transpose())
        setp(fig3, xticks=[], yticks=[])

        fig4.cla()
        fig4.plot(NO_max_t, 'b--')
        fig4.plot(NO_min_t, 'b--')
        fig4.fill_between(range(len(NO_mean_t)),add(NO_mean_t,NO_std_t),subtract(NO_mean_t,NO_std_t), facecolor='blue', alpha=0.5)
        fig4.set_ylabel("NO")
        fig4.set_xlabel("t")

    #### DUMP NOconc for plotting illustrative figure ###
    #NOfile = open('NOconc_dump/NOconc_'+str(len(NO_mean_t))+'.pickle','w')
    #pickle.dump(NO_concs.transpose(), NOfile,2)
    #NOfile.close()


def show_Vt(show):
    Vt_max_t.append(numpy.max(all_cells.Vt))
    Vt_min_t.append(numpy.min(all_cells.Vt))
    if not sim_properties['single_group']:
        for i in xrange(num_input_groups):
            if i ==0:
                VtA_std.append(numpy.std(all_cells.Vt[i*input_group_pop:(i+1)*input_group_pop]))
                VtA_mean.append(numpy.mean(all_cells.Vt[i*input_group_pop:(i+1)*input_group_pop]))
            else:
                VtB_std.append(numpy.std(all_cells.Vt[i*input_group_pop:(i+1)*input_group_pop]))
                VtB_mean.append(numpy.mean(all_cells.Vt[i*input_group_pop:(i+1)*input_group_pop]))
    Vt_std.append(std(all_cells.Vt[num_input_groups*input_group_pop:]))
    Vt_mean.append(mean(all_cells.Vt[num_input_groups*input_group_pop:]))
    Vt_samples.append(all_cells.Vt[0:sim_properties['sample_neurons_N']])

    if show:
        fig6.cla()
        fig6.hold(True)
        im = fig6.scatter(all_cells.position[:, 0], all_cells.position[:, 1], c=all_cells.Vt)
        im.set_clim(Vt_min_t[-1:][0], Vt_max_t[-1:][0])
        fig6.set_ylabel("Vt")
        #im.set_clim(min(NO_all.values[-1]), max(NO_all.values[-1]))
        setp(fig6, xticks=[], yticks=[])
        fig7.cla()
        fig7.hold(True)
        if not sim_properties['single_group']:
            fig7.errorbar(multiply(range(len(VtA_mean)),sim_step),VtA_mean,VtA_std,label='A')
            if sim_properties['num_input_groups'] > 1:
                fig7.errorbar(multiply(range(len(VtB_mean)),sim_step),VtB_mean,VtB_std,label='B')
        fig7.errorbar(multiply(range(len(Vt_mean)),sim_step),Vt_mean,Vt_std,label='rest')

        fig7.legend(loc="upper left", prop={'size':10})# fig5.legend()
        fig7.set_ylabel("max/min Vt")
        fig7.set_xlabel("tstep")

        fig12.cla()
        fig12.hold(True)
     #   bins = numpy.linspace(Vt_min, Vt_max, 50)
        colors = ["r","g","b"]
        if not sim_properties['single_group']:
            if num_input_groups > 1:
                fig12.hist([all_cells.Vt[num_input_groups*input_group_pop:],all_cells.Vt[input_group_pop:2*input_group_pop],all_cells.Vt[0:input_group_pop]],bins=50,
                        color=colors,normed=True, histtype='barstacked')
            else:
                fig12.hist([all_cells.Vt[num_input_groups*input_group_pop:],all_cells.Vt[0:input_group_pop]],bins=50,
                        color=pylab.colors(),normed=True, histtype='barstacked')
        else:
            try:
                fig12.hist([all_cells.Vt],bins=50,normed=True)
            except:
                pass

def show_dist_Vm_sample():
    fig13.cla()
    fig13.hold(True)
    fig13.plot(Vm_samplecell.values)
    fig13.set_xlabel("mV")

def modulate_neurons():
    feed_NOgrid_to_states()
    if sim_properties['mod_targetNOconc']:
        NOcells = all_cells.NO
    else:
        NOcells = numpy.minimum(all_cells.NO,maxthresh_NO)
    if sim_properties['global_NOreadout']:
        NOcells = numpy.ones(n_cells)*numpy.mean(NOcells)
    if sim_properties['variable_targets']:
        if not sim_properties['mod_targetNOconc']:
            if sim_properties['modulate_only_input_groups']:
                all_cells.Vt[0:num_input_groups*input_group_pop] = all_cells.Vt[0:num_input_groups*input_group_pop] + (1./mod_tau)*((Vt_min+(NOcells[0:num_input_groups*input_group_pop])*(1./maxthresholds_NO[0:num_input_groups*input_group_pop])*(Vt_max-Vt_min))-all_cells.Vt[0:num_input_groups*input_group_pop])
            else:
                all_cells.Vt[:] = all_cells.Vt[:] + (1./mod_tau)*((Vt_min+(NOcells[:])*(1./maxthresholds_NO[:])*(Vt_max-Vt_min))-all_cells.Vt[:])
        else:
            if sim_properties['modulate_only_input_groups']:
                all_cells.Vt[0:num_input_groups*input_group_pop] = all_cells.Vt[0:num_input_groups*input_group_pop] + 0.5*mV*(1./mod_tau)*((NOcells[0:num_input_groups*input_group_pop]-NO_targets[0:num_input_groups*input_group_pop])/NO_targets[0:num_input_groups*input_group_pop])
            else:
                all_cells.Vt[:] = all_cells.Vt[:] + 0.5*mV*(1./mod_tau)*((NOcells[:]-NO_targets[:])/NO_targets[:])
    else:
        if not sim_properties['mod_targetNOconc']:
            if sim_properties['modulate_only_input_groups']:
                all_cells.Vt[0:num_input_groups*input_group_pop] = all_cells.Vt[0:num_input_groups*input_group_pop] + (1./mod_tau)*((Vt_min+(NOcells[0:num_input_groups*input_group_pop])*(1./maxthresh_NO)*(Vt_max-Vt_min))-all_cells.Vt[0:num_input_groups*input_group_pop])
            else:
                all_cells.Vt[:] = all_cells.Vt[:] + (1./mod_tau)*((Vt_min+(NOcells[:])*(1./maxthresh_NO)*(Vt_max-Vt_min))-all_cells.Vt[:])
        else:
            if sim_properties['modulate_only_input_groups']:
                all_cells.Vt[0:num_input_groups*input_group_pop] = all_cells.Vt[0:num_input_groups*input_group_pop] + 0.5*mV*(1./mod_tau)*((NOcells[0:num_input_groups*input_group_pop]-target_NO)/target_NO)
            else:
                all_cells.Vt[:] = all_cells.Vt[:] + 0.5*mV*(1./mod_tau)*((NOcells[:]-target_NO)/target_NO)

    # bound the values of Vt to be < max, > min
    if sim_properties['bounded_Vt']:
        all_cells.Vt.clip(Vt_min,Vt_max,out=all_cells.Vt)


def interactive_step(time_step, init, plot_hist, show):
    global max_initrate
    global max_std
    global max_cov

    fig1.cla()
    if not sim_properties['single_group'] and show:
        raster_plot(raster_A,raster_B,raster_rest,showgrouplines=True,spacebetweengroups=0.01)
    elif show:
        raster_plot(raster_rest,showgrouplines=True,spacebetweengroups=0.01)

    rates = multiply(s_all_short.count,1000.*msecond/(sim_step))

    if not sim_properties['single_group']:
        for i in xrange(num_input_groups):
            if i ==0:
                ratesA_min.append(min(rates[i*input_group_pop:(i+1)*input_group_pop]))
                ratesA_max.append(max(rates[i*input_group_pop:(i+1)*input_group_pop]))
                ratesA_std.append(numpy.std(rates[i*input_group_pop:(i+1)*input_group_pop]))
                ratesA_mean.append(numpy.mean(rates[i*input_group_pop:(i+1)*input_group_pop]))
            else:
                ratesB_min.append(min(rates[i*input_group_pop:(i+1)*input_group_pop]))
                ratesB_max.append(max(rates[i*input_group_pop:(i+1)*input_group_pop]))
                ratesB_std.append(numpy.std(rates[i*input_group_pop:(i+1)*input_group_pop]))
                ratesB_mean.append(numpy.mean(rates[i*input_group_pop:(i+1)*input_group_pop]))
    rates_min.append(min(rates[num_input_groups*input_group_pop:]))
    rates_max.append(max(rates[num_input_groups*input_group_pop:]))
    rates_std.append(std(rates[num_input_groups*input_group_pop:]))
    rates_mean.append(mean(rates[num_input_groups*input_group_pop:]))
    rate_samples_exc.append(rates[0:sim_properties['sample_neurons_N']])
    rate_samples_inh.append(rates[n_exc:n_exc+sim_properties['sample_neurons_N']])
    rates_kurt_bias.append(scipy.stats.kurtosis(rates))
    rates_kurt_nobias.append(scipy.stats.kurtosis(rates,bias=False))

    print 'mean,std rates ', rates_mean[-1],rates_std[-1]

    #### DUMP rates for plotting illustrative figure ###
    #ratefile = open('rate_dump/rates_'+str(len(rates_mean))+'.pickle','w')
    #pickle.dump(rates.transpose(), ratefile,2)
    #ratefile.close()

    global CV_hist

    CV_all = []
    for i in xrange(n_cells):
        temp = CV(SM_all[i])
        if temp > 0 and not isnan(temp):
            CV_all.append(temp)
    CV_samples_exc.append(CV_all[0:sim_properties['sample_neurons_N']])
    CV_samples_inh.append(CV_all[n_exc:n_exc+sim_properties['sample_neurons_N']])
    CV_mean.append(mean(CV_all))
    CV_hist = CV_all

    rates = multiply(s_all.count,1000.*msecond/(sim_step*(time_step+1)))
    ## plot initial hist
    if init and plot_hist and show:
        fig8.cla()
        fig8.hold(True)
        max_initrate = max(rates)
        bins = numpy.linspace(0, max(rates), 50)
        colors = ["r","g","b"]

        if not sim_properties['single_group']:
            if sim_properties['num_input_groups'] > 1:
                fig8.hist([rates[num_input_groups*input_group_pop:],rates[input_group_pop:2*input_group_pop],rates[0:input_group_pop]],bins,
                        color=colors,normed=True, histtype='barstacked')
            else:
                fig8.hist([rates[0:input_group_pop],rates[num_input_groups*input_group_pop:]],bins,
                        color=pylab.colors(),normed=True,
                        histtype='barstacked',alpha=0.3)
                fig8.hist([rates],color='black',bins=bins,normed=True,histtype='step',alpha=1.0)
        else:
            fig8.hist([rates],color='r',bins=bins,normed=True)

        fig11.cla()
        fig11.hold(True)
        bins = numpy.linspace(Vt_min, Vt_max, 50)
        colors = ["r","g","b"]

        try:
            fig11.hist(all_cells.Vt,color='r',bins=50,normed=True)
        except:
            pass

        max_std = rates_std[0]
        max_cov = numpy.divide(rates_std[0],rates_mean[0])
    elif plot_hist and show:
        fig9.cla()
        fig9.hold(True)
        bins = numpy.linspace(0,max(rates), 50)
        colors = ["r","g","b"]

        if not sim_properties['single_group']:
            if sim_properties['num_input_groups'] > 1:
                fig9.hist([rates[num_input_groups*input_group_pop:],rates[input_group_pop:2*input_group_pop],rates[0:input_group_pop]],bins,
                    color=colors,normed=True, histtype='barstacked')
            else:
                fig9.hist([rates[0:input_group_pop],rates[num_input_groups*input_group_pop:]],bins,
                        color=pylab.colors(),normed=True,
                        histtype='barstacked',alpha=0.3)
                fig9.hist([rates],color='black',bins=bins,normed=True,histtype='step',alpha=1.0)
        else:
            fig9.hist([rates],bins=bins,normed=True)
    if show:
        fig5.cla()
        fig5.hold(True)
        if not sim_properties['single_group']:
            fig5.errorbar(multiply(range(len(ratesA_mean)),sim_step),ratesA_mean,ratesA_std,label='A')
            if sim_properties['num_input_groups'] > 1:
                fig5.errorbar(multiply(range(len(ratesB_mean)),sim_step),ratesB_mean,ratesB_std,label='B')
        fig5.errorbar(multiply(range(len(rates_mean)),sim_step),rates_mean,rates_std,label='rest')

        fig5.set_ylabel("rates")
        fig5.set_xlabel("t")
        fig5.legend(loc="upper left", bbox_to_anchor=(0.5,1), prop={'size':10})# fig5.legend()

        fig10.cla()
        fig10.hold(True)
        if not sim_properties['single_group']:
            fig10.plot(multiply(range(len(ratesA_std)),sim_step),numpy.divide(ratesA_std,ratesA_mean))
            if sim_properties['num_input_groups'] > 1:
                fig10.plot(multiply(range(len(ratesB_std)),sim_step),numpy.divide(ratesB_std,ratesB_mean))
            fig10.plot(multiply(range(len(rates_std)),sim_step),numpy.divide(rates_std,rates_mean))
        else:
            fig10.plot(multiply(range(len(rates_std)),sim_step),numpy.divide(numpy.divide(rates_std,rates_mean),max_cov),'b')
            fig10.plot(multiply(range(len(rates_std)),sim_step),numpy.divide(rates_std,max_std),'b--')

        show_dist_Vm_sample()


    show_NOgrid(show)
    if sim_properties['modulating']:
        show_Vt(show)

    s_all_short.count = numpy.zeros(n_cells)
    SM_all.reinit()

    if sim_properties['dump_each_simstep'] and int(time_step)%int(sim_properties['N_simstep_dump'])==0:
        results = {
                    'NO_concs': NO_concs,
                    'rate_samples_e': rate_samples_exc.pop(),
                    'rate_samples_i': rate_samples_inh.pop(),
                    'CV_samples_e': CV_samples_exc.pop(),
                    'CV_samples_i': CV_samples_inh.pop(),
                    'NO_postmod': all_cells.NO,
                    'thresholds': all_cells.Vt
                    #'NO_all': NO_all
                }
        if sim_properties['modulating']:
            results['Vt_samples'] = Vt_samples.pop()

        if time_step==0:
            print 'saving initial dump'
            results['sim_pars'] = sim_properties
            results['NO_threshold'] = maxthresh_NO
            results['NO_target'] = target_NO
            results['cell_positions'] = all_cells.position
            results['NO_postmod'] = all_cells.NO
            results['weights_dist'] =   Ce.W.todense()
            results['Ci_W'] = Ci.W.todense()
            results['nNOS_dist'] = all_cells.nNOS

            results['inputs_original'] = Cext.rate

        if sim_properties['variable_targets']:
            try:
                results['NO_targets'] = NO_targets
            except:
                pass

        if not sim_properties['single_group']:
            results['ratesA_mean'] = ratesA_mean
            results['ratesB_mean'] = ratesB_mean
            results['ratesA_std'] = ratesA_std
            results['ratesA_max'] = ratesA_max
            results['ratesA_min'] = ratesA_min
            results['ratesB_std'] = ratesB_std
            results['ratesB_min'] = ratesB_min
            results['ratesB_max'] = ratesB_max
            results['rates_min'] = rates_min
            results['rates_max'] = rates_max
            results['rates_mean'] = rates_mean
            results['rates_std'] = rates_std

        logfile = sim_properties['logfile']
        common.save_pickle_safe(logfile+'_dump_',results)


def main(dict_pars):
    sim_properties.update(dict_pars)


    ### reload network state from a previous run if below is true
    if sim_properties['reload_sim']:
        reload_results = common.load_pickle(sim_properties['reload_file'])
        sim_properties['reload_results'] = reload_results
        sim_properties['maxthresh_NO_auto'] = False
    elif sim_properties['reload_sim_from_dir']:
        print 'loading reloadable sim_file from ', glob(sim_properties['reload_dir']+"{0}{1}".format('results',sim_properties['sim_id'])+'*.pickle')[0]
        reload_results = common.load_pickle(glob(sim_properties['reload_dir']+"{0}{1}".format('results',sim_properties['sim_id'])+'*.pickle')[0])

        sim_properties['reload_results'] = reload_results
        sim_properties['maxthresh_NO_auto'] = False
        sim_properties['reload_sim'] = True

    #import matplotlib
    #if sim_properties['interactive_plotting']:
    #    matplotlib.use('TkAgg') # You may need to experiment, try WXAgg, GTKAgg, QTAgg, TkAgg

    import pylab

    now = datetime.datetime.now()
    logdir = '/Users/yann/NO/simresults/'+sys.argv[0].rpartition('.')[0]+'/' +now.strftime("%Y%m%d%H%M")

    ### for remote/ssh sessions ; CHANGEME - this is my username... need to upate with your ssh details
    ### Also, writes to /tmp file by default, so should probably change that
    if os.environ['USER'] == 's1144270':
        sim_properties['remote'] = True
        #logdir = '/afs/inf.ed.ac.uk/user/s11/s1144270/repos/NO/simresults/'+sys.argv[0].rpartition('.')[0]+'/' +now.strftime("%Y%m%d%H%M")
        logdir = '/tmp/' + now.strftime("%Y%m%d%H%M")

    ### We are setting the global timestep of the simulation
    simulation_clock=Clock(dt=0.1*ms)
    Clock(0.1 * ms)

    ### Cell parameters ###
    tau_m   = sim_properties['tau_m'] # Membrane time constant
    c_m     =  sim_properties['c_m'] # Capacitance
    tau_exc =  sim_properties['tau_exc']
    tau_inh =  sim_properties['tau_inh']
    tau_ref = sim_properties['tau_ref']
    El      =sim_properties['El']
    Ee      = sim_properties['Ee']
    Ei      = sim_properties['Ei']
    Vt_init =sim_properties['Vt_init']
    Vr      =sim_properties['Vr']
    Ca_tau  =sim_properties['Ca_tau']
    nNOS_tau=sim_properties['nNOS_tau']
    ### Hill function (nNOS =  hillfunction(Ca)) coefficients ###
    hill_K  =sim_properties['hill_K']
    hill_n  =sim_properties['hill_n']
    ## range for modulation of threshold
    global Vt_min, Vt_max, mod_tau
    Vt_min  =sim_properties['Vt_min']
    Vt_max  =sim_properties['Vt_max']
    Vt_init_nomod = sim_properties['Vt_init_nomod']
    mod_tau = sim_properties['mod_tau']

    ### Equation for a Conductance-based IAF ####

    if sim_properties['ind_ext_input']:
        tau_OU = sim_properties['tau_OU']
        sigma = sim_properties['sigma_OU']
        nu = sim_properties['nu_ext']
        eqs_iaf = Equations('''
        dv/dt  = (El-v)/tau_m + (ge*(Ee-v)+gi*(Ei-v))/c_m + nu/tau_OU + sigma*xi/tau_OU**.5 : volt
        dge/dt = -ge*(1./tau_exc) : uS
        dgi/dt = -gi*(1./tau_inh) : uS
        ''')
    else:
        eqs_iaf = Equations('''
        dv/dt  = (El-v)/tau_m + (ge*(Ee-v)+gi*(Ei-v))/c_m : volt
        dge/dt = -ge*(1./tau_exc) : uS
        dgi/dt = -gi*(1./tau_inh) : uS
        ''')

    eqs_Ca = """
    dCa/dt =  -Ca*(1./Ca_tau) : 1
    """

    eqs_nNOS = """
    dnNOS/dt = ((Ca**hill_n)/(Ca**hill_n+hill_K**hill_n))*(1./nNOS_tau) - nNOS*(1./nNOS_tau) : 1
    """

    eqs_nNOS_ampa = """
    dnNOS/dt = (((ge/(1*uS))**hill_n)/((ge/(1*uS))**hill_n+hill_K**hill_n))*(1./nNOS_tau) - nNOS*(1./nNOS_tau) : 1
    """

    if sim_properties['nNOS_ampa_activated']:
        eqs = eqs_iaf + eqs_nNOS_ampa
    else:
        eqs = eqs_iaf + eqs_Ca + eqs_nNOS

    global n_cells, n_exc, size
    n_cells       = sim_properties['n_cells']
    n_exc         = int(sim_properties['n_exc_ratio']*n_cells)
    size          = sim_properties['size']                      # Size of the network - in mm (roughly... it is a torus)

    if sim_properties['use_C_over_epsilon']:
        epsilon = min(float(float(sim_properties['C'])/n_cells),1.0)
        print 'epsilon: ', epsilon
    else:
        epsilon  = sim_properties['epsilon']                    # Probability density
    g_exc         = sim_properties['g_exc']
    g_inh         = sim_properties['g_inh']
    g_cov         = sim_properties['g_cov']
    g_ext         = sim_properties['g_ext']
    ext_rate      = sim_properties['ext_rate']
    ext_mean      = sim_properties['ext_mean']
    ext_std       = sim_properties['ext_std']
    max_distance  = size * mm/numpy.sqrt(2) # Since this is a torus
    max_delay     = max_distance/sim_properties['velocity']   # Needed for the connectors
    velocity = sim_properties['velocity']
    dt_diffusion  = sim_properties['dt_diffusion']
    NO_diff       = sim_properties['NO_diff']
    dt_modulation = 5*ms

    global target_NO, maxthresh_NO
    target_NO = sim_properties['target_NO']
    maxthresh_NO = sim_properties['maxthresh_NO']         #value for NO concentration at which Vt -> Vt_max

    global num_input_groups, input_group_pop, input_group_size
    num_input_groups = sim_properties['num_input_groups']
    input_group_size = sim_properties['input_group_size']
    input_group_pop = sim_properties['input_group_pop']

    input_rates = sim_properties['input_rates']

    size_display_group = sim_properties['size_display_group']

    ##### DON'T DEFINE ANY SIMULATION PARS AFTER THIS

    global ds, dx2, dy2, dt_diff
    ds = sim_properties['ds']
    dx2 = ds**2
    dy2 = ds**2
    #dt_diff = dx2*dy2/( 2*NO_diff*(dx2+dy2) )
    dt_diff = dt_diffusion/(1*ms)

    print 'DT DIFF , max stable', dx2*dy2/( 2*NO_diff*(dx2+dy2) )
    print 'DT DIFF, using ;', dt_diff

    if sim_properties['log']:
        if os.environ['USER'] == 's1144270':
            os.system('mkdir ~/repos/NO/simresults/'+sys.argv[0].rpartition('.')[0])
        else:
            os.system('mkdir /Users/yann/NO/simresults/'+sys.argv[0].rpartition('.')[0])
        os.system('mkdir ' + logdir )
        os.system('cp ' + sys.argv[0] + ' ' + logdir + '/')



    ### We create the cells and generate random positons in [0, size]x[0, size]
    global all_cells, exc_cells, inh_cells
    all_cells          = NeuronGroup(n_cells, model=eqs, threshold=threshold_fun, reset=Vr, refractory=tau_ref, clock=simulation_clock)
    all_cells.position = ds*size*generate_pos()
    exc_cells          = all_cells[0:n_exc]
    inh_cells          = all_cells[n_exc:n_cells]

    ### We initialize v values slightly above Vt, to have initial spikes
    all_cells.v        = El + 1.1*numpy.random.rand(n_cells) * (Vt_init - El)
    ### initialize threshold to be uniform at first
    if sim_properties['modulating']:
        all_cells.Vt = Vt_init*numpy.ones(n_cells)
    else:
        all_cells.Vt = Vt_init_nomod*numpy.ones(n_cells)

    if sim_properties['mod_targetNOconc'] and not sim_properties['maxthresh_NO_auto']:
        all_cells.NO = numpy.ones(n_cells)*target_NO
    else:
        all_cells.NO = numpy.zeros(n_cells)

    input_groups = []
    for i in xrange(num_input_groups):
        input_group = all_cells[i*input_group_pop:((i+1)*input_group_pop)]
        input_groups.append(input_group)


    ### Create 2-d array containing concentrations of NO at each x,y position
    global NO_concs
    if sim_properties['mod_targetNOconc'] and not sim_properties['maxthresh_NO_auto']:
        NO_concs = target_NO*numpy.ones((int(size/ds),int(size/ds)))
    else:
        NO_concs = numpy.zeros((int(size/ds),int(size/ds)))

    global NO_rates, NO_rate_max_t, NO_rate_min_t, NO_max_t, NO_min_t, NO_mean_t, NO_std_t, Vt_max_t, Vt_min_t, Vt_mean, Vt_std, VtA_mean, VtA_std, VtB_mean, VtB_std
    NO_rates = numpy.zeros((int(size/ds),int(size/ds)))
    NO_max_t = []
    NO_min_t = []
    NO_mean_t = []
    NO_std_t =  []
    NO_rate_max_t = []
    NO_rate_min_t = []
    Vt_max_t = []
    Vt_min_t = []
    Vt_mean = []
    Vt_std = []
    VtA_mean = []
    VtA_std = []
    VtB_mean = []
    VtB_std = []

    global ratesA_mean,ratesA_min,ratesA_max,ratesA_std, ratesB_mean,ratesB_min,ratesB_max,ratesB_std, rates_min, rates_max, rates_mean, rates_std, rates_kurt_nobias, rates_kurt_bias, CV_mean, CV_hist

    ratesA_max = []
    ratesA_min = []
    ratesA_mean = []
    ratesA_std =[]
    ratesB_max = []
    ratesB_min = []
    ratesB_mean = []
    ratesB_std =[]
    rates_max = []
    rates_min = []
    rates_mean = []
    rates_std =[]
    rates_kurt_bias = []
    rates_kurt_nobias = []
    CV_mean = []
    CV_hist = []

    global VmB_mean, VmB_std
    VmB_mean = []
    VmB_std = []


    diffusion_clock = Clock(dt=dt_diffusion)
    #diffusion_clock = Clock(dt=dt_diff*ms)  # setting clock to actual dt used in diffusion eq, which is minimum allowable dt for NO_diff and ds parameters
    modulation_clock = Clock(dt=dt_modulation)

    @network_operation(diffusion_clock, when='end')
    def do_diffusion(diffusion_clock):
        if(sim_properties['diffusing']):
            feed_nNOSstates_to_grid()
            diffuse_NOgrid()
        elif(sim_properties['NO_monitor']):
            feed_nNOSstates_to_grid()
            nodiffuse_NOgrid()

    @network_operation(modulation_clock, when='start')
    def do_modulation(modulation_clock):
        if sim_properties['NO_mod'] and sim_properties['modulating']:
            modulate_neurons()


    print "Building network..."


    C_inputs = []
    global Ce, Ci, Cext

    if sim_properties['reload_sim']:
        all_cells.position = sim_properties['reload_results']['cell_positions']
        try:
            all_cells.Vt = sim_properties['reload_results']['thresholds']
        except:
            try:
                all_cells.Vt = sim_properties['reload_results']['Vt_hist']
            except:
                print 'error: no thresholds to reload from'
                exit(1)
        NO_concs = sim_properties['reload_results']['NO_concs'].copy()
        all_cells.nNOS = sim_properties['reload_results']['nNOS_dist']

        Ce = Connection(exc_cells, all_cells, 'ge',max_delay=max_delay,
                                       delay     =0.1*ms)
        Ci = Connection(inh_cells, all_cells, 'gi',max_delay=max_delay,
                               delay     =0.1*ms)

        Ce.connect_from_sparse(W=scipy.sparse.csr_matrix(sim_properties['reload_results']['weights_dist']))
        Ci.connect(W=scipy.sparse.csr_matrix(sim_properties['reload_results']['Ci_W']))
        #print 'MEAN CE_E_W ', pylab.mean(Ce.todense().flatten())

    elif sim_properties['weight_dist'] == 'uniform':
        Ce = Connection(exc_cells, all_cells, 'ge', weight=g_exc, max_delay=max_delay,
                sparseness=epsilon,
                delay     =0.1*ms)
        Ci = Connection(inh_cells, all_cells, 'gi', weight=g_inh, max_delay=max_delay,
                sparseness=epsilon,
                delay     =0.1*ms)
    elif sim_properties['weight_dist'] == 'lognormal':
        Ce = Connection(exc_cells, all_cells, 'ge',
                weight=lambda:numpy.random.lognormal(numpy.log(g_exc),g_cov)*siemens, max_delay=max_delay,
                sparseness=epsilon,
                delay     =0.1*ms)
        Ci = Connection(inh_cells, all_cells, 'gi',
                weight=lambda:numpy.random.lognormal(numpy.log(g_inh),g_cov)*siemens, max_delay=max_delay,
                sparseness=epsilon,
                delay     =0.1*ms)
    elif sim_properties['weight_dist'] == 'seperate_groups_uniform':
        Ce_inputs = []
        for i in xrange(num_input_groups):
            Ce_input = Connection(all_cells[i*input_group_pop:((i+1)*input_group_pop)], all_cells[i*input_group_pop:((i+1)*input_group_pop)], 'ge', weight=g_exc, max_delay=max_delay,
                    sparseness=epsilon,
                    delay     =0.1*ms)
            Ce_inputs.append(Ce_input)
        Ce = Connection(exc_cells[((num_input_groups)*input_group_pop):],
                all_cells[((num_input_groups)*input_group_pop):], 'ge',
                weight=g_exc, max_delay=max_delay,
                sparseness=epsilon,
                delay     =0.1*ms)
        Ci = Connection(inh_cells, all_cells, 'gi', weight=g_inh, max_delay=max_delay,
                sparseness=epsilon,
                delay     =0.1*ms)

    elif sim_properties['weight_dist'] == 'spatial_s_lat':
        print 'N', len(exc_cells), len(all_cells)
        Ce = Connection(exc_cells, all_cells, 'ge', weight=g_exc, max_delay=max_delay,
                    sparseness=lambda i, j : probas(i, j, exc_cells.position, all_cells.position,epsilon),
                    delay     =lambda i, j : delays(i, j, exc_cells.position, all_cells.position,velocity))
        Ci = Connection(inh_cells, all_cells, 'gi', weight=g_inh, max_delay=max_delay,
                    sparseness=lambda i, j : probas(i, j, inh_cells.position, all_cells.position,epsilon),
                    delay     =lambda i, j : delays(i, j, inh_cells.position, all_cells.position,velocity))


    for i in xrange(num_input_groups):
        Cin = PoissonInput(all_cells[i*input_group_pop:((i+1)*input_group_pop)],1,input_rates[i],g_ext)
        C_inputs.append(Cin)

    if sim_properties['main_group_input_type'] == 'uniform':
        Cext = PoissonInput(all_cells[input_group_pop*num_input_groups:],1,ext_rate,g_ext)
    elif sim_properties['main_group_input_type'] == 'normal':
        norm = numpy.random.normal(ext_mean,ext_std,n_cells*2)
        Cext = PoissonInput(all_cells[input_group_pop*num_input_groups:],1,1*Hz*norm[norm>0][:(n_cells-input_group_pop*num_input_groups)],g_ext)
    elif sim_properties['main_group_input_type'] == 'lognormal':
        Cext = PoissonInput(all_cells[input_group_pop*num_input_groups:],1,1*Hz*(numpy.random.lognormal(numpy.log(ext_mean),ext_std,n_cells)),g_ext)

    if sim_properties['reload_sim']:
        Cext.rate = sim_properties['reload_results']['inputs_original']

    print "Cext_std, diffusing, input_type = ", ext_std, sim_properties['diffusing'], sim_properties['main_group_input_type']

    print "--> mean probability from excitatory synapses:", Ce.W.getnnz()/float(n_exc*n_cells) * 100, "%"
    print "--> mean probability from inhibitory synapses:", Ci.W.getnnz()/float((n_cells - n_exc)*n_cells) * 100, "%"


    print "Setting the recorders..."
    #v_exc   = RecentStateMonitor(exc_cells, 'v', record=True, clock=simulation_clock)
    global Ca_all, nNOS_all, NO_all, s_all_short, s_all, SM_ca, SM_all
    #Ca_all   = RecentStateMonitor(all_cells, 'Ca', record=True, clock=simulation_clock)
    #NO_all   = StateMonitor(all_cells, 'NO', record=True, clock=modulation_clock)
    #nNOS_all   = RecentStateMonitor(all_cells, 'nNOS', record=True, clock=simulation_clock)
    s_all   = SpikeCounter(all_cells)
    s_all_short = SpikeCounter(all_cells)
    ### Spike Monitor to trigger postsynaptic Ca influx after each spike
    if sim_properties['nNOS_inhibitory_only'] and not sim_properties['nNOS_ampa_activated']: ## THIS IS DONE INEFFICIENTLY: since nNOS and Ca values are still tracked for
    ## excitatory cells: best would be to have different eqs for exc cells
        SM_ca = SpikeMonitor(inh_cells, function=ca_influx,record=False)    #try doing this with e.g. only inhibitory neurons: i.e. only interneurons express nNOS
    elif not sim_properties['nNOS_ampa_activated']:
        SM_ca = SpikeMonitor(all_cells, function=ca_influx,record=False)

    SM_all = SpikeMonitor(all_cells)
    SM_sample_exc = SpikeMonitor(exc_cells[0:sim_properties['sample_neurons_N']])
    SM_sample_inh = SpikeMonitor(inh_cells[0:sim_properties['sample_neurons_N']])
    if sim_properties['record_pop_rate']:
        rates_exc = PopulationRateMonitor(exc_cells)
        rates_inh = PopulationRateMonitor(inh_cells)

    global raster_rest
    if not sim_properties['single_group']:
        global raster_A, raster_B
        raster_A = SpikeMonitor(all_cells[0:size_display_group])
        raster_B = SpikeMonitor(all_cells[input_group_pop:input_group_pop+size_display_group])
        raster_rest = SpikeMonitor(all_cells[input_group_pop*num_input_groups:(input_group_pop*num_input_groups)+size_display_group])
    else:
        print 'OVERWRITING num_input_groups: setting to 0'
        num_input_groups = 0
        raster_rest = SpikeMonitor(all_cells[:size_display_group*3])

    global Vm_B, Vm_samplecell
    Vm_B = RecentStateMonitor(all_cells[input_group_pop:input_group_pop+size_display_group],'v', record=True, clock=simulation_clock)
    Vm_samplecell = RecentStateMonitor(all_cells[0],'v', record=True, clock=simulation_clock)

    global Vt_samples, rate_samples_exc, rate_samples_inh, CV_samples_exc, CV_samples_inh
    Vt_samples = []
    rate_samples_exc = []
    rate_samples_inh = []
    CV_samples_exc = []
    CV_samples_inh = []

    if sim_properties['interactive_plotting']:
        ion() # To enter the interactive mode


    print "Initializing the plots..."


    global fig3, fig4, fig5, fig6, fig7, fig8, fig9, fig10, fig11, fig12, fig13
    figure_container = pylab.figure(figsize=(12,9))

    fig3    = pylab.subplot(3,4,3)
    im      = fig3.scatter(all_cells.position[:, 0], all_cells.position[:, 1], c=[0]*n_cells)
    fig4    = pylab.subplot(3,4,4)
    fig4.set_ylabel("NO")

    fig5    = pylab.subplot(3,4,2)
    fig5.set_ylabel("nNOS")

    fig6  = pylab.subplot(3,4,5)
    im = fig6.scatter(all_cells.position[:, 0], all_cells.position[:, 1],c=[0]*n_cells)

    fig7  = pylab.subplot(3,4,6)
    fig7.set_ylabel("max/min Vt")

    fig8 = pylab.subplot(3,4,7)
    fig8.set_ylabel("Vm_B")

    fig9 = pylab.subplot(3,4,8)
    fig9.set_ylabel("rate hist")

    fig10 = pylab.subplot(3,4,9)
    fig10.set_ylabel("std rates")

    fig11 = pylab.subplot(3,4,10)
    fig11.set_ylabel("Vt_settle hist")

    fig12 = pylab.subplot(3,4,11)
    fig12.set_ylabel("Vt hist")

    fig13 = pylab.subplot(3,4,12)

    manager = pylab.get_current_fig_manager()

    NOsettletime = sim_properties['NOsettletime']
    autothresh_time = sim_properties['autothresh_time']
    global sim_step
    sim_step = sim_properties['sim_step']
    sampletime = sim_properties['sampletime']
    settletime = sim_properties['settletime']

    print 'Nsteps:', int((NOsettletime/sim_step))

    print "Running network ..."

    tempd = sim_properties['diffusing']
    tempm = sim_properties['modulating']
    tempmn = sim_properties['NO_monitor']

    #diffusing = False
    sim_properties['modulating'] = False
    sim_properties['NO_monitor'] = False

    run(settletime)

    sim_properties['NO_monitor'] = tempmn
    #diffusing = tempd
   #sim_properties['modulating'] = tempm
    s_all.count = numpy.zeros(n_cells)
    s_all_short.count = numpy.zeros(n_cells)
 #   all_cells.nNOS = numpy.zeros(n_cells)

    print "Finished settling network ..."

    global fig1
    fig1 = subplot(3,4,1)

    global max_std, max_cov, max_initrate
    max_std = 0
    max_cov = 0
    max_initrate = 0


    if sim_properties['do_discrim_task']:
        print 'discrim task type (difference): ', sim_properties['discrim_stim_difference']

   ### automatically computes maxthresh_NO/NO_target depending on natural NO
   ## concentration obversved in network at basal state
    if sim_properties['maxthresh_NO_auto']:
        for C_input in C_inputs:
            C_input.rate = ext_rate

        Cext_rate_temp = Cext.rate
        if sim_properties['main_group_input_type'] == 'normal' or sim_properties['main_group_input_type'] == 'lognormal':
            ### CAN'T DO THIS IF WANT PROPER NO_PREMOD
            if sim_properties['variable_targets'] and not sim_properties['diffusing']:
                print 'NORMAL INPUT PREMOD FOR TARGET_VAR'
                norm_premod = numpy.random.normal(sim_properties['var_targets_premod_mean'],sim_properties['var_targets_premod_std'],n_cells*2)
                Cext.rate = 1*Hz*norm_premod[norm_premod>0][:(n_cells-input_group_pop*num_input_groups)]
            else:
                Cext.rate = ext_rate*numpy.ones(n_cells)
        elif sim_properties['main_group_input_type'] == 'uniform' and sim_properties['variable_targets'] and not sim_properties['diffusing']:
            print 'NORMAL INPUT PREMOD FOR TARGET_VAR'
            norm_premod = numpy.random.normal(sim_properties['var_targets_premod_mean'],sim_properties['var_targets_premod_std'],n_cells*2)
            Cext.rate = 1*Hz*norm_premod[norm_premod>0][:(n_cells-input_group_pop*num_input_groups)]

        for time in xrange(int(pylab.ceil(autothresh_time/sim_step))):
            run(sim_step)

            interactive_step(time,True,False,sim_properties['interactive_plotting'])

            try:
                print "autothresh progress:", float(time/float((autothresh_time/sim_step)))
            except: pass
            if sim_properties['interactive_plotting']:
                manager.canvas.draw()

        feed_NOgrid_to_states()
        mean_NO = numpy.mean(all_cells.NO)
        if sim_properties['scale_to_mean']:
            if sim_properties['scale_to_rate']:
                feed_NOgrid_to_states()
                rates_autothresh = multiply(s_all.count,1000.*msecond/sampletime)
                mean_rate = numpy.mean(rates_autothresh)
                meanrate_idx = []
                for i in xrange(sim_properties['Nsample_rate_scale']):
                    meanrate_idx.append(numpy.argmin(numpy.abs(rates_autothresh-mean_rate)))
                    rates_autothresh[meanrate_idx[-1]]=9999.0
                target_NO = numpy.mean(all_cells.NO[meanrate_idx])
                maxthresh_NO = numpy.mean(all_cells.NO[meanrate_idx])
            else:
                target_NO = numpy.mean(all_cells.NO)*sim_properties['NOthresh_scaling']
                maxthresh_NO = numpy.mean(all_cells.NO)*sim_properties['NOthresh_scaling']
        else:
            if sim_properties['scale_to_rate']:
                feed_NOgrid_to_states()
                rates_autothresh = multiply(s_all.count,1000.*msecond/sampletime)
                mean_rate = numpy.mean(rates_autothresh)
                meanrate_idx = []
                for i in xrange(sim_properties['Nsample_rate_scale']):
                    meanrate_idx.append(numpy.argmin(numpy.abs(rates_autothresh-mean_rate)))
                    rates_autothresh[meanrate_idx[-1]]=9999.0
                target_NO = numpy.max(all_cells.NO[meanrate_idx])
                maxthresh_NO = numpy.max(all_cells.NO[meanrate_idx])
            else:
                target_NO = numpy.max(all_cells.NO)*sim_properties['NOthresh_scaling']
                maxthresh_NO = numpy.max(all_cells.NO)*sim_properties['NOthresh_scaling']

        print 'mean NO :', mean_NO

        if sim_properties['mod_targetNOconc']:
            print 'target_NO :', target_NO
        else:
            print 'maxthresh_NO :', maxthresh_NO

        for i in xrange(len(C_inputs)):
            C_inputs[i].rate = input_rates[i]

        if sim_properties['main_group_input_type'] == 'normal' or sim_properties['main_group_input_type'] == 'lognormal':
            Cext.rate = Cext_rate_temp

        rates_premod = multiply(s_all.count,1000.*msecond/autothresh_time)
        NO_premod = all_cells.NO
        s_all.count = numpy.zeros(n_cells)

        if sim_properties['variable_targets']:
            ## assigning targets from NO_premod distribution
            if sim_properties['mod_targetNOconc']:
                global NO_targets
                #NO_targets = NO_premod + target_NO*sim_properties['target_variance']*(numpy.random.rand(n_cells)-0.5)
                #### TRYING TO USE PREMOD RATE DISTRIBUTIONS INSTEAD
                print 'using targets from NO_premod'
                NO_targets = pylab.mean(NO_premod)*((rates_premod)/pylab.mean(rates_premod))
            else:
                global maxthresholds_NO
                maxthresholds_NO = NO_premod + maxthresh_NO*sim_properties['target_variance']*(numpy.random.rand(n_cells)-0.5)



    sim_properties['modulating'] = tempm


    if sim_properties['reload_sim']:
        if sim_properties['variable_targets']:
            global NO_targets
            print 'using relaoded NO targets'
            NO_targets = sim_properties['reload_results']['NO_targets']
        else:
            if sim_properties['mod_targetNOconc']:
                target_NO = sim_properties['reload_results']['NO_target']


    if sim_properties['init_homogenous_inputs']:
        for C_input in C_inputs:
            C_input.rate = ext_rate
        if not sim_properties['maxthresh_NO_auto']:
            for time in xrange(int((NOsettletime/sim_step))):
                run(sim_step)

                interactive_step(time,True,False,sim_properties['interactive_plotting'])

                try:
                    print "init_stage_1 progress:", float(time/float(((NOsettletime)/sim_step)))
                except: pass
                if sim_properties['interactive_plotting']:
                    manager.canvas.draw()

            s_all.count = numpy.zeros(n_cells)

        for time in xrange(int((sampletime/sim_step))):
            run(sim_step)

            interactive_step(time,True,True,sim_properties['interactive_plotting'])

            try:
                print "init_stage_2 progress:", float(time/float(((sampletime)/sim_step)))
            except: pass
            if sim_properties['interactive_plotting']:
                manager.canvas.draw()

        rates_init = multiply(s_all.count,1000.*msecond/sampletime)
        Vt_init_hist = all_cells.Vt.copy()
        s_all.count = numpy.zeros(n_cells)

        for i in xrange(len(C_inputs)):
            C_inputs[i].rate = input_rates[i]


    Cext_rate_original = Cext.rate


    for time in xrange(int((NOsettletime/sim_step))):
        if sim_properties['subgroup_presentations_during_homeostasis']:
            print 'new subgroup pattern during HIP'
            y_1 = Cext_rate_original
            subgroup_pattern = []
            for sg in xrange(int(n_cells/sim_properties['discrim_pattern_size'])):
                subgroup_pattern = subgroup_pattern + [(numpy.random.normal(0,ext_std))]*sim_properties['discrim_pattern_size']
            y_2 = numpy.add(y_1,subgroup_pattern).clip(min=0.1)

            Cext.rate = 1*Hz*numpy.array(list(y_2))

        run(sim_step)

        Cext.rate = Cext_rate_original

        interactive_step(time,False,False,sim_properties['interactive_plotting'])

        try:
            print "pert_stage_1 progress:", float(time/float(((NOsettletime)/sim_step)))
        except: pass
        if sim_properties['interactive_plotting']:
                manager.canvas.draw()

    s_all.count = numpy.zeros(n_cells)

    for time in xrange(int((sampletime/sim_step))):
        if sim_properties['subgroup_presentations_during_homeostasis']:
            print 'new subgroup pattern during HIP'
            y_1 = Cext_rate_original
            subgroup_pattern = []
            for sg in xrange(int(n_cells/sim_properties['discrim_pattern_size'])):
                subgroup_pattern = subgroup_pattern + [(numpy.random.normal(0,ext_std))]*sim_properties['discrim_pattern_size']
            y_2 = numpy.add(y_1,subgroup_pattern).clip(min=0.1)

            Cext.rate = 1*Hz*numpy.array(list(y_2))

        run(sim_step)

        Cext.rate = Cext_rate_original

        if time == range(int((sampletime/sim_step)))[-1]:
            interactive_step(time,False,True,True)
        else:
            interactive_step(time,False,True,sim_properties['interactive_plotting'])
        try:
            print "pert_stage_2 progress:", float(time/float(((sampletime)/sim_step)))
        except: pass
        if sim_properties['interactive_plotting']:
                manager.canvas.draw()


    print 'finished main sim, doing options'

    rates_pert = multiply(s_all.count,1000.*msecond/sampletime)
    Vt_pert = all_cells.Vt

    if sim_properties['reconfigure_inputs']:
        s_all.count = numpy.zeros(n_cells)

        Cext_rate_original = Cext.rate
        if sim_properties['main_group_input_type'] == 'normal':
            norm = numpy.random.normal(ext_mean,ext_std,n_cells*2)
            Cext.rate = 1*Hz*norm[norm>0][:n_cells]
        elif sim_properties['main_group_input_type'] == 'lognormal':
            Cext.rate = 1*Hz*abs(numpy.random.lognormal(numpy.log(ext_mean),ext_std,n_cells))
        Cext_rate_diff = Cext_rate_original - Cext.rate
        Cext_rate_reconf = Cext.rate

        sim_properties['modulating'] = False
        tempmon = sim_properties['NO_monitor']
        sim_properties['NO_monitor'] = False

        sim_properties['diffusing'] = False

        for time in xrange(int(sampletime/sim_step)):
            run(sim_step)
            interactive_step(time,False,True,sim_properties['interactive_plotting'])
            try:
                print "reconfigured inputs progress:", float(time/float(sampletime/sim_step))
            except: pass
            if sim_properties['interactive_plotting']:
                    manager.canvas.draw()

        rates_reconf = multiply(s_all.count,1000*msecond/sampletime)
        sim_properties['modulating'] = tempm

        sim_properties['diffusing'] = tempd
        sim_properties['NO_monitor'] = tempmon

    if sim_properties['get_transfer_funs']:
        if sim_properties['main_group_input_type'] == 'normal':
            norm = numpy.random.normal(ext_mean,ext_std,n_cells*2)
            Cext.rate = 1*Hz*norm[norm>0][:n_cells]
            assert len(Cext.rate)==n_cells, "not enough rate inputs above zero"
        elif sim_properties['main_group_input_type'] == 'lognormal':
            Cext.rate = 1*Hz*abs(numpy.random.lognormal(numpy.log(ext_mean),ext_std,n_cells))

        Cext_rate_transf = Cext.rate

        sim_properties['modulating'] = False
        sim_properties['diffusing'] = False
        tempmon = sim_properties['NO_monitor']
        sim_properties['NO_monitor'] = False

        ### COMMENTING BELOW LINE FOR REMOTE CAPABILITY, USE IF WANT TRANSFER FUNCTIONS!!!
        #transfer_funs = {i:[] for i in range(n_cells)}
        ## replacing above line, for python 2.6 compability
        print 'WARNING: tranfser functions probably wont work; uncomment above line and run on abu/mbp to use tranf_fun'

        transfer_funs = dict.fromkeys(range(n_cells),[])


        for i in xrange(sim_properties['N_transf_presentations']):
            print 'presenting transfer samples: ', float(i)/sim_properties['N_transf_presentations']
            s_all.count = numpy.zeros(n_cells)
            if not i == 0:
                Cext.rate = 1*Hz*numpy.array((list(Cext_rate_transf[-i:])+list(Cext_rate_transf[:-i])))
            run(sim_properties['t_transf_presentation'])
            rates_transf = multiply(s_all.count,1000*msecond/sim_properties['t_transf_presentation'])
            #temp_transf = []
            for j in xrange(n_cells):
                transfer_funs[j].append((Cext.rate[j],rates_transf[j]))
                #temp = transfer_funs[j]
                #temp.append((Cext.rate[j],rates_transf[j]))

            #transfer_funs[i] = temp_transf

        sim_properties['modulating'] = tempm

        sim_properties['diffusing'] = tempd
        sim_properties['NO_monitor'] = tempmon

    if sim_properties['restore_inputs']:
        s_all.count = numpy.zeros(n_cells)
        for C_input in C_inputs:
            C_input.rate = ext_rate
        for time in xrange(int((NOsettletime/sim_step))):
            run(sim_step)

            interactive_step(time,True,False,sim_properties['interactive_plotting'])

            try:
                print "restore_stage_1 progress:", float(time/float(((NOsettletime)/sim_step)))
            except: pass
            if sim_properties['interactive_plotting']:
                manager.canvas.draw()

        s_all.count = numpy.zeros(n_cells)

        for time in xrange(int((sampletime/sim_step))):
            run(sim_step)

            interactive_step(time,True,True,sim_properties['interactive_plotting'])

            try:
                print "restore_stage_2 progress:", float(time/float(((sampletime)/sim_step)))
            except: pass
            if sim_properties['interactive_plotting']:
                manager.canvas.draw()

        rates_restore = multiply(s_all.count,1000.*msecond/sampletime)
        Vt_restore_hist = all_cells.Vt.copy()

    if sim_properties['do_discrim_task']:
        sim_properties['modulating'] = False
        sim_properties['diffusing'] = False
        tempmon = sim_properties['NO_monitor']
        sim_properties['NO_monitor'] = False

        responses_Nstim_1 = []
        responses_Nstim_2 = []
        patterns_Nstim_1 = []
        patterns_Nstim_2 = []

        Cext_rate_original = Cext.rate


        for stim_n in xrange(sim_properties['N_discrim_stim']):
            norm = numpy.random.normal(ext_mean,ext_std,n_cells*2)

            ## SINUDOIDAL INPUT
            #x = pylab.arange(0,1,(float(1)/n_cells))*4*3.14
            #y_1 = 1*Hz*ext_mean*(abs(pylab.sin(x)+1))
            #y_2 = 1*Hz*ext_mean**(abs(pylab.cos(x)+1))

            ## RANDOM INPUTS, DIFFERENT MEANS
            #y_1 = numpy.random.normal(ext_mean*(1.0+sim_properties['discrim_stim_difference']),ext_std,n_cells*2)
            #y_1 = y_1[y_1>0][:n_cells]

            #y_2 = numpy.random.normal(ext_mean*(1.0-sim_properties['discrim_stim_difference']),ext_std,n_cells*2)
            #y_2 = y_2[y_2>0][:n_cells]

            ## RANDOM FIRST INPUT, ADD NOISE TO FIRST TO GET SECOND
#            y_1 = numpy.random.normal(ext_mean,ext_std,n_cells*2)
#            y_1 = y_1[y_1>0][:n_cells]

            if sim_properties['reconfigure_inputs']:
                print 'cant use original rates as pattern'
            y_1 = Cext_rate_original

            if sim_properties['discrim_stim_difference'] == 'independent':
                y_2 = numpy.random.normal(ext_mean,ext_std,n_cells*2)
                y_2 = y_2[y_2>0][:n_cells]
            else:
                ## pattern_size < n_cells
                pattern_idx = numpy.random.randint(0,n_cells,sim_properties['discrim_pattern_size'])
                temp_pattern = numpy.zeros(n_cells)
                temp_pattern[pattern_idx] = scipy.random.normal(0,sim_properties['discrim_stim_difference']*ext_std,sim_properties['discrim_pattern_size'])
                y_2 = numpy.abs(numpy.array(y_1) + temp_pattern)

                ## pattern_size = n_cells
                if sim_properties['discrim_pattern_size'] == n_cells:
                    y_2 = numpy.abs(numpy.array(y_1) + numpy.random.normal(0,sim_properties['discrim_stim_difference']*ext_std,n_cells))
                #y_2 = numpy.abs(numpy.array(y_1) + numpy.random.normal(0,sim_properties['discrim_stim_difference']*ext_std,n_cells))


            #y_2 = numpy.random.normal(ext_mean*(1.0-sim_properties['discrim_stim_difference']),ext_std,n_cells*2)
            #y_2 = y_2[y_2>0][:n_cells]


            pattern_1 = list(y_1)
            #pattern_2 = list(y[n_cells/2:])+list(y[:n_cells/2])
            pattern_2 = list(y_2)
            responses_1 = []
            responses_2 = []


            print 'presenting stim_N: ', stim_n

            for i in xrange(sim_properties['N_discrim_presentations']):
                #print 'presenting discrim patterns trial #: ', float(i)/sim_properties['N_discrim_presentations']

                s_all.count = numpy.zeros(n_cells)

                Cext.rate = 1*Hz*numpy.array(pattern_1)
                run(sim_properties['t_discrim_presentation'])
                rates_pres = multiply(s_all.count,1000*msecond/sim_properties['t_discrim_presentation'])
                responses_1.append(rates_pres)


                s_all.count = numpy.zeros(n_cells)

                Cext.rate = 1*Hz*numpy.array(pattern_2)
                run(sim_properties['t_discrim_presentation'])
                rates_pres = multiply(s_all.count,1000*msecond/sim_properties['t_discrim_presentation'])
                responses_2.append(rates_pres)


            patterns_Nstim_1.append(pattern_1)
            patterns_Nstim_2.append(pattern_2)
            responses_Nstim_1.append(responses_1)
            responses_Nstim_2.append(responses_2)

        sim_properties['modulating'] = tempm

        sim_properties['diffusing'] = tempd
        sim_properties['NO_monitor'] = tempmon

    if sim_properties['get_pattern_responses']:
        sim_properties['modulating'] = False
        sim_properties['diffusing'] = False
        tempmon = sim_properties['NO_monitor']
        sim_properties['NO_monitor'] = False

        Cext_rate_original = Cext.rate

        if sim_properties['reconfigure_inputs']:
                print 'cant use original rates as pattern'
        y_1 = Cext_rate_original

        pattern_responses = []
        pattern_stims = []

        for i in xrange(sim_properties['N_discrim_stim']):
            print 'presenting pattern : ', float(i)/sim_properties['N_discrim_stim']

            if sim_properties['discrim_stim_difference'] == 'independent':
                    y_2 = numpy.random.normal(ext_mean,ext_std,n_cells*2)
                    y_2 = y_2[y_2>0][:n_cells]

            elif sim_properties['discrim_stim_difference'] == 'uniform':
                    print 'UNIFORM PATTERNS'
                    y_2 = numpy.random.uniform(ext_mean-ext_std*sim_properties['uniform_pattern_range'],ext_mean+ext_std*sim_properties['uniform_pattern_range'],n_cells)

            elif sim_properties['discrim_stim_difference'] == 'subgroup_pattern':
                    print 'SUBGROUP PATTERNS'
                    #y_2 = numpy.random.normal(ext_mean,ext_std,n_cells*2)
                    #temp_pattern = y_2[y_2>0][:n_cells]

                    ## uniform, randomly drawn base patten
                    #temp_pattern = [abs(numpy.random.normal(ext_mean,ext_std))]*n_cells

                    #learned_group = int(sim_properties['discrim_pattern_size'])*[sim_properties['subgroup_pattern_mag']] + int(n_cells-sim_properties['discrim_pattern_size'])*[0]

                    if not sim_properties['pattern_provided']:
                        subgroup_pattern = [(numpy.random.normal(0,ext_std))]*sim_properties['discrim_pattern_size'] + int(n_cells-sim_properties['discrim_pattern_size'])*[0]
                    else:
                        print 'provided pattern :',sim_properties['pattern_passed'][i]
                        subgroup_pattern = [sim_properties['pattern_passed'][i]]*sim_properties['discrim_pattern_size'] + int(n_cells-sim_properties['discrim_pattern_size'])*[0]

                    y_2 = numpy.add(y_1,subgroup_pattern).clip(min=0.1)
            else:
                ## pattern_size < n_cells
                pattern_idx = numpy.random.randint(0,n_cells,sim_properties['discrim_pattern_size'])
                temp_pattern = numpy.zeros(n_cells)
                temp_pattern[pattern_idx] = scipy.random.normal(0,sim_properties['discrim_stim_difference']*ext_std,sim_properties['discrim_pattern_size'])
                y_2 = numpy.abs(numpy.array(y_1) + temp_pattern)

                ## pattern_size = n_cells
                if sim_properties['discrim_pattern_size'] == n_cells:
                    y_2 = numpy.abs(numpy.array(y_1) + numpy.random.normal(0,sim_properties['discrim_stim_difference']*ext_std,n_cells))
                #y_2 = numpy.abs(numpy.array(y_1) + numpy.random.normal(0,sim_properties['discrim_stim_difference']*ext_std,n_cells))

                #y_2 = numpy.random.normal(ext_mean*(1.0-sim_properties['discrim_stim_difference']),ext_std,n_cells*2)
                #y_2 = y_2[y_2>0][:n_cells]

            pattern_stim = list(y_2)
            pattern_stims.append(pattern_stim)

            s_all.count = numpy.zeros(n_cells)

            Cext.rate = 1*Hz*numpy.array(pattern_stim)
            run(sim_properties['t_discrim_presentation'])
            rates_pres = multiply(s_all.count,1000*msecond/sim_properties['t_discrim_presentation'])
            pattern_responses.append(rates_pres)

        sim_properties['modulating'] = tempm

        sim_properties['diffusing'] = tempd
        sim_properties['NO_monitor'] = tempmon

    if sim_properties['present_new_input_stats']:
        s_all.count = numpy.zeros(n_cells)

        Cext_rate_original = Cext.rate

        if sim_properties['scale_current_input_stats']:
            Cext_scaled = Cext.rate
            temp_mean = pylab.mean(Cext_scaled)
            temp_scale = sim_properties['new_ext_std']/pylab.std(Cext_scaled)
            for i in xrange(len(Cext_rate_original)):
                Cext_scaled[i] = sim_properties['new_ext_mean'] + (Cext_scaled[i]-temp_mean)*temp_scale
            Cext_scaled = pylab.clip(Cext_scaled,0.0,numpy.max(Cext_scaled))
            Cext.rate = Cext_scaled
        else:
            if sim_properties['new_main_group_input_type'] == 'uniform':
                Cext = PoissonInput(all_cells[input_group_pop*num_input_groups:],1,sim_properties['new_ext_rate'],g_ext)
            if sim_properties['new_main_group_input_type'] == 'normal':
                norm = numpy.random.normal(sim_properties['new_ext_mean'],sim_properties['new_ext_std'],n_cells*2)
                Cext.rate = 1*Hz*norm[norm>0][:n_cells]
            elif sim_properties['new_main_group_input_type'] == 'lognormal':
                Cext.rate = 1*Hz*abs(numpy.random.lognormal(numpy.log(sim_properties['new_ext_mean']),sim_properties['new_ext_std'],n_cells))

        Cext_rate_new_input_stats_diff = Cext_rate_original - Cext.rate
        Cext_rate_new_input_stats = Cext.rate

        sim_properties['modulating'] = False
        tempmon = sim_properties['NO_monitor']
        sim_properties['NO_monitor'] = False

        sim_properties['diffusing'] = False

        for time in xrange(int(sampletime/sim_step)):
            run(sim_step)
            interactive_step(time,False,True,sim_properties['interactive_plotting'])
            try:
                print "new inputs progress:", float(time/float(sampletime/sim_step))
            except: pass
            if sim_properties['interactive_plotting']:
                    manager.canvas.draw()

        rates_new_input_stats = multiply(s_all.count,1000*msecond/sampletime)
        sim_properties['modulating'] = tempm

        sim_properties['diffusing'] = tempd
        sim_properties['NO_monitor'] = tempmon

    if sim_properties['get_network_dynamic_range']:
        print 'testing dynamic range'
        s_all.count = numpy.zeros(n_cells)

        sim_properties['modulating'] = False
        tempmon = sim_properties['NO_monitor']
        sim_properties['NO_monitor'] = False

        sim_properties['diffusing'] = False

        network_input_range = sim_properties['network_input_range']
        dynamic_range_sampletime = sim_properties['dynamic_range_sampletime']
        mean_network_response = []

        for input_i in network_input_range:
            print 'testing input ', input_i
            Cext = PoissonInput(all_cells[input_group_pop*num_input_groups:],1,input_i,g_ext)
            for time in xrange(int(dynamic_range_sampletime/sim_step)):
                run(sim_step)
                interactive_step(time,False,True,sim_properties['interactive_plotting'])
                if sim_properties['interactive_plotting']:
                        manager.canvas.draw()
            mean_network_response.append(multiply(s_all.count,1000*msecond/dynamic_range_sampletime))
            s_all.count = numpy.zeros(n_cells)

        sim_properties['modulating'] = tempm

        sim_properties['diffusing'] = tempd
        sim_properties['NO_monitor'] = tempmon

        if sim_properties['just_save_dynamic_range']:
            results = mean_network_response
            #resultfile = open(logdir+'/network_dynamic_range.pickle','w')
            common.save_pickle_safe(logdir+'_network_dynamic_range.pickle',results)
            #pickle.dump(results, resultfile)
            #resultfile.close()
            return results


    if sim_properties['do_orientation_tasks']:
        sim_properties['modulating'] = False
        sim_properties['diffusing'] = False
        tempmon = sim_properties['NO_monitor']
        sim_properties['NO_monitor'] = False

        Cext_rate_original = Cext.rate

        if sim_properties['reconfigure_inputs']:
                print 'cant use original rates as pattern'
        y_1 = Cext_rate_original

        orientation_responses = []
        orientation_stims = []

        pref_orientations = numpy.arange(0,1,1.0/n_exc)  ## just mulitply by 2pi to get actual orientation

        presented_orientations = numpy.random.rand(sim_properties['N_discrim_stim'])

        for i in xrange(sim_properties['N_discrim_stim']):
            print 'presenting orientation #, degree ', float(i)/sim_properties['N_discrim_stim'], presented_orientations[i]

            y_2 = numpy.zeros(n_exc)
            for j in xrange(n_exc):
                y_2[j] = scipy.stats.norm.pdf(numpy.min([numpy.abs(presented_orientations[i]-pref_orientations[j]),1.0+presented_orientations[i]-pref_orientations[j]]),0.0,sim_properties['orientation_width'])
            #print 'orientation stim #-2, ', y_2
            y_2*=sim_properties['subgroup_pattern_mag']
            y_2 += sim_properties['no_orient_pref_bias']

            #print 'orientation stim #-1', y_2

            y_2 = y_2.clip(min=0.1)

            pattern_stim = list(y_2)+list(y_1)[n_exc:]

            orientation_stims.append(pattern_stim)

            #print 'orientation stim ', pattern_stim

            s_all.count = numpy.zeros(n_cells)

            Cext.rate = 1*Hz*numpy.array(pattern_stim)
            run(sim_properties['t_discrim_presentation'])
            rates_pres = multiply(s_all.count,1000*msecond/sim_properties['t_discrim_presentation'])
            orientation_responses.append(rates_pres)

        sim_properties['modulating'] = tempm

        sim_properties['diffusing'] = tempd
        sim_properties['NO_monitor'] = tempmon


    if sim_properties['interactive_plotting']:
        ioff() # To leave the interactive mode
        show() # and wait for user to close the window before shutting down

    if sim_properties['savefig_id']:
        try:
            figure_container.savefig(logdir+'/plots.png')

            #rast_fig = pylab.figure(figsize=(10,10))
            #if not sim_properties['single_group'] and show:
            #    raster_plot(raster_A,raster_B,raster_rest,showgrouplines=True,spacebetweengroups=0.01)
            #elif show:
            #    raster_plot(raster_rest,showgrouplines=True,spacebetweengroups=0.01)

            #fig3.savefig(logdir + '/NOconc_map.png')
    #        rast_fig.savefig(logdir + '/raster.png')
    #        fig6.savefig(logdir + '/Vt_map.png')
        except:
            if sim_properties['remote']:
                os.system('mkdir /tmp/' + now.strftime("%Y%m%d%H%M"))
                logdir = '/tmp/' + now.strftime("%Y%m%d%H%M")
                figure_container.savefig(logdir + '/plots.png')
    if sim_properties['savefig_id']:
        figure_container.savefig(sim_properties['logdir']+ 'plot' +
            sim_properties['sim_id']  + '.png')
        figNO = pylab.figure()
        pylab.pcolor(NO_concs.transpose())
        figNO.savefig(sim_properties['logdir']+ 'NOconc_map.png' +
            sim_properties['sim_id'] + '.png')
        rast_fig = pylab.figure(figsize=(10,10))
        if not sim_properties['single_group']:
            raster_plot(raster_A,raster_B,raster_rest,showgrouplines=True,spacebetweengroups=0.01)
        elif show:
            raster_plot(raster_rest,showgrouplines=True,spacebetweengroups=0.01)
            xlims =  pylab.xlim()
            try:
                pylab.xlim(xmin=xlims[1]*.7-10000,xmax=xlims[1]*.7)
            except:
                pass
        rast_fig.savefig(sim_properties['logdir']+ 'raster.png' +
            sim_properties['sim_id'] + '.png')
        figVt = pylab.figure()
        pylab.scatter(all_cells.position[:, 0], all_cells.position[:, 1], c=all_cells.Vt)
        figVt.savefig(sim_properties['logdir']+ 'Vt_map.png' +
            sim_properties['sim_id'] + '.png')

    autocorr_e = []
    autocorr_i = []
    for i in xrange(sim_properties['sample_neurons_N']):
        if len(SM_sample_exc[i])>1:
            autocorr_e.append(autocorrelogram(SM_sample_exc[i],width=100*msecond))
        #if len(SM_sample_inh[i])>1:
        #    autocorr_i.append(autocorrelogram(SM_sample_inh[i],width=100*msecond))


    if sim_properties['dump_each_simstep']:
        print 'saving final dump'
        results = {
                'sim_pars' : sim_properties,
                'NO_threshold': maxthresh_NO,
                'NO_target': target_NO,
                'NO_concs': NO_concs,
                'cell_positions': all_cells.position,
                'CV_hist': CV_hist,
                'weights_dist': Ce.W.todense(),
                'NO_postmod': all_cells.NO,
                'Ci_W': Ci.W.todense(),
                'nNOS_dist': all_cells.nNOS,
                'thresholds': all_cells.Vt

        }
        if sim_properties['maxthresh_NO_auto']:
            results['rates_premod'] = rates_premod
            results['NO_premod'] = NO_premod
        if sim_properties['init_homogenous_inputs']:
            if sim_properties['restore_inputs']:
                results['rates_restore'] = rates_restore
                results['Vt_restore'] = Vt_restore_hist
            results['rates_init'] = rates_init
            results['Vt_init'] = Vt_init_hist
            results['rates_pert'] = rates_pert
            results['Vt_pert'] = Vt_pert
        else:
            results['rates_hist'] = rates_pert
            results['Vt_hist'] = Vt_pert

        if sim_properties['variable_targets']:
            results['NO_targets'] = NO_targets

        if sim_properties['reconfigure_inputs']:
            results['inputs_original'] = Cext_rate_original
            results['inputs_diff'] = Cext_rate_diff
            results['inputs_reconf'] = Cext_rate_reconf
            results['rates_reconf'] = rates_reconf
        else:
            results['inputs_original'] = Cext.rate

        if sim_properties['get_network_dynamic_range']:
            results['network_dynamic_range'] = mean_network_response
        if sim_properties['save_raster_spikes']:
            results['raster_spikes'] = raster_rest.it
        if sim_properties['present_new_input_stats']:
            results['rates_new_input_stats'] = rates_new_input_stats
            results['inputs_new_input_stats'] = Cext_rate_new_input_stats
            results['inputs_original'] = Cext_rate_original
            results['inputs_diff_new_inputs_stats'] = Cext_rate_new_input_stats_diff



    else:
        results = {
                    'NO_target': target_NO,
                    'CV_hist': CV_hist,
                    'weights_dist': Ce.W.todense(),
                    'NO_postmod': all_cells.NO,
                    'Ci_W': Ci.W.todense(),
                    'nNOS_dist': all_cells.nNOS,
                    'thresholds': all_cells.Vt,
                    'sim_pars' : sim_properties,
                    'ratesA_max': ratesA_max,
                    'ratesA_min': ratesA_min,
                    'ratesA_std': ratesA_std,
                    'ratesA_mean': ratesA_mean,
                    'ratesB_max': ratesB_max,
                    'ratesB_min': ratesB_min,
                    'ratesB_std': ratesB_std,
                    'ratesB_mean': ratesB_mean,
                    'rates_max': rates_max,
                    'rates_min': rates_min,
                    'rates_std': rates_std,
                    'rates_mean': rates_mean,
                    'NO_max_t': NO_max_t,
                    'NO_min_t': NO_min_t,
                    'NO_mean_t': NO_mean_t,
                    'NO_std_t': NO_std_t,
                    'Vt_mean': Vt_mean,
                    'Vt_max_t': Vt_max_t,
                    'VtA_mean': VtA_mean,
                    'VtA_std': VtA_std,
                    'VtB_mean': VtB_mean,
                    'VtB_std': VtB_std,
                    'Vt_std': Vt_std,
                    'Vt_min_t': Vt_min_t,
                    'NO_threshold': maxthresh_NO,
                    'NO_target': target_NO,
                    'NO_concs': NO_concs,
                    'cell_positions': all_cells.position,
                    'Vt_samples': Vt_samples,
                    'rate_samples_e': rate_samples_exc,
                    'rate_samples_i': rate_samples_inh,
                    'CV_mean': CV_mean,
                    'CV_samples_e': CV_samples_exc,
                    'CV_samples_i': CV_samples_inh,
                    'autocorr_e': autocorr_e,
                    'autocorr_i': autocorr_i,
                    'NO_conc': NO_concs.transpose(),
                    #'transfer_funs': transfer_funs
                    #'NO_all': NO_all
                }
        if sim_properties['record_pop_rate']:
            results['rates_exc'] = rates_exc.rate
            results['rates_inh'] = rates_inh.rate
        if sim_properties['maxthresh_NO_auto']:
            results['rates_premod'] = rates_premod
            results['NO_premod'] = NO_premod
        if sim_properties['init_homogenous_inputs']:
            if sim_properties['restore_inputs']:
                results['rates_restore'] = rates_restore
                results['Vt_restore'] = Vt_restore_hist
            results['rates_init'] = rates_init
            results['Vt_init'] = Vt_init_hist
            results['rates_pert'] = rates_pert
            results['Vt_pert'] = Vt_pert
        else:
            results['rates_hist'] = rates_pert
            results['Vt_hist'] = Vt_pert

        if sim_properties['variable_targets']:
            results['NO_targets'] = NO_targets,

        if sim_properties['reconfigure_inputs']:
            results['inputs_original'] = Cext_rate_original
            results['inputs_diff'] = Cext_rate_diff
            results['inputs_reconf'] = Cext_rate_reconf
            results['rates_reconf'] = rates_reconf
        if sim_properties['get_transfer_funs']:
            results['transfer_funs'] = transfer_funs
            results['inputs_transf'] = Cext_rate_transf
        if sim_properties['do_discrim_task']:
            #results['pattern_1'] = pattern_1
            #results['pattern_2'] = pattern_2
            #results['responses_1'] = responses_1
            #results['responses_2'] = responses_2
            results['patterns_Nstim_1'] = patterns_Nstim_1
            results['patterns_Nstim_2'] = patterns_Nstim_2
            results['responses_Nstim_1'] = responses_Nstim_1
            results['responses_Nstim_2'] = responses_Nstim_2
        if sim_properties['get_pattern_responses']:
            results['pattern_responses'] = pattern_responses
            results['pattern_stims'] = pattern_stims
            results['inputs_original'] = Cext_rate_original
        if sim_properties['get_network_dynamic_range']:
            results['network_dynamic_range'] = mean_network_response
        if sim_properties['save_raster_spikes']:
            results['raster_spikes'] = raster_rest.it
        if sim_properties['present_new_input_stats']:
            results['rates_new_input_stats'] = rates_new_input_stats
            results['inputs_new_input_stats'] = Cext_rate_new_input_stats
            results['inputs_original'] = Cext_rate_original
            results['inputs_diff_new_inputs_stats'] = Cext_rate_new_input_stats_diff
        if sim_properties['do_orientation_tasks']:
            results['presented_orientations'] = presented_orientations
            results['orientation_stims'] = orientation_stims
            results['orientation_responses'] = orientation_responses
            results['preferred_orientations'] = pref_orientations

    if sim_properties['log']:
        resultfile = open(logdir+'/results.pickle','w')
        pickle.dump(results, resultfile)
        resultfile.close()

    return results


if __name__ == '__main__':
    main({})
