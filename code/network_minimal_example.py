## This script simulates NO diffusion over a randomly connected network of excitatory
## and inhibitory integrate-and-fire neurons. The concenctration of NO is used as
## homeostatic signal by a mechanism which adjusts the firing threshold of each neuron
## depending on how the NO concentration at that point differs from a target concentration

### Copyright (C) 2014  Yann Sweeney ; yann.sweeney@ed.ac.uk

###    This program is free software: you can redistribute it and/or modify
###    it under the terms of the GNU General Public License as published by
###    the Free Software Foundation, either version 3 of the License, or
###    (at your option) any later version.

###    This program is distributed in the hope that it will be useful,
###    but WITHOUT ANY WARRANTY; without even the implied warranty of
###    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
###    GNU General Public License for more details.

###    You should have received a copy of the GNU General Public License
###    along with this program.  If not, see <http://www.gnu.org/licenses/>


## This script, meant primarily as a demo, runs the network for 50s with diffusive homeostasis
## active and random, Gaussian-distributed external Poisson spiking input. The mean firing rate
## over time and final distribution of firing rates is plotted. Feel free to change defined parameters.

## The below code can easily be incorporated into an existing brian model,
## or extended to contain other features.

## apologies for the use of global variables...


from brian import *
import pylab
import numpy
from random import shuffle
from itertools import product

#### returns random unique x,y positions for each cell
def random_pos(N,L):
    allpairs = list(product(*[xrange(int(L/ds))]*2))
    shuffle(allpairs)
    return numpy.array(allpairs[:N])
#### WARNING - unique positions not guaranteed here

def generate_pos():
    positions = numpy.zeros(shape=(n_cells,2))

    positions = random_pos(n_cells,size)

    positions[positions>=size/ds]=positions[positions>=size/ds]-size/ds
    positions[positions<0]=positions[positions<0]+size/ds

    return positions

def threshold_fun(states):
    return states >= all_cells.Vt

## this is called every timestep with a list of all neurons that spiked
def ca_influx(spikes):
     all_cells.Ca[spikes] = all_cells.Ca[spikes] + 1.0

def feed_nNOSstates_to_grid():
    NO_rates[(all_cells.position[:]/(size*ds))[:].astype(int)[:,0],(all_cells.position[:]/(size*ds))[:].astype(int)[:,1]] = all_cells.nNOS

def feed_NOgrid_to_states():
    all_cells.NO = NO_concs[(all_cells.position[:]/(size*ds))[:].astype(int)[:,0],(all_cells.position[:]/(size*ds))[:].astype(int)[:,1]]
    #print numpy.mean(all_cells.NO)

def diffuse_NOgrid(): # taken from Philippedes et al. (2000, J Neurosci). Also code snippet from http://goo.gl/18lc8
    #	This function uses a numpy expression to evaluate the derivatives in the Laplacian, and calculates u[i,j] based on ui[i,j].
    global NO_concs

    NO_diff =  1.0e-6
    NOdecay_rate = 1.0e-4

    ui = NO_concs.copy()
    nx = int(size/ds)
    ny = int(size/ds)
    leftX = arange(nx)-1
    rightX = (arange(nx)+1)%nx
    leftY = arange(ny)-1
    rightY = (arange(ny)+1)%ny

    NO_concs = ui + dt_diff*(NO_diff*((ui[rightX,:]-2*ui + ui[leftX,:])/dx2 + (ui[:,rightY] - 2*ui
        + ui[:,leftY])/dy2)+ NO_rates - NOdecay_rate*ui)

    # make sure all values for NO_conc >= 0
    NO_concs.clip(min=0,out=NO_concs)


def nodiffuse_NOgrid():
    #	This function uses a numpy expression to evaluate the derivatives in the Laplacian, and calculates u[i,j] based on ui[i,j].
    global NO_concs

    NOdecay_rate = 1.0e-4

    ui = NO_concs.copy()
    NO_concs = ui + dt_diff*(NO_rates - NOdecay_rate*ui)

    # make sure all values for NO_conc >= 0
    NO_concs.clip(min=0,out=NO_concs)


def modulate_neurons(mod_tau,target_NO):
    feed_NOgrid_to_states()

    NOcells = all_cells.NO

    all_cells.Vt[:] = all_cells.Vt[:] + 0.5*mV*(1./mod_tau)*((NOcells[:]-target_NO)/target_NO)


def run_step(time_step):
    rates = multiply(s_all_short.count,1000.*msecond/(sim_step))

    rates_min.append(min(rates))
    rates_max.append(max(rates))
    rates_std.append(std(rates))
    rates_mean.append(mean(rates))

    print 'mean,std rates ', rates_mean[-1],rates_std[-1]

    s_all_short.count = numpy.zeros(n_cells)


def main():

    ### We are setting the global timestep of the simulation
    simulation_clock=Clock(dt=0.1*ms)
    Clock(0.1 * ms)

    ### Cell parameters ###
    tau_m   = 20. * ms # Membrane time constant
    c_m     =  0.2 * nF # Capacitance
    tau_exc =  3. * ms
    tau_inh =  7. * ms
    tau_ref = 5. * ms
    El      = -80 * mV
    Ee      = 0. * mV
    Ei      = -70.* mV
    Vt_init = -50 * mV
    Vr      = -60 * mV
    Ca_tau  = 10 * msecond #tau for Ca decay
    nNOS_tau= 100 * msecond            #tau for nNOS activation from Ca conc.
    ### Hill function (nNOS =  hillfunction(Ca)) coefficients ###
    hill_K  = 1.0
    hill_n  = 3.0
    mod_tau = 1500

    modulating = True
    diffusing = True #set =False if you want to look at non-diffusive HIP (may have to change target_NO to a suitable value if so)


    ### Equation for a Conductance-based IAF ####
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

    eqs = eqs_iaf + eqs_Ca + eqs_nNOS

    global n_cells, n_exc, size
    n_cells       = 1000
    n_exc         = int(0.8*n_cells)
    size          = 1.0 # Size of the network - in mm (roughly... it is a torus)

    epsilon = min(float(float(50.0)/n_cells),1.0)
    g_exc         = 5.5*nS
    g_inh         = 64.0 * nS
    g_ext         = 80 * nS
    ext_mean      = 10
    ext_std       = 10
    max_distance  = size * mm/numpy.sqrt(2) # Since this is a torus
    max_delay     = max_distance/(0.3 * mm/ms)  # Needed for the connectors
    dt_diffusion  = 1*ms
    NO_diff       = 1e-6
    dt_modulation = 5*ms

    global target_NO
    target_NO = 0.50    #this is the homeostatic target of NO concentration (change to get different 'target' firing rates)

    NOsettletime = 50000*ms
    global sim_step
    sim_step = 1000*ms
    sampletime = 5000*ms
    settletime = 500*ms

  ##### DON'T DEFINE ANY SIMULATION PARS AFTER THIS

    global ds, dx2, dy2, dt_diff
    ds = 0.002
    dx2 = ds**2
    dy2 = ds**2
    dt_diff = dt_diffusion/(1*ms)

    print 'DT DIFF , max stable', dx2*dy2/( 2*NO_diff*(dx2+dy2) )
    print 'DT DIFF, using ;', dt_diff


    ### We create the cells and generate random positons in [0, size]x[0, size]
    global all_cells, exc_cells, inh_cells
    all_cells          = NeuronGroup(n_cells, model=eqs, threshold=threshold_fun, reset=Vr, refractory=tau_ref, clock=simulation_clock)
    all_cells.position = ds*size*generate_pos()
    exc_cells          = all_cells[0:n_exc]
    inh_cells          = all_cells[n_exc:n_cells]

    ### We initialize v values slightly above Vt, to have initial spikes
    all_cells.v        = El + 1.1*numpy.random.rand(n_cells) * (Vt_init - El)

    ### initialize threshold to be uniform at first
    all_cells.Vt = Vt_init*numpy.ones(n_cells)

    #all_cells.NO = numpy.ones(n_cells)*target_NO
    all_cells.NO = numpy.zeros(n_cells)

    ### Create 2-d array containing concentrations of NO at each x,y position
    global NO_concs
    #NO_concs = target_NO*numpy.ones((int(size/ds),int(size/ds)))
    NO_concs = numpy.zeros((int(size/ds),int(size/ds)))

    global NO_rates
    NO_rates = numpy.zeros((int(size/ds),int(size/ds)))

    global rates_min, rates_max, rates_mean, rates_std

    rates_max = []
    rates_min = []
    rates_mean = []
    rates_std =[]

    diffusion_clock = Clock(dt=dt_diffusion)
    modulation_clock = Clock(dt=dt_modulation)

    @network_operation(diffusion_clock, when='end')
    def do_diffusion(diffusion_clock):
        if diffusing:
            feed_nNOSstates_to_grid()
            diffuse_NOgrid()
        else:
            feed_nNOSstates_to_grid()
            nodiffuse_NOgrid()

    @network_operation(modulation_clock, when='start')
    def do_modulation(modulation_clock):
        if modulating:
            modulate_neurons(mod_tau,target_NO)


    print "Building network..."

    Ce = Connection(exc_cells, all_cells, 'ge', weight=g_exc, max_delay=max_delay,
                sparseness=epsilon,
                delay     =0.1*ms)
    Ci = Connection(inh_cells, all_cells, 'gi', weight=g_inh, max_delay=max_delay,
                sparseness=epsilon,
                delay     =0.1*ms)

    norm = numpy.random.normal(ext_mean,ext_std,n_cells*2)
    Cext = PoissonInput(all_cells,1,1*Hz*norm[norm>0][:n_cells],g_ext)

    print "Setting the recorders..."
    global s_all_short, s_all
    s_all   = SpikeCounter(all_cells)
    s_all_short = SpikeCounter(all_cells)

    ### Spike Monitor to trigger postsynaptic Ca influx after each spike
    SM_ca = SpikeMonitor(all_cells, function=ca_influx,record=False)

    print 'Nsteps:', int((NOsettletime/sim_step))

    print "settling network ..."

    run(settletime)

    print "Finished settling network ..."


    for time in xrange(int((NOsettletime/sim_step))):
        run(sim_step)

        run_step(time)

        try:
            print "simulation progress:", float(time/float(((NOsettletime)/sim_step)))
        except: pass

        s_all.count = numpy.zeros(n_cells)

    for time in xrange(int((sampletime/sim_step))):
        run(sim_step)

        run_step(time)

        try:
            print "sampling progress:", float(time/float(((sampletime)/sim_step)))
        except: pass

    rates_final = multiply(s_all.count,1000.*msecond/sampletime)
    Vt_final = all_cells.Vt.copy()
    s_all.count = numpy.zeros(n_cells)

    #print 'rates_pert ', rates_final
    #print 'Vt_pert ', Vt_final
    #print 'NO_concs ', all_cells.NO

    pylab.figure()
    pylab.plot(rates_mean)
    pylab.ylabel('mean firing rate (Hz)')
    pylab.xlabel('timestep')
    pylab.show()

    pylab.figure()
    pylab.hist(rates_final,20)
    pylab.xlabel('firing rate (Hz)')
    pylab.ylabel('count')
    pylab.show()

if __name__ == '__main__':
    main()
