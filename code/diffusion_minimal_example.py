## This script simulates NO diffusion over a 2D grid with periodic boundary conditions,
## with randomly positioned neurons emitting NO at random rates.
## A plot showing NO concentrations is shown every 100 timesteps

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


## The below code can easily be implemented in an existing
## brian model by defining a clock and associating a network
## operation with this clock. This operation then calls
## diffuse_NOgrid(), which is the code in the for loop below.
## The parameter NO_rates can be altered by other processes, and
## passed to this function. NO_rates describes the rate of production
## of NO for each 'neuron' on the 2D grid.
# @network_operation(diffusion_clock, when='end')
#    def do_diffusion(diffusion_clock):
#        diffuse_NOgrid()


## method from Philippedes et al. (2000, J Neurosci). Also a code snippet from http://goo.gl/18lc8


import numpy
from numpy import arange as arange

import pylab

NO_diff       = 0.1   # diffusion constant for NO
NOdecay_rate  = 0.90  # decay rate for NO
dt_diff = 10 # *ms (use ms units if using brian)

ds = 0.005
size = 1.0
time = 1000

dx2 = ds**2
dy2 = ds**2
dt_diff = dx2*dy2/( 2*NO_diff*(dx2+dy2) )

NO_concs = numpy.zeros((int(size/ds),int(size/ds)))

## populates grid with neurons at a probability of 0.02
neuron_present = numpy.random.binomial(2,0.02,(int(size/ds),int(size/ds)))

## assigns random NO production rate to each neuron; 0<rate<10
NO_rates =  numpy.random.randint(0,10,(int(size/ds),int(size/ds)))*neuron_present

NO_mean_t = []

for t in xrange(time):
    nx = int(size/ds)
    ny = int(size/ds)
    leftX = arange(nx)-1
    rightX = (arange(nx)+1)%nx
    leftY = arange(ny)-1
    rightY = (arange(ny)+1)%ny
    NO_concs = NO_concs + dt_diff*(NO_diff*((NO_concs[rightX,:]-2*NO_concs + NO_concs[leftX,:])/dx2 + (NO_concs[:,rightY] - 2*NO_concs
        + NO_concs[:,leftY])/dy2)+ NO_rates - NOdecay_rate*NO_concs)

    # make sure all values for NO_conc >= 0
    NO_concs.clip(min=0,out=NO_concs)
    NO_mean_t.append(numpy.mean(NO_concs))

    if t%100 == 0:
        print 't = ', t
        pylab.pcolor(NO_concs)
        pylab.show()

