# Yann Sweeney, Jeanette Hellgren Kotaleski, Matthias H. Hennig
# A diffusive homeostatic signal maintains neural heterogeneity and responsiveness in cortical networks

To see the code used to implement the numerical diffusion method, see code/diffusion_minimal_example.py . Brief instructions on how to include this in an existing brian model are contained within. 

To generate the data for each figure run the run_simulations_figX.py script from command line using
'python run_simulations_figX.py' , where X is the figure number. A lot of text output will be created showing the progress. You will need to be in the same directory as the 'network_simulation.py' file, and other .py files.  

Python with numpy and the brian simulator (briansimulator.org) will need to be installed. These scripts were tested on a machine with OSX 10.10, Python v2.7.10, numpy v1.9.2, and brian v1.4.1

As the simulations can take quite a long time to run (several hours on a multicore machine) and the resulting datafiles are quite large (~1 GB total), sample data generated from the code are provided for each figure and may be downloaded from figshare.

Figure 1: http://dx.doi.org/10.6084/m9.figshare.1454450

Figure 4: http://dx.doi.org/10.6084/m9.figshare.1454455 

Figure 5: http://dx.doi.org/10.6084/m9.figshare.1454469

For any queries, comments or request please do not hesitate to contact me at yann.sweeney@ed.ac.uk. 
