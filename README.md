There are three files in this repository. 
mainfreq_CLM.m is a Matlab function that evaluates the near-inertial wave minimum frequency expression for an idealized Gaussian shielded vortex.
raytrace_CLM_ode23.m is a Matlab function that performs ray-tracing for a given wavevector, position and frequency, in an idealized Gaussian shielded vortex. It uses ODE23 solver.
maincode.m contains all the numerical values (that can be altered easily) and links to the two functions above. Finally, it generates plots of the rays and minimum frequency expression.
