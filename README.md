# CHAOTIC BEHAVIOR IN A PLANE PENDULUM WITH A SINUSOIDAL DRIVING TORQUE AT THE PIVOT POINT

In this work, a plane pendulum driven by a sinusoidal driving torque at the pivot point was simulated. Even though for the case of small θ (<10°) and no external torque, the motion can be simplified into a simple pendulum, the nonlinear term sinθ will become significant as θ increases and transit the system into a chaotic one. 

To solve the second order nonlinear ordinary differential equations, numerical methods including Euler’s Method and forth order Runge-Kutta method were implemented via MATLAB. Two methods were compared in simple harmonic oscillator to determine the accuracy of method with same step size. Through converting the iterative methods into programs, simulation of the motion was done. Some features of the chaotic dynamics (violation of strong principle of causality) were observed and discussed.

[__motion_sim.m__](motion_sim.m) contains all functions for ODE solver, graph plotting, and animation simulation. 
__PPT_for_the_project.pdf__ and [__Simulation_on_the_chaotic_behaviours_of_pendulum_motion.pdf__](Simulation_on_the_chaotic_behaviours_of_pendulum_motion.pdf) contain the detailed report on the simulation results, data and error analysis.
