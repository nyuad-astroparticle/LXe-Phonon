# LXe-Phonon
**Theoretical Study on the possibility of acoustic detection of particles through noble liquids**

Single particle detection in cryogenic and noble liqiuds has become prominent in recent years due to the seach of dark matter. In this work we examine the possibility of augementing our detection capabilities by adding an additional detection channel through the acoustic detection of single particles. We examine multiple methodologies of energy deposition and present theoretical estimates for the peak pressure over time for multiple *Minimum Ionising Particles* on different noble liquids.

This repository contains the relevant simulations, code and background for said project. Specifically, we developed a framework in C++ to calculate single particle energy deposition using finite element methods, as well as multiple python notebooks to examine the properties of certain analytic solutions.

# Overview
The contents are as follows.

1. **[Phonon Background Simulation in LXe](https://github.com/PanosEconomou/LXe-Phonon/tree/master/Simulation)**  
  Simulates the thermal noise of phonons in Liquid Xenon through statistical arguments, by modelling phonons as sound bozons.
2. **[Analytic Solution of Viscous Wave Equation](https://github.com/PanosEconomou/LXe-Phonon/tree/master/Simulation/Fluid_Wave)**  
  Solving the Wave equation with a viscous damping term introduced from allowing viscosity in Navier Stokes Equations
3. **[Particle Energy Deposition Simulation](https://github.com/PanosEconomou/LXe-Phonon/tree/master/Simulation/Particle_Sound)**  
  Computational model for calculating the energy deposition of MIPs in a noble liquid through finite element methods.
  

## Notes
The theoretical paper can be fonud [here](panos.pw "No link yet : (")

Here is an image of the proposed detector

![Single Particle Phonon Detector](https://github.com/PanosEconomou/LXe-Phonon/blob/master/Simulation/Detector/detector%20ortho.png)
