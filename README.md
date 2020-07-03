# tfm-laminarflow

The Matlab code is associated with the publications: 

M. Schaefer, W. Wicke, L. Brand, R. Rabenstein and R. Schober, “Transfer Function Models for Diffusion and Laminar Flow in Cylindrical Systems”, submitted to IEEE Transactions on Molecular, Biological and Multi-scale Communications, 2020, [online]: …

Further information regarding the implementation can be found on: 
maximilianschaefer.org/publication/tfm-laminar

Author:         Dr.-Ing Maximilian Schaefer
Affiliation:    University of Erlangen-Nuremberg
email:          max.schaefer@fau.de
website:        maximilianschaefer.org
Created:        April-July 2020
Last update:    29.06.2020 

The Matlab code is updated and extended regularly. 

# How to use? 

- Use the "startup" script to set all necessary paths 
- Parameters are defined in "main.m" 


# known issues

- Depending on the value of alpha and the type of particle release a lot of 
eigenvalues may be necessary. Even if the complete Code exploits the block 
partitioning of matrices, large matrices may occur. This can lead to a high 
computational load for the pre-calculation of the state matrix. 


# planned extensions

- Finish framework for the incorporation of external sources 
    - add functions for different temporal pulse shaping
    - add functions for different spatial pulse shapting (can be taken from 
initial conditions)

- add first order degradation reaction (straightforward by scaling of the
open loop state matrix state.As) 

- add more complex boundary conditions (permeable membrane) in radial 
direction

- Solve the issue with reflections at the boundaries in z-direction by the 
design of a perfectly matched layer

- Incorporate additional receiver models

