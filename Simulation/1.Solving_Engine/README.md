# Strongly Damped Wave Simulation Engine

## Outline

A simulation engine written from scratch in C++ designed to be used for numerically evaluating the strongly damped wave equation derived in earlier parts of this project.

The methods that are attempted are finite differencing and finite element. The method with the best results will be numerically eveluated and parallelised in order to be run in NYU's high performance computing facility **Dalma**.

## Implementation

To compile the simulation you need CMAKE. With a configured CMAKE download this repository and run

    $ ccmake .
    $ cmake --build .
    $ cmake --install .

Then magik will happen. :_ ) 

## Contents

Here we outline the various contents of our simulation.