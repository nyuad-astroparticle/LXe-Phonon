#ifndef SIMULATION_H
#define SIMULATION_H

#include <iostream>

#include <linearc.h>

// This class is in charge of the simulation
// Everything that is happening is managed by this class
class Simulation{

    private:
    // The discretised step lengths
    double dt;
    double dx;
    double dz;

    // The relevant lengths of the discretisation 
    int Nx;
    int Nz;
    int N;

    // The length variables 
    // (assuming that the sim starts ar r = 0, z = 0 and renders until r = R_max, z = Z_max)
    double R_max;
    double Z_max;

    // Pointers to arrays that store the pressure elements
    Array* P_prev;      // Stores the previous step (P^j-1)
    Array* P_curr;      // Stores the current  step (P^j)
    
    // Function pointer to forcing function dE/dt
    // We use E instead of dE/dt for obvious reasons here.
    double (*E)(double,double);

    public:
    // Constructors
    Simulation();
    
    // Function that performs one step of the simulation in time
    void step();

    // Function that runs for an arbitrary time T, or number of steps
    void run(double);
    void run(int);

    // Mutators
    double get_dt();
    double get_dx();
    double get_dz();

    int get_Nx();
    int get_Nz();
    int get_N();

    double get_R_max();
    double get_Z_max();

    Array* get_P_prev();
    Array* get_P_curr();


    void set_dt(double);
    void set_dx(double);
    void set_dz(double);

    void set_Nx(int);
    void set_Nz(int);
    void set_N(int);
 
    void set_R_max(double);
    void set_Z_max(double);

    void set_P_prev(Array*);
    void set_P_curr(Array*);



};

#endif /* End of SIMULATION_H */