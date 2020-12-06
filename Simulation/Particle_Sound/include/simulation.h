#ifndef SIMULATION_H
#define SIMULATION_H

#include <iostream>

#include <fluid.h>
#include <linearc.h>


/////////////////////////////////////////////////////
// Simulation Class /////////////////////////////////
/////////////////////////////////////////////////////

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

    // The length variables 
    // (assuming that the sim starts ar r = 0, z = 0 and renders until r = R_max, z = Z_max)
    double R_max;
    double Z_max;

    // Pointers to arrays that store the pressure elements
    Array* P_prev;      // Stores the previous step (P^j-1)
    Array* P_curr;      // Stores the current  step (P^j)
    Array* P_next;      // Stores the next     step (R^j+1)

    // Relevant Tridiagonal element Arrays
    // Since the matrix is block tridiagonal, instead of storing 
    // a bunch of matrices that have zero elements in them, we
    // Are storing the diagonals and the left hand side only
    Matrix* A;  // Upper Diagonal
    Matrix* B;  // Lower Diagonal
    Matrix* D;  // Center Diagonal
    Array* lhs; // Left hand side array
    
    // Function pointer to forcing function dE/dt
    // We use E instead of dE/dt for obvious reasons here.
    double (*E)(double,double);

    // Pointer to fluid object
    Fluid* fluid;



    public:
    // Constructors
    Simulation();
    Simulation(double,double,double,double,double,Fluid*,double (*)(double,double));
    Simulation(double,double,double,double,Fluid*,double (*)(double,double));
    Simulation(double,double,double,int,int,Fluid*,double (*)(double,double));
    Simulation(double,double,int,int,Fluid*,double (*)(double,double));

    // Function that performs one step of the simulation in time
    void step();

    // Function that runs for an arbitrary time T, or number of steps
    void run(double);
    void run(int);

    // Set's up the simulation, after this is called, the sim is ready to run
    // Called in constructor by default.
    void setup();

    // Takes as input a linear vector Nx x Nz and decomposes it to array vectors for tridiag solving.
    void convert_to_Array(Array*,Array*);

    // Set intial conditions for P by passing P_prev and P_curr
    void set_initial(Array&,Array&);

    // Creates the block tridiagonal matrix to solve
    void create_tridiag(Matrix*,Matrix*,Matrix*);
    
    // Creates the left hand side of the matrix equation as a list of vectors
    void create_lhs(Array*);

    // Mutators
    double get_dt();
    double get_dx();
    double get_dz();

    int get_Nx();
    int get_Nz();

    double get_R_max();
    double get_Z_max();

    Array* get_P_prev();
    Array* get_P_curr();
    Array* get_P_next();

    Fluid* get_fluid();


    void set_dt(double);
    void set_dx(double);
    void set_dz(double);

    void set_Nx(int);
    void set_Nz(int);
 
    void set_R_max(double);
    void set_Z_max(double);

    void set_fluid(Fluid*);

    // Printing functions
    void print();
    void print(std::ostream&);

};

#endif /* End of SIMULATION_H */