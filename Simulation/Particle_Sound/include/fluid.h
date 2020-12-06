#ifndef FLUID_H
#define FLUID_H

/////////////////////////////////////////////////////
// Fluid Class //////////////////////////////////////
/////////////////////////////////////////////////////

// Class that represents a fluid. Contains all the 
// usefull constants, nondimensionalisation parameters, etc.

class Fluid{

    private:
    // Relevant Variables
    double rest_density;
    double viscocity;
    double compressibility;
    double coefficient_thermal_expansion;
    double specific_heat_pressure;


    public:
    // Constructors
    Fluid();
    Fluid(const char*);                          // Initializes from file
    Fluid(double,double,double,double,double);   // Initializes throuhg varibles

    // Mutators
    double get_rest_density();
    double r0();

    double get_viscocity();
    double mu();

    double get_compressibility();
    double K();

    double get_coefficient_thermal_expansion();
    double beta();

    double get_specific_heat_pressure();
    double Cp();
};


#endif /* End of FLUID_H*/