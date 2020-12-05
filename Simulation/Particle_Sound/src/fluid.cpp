#include <fluid.h>
#include <fstream>
#include <iostream>

using namespace std;

// Default Constructor for fluid class
Fluid::Fluid(){
    rest_density = 0;
    viscocity = 0;
    compressibility = 0;
    coefficient_thermal_expansion = 0;
    specific_heat_pressure = 0;
}

// Constructor
Fluid::Fluid(double rest_density, double viscocity, double compressibility, double coefficient_thermal_expansion, double specific_heat_pressure){
    this->rest_density = rest_density;
    this->viscocity = viscocity;
    this->compressibility = compressibility;
    this->coefficient_thermal_expansion = coefficient_thermal_expansion;
    this->specific_heat_pressure = specific_heat_pressure;
}

// Constructor from a file
Fluid::Fluid(const char* filename){

    // Open file
    ifstream file;
    file.open(filename);

    // If the file fails load the default constructor
    if (file.fail()){
        cout << "Failed to open the file, going to default constructor" << endl;
        rest_density = 0;
        viscocity = 0;
        compressibility = 0;
        coefficient_thermal_expansion = 0;
        specific_heat_pressure = 0;

        return;
    }

    // Now read the file
    char c;
    int i = 0;
    double x[5];
    while (!file.eof()) if (file.get() == ':') file >> x[i++];

    // Close the file after reading
    file.close();

    // Assign the variables
    rest_density = x[0];
    viscocity = x[1];
    compressibility = x[2];
    coefficient_thermal_expansion = x[3];
    specific_heat_pressure = x[4];
}


// Mutators
double Fluid::get_rest_density(){
    return rest_density;
}

double Fluid::r0(){
    return get_rest_density();
}


double Fluid::get_viscocity(){
    return rest_density;
}

double Fluid::eta(){
    return get_viscocity();
}


double Fluid::get_compressibility(){
    return rest_density;
}

double Fluid::K(){
    return get_compressibility();
}


double Fluid::get_coefficient_thermal_expansion(){
    return rest_density;
}

double Fluid::beta(){
    return get_coefficient_thermal_expansion();
}



double Fluid::get_specific_heat_pressure(){
    return rest_density;
}

double Fluid::Cp(){
    return get_specific_heat_pressure();
}
