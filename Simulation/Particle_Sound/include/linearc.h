#ifndef LINEARC_H
#define LINEARC_H

#include <iostream>

/////////////////////////////////////////////////////
// Array Class //////////////////////////////////////
/////////////////////////////////////////////////////

// A class that represents an array the way we want 
// to, with properly overloaded operators, etc.
class Array{
    private:
    double* ptr;   // Array pointer
    int SIZE;   // Size variable :)

    public:
    // Constructors
    Array();
    Array(int);
    Array(double*, int);

    // Overloading [] to allow for periodicity in indexing
    double& operator[](int);

    // Overloading other operators to allow for some elementary manipulation
    Array operator+(Array);
    Array operator-(Array);
    Array operator*(double);
    Array operator/(double);
    Array operator*(Array);

    // Returns the size of the array
    int size();

    // Prints the contents of the array
    void print();
    void print(bool);
    void print(std::ostream&);
    void print(bool,std::ostream&);

    // Destructor
    ~Array();
};

Array operator*(double,Array);



/////////////////////////////////////////////////////
// Matrix Class /////////////////////////////////////
/////////////////////////////////////////////////////

class Matrix{
    private:
    double** ptr;
    int ROW_SIZE;
    int COL_SIZE;

    public:
    // Constructors
    Matrix();
    Matrix(int,int);
    Matrix(double**,int,int);

    // Overload some opeartors
    Matrix operator*(Matrix);
    Array operator*(Array);
    Matrix operator+(Matrix);
    Matrix operator-(Matrix);
    Matrix operator*(double);
    Matrix operator/(double);
    Matrix operator+=(Matrix);
    Matrix operator-=(Matrix);
    double*& operator[](int);

    // true if NxN
    bool is_square();

    // Get the determinant
    double det();
    static double det(Matrix);

    // Get the inverse (the slow way)
    // Matrix inv();

    //Print the matrix
    void print();
    void print(std::ostream&);

    // Mutators
    int row_size();
    int col_size();

    // Destructor
    ~Matrix();
};

Matrix operator*(double,Matrix);




/////////////////////////////////////////////////////
// LU Decomposition /////////////////////////////////
/////////////////////////////////////////////////////
// Uses Crout's Algorithm to compute LU decomposition in O(n^3) for nxn
Matrix lu_dcmp(Matrix);            // Returns the decomposed version (takes O(n^2)) to create one
void lu_dcmp_replace(Matrix&);      // replaces the current matrix with it's decomposed vesion
void lu_dcmp(Matrix&,Matrix&);      // replaces the second matrix provided

/////////////////////////////////////////////////////
// Various Solvers //////////////////////////////////
/////////////////////////////////////////////////////

// Solve equation A X = B where X,B are either matrices or vectors, using LU decomposition.
Array solve(Matrix,Array);                // bool false
Array solve(Matrix,Array,bool);           // Set bool to True if M is already lu-decomposed
void solve(Matrix&,Array&,Array&);          // bool false | Replace the second array
void solve(Matrix&,Array&,Array&,bool);     // Set bool to True if M is already lu-decomposed | Replace the second array

Matrix solve(Matrix,Matrix);              // bool false
Matrix solve(Matrix,Matrix,bool);         // Set bool to True if M is already lu-decomposed
void solve(Matrix&,Matrix&,Matrix&);        // bool false | Replace the third matrix
void solve(Matrix&,Matrix&,Matrix&,bool);   // Set bool to True if M is already lu-decomposed | Replace the third matrix


/////////////////////////////////////////////////////
// Tridiagonal Solver ///////////////////////////////
/////////////////////////////////////////////////////

// Solves a tridiagonal system of equations in O(N) time
Array tridiag(Array&,Array&,Array&,Array&);
Array tridiag(Matrix&,Array&);


/////////////////////////////////////////////////////
// Block Tridiagonal Solver /////////////////////////
/////////////////////////////////////////////////////

// Solves a tridiagonal system of equations in O(N M^2) time
Array* block_tridiag(Matrix*,int,Matrix*,int,Matrix*,int,Array*,int);
Array block_tridiag(Matrix,Array,int,int);
Array block_tridiag(Matrix,Array);

#endif /* End of LINEARC_H */