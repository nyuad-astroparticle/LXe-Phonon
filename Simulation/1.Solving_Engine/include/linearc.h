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
    int SIZE;      // Size variable :)

    public:
    // Constructors
    Array();
    Array(int);
    Array(double*, int);
    void init(int);

    // Overloading [] to allow for periodicity in indexing
    double& operator[](int);
    double& operator[](int) const;

    // Overloading other operators to allow for some elementary manipulation
    Array operator+(const Array&) const;
    Array operator-(const Array&) const;
    Array operator*(const double&) const;
    Array operator/(const double&) const;
    Array operator*(const Array&) const;
    Array& operator=(const Array&);

    // Returns the size of the array
    int size() const;

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
    void init(int,int);

    // Overload some opeartors
    Matrix operator*(const Matrix&) const;
    Array operator*(const Array&) const;
    Matrix operator+(const Matrix&) const;
    Matrix operator-(const Matrix&) const;
    Matrix operator*(const double&) const;
    Matrix operator/(const double&) const;
    Matrix& operator+=(const Matrix&);
    Matrix& operator-=(const Matrix&);
    Matrix& operator=(const Matrix&);
    double*& operator[](int);
    double*& operator[](int) const;

    // true if NxN
    bool is_square() const;

    // Get the determinant
    double det();
    static double det(Matrix);

    // Get the inverse (the slow way)
    // Matrix inv();

    //Print the matrix
    void print();
    void print(std::ostream&);

    // Mutators
    int row_size() const;
    int col_size() const;
    double** get_ptr() const;

    // Destructor
    ~Matrix();
};

Matrix operator*(double,Matrix);




/////////////////////////////////////////////////////
// LU Decomposition /////////////////////////////////
/////////////////////////////////////////////////////
// Uses Crout's Algorithm to compute LU decomposition in O(n^3) for nxn
Matrix lu_dcmp(const Matrix&);          // Returns the decomposed version (takes O(n^2)) to create one
void lu_dcmp_replace(Matrix&);          // replaces the current matrix with it's decomposed vesion
void lu_dcmp(const Matrix&,Matrix&);    // replaces the second matrix provided

/////////////////////////////////////////////////////
// Various Solvers //////////////////////////////////
/////////////////////////////////////////////////////

// Solve equation A X = B where X,B are either matrices or vectors, using LU decomposition.
Array solve(const Matrix&,const Array&);                // bool false
Array solve(const Matrix&,const Array&,bool);           // Set bool to True if M is already lu-decomposed
void solve(const Matrix&,const Array&,Array&);          // bool false | Replace the second array
void solve(const Matrix&,const Array&,Array&,bool);     // Set bool to True if M is already lu-decomposed | Replace the second array

Matrix solve(Matrix&,Matrix&);                          // bool false
Matrix solve(Matrix&,Matrix&,bool);                     // Set bool to True if M is already lu-decomposed
void solve(Matrix&,Matrix&,Matrix&);                    // bool false | Replace the third matrix
void solve(Matrix&,Matrix&,Matrix&,bool);               // Set bool to True if M is already lu-decomposed | Replace the third matrix


/////////////////////////////////////////////////////
// Tridiagonal Solver ///////////////////////////////
/////////////////////////////////////////////////////

// Solves a tridiagonal system of equations in O(N) time
Array tridiag(const Array&,const Array&,const Array&,const Array&);
void tridiag(const Array&,const Array&,const Array&,const Array&, Array&);
Array tridiag(const Matrix&,const Array&);
void tridiag(const Matrix&,const Array&, Array&);


/////////////////////////////////////////////////////
// Block Tridiagonal Solver /////////////////////////
/////////////////////////////////////////////////////

// Solves a tridiagonal system of equations in O(N M^2) time
Array* block_tridiag(const Matrix*,int,const Matrix*,int,const Matrix*,int,const Array*,int);
void block_tridiag(const Matrix*,int,const Matrix*,int,const Matrix*,int,const Array*,int,Array*); // Same thing with replacement
Array block_tridiag(const Matrix&,const Array&,int,int);
Array block_tridiag(const Matrix&,const Array&);

#endif /* End of LINEARC_H */