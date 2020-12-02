#include <iostream>

#include <linearc.h>

using namespace std;

/////////////////////////////////////////////////////
// Array Class //////////////////////////////////////
/////////////////////////////////////////////////////
Array::Array(){
    this->ptr = nullptr;
    this->SIZE = 0;
}

Array::Array(int SIZE){
    double* arr = new double[SIZE];
    for (int i=0;i<SIZE;i++) arr[i] = 0;
    this->ptr = arr;
    this->SIZE = SIZE;
}

Array::Array(double* ptr,int SIZE){
    this->ptr = ptr;
    this->SIZE = SIZE;
}

/////////////////////////////////////////////////////
// CLASS FUNCTIONS //////////////////////////////////

int Array::size(){
    return this->SIZE;
}

void Array::print(bool LINEAR, ostream& stream){

    if (LINEAR){
        stream << "[ ";
        for (int i=0;i<this->size();i++){
            stream << this->ptr[i];
            if (i!=this->size()-1) stream << ", ";
        }
        stream << " ]" << endl;
    }

    else {
        for (int i=0;i<this->size();i++){
            stream << this->ptr[i] << endl;
        }
    }
}

void Array::print(bool Linear){
    print(Linear,cout);
}

void Array::print(ostream& stream){
    print(true,stream);
}

void Array::print(){
    print(true,cout);
}

/////////////////////////////////////////////////////
// Operator OVERLOADING ////////////////////////////

// This implements periodic boundaries on the array
double& Array::operator[](int i){
    if (i<this->size() && i >= 0) return this->ptr[i];
    else if (i>=this->size()) return this->ptr[i%this->size()];
    else return this->ptr[this->size()+i%(this->size())];
}

// Multiplication by a constant
Array Array::operator*(double C){
    Array A(new double[this->SIZE],this->SIZE);
    for (int i=0; i<this->size(); i++) A[i] = C*ptr[i];

    return A;
}

// Allowing for reverse multiplication of arrays with scalars
Array operator*(double C,Array A){
    return A*C;
}


// Dot product
Array Array::operator*(Array arr){
    // Check if they have the same length
    if (this->size() != arr.size()) throw "Dot product of non-equal arrays";
    
    Array A(new double[this->SIZE],this->SIZE);
    for (int i=0; i<this->size(); i++) A[i] = ptr[i]*arr[i];

    return A;
}

// Division by a constant
Array Array::operator/(double C){
    if (!C) throw "Division of array by zero";
    return *(this)*(1/C);
}


// Array Addition
Array Array::operator+(Array arr){
    // Check if they have the same length
    if (this->size() != arr.size()) throw "Addition of non-equal arrays";
    
    Array A(new double[this->SIZE],this->SIZE);
    for (int i=0; i<this->size(); i++) A[i] = ptr[i] + arr[i];

    return A;
}

// Array subtraction
Array Array::operator-(Array arr){
    return this->operator+(arr*(-1));
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////
// Matrix Class /////////////////////////////////////
/////////////////////////////////////////////////////

// Constructors
Matrix::Matrix(){
    this->ptr = nullptr;
    this->ROW_SIZE = 0;
    this->COL_SIZE = 0;
}

// Initializes a Matrix (ROW,COL) to zeroes
Matrix::Matrix(int ROW_SIZE, int COL_SIZE){
    double** mat = new double*[ROW_SIZE];
    for (int i=0;i<ROW_SIZE;i++){
        mat[i] = new double[COL_SIZE];
        for (int j=0;j<COL_SIZE;j++){
            mat[i][j] = 0;
        }
    }

    this->ptr = mat;
    this->ROW_SIZE = ROW_SIZE;
    this->COL_SIZE = COL_SIZE;
}

Matrix::Matrix(double** ptr,int ROW_SIZE, int COL_SIZE){
    this->ptr = ptr;
    this->ROW_SIZE = ROW_SIZE;
    this->COL_SIZE = COL_SIZE;
}

/////////////////////////////////////////////////////
// CLASS FUNCTIONS //////////////////////////////////

// To get the determinant of a square matrix
double Matrix::det(){
    if (!is_square()) throw "Attempted determinant of non-square matrix";
    Matrix M = lu_dcmp(*this);
    double dd = 1;
    for (int i=0;i<M.col_size();i++){
        dd *= M[i][i];
    }

    return dd;
}

// True if it's a square matrix
bool Matrix::is_square(){
    return ROW_SIZE==COL_SIZE;
}


// Prints the matrix through the output stream given 
// i.e. M.print(cout) etc.
void Matrix::print(std::ostream& stream){
    for (int i=0;i<ROW_SIZE;i++){
        for (int j=0;j<COL_SIZE;j++){
            stream << this->ptr[i][j]; 
            if (j<COL_SIZE-1) stream << ", ";
        }
        stream << "\n";
    }
}

// Prints the matrix in the console using cout
void Matrix::print(){
    print(cout);
}

// Returns the row size
int Matrix::row_size(){
    return ROW_SIZE;
}

// Returns the column size
int Matrix::col_size(){
    return COL_SIZE;
}

/////////////////////////////////////////////////////
// Operator OVERLOADING /////////////////////////////
double*& Matrix::operator[](int i){
    return ptr[i];
}

// Multiplies matrix by a double
Matrix Matrix::operator*(double X){

    // First create a new matrix object identical to M
    double** ptr = new double*[this->row_size()];
    Matrix mat(ptr,this->row_size(),this->col_size());

    // Copy M to mat 
    for (int i=0;i<mat.row_size();i++){
        ptr[i] = new double[mat.col_size()];
        for (int j=0; j < mat.col_size();j++){
            mat[i][j] = this->ptr[i][j]*X;
        }
    }

    return mat;
}

// Define the orher way for multiplication by a scalar
Matrix operator*(double C, Matrix M){
    return M*C;
}

// Divides matrix by a double
Matrix Matrix::operator/(double X){
    if (!X) throw "Matrix division by 0";
    return this->operator*(1/X);
}

// Adds a matrix to a matrix
Matrix Matrix::operator+(Matrix M){

    if (!((M.col_size() == col_size()) && (M.row_size() == row_size()))){
        throw "Adding matrices with unequal dimensions";
    }

    // First create a new matrix object identical to M
    double** ptr = new double*[this->row_size()];
    Matrix mat(ptr,this->row_size(),this->col_size());

    // Do the addition
    for (int i=0;i<mat.row_size();i++){
        ptr[i] = new double[mat.col_size()];
        for (int j=0; j < mat.col_size();j++){
            mat[i][j] = this->ptr[i][j]+M[i][j];
        }
    }

    return mat;
}

// Subtracting matrices
Matrix Matrix::operator-(Matrix M){
    return this->operator+(M*(-1));
}


// Matrix multiplication
Matrix Matrix::operator*(Matrix M){
    if (!((this->col_size() == M.row_size()))){
        throw "Multiplying incompatible matrices";
    }

    // First create a new matrix object identical to M
    double** ptr = new double*[this->row_size()];
    Matrix mat(ptr,this->row_size(),M.col_size());

    // Do the product
    for (int i=0;i<mat.row_size();i++){
        ptr[i] = new double[mat.col_size()];
        for (int j=0; j < mat.col_size();j++){
            double sum = 0;
            for (int n=0;n<this->col_size();n++) sum += this->ptr[i][n]*M[n][j];
            mat[i][j] = sum;
        }
    }

    return mat;
}

// Matrix multiplication with vector
Array Matrix::operator*(Array M){
    if (!((this->col_size() == M.size()))){
        throw "Multiplying incompatible matrices";
    }

    // First create a new matrix object identical to M
    double* ptr = new double[this->row_size()];
    Array mat(ptr,this->row_size());

    // Do the product
    for (int i=0;i<mat.size();i++){
        double sum = 0;
        for (int n=0;n<this->col_size();n++) sum += this->ptr[i][n]*M[n];
        mat[i] = sum;
    }

    return mat;
}


// Adds one matrix to another
Matrix Matrix::operator+=(Matrix M){
    if (!((M.col_size() == col_size()) && (M.row_size() == row_size()))){
        throw "Adding matrices with unequal dimensions";
    }

    for (int i=0;i<this->row_size();i++){
        for (int j=0; j < this->col_size();j++){
            this->ptr[i][j]+=M[i][j];
        }
    }

    return *this;
}


// Subtracts one matrix from another
Matrix Matrix::operator-=(Matrix M){
    if (!((M.col_size() == col_size()) && (M.row_size() == row_size()))){
        throw "Adding matrices with unequal dimensions";
    }

    for (int i=0;i<this->row_size();i++){
        for (int j=0; j < this->col_size();j++){
            this->ptr[i][j]-=M[i][j];
        }
    }

    return *this;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////
// LU Decomposition /////////////////////////////////
/////////////////////////////////////////////////////
// Uses Crout's Algorithm to compute LU decomposition in O(n^3) for nxn
/*
   The returned matrix is a factorisation of two matrices: M = AB
   but placed on the same matrix, like so:

        |1  0  0  0|
   A =  |a  1  0  0|
        |a  a  1  0|
        |a  a  a  1|

        |b  b  b  b|
   B =  |0  b  b  b|
        |0  0  b  b|
        |0  0  0  b|
 
   So M becomes:

        |b  b  b  b|
   M =  |a  b  b  b|
        |a  a  b  b|
        |a  a  a  b|
*/
void lu_dcmp_replace(Matrix& M){

    for(int j=0; j<M.col_size(); j++){

        // Forward
        for(int i=0; i<=j; i++){
            double sum = 0;
            for (int n=0;n<i;n++){
                sum += M[i][n] * M[n][j];
            }

            M[i][j] -= sum;
        }

        // Backward
        for(int i=j+1; i<M.col_size(); i++){
            double sum = 0;
            for (int n=0;n<j;n++){
                sum += M[i][n] * M[n][j];
            }

            M[i][j] = (M[i][j] - sum)/(M[j][j]);
        }
    }
}


//Decomposes M into a new matrix without replacing M
Matrix lu_dcmp(Matrix M){
    // First create a new matrix object identical to M
    double** ptr = new double*[M.row_size()];
    Matrix mat(ptr,M.row_size(),M.col_size());

    // Copy M to mat 
    for (int i=0;i<mat.row_size();i++){
        ptr[i] = new double[M.col_size()];
        for (int j=0; j < mat.col_size();j++){
            mat[i][j] = M[i][j];
        }
    }

    // Replace mat with it's decomposition
    lu_dcmp_replace(mat);
    return mat;
}

//Decomposes M into a new matrix and places it to N
void lu_dcmp(Matrix& M, Matrix& N){
    for(int j=0; j<M.col_size(); j++){

        // Forward
        for(int i=0; i<=j; i++){
            double sum = 0;
            for (int n=0;n<i;n++){
                sum += M[i][n] * M[n][j];
            }

            N[i][j] = M[i][j] - sum;
        }

        // Backward
        for(int i=j+1; i<M.col_size(); i++){
            double sum = 0;
            for (int n=0;n<j;n++){
                sum += M[i][n] * M[n][j];
            }

            N[i][j] = (N[i][j] - sum)/(M[j][j]);
        }
    }
}





//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////
// Various Solvers //////////////////////////////////
/////////////////////////////////////////////////////

// Solve equation A X = B where X,B are either matrices or vectors, using LU decomposition.
Array solve(Matrix& A,Array& b,bool decomposed){
    
    // If the matrix is not in LU decompose it to LU
    Matrix M = !decomposed ? lu_dcmp(A) : A; 
    
    // Create the output vector
    Array x(b.size());

    // Forward Substitution
    for(int i=0;i<x.size();i++){
        double sum = 0;
        for (int j=0;j<i-1;j++) sum+=M[i][j]*x[j];
        x[i] = b[i] - sum;
    }

    //backward subsitution
    for(int i=x.size()-1;i>=0;i--){
        double sum = 0;
        for (int j=i+1;j<x.size();j++) sum+=M[i][j]*x[j];
        x[i] = (x[i]-sum)/M[i][i];
    }

    return x;
}

// Assuming the matrix is not decomposed
Array solve(Matrix& A,Array& b){
    return solve(A,b,false);
}


// Solve equation A X = B where X,B are either matrices or vectors, using LU decomposition.
void solve(Matrix& A,Array& b, Array& x,bool decomposed){
    
    // If the matrix is not in LU decompose it to LU
    Matrix M = !decomposed ? lu_dcmp(A) : A; 

    // Forward Substitution
    for(int i=0;i<x.size();i++){
        double sum = 0;
        for (int j=0;j<i-1;j++) sum+=M[i][j]*x[j];
        x[i] = b[i] - sum;
    }

    //backward subsitution
    for(int i=x.size()-1;i>=0;i--){
        double sum = 0;
        for (int j=i+1;j<x.size();j++) sum+=M[i][j]*x[j];
        x[i] = (x[i]-sum)/M[i][i];
    }
}

// Assuming hte matrix is not decomposed
void solve(Matrix& A,Array& b, Array& x){
    solve(A,b,x,false);
}


// To do the same thing but with X, B matrices:
Matrix solve(Matrix& A, Matrix& b, bool decomposed){
    
    // If the matrix is not in LU decompose it to LU
    Matrix M = !decomposed ? lu_dcmp(A) : A; 
    
    //Create the output matrix
    Matrix x(M.row_size(),b.col_size()); 

    //Solve the equivalent vecotr equtions
    Array xx(M.row_size());
    for (int i=0; i < b.col_size(); i++){
        for (int j=0; j < xx.size(); j++) xx[j] = b[j][i];
        Array tmp = solve(M,xx,true);
        for (int j=0; j < xx.size(); j++)x[j][i] = tmp[j];
    }

    return x;
}


// To do the same thing but with X, B matrices assuming decomposed = false:
Matrix solve(Matrix& A, Matrix& b){
    return solve(A,b,false);
}


// To do the same thing but with X, B matrices:
void solve(Matrix& A, Matrix& b, Matrix& x, bool decomposed){
    
    // If the matrix is not in LU decompose it to LU
    Matrix M = !decomposed ? lu_dcmp(A) : A; 

    //Solve the equivalent vecotr equtions
    Array xx(M.row_size());
    for (int i=0; i < b.col_size(); i++){
        for (int j=0; j < xx.size(); j++) xx[j] = b[j][i];
        Array tmp = solve(M,xx,true);
        for (int j=0; j < xx.size(); j++)x[j][i] = tmp[j];
    }
}


// To do the same thing but with X, B matrices assuming decomposed = false:
void solve(Matrix& A, Matrix& b, Matrix& x){
    solve(A,b,x,false);
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////
// Tridiagonal Solver ///////////////////////////////
/////////////////////////////////////////////////////

// Solves a tridiagonal system of equations in O(N) time
// Specifically we solve equation Mx = b where:
//      - A is the upper diagonal
//      - B is the lower diagonal
//      - D is the diagonal 
//      - b is the equal vector
Array tridiag(Array& A, Array& B, Array& D, Array& b){

    // We do two passes, a forward and a backward pass
    // Forward pass
    double* gamma = new double[A.size()];
    double* beta  = new double[D.size()];

    if(!D[0]) throw "Tridiag division by 0";
    gamma[0] = A[0]/D[0];
    beta[0]  = b[0]/D[0];
    for (int i=1; i<D.size(); i++){
        if(!(D[i] - B[i-1]*gamma[i-1])) throw "Tridiag division by 0";

        if (i<A.size()) gamma[i] = A[i]/(D[i] - B[i-1]*gamma[i-1]);
        beta[i] = (b[i]-beta[i-1]*B[i-1])/(D[i] - B[i-1]*gamma[i-1]);
}

    // Backward Pass
    Array x(D.size());
    x[x.size()-1] = beta[x.size()-1];
    for (int i=x.size()-2; i>=0;i--){
        x[i] = beta[i] - gamma[i] * x[i+1];
    }

    // Return the final vector
    return x;
}

// Solves tridiagonal matrix equation Mx = b in O(n) 
Array tridiag(Matrix& M, Array& b){
    // Check for the matrix being square
    if(!M.is_square()) throw "Tridiagonal solving of non-square matrix";
    
    // obtain the arrays for the tridiagonal matrix
    Array D(M.row_size()), A(M.row_size()-1), B(M.row_size()-1);
    for (int i=0;i<M.row_size();i++){
        if(i<A.size()){
            A[i] = M[i][i+1];
            B[i] = M[i+1][i];
        } 
        D[i] = M[i][i];
    }

    return tridiag(A,B,D,b);
}






//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////
// Block Tridiagonal Solver /////////////////////////
/////////////////////////////////////////////////////

// Solves a block tridiagonal system
// We are slightly modifying thomas algorithm
Array block_tridiag(Matrix* A,int A_SIZE,Matrix* B,int B_SIZE,Matrix* D,int D_SIZE,Array* b,int b_SIZE){

    Matrix* Gamma = new Matrix[A_SIZE];
    Array* Beta = new Array[D_SIZE];

    // Decompose the original matrix
    Matrix lu = lu_dcmp(D[0]);
    Gamma[0] = solve(lu,A[0],true);
    Beta[0] = solve(lu,b[0],true);

    // Forward propagation
    for(int i=1; i<D_SIZE; i++){
        
        // get the decomposition of the divisor
        lu = lu_dcmp(D[i] - B[i-1]*Gamma[i-1]);
        if (i < A_SIZE)Gamma[i] = solve(lu,A[i],true);
        Beta[i] = 

        // if(!(D[i] - B[i-1]*gamma[i-1])) throw "Tridiag division by 0";

        // if (i<A.size()) gamma[i] = A[i]/(D[i] - B[i-1]*gamma[i-1]);
        // beta[i] = (b[i]-beta[i-1]*B[i-1])/(D[i] - B[i-1]*gamma[i-1]);

    }

}