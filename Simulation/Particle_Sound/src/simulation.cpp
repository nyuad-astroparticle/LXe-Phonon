#include <simulation.h>

/////////////////////////////////////////////////////
// Simulation Class /////////////////////////////////
/////////////////////////////////////////////////////


/////////////////////////////////////////////////////
// Constructors /////////////////////////////////////





/////////////////////////////////////////////////////
// Class Functions //////////////////////////////////

// Creates the tridiagonal matrix to be solved in each iteration.
// Assumes that the arrays A, B, D, and b exist and are initialised
// the partition happens every Z_max cells for R_max repetitions
void Simulation::create_tridiag(Matrix* A, Matrix* B, Matrix* D){
    
    // For all the submatrices
    for(int n=0; n<R_max; n++){

        // Initialize the D matrix
        for (int i=0; i<Z_max; i++){
            D[n][i][i]   = -4*fluid->mu()/(fluid->K()*dx*dx*dt) + fluid->r0()/(fluid->K()*dt*dt);
            D[n][i][i+1] = fluid->mu()/(fluid->K()*dx*dx*dt);
            D[n][i][i-1] = fluid->mu()/(fluid->K()*dx*dx*dt);
        }

        // Initialize the A and B matrices
        for (int i=0; i<Z_max-1; i++){
            double r = n*dx;
            A[n][i][i]  = fluid->mu()/(fluid->K()*dx*dx*dt) + fluid->mu()/(2*r*fluid->K()*dx*dt);
            
            r = (n+1)*dx;
            B[n][i][i]  = fluid->mu()/(fluid->K()*dx*dx*dt) - fluid->mu()/(2*r*fluid->K()*dx*dt);
        }
    }

}


// Creates an array for the left hand side to be used for tridiagonal solving
// Assumes that the array b exists and has already the correct dimensions
void Simulation::create_lhs(Array* b){
    for(int n=0; n<R_max; n++){
        for(int i=0; i<Z_max; i++){
            double r = n*dx;
            // For the cells that are not in the boundaries
            if ((n+1 < R_max && n-1 >= 0) && (i+1<Z_max && i-1>=0)){
                b[n][i] = 
                (1 + fluid->mu()/(fluid->K() * dt))/(dx*dx)*(
                    P_curr[n][i+1] - 2*P_curr[n][i] + P_curr[n][i-1] + 
                    P_curr[n+1][i] - 2*P_curr[n][i] + P_curr[n-1][i]
                ) + 1/(2*r*dx)* (P_curr[n+1][i] - P_curr[n-1][i]) +
                fluid->r0()/(fluid->K()*dt*dt)* (2*P_curr[n][i] - P_prev[n][i]) - 
                fluid->beta()/fluid->Cp() * E(r,i*dx);
            }


            // Now we are going to implement a box because it is easy
            // But it actually may be better to change this in the future
            // Specifically we are implementing fully absorbent boundary conditions
            else if(n+1 >= R_max){
                if (i+1 >= Z_max){

                    b[n][i] = 
                    (1 + fluid->mu()/(fluid->K() * dt))/(dx*dx)*(
                        - 2*P_curr[n][i] + P_curr[n][i-1] + 
                        - 2*P_curr[n][i] + P_curr[n-1][i]
                    ) + 1/(2*r*dx)* (- P_curr[n-1][i]) +
                    fluid->r0()/(fluid->K()*dt*dt)* (2*P_curr[n][i] - P_prev[n][i]) - 
                    fluid->beta()/fluid->Cp() * E(r,i*dx);

                }

                else if (i-1<0){

                    b[n][i] = 
                    (1 + fluid->mu()/(fluid->K() * dt))/(dx*dx)*(
                        P_curr[n][i+1] - 2*P_curr[n][i] + 
                        - 2*P_curr[n][i] + P_curr[n-1][i]
                    ) + 1/(2*r*dx)* (- P_curr[n-1][i]) +
                    fluid->r0()/(fluid->K()*dt*dt)* (2*P_curr[n][i] - P_prev[n][i]) - 
                    fluid->beta()/fluid->Cp() * E(r,i*dx);

                }
            }

            else if(n-1 < 0){
                if(i+1 >= Z_max){
                    b[n][i] = 
                    (1 + fluid->mu()/(fluid->K() * dt))/(dx*dx)*(
                        - 2*P_curr[n][i] + P_curr[n][i-1] + 
                        P_curr[n+1][i] - 2*P_curr[n][i]
                    ) + 1/(2*r*dx)* (P_curr[n+1][i]) +
                    fluid->r0()/(fluid->K()*dt*dt)* (2*P_curr[n][i] - P_prev[n][i]) - 
                    fluid->beta()/fluid->Cp() * E(r,i*dx);
                }

                else if(i-1 < 0){
                    b[n][i] = 
                    (1 + fluid->mu()/(fluid->K() * dt))/(dx*dx)*(
                        P_curr[n][i+1] - 2*P_curr[n][i] + 
                        P_curr[n+1][i] - 2*P_curr[n][i]
                    ) + 1/(2*r*dx)* (P_curr[n+1][i]) +
                    fluid->r0()/(fluid->K()*dt*dt)* (2*P_curr[n][i] - P_prev[n][i]) - 
                    fluid->beta()/fluid->Cp() * E(r,i*dx);
                }
            }
        }
    }
}



/////////////////////////////////////////////////////
// Mutators /////////////////////////////////////////
