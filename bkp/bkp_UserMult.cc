//Dummy routine for mat-vec multiply (for now)

#include <petscmat.h>
#include <petscvec.h>
#include "UserMult.hh"

/*
    Constructor
*/
UserMult::UserMult()
{
}

void UserMult::mult(Mat mat, Vec x, Vec y){
    
    void *ptr;
    MatShellGetContext(mat, &ptr);	

    PetscErrorCode ierr;
    PetscInt size;
	
    //Get the size of the vector
    ierr = VecGetSize(x, &size);
    int init = size;

    //Allocate indArray and temp
    PetscInt indArray[init];
    PetscScalar temp[init];
    
    //Sets indArray to 1 to size
    for(int index=0; index<init; index++){
    	indArray[index] = index;
    }    

    //Gets the values of x and puts them in temp[]
    ierr = VecGetValues(x, init, indArray, temp);
    
    //Alters temp
    for(int index=0; index<init; index++){
    	temp[index] = temp[index];
    }  

    //Takes the values of temp[] and puts them in y
    ierr = VecSetValues(y, init, indArray, temp, INSERT_VALUES);



//return 0;
}

