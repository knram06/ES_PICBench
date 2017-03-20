
static char help[] = "Solves a tridiagonal linear system with KSP.\n\n";

/*T
   Concepts: KSP^solving a system of linear equations
   Processors: 1
T*/

/*
  Include "petscksp.h" so that we can use KSP solvers.  Note that this file
  automatically includes:
     petscsys.h       - base PETSc routines   petscvec.h - vectors
     petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners

  Note:  The corresponding parallel example is ex23.c
*/
#include <petscksp.h>

// define the matrix-free mat mult operation
// y = mat*x
#undef __FUNCT__
#define __FUNCT__ "userMult"
PetscErrorCode userMult(Mat mat, Vec x, Vec y)
{
   PetscErrorCode ierr;
   // get the vec size
   PetscInt n, i;
   VecGetSize(x, &n);

   PetscScalar *out = malloc(sizeof(PetscScalar)*n);
   PetscInt *indices = malloc(sizeof(PetscInt)*n);

   // get elements of Vec x
   PetscScalar *valX;
   VecGetArray(x, &valX);

   // loop through petsc scalar array and fill in values based on Mat values
   // fill in first element
   out[0] = 2*valX[0] - 1*valX[1];
   indices[0] = 0;
   // fill in remaining elements
   for(i = 1; i < n-1; i++)
   {
       out[i] = -valX[i-1] + 2*valX[i] - valX[i+1];
       indices[i] = i;
   }
   // fill in last element
   out[n-1] = -valX[n-2] + 2*valX[n-1];
   indices[n-1] = n-1;

   // set this to output array
   ierr = VecSetValues(y, n, indices, out, INSERT_VALUES); CHKERRQ(ierr);
   ierr = VecAssemblyBegin(y); CHKERRQ(ierr);
   ierr = VecAssemblyEnd(y); CHKERRQ(ierr);

   // restore values to x
   ierr = VecRestoreArray(x, &valX); CHKERRQ(ierr);

   free(indices);
   free(out);
   return ierr;
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  Vec            x, b, u;      /* approx solution, RHS, exact solution */
  Mat            A;            /* linear system matrix */
  KSP            ksp;         /* linear solver context */
  //PC             pc;           /* preconditioner context */
  //PetscReal      norm;  /* norm of solution error */
  PetscErrorCode ierr;
  PetscInt       n = 10;  //i, col[3],its;
  PetscScalar    neg_one      = -1.0,one = 1.0; //value[3];
  //PetscBool      nonzeroguess = PETSC_FALSE;

  PetscInitialize(&argc,&args,(char*)0,help);

  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
  ierr = VecSetSizes(x,PETSC_DECIDE,n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&b);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&u);CHKERRQ(ierr);

  ierr = MatCreateShell(PETSC_COMM_WORLD, n, n, PETSC_DETERMINE, PETSC_DETERMINE, NULL, &A);CHKERRQ(ierr);
  /*
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetType(A, MATSHELL);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);
  */

  // test MatShellSetOp
  ierr = MatShellSetOperation(A, MATOP_MULT, (void (*)(void))userMult);

  /*
  value[0] = -1.0; value[1] = 2.0; value[2] = -1.0;
  for (i=1; i<n-1; i++) {
    col[0] = i-1; col[1] = i; col[2] = i+1;
    ierr   = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
  }
  i    = n - 1; col[0] = n - 2; col[1] = n - 1;
  ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
  i    = 0; col[0] = 0; col[1] = 1; value[0] = 2.0; value[1] = -1.0;
  ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  */

  ierr = VecSet(u,one);CHKERRQ(ierr);
  ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  ierr = MatMult(A,u,b);CHKERRQ(ierr);
  ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);

  /*
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);

  //ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  //ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp,1.e-5,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);

  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

  if (nonzeroguess) {
    PetscScalar p = .5;
    ierr = VecSet(x,p);CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);CHKERRQ(ierr);
  }

  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);

  ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  ierr = VecAXPY(x,neg_one,u);CHKERRQ(ierr);
  ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g, Iterations %D\n",(double)norm,its);CHKERRQ(ierr);
  */

  ierr = VecDestroy(&x);CHKERRQ(ierr); ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr); ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);

  ierr = PetscFinalize();
  return 0;
}
