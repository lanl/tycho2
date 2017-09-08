
/* nonlinear_krylov_accelerator.h */

#ifndef operator_nonlinear_krylov_accelerator_h
#define operator_nonlinear_krylov_accelerator_h


typedef struct nka_state *NKA;
NKA nka_create (int, int, double, double (*dp)(int, double *, double *));
void nka_destroy (NKA);
void nka_correction (NKA, double *);
void nka_restart (NKA);
void nka_relax (NKA);


#endif   /* operator_nonlinear_krylov_accelerator_h */
