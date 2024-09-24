/**************************************************************************//**
 * @file fitting.h
 * @brief Header file for functions for fitting routines for the model.
 *****************************************************************************/

#ifndef FITTING_H_
#define FITTING_H_

int load_fitting_data_n(char *);
void load_fitting_data(char *, fitting_data_struct *, int);
void check_fitting_data(fitting_data_struct *, int );
int fit_data(int , int , fitting_data_struct *, patch_struct *, int );
double perform_target_fit(fitting_data_struct *, double );

#endif /* FITTING_H_ */
