/**************************************************************************//**
 * @file input.h
 * @brief Header file related to reading data and parameters
*****************************************************************************/

#ifndef INPUT_H_
#define INPUT_H_

#include "structures.h"

void read_param(char *, parameters **, int, patch_struct *);
void read_patch_info(char *, patch_struct *);
void read_demographic_params(char *, parameters *, int);
void read_hiv_params(char *, parameters *, int, int);
void read_partnership_params(char *, parameters *, int);
void read_time_params(char *, parameters *, int, int);
void read_cascade_params(char *, parameters *, int);
void read_chips_uptake_params(char *, parameters *);
void copy_chips_params(parameters **, int );
void read_pc0_enrolment_params(char *, int , parameters *, int , int );
void copy_pc_params( parameters **, int );
void read_pc_future_params(char *, parameters *, int );
void read_initial_params(char *, parameters *, int);
long get_python_seed(char *);
#endif /* INPUT_H_ */
