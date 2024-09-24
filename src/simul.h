/**************************************************************************//**
 * @file simul.h
 * @brief Header file for functions associated with the main simulation 
 * processes
*****************************************************************************/

#ifndef SIMUL_H_
#define SIMUL_H_

#include "structures.h"

void carry_out_partnership_processes_by_time_step(int , int , patch_struct *, all_partnerships *, 
    output_struct *, debug_struct *, file_struct *);
int carry_out_processes(int, fitting_data_struct *, patch_struct *, all_partnerships *, 
    output_struct *, int, int, debug_struct *, file_struct *, int);
int carry_out_processes_by_patch_by_time_step(int , int , fitting_data_struct *, patch_struct *, 
    int , all_partnerships *, output_struct *, int, int, debug_struct *, file_struct *, int);

#endif /* SIMUL_H_ */
