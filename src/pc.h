/**************************************************************************//**
 * @file pc.h
 * @brief Header file for functions for simulating the PC sample
 * @details This code is not currently used in the model
*****************************************************************************/

#ifndef PC_H_
#define PC_H_

#include "constants.h"
#include "structures.h"
#include "utilities.h"
#include "output.h"
#include "debug.h"

int get_PC_HIV_stratum(individual *);

void remove_extras_from_timestep_recruitment(patch_struct *, int , int , int , int , 
    int , int , int );
void create_popart_pc_sample(patch_struct *, age_list_struct *, PC_sample_struct *, 
    parameters *, int , int );
void carry_out_PC_enrolment_per_timestep(int , int , patch_struct *, int , int );
void PC_enroll_person(individual *, patch_struct *, double , int , int , int , int , int );
void PC_next_cohort_round(patch_struct *, int , int );
void carry_out_PC_visits_per_timestep(int , int , patch_struct *, int , int , int );
void PC_visit_person(individual *, patch_struct *,  double , int , int , int );
void reset_visit_counter(patch_struct *, int );
int get_pc_round(int ,int , patch_struct *, int );
#endif /* PC_H_ */
