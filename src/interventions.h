/**************************************************************************//**
 * @file interventions.h
 * @brief Header file for functions related to PopART/CHiPs interventions
 *****************************************************************************/

#ifndef INTERVENTIONS_H_
#define INTERVENTIONS_H_

#include "structures.h"

void create_popart_chips_samples(age_list_struct *, chips_sample_struct *,
    parameters *, int, int);
void schedule_chips_visits(chips_sample_struct *, parameters *, int);
void carry_out_chips_visits_per_timestep(int, int , patch_struct *,
    int , int , debug_struct *, output_struct *);
void chips_visit_person(individual *, cumulative_outputs_struct *, calendar_outputs_struct *, 
    double ,individual ***, long *, long *, individual ***, long *, long *, parameters *, 
    individual ***, long *, long *, patch_struct *, int , int , debug_struct *, output_struct *,
    int, int);
void draw_if_VMMC(individual *, parameters *, individual ***, long *, long *, double , int );
void schedule_vmmc(individual *,  parameters *, individual ***, long *, long *, double , int );
void schedule_vmmc_healing(individual *, parameters *, individual ***, long *, long *, double );
void schedule_generic_vmmc_event(individual *, parameters *, individual ***, long *, 
    long *, double, double );
void carry_out_VMMC_events_per_timestep(int , double , patch_struct *, int );

#endif /* INTERVENTIONS_H_ */
