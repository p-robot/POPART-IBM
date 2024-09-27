/**************************************************************************//**
 * @file partnership.h
 * @brief Header file for functions for partnership formation and dissolution
*****************************************************************************/

#ifndef PARTNERSHIP_H_
#define PARTNERSHIP_H_

#include "structures.h"

void new_partnership(individual* , individual* , int , all_partnerships *, parameters *, 
    debug_struct *, file_struct *);
int time_to_partnership_dissolution(parameters *, int r_m, int r_f, int p_m, int p_f);
void breakup(double, partnership*, all_partnerships *);
void update_list_available_partners_breakup(double , partnership* , population_partners*, 
    population_size_all_patches *);
void add_susceptible_to_list_serodiscordant_partnership(individual* , individual** , long *);
void remove_susceptible_from_list_serodiscordant_partnership(individual* , individual** , long *);
void update_list_susceptibles_in_serodiscordant_partnerships_breakup(partnership* , 
    individual** , long *);
void draw_nb_new_partnerships(patch_struct *, parameters *, int, int);
void draw_n_new_partnerships(int , long, parameters *, int , int , int , int , int *,
        all_partnerships *, patch_struct *, int , int , debug_struct *, file_struct *);
void draw_new_partnerships(int , all_partnerships *, patch_struct *, parameters *, int , int , 
    debug_struct *, file_struct *);

#endif /* PARTNERSHIP_H_ */
