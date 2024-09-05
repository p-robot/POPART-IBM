/**************************************************************************//**
 * @file checks.h
 * @brief Header file for functions that check the validity of the model.
*****************************************************************************/

#ifndef CHECKS_H_
#define CHECKS_H_

#include "structures.h"

void check_partnership_formation(all_partnerships *, parameters *, debug_struct *, file_struct *);
void check_partnership_formation_and_HIV_acquisition(patch_struct *, int , all_partnerships *,
    output_struct *, debug_struct *, file_struct *);
void check_partnership_dissolution(all_partnerships *, parameters *, debug_struct *, file_struct *);
void make_fake_population(population_size *, stratified_population_size *);
void check_draw_number_partnership(patch_struct *, int);
void check_available_partnerships(population_partners *, population_size_all_patches *);
void check_males_females(stratified_population_size *,individual *);
void print_dob(stratified_population_size *,individual *);
void validate_ages_based_on_age_group(age_list_struct *, int , double );

#endif /* CHECKS_H_ */
