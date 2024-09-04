/**************************************************************************//**
 * @file init.h
 * @brief Header file related to initializing the simulated population
 *****************************************************************************/

#ifndef INIT_H_
#define INIT_H_

#include "demographics.h"
#include "constants.h"
#include "structures.h"
#include "demographics.h"

void get_initial_population_distribution(population_size *, parameters *);
int set_max_n_partners(int , int, int, parameters *);
double make_DoB(int , double, int *);
void set_population_count_zero(population_size*);
void set_population_count_one_year_zero(population_size_one_year_age *n_population);
void set_population_count_stratified(stratified_population_size*, population_size*);
void set_population_count_stratified_zero(stratified_population_size*);
void initialize_child_population(parameters *, child_population_struct *, stratified_population_size *, int, age_list_struct *);
void set_up_population(int, patch_struct *, population *);
void init_available_partnerships(int , patch_struct *, all_partnerships *,population *);
void init_cumulative_counters(cumulative_outputs_struct *);
void init_calendar_counters(calendar_outputs_struct *);
void initialise_debug_variables(debug_struct *);

#endif /* INIT_H_ */
