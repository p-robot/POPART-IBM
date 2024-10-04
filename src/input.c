/**************************************************************************//**
 * @file input.c
 * @brief Functions for reading data and parameters
 * @details These functions read data and parameters from file and populate
 * relevant structures in the IBM.
*****************************************************************************/

#include "input.h"
#include "constants.h"
#include "utilities.h"

/**************************************************************************//**
 * @brief Call all other functions that read in parameters from file
 * 
 * @details Nothing, other functions for reading parameters are called and 
 * different structures are populated.
 * 
 * @param file_directory Name of dir where parameter files 
 *   ("param_processed*.csv") are stored.
 * @param param Pointer to a parameters struct, in which to store parameter 
 *   values for each simulation run.
 * @param n_runs Number of runs of the simulation to be performed
 * @param patch The (mainly empty) @ref patch_struct structure
 *****************************************************************************/

void read_param(char *file_directory, parameters **param, int n_runs, patch_struct *patch){
    int p;
    char patch_tag[LONGSTRINGLENGTH];
    char patch_number[10];
    
    /* Read in the information for each patch that doesn't change over time - community id, arm.*/
    read_patch_info(file_directory, patch);
    
    for(p = 0; p < NPATCHES; p++){
        sprintf(patch_number, "%i", p);
        strncpy(patch_tag, file_directory, LONGSTRINGLENGTH);
        
        /* Adds a / or \ as needed if working in directory other than current local dir. */
        add_slash(patch_tag);
        strcat(patch_tag, "param_processed_patch");
        strcat(patch_tag, patch_number);
        strcat(patch_tag, "_");

        read_demographic_params(patch_tag, param[p], n_runs);
        read_hiv_params(patch_tag, param[p], n_runs, p);
        read_partnership_params(patch_tag, param[p], n_runs);
        read_time_params(patch_tag, param[p], n_runs, p);
        read_cascade_params(patch_tag, param[p], n_runs);
        
        /* Note that this MUST be called after read_cascade_params(). */
        read_chips_uptake_params(patch_tag, param[p]);
        read_pc0_enrolment_params(patch_tag, patch[p].community_id, param[p], n_runs, p);
        read_pc_future_params(patch_tag, param[p], n_runs);
        
        /* Read in the parameters related to initial conditions. */
        read_initial_params(patch_tag, param[p], n_runs);
    }
    
    // Calling outside the patch loop to avoid a valgrind error where patch[1] is not initialised.
    copy_chips_params(param, n_runs);
    copy_pc_params(param, n_runs);
    return;
}


/**************************************************************************//**
 * @brief Read in the community id and the arm for each patch
 * 
 * @details
 * This function reads the file `param_processed_patchinfo.txt` within the 
 * folder `file_directory` and saved the community ID and trial arm within 
 * that file to the attributes `community_id` and and `trial_arm` of the 
 * patch structure.
 * 
 * Note that trial_arm can be overwritten later on if we are in a counterfactual 
 * scenario (`is_counterfactual = 1`).
 * 
 * Nothing is returned from this function, values are read from file and 
 * stored within attributes of the patch structure.
 * 
 * @param file_directory Name of directory where parameter files 
 *  (`param_processed*.csv`) are stored.
 * @param patch The patch structure (a @ref patch_struct)
 *****************************************************************************/

void read_patch_info(char *file_directory, patch_struct *patch){
    FILE *patchinfo_file;
    char patchinfo_file_name[LONGSTRINGLENGTH];
    int p, checkreadok;
    double temp_int; /* Used to convert float/double to int. */

    // Add path before file name
    strncpy(patchinfo_file_name, file_directory, LONGSTRINGLENGTH);
    
    /* Add a / or \ as needed if working in directory other than current local dir. */
    add_slash(patchinfo_file_name);
    strcat(patchinfo_file_name, "param_processed_patchinfo.txt");

    // Open parameter file
    if ((patchinfo_file = fopen(patchinfo_file_name, "r")) == NULL){
        printf("Cannot open: %s", patchinfo_file_name);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }else{
        if(VERBOSE_OUTPUT == 1){
            printf("Patch info read from: %s:\n",patchinfo_file_name);
        }
    }
    
    // For each patch, read the community ID and store in the patch structure
    for(p = 0; p < NPATCHES; p++){
        checkreadok = fscanf(patchinfo_file, "%lf", &temp_int);
        (patch + p)->community_id = (int) floor(temp_int);
        check_if_cannot_read_param(checkreadok, "patch[p]->community_id");
    }

    // For each patch, read trial arm and store in the patch structure
    for(p = 0; p < NPATCHES; p++){
        checkreadok = fscanf(patchinfo_file, "%lf", &temp_int);
        (patch + p)->trial_arm = (int) floor(temp_int);
        check_if_cannot_read_param(checkreadok, "param_local->trial_arm");
    }
    // Close patch info file
    fclose(patchinfo_file);
    return;
}


/**************************************************************************//**
 * @brief Read parameters related to demographics; copy mortality and 
 * fertility parameters across to all simulation runs.
 * 
 * @details
 * Values from the file param_processed_demographics.csv are read into the 
 * IBM and stored in the allrunparameters structure (an array of 
 * @ref parameters structures, one for each run).  Demographic parameters 
 * may be different across the different simulation runs.  Fertility and 
 * mortality parameters are assumed to be the same across simulation runs 
 * so these parameters are read into the first element of the allrunparameters
 * array (an array of @ref parameters structs) and then copied across to 
 * all other array elements.
 * 
 * The fertility and mortality parameters refer to several macros that are
 * defined within @ref constants.h.  These are:\n
 * `N_UNPD_TIMEPOINTS`\n
 *     the number of time periods for which fertility data is given by 
 *      the UNPD (30)
 * `N_AGE_UNPD_FERTILITY`\n
 *     the number of age groups in which fertility data is specified by 
 *      the UNPD (7)
 * `N_AGE_UNPD_MORTALITY`\n
 *     the number of age groups in which mortality data is specified by 
 *      the UNPD (17)
 * 
 * @param patch_tag pointer to directory and file prefix for 
 *  patch-specific file.  For instance, this will be:\n
 * `./data/SAMPLED_PARAMETERS/PARAMS_COMMUNITY5/param_processed_patch1_`\n
 * and is concatenated with `demographics.csv` or `fertility.csv` to 
 * create the full filenames.
 * 
 * @param allrunparameters pointer to an array of @ref parameters structures, 
 * one for each simulation run (`n_runs`)
 * @param n_runs Number of runs in the simulation
 *****************************************************************************/

void read_demographic_params(char *patch_tag, parameters *allrunparameters, int n_runs){
    FILE * param_file;
    char param_file_name[LONGSTRINGLENGTH];
    int i_run, a_unpd, g, y;
    int checkreadok;
    // This is a temporary var so as not to keep writing allparameters+i_run
    // (or equivalently &allparameters[i_run]). 
    parameters *param_local;

    //  Read in demographic data (apart from mortality/fertility)
    // Add path before file name.
    
    strncpy(param_file_name, patch_tag, LONGSTRINGLENGTH);
    strcat(param_file_name, "demographics.csv");

    // Open parameter file; print error if file not found.
    if((param_file = fopen(param_file_name, "r")) == NULL){
        printf("Cannot open %s", param_file_name);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }else{
        if(VERBOSE_OUTPUT == 1){
            printf("Demographics parameters read from: %s:\n", param_file_name);
        }
    }
    
    // Throw away first line of the file (the header line)
    fscanf(param_file, "%*[^\n]\n");
    
    // Read parameters from each line (i_run) of the file
    for(i_run = 0; i_run < n_runs; i_run++){
        param_local = allrunparameters + i_run;
        checkreadok = fscanf(param_file, "%lg", &(param_local->sex_ratio));
        check_if_cannot_read_param(checkreadok, "param_local->sex_ratio");
    }
    fclose(param_file);     // Closing demographics parameter file
    
    //  Mortality data
    // Add path before file name.
    strncpy(param_file_name, patch_tag, LONGSTRINGLENGTH);
    strcat(param_file_name, "mortality.csv");
    
    // Open parameter file; print error if file not found.
    if((param_file = fopen(param_file_name, "r")) == NULL){
        printf("Cannot open %s", param_file_name);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }else{
        if(VERBOSE_OUTPUT == 1){
            printf("Demographic mortality parameters read from: %s:\n", param_file_name);
        }
    }
    
    // Read mortality parameters
    // Store outputs in the i_run=0 parameter struct; later copy this to every other run.
    i_run = 0;
    param_local = allrunparameters + i_run;
    for(g = 0; g < N_GENDER; g++){
        for(a_unpd = 0; a_unpd < N_AGE_UNPD_MORTALITY; a_unpd++){
            checkreadok = fscanf(param_file, "%lg", 
                &(param_local->mortality_rate_by_gender_age_intercept[g][a_unpd]));
            check_if_cannot_read_param(checkreadok, 
                "param_local->mortality_rate_by_gender_age_intercept[g][a_unpd]");
        }
        for(a_unpd = 0; a_unpd < N_AGE_UNPD_MORTALITY; a_unpd++){
            param_local = allrunparameters + i_run;
            checkreadok = fscanf(param_file, "%lg", 
                &(param_local->mortality_rate_by_gender_age_slope[g][a_unpd]));
            check_if_cannot_read_param(checkreadok, 
                "param_local->mortality_rate_by_gender_age_slope[g][a_unpd]");
        }
    }
    // Close mortality parameter file
    fclose(param_file);
    
    // Copy parameters from the first simulation run across to all simulation runs
    // allrunparameters stores a copy of parameters used for each simulation run
    for(g = 0; g < N_GENDER; g++){
        for(a_unpd = 0; a_unpd < N_AGE_UNPD_MORTALITY; a_unpd++){
            for(i_run = 1; i_run < n_runs; i_run++){
                
                allrunparameters[i_run].mortality_rate_by_gender_age_intercept[g][a_unpd] =
                    allrunparameters[0].mortality_rate_by_gender_age_intercept[g][a_unpd];
                
                allrunparameters[i_run].mortality_rate_by_gender_age_slope[g][a_unpd] =
                    allrunparameters[0].mortality_rate_by_gender_age_slope[g][a_unpd];
            }
        }
    }

    // Fertility data
    // Add path before file name
    strncpy(param_file_name, patch_tag, LONGSTRINGLENGTH);
    strcat(param_file_name, "fertility.csv");

    // Open parameter file; print error if file not found
    if((param_file = fopen(param_file_name, "r")) == NULL){
        printf("Cannot open %s", param_file_name);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }else{
        if(VERBOSE_OUTPUT == 1){
            printf("Demographic mortality parameters read from: %s:\n",param_file_name);
        }
    }
    
    // Read parameters; throw error if not read correctly
    i_run = 0;
    param_local = allrunparameters + i_run;
    for(y = 0; y < N_UNPD_TIMEPOINTS; y++){
        for(a_unpd = 0; a_unpd < N_AGE_UNPD_FERTILITY; a_unpd++){
            
            checkreadok = fscanf(param_file, "%lg", 
                &(param_local->fertility_rate_by_age[a_unpd][y]));
            check_if_cannot_read_param(checkreadok,
                "param_local->fertility_rate_by_age[a_unpd][y]");
        }
    }
    // Close mortality parameter file
    fclose(param_file);
    
    // Copy parameters from the first simulation run across to all simulation runs
    // allrunparameters stores a copy of parameters used for each simulation run
    for(y = 0; y < N_UNPD_TIMEPOINTS; y++){
        for(a_unpd = 0; a_unpd < N_AGE_UNPD_FERTILITY; a_unpd++){
            for(i_run = 1; i_run < n_runs; i_run++){
                
                allrunparameters[i_run].fertility_rate_by_age[a_unpd][y] =
                    allrunparameters[0].fertility_rate_by_age[a_unpd][y];
            }
        }
    }
    return;
}


/**************************************************************************//**
 * @brief Read parameters related to HIV
 * 
 * @details
 * This function also does some small error checking.
 * 
 * For a given true CD4, p_misclassify_cd4[j] stores the probability that 
 * the measured cd4 cat is j.  Note that icd4 is the true cd4, and the other
 * index is the measured cd4.
 * 
 * @param patch_tag Pointer to a character array for denoting the patch
 * @param allrunparameters pointer to an array of 'parameters' structures, 
one for each simulation run (n_runs)
 * @param n_runs Number of runs of the simulation
 * @param p Patch number
 *****************************************************************************/

void read_hiv_params(char *patch_tag, parameters *allrunparameters, int n_runs, int p){
    FILE * param_file;
    char param_file_name[LONGSTRINGLENGTH];
    int spvl, icd4, checkreadok, i_run;
    double p_misclassify_cd4[NCD4], temp;
    
    // Temp variable to avoid having to write allparameters+i_run 
    //(or equivalently &allparameters[i_run]).
    parameters *param_local;

    strncpy(param_file_name, patch_tag, LONGSTRINGLENGTH);
    strcat(param_file_name, "HIV.csv");

    // Open parameter file
    if((param_file = fopen(param_file_name, "r")) == NULL){
        printf("Cannot open %s", param_file_name);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }else{
        if(VERBOSE_OUTPUT == 1){
            printf("HIV parameters read from: %s:\n", param_file_name);
        }
    }
    // Throw away the first line of the parameter file (the header line)
    fscanf(param_file, "%*[^\n]\n");

    // Read parameters from each line (i_run) of the file
    for(i_run = 0; i_run < n_runs; i_run++){
        param_local = allrunparameters + i_run;

        checkreadok = fscanf(param_file,"%lg", &(param_local->p_child_circ));
        check_if_cannot_read_param(checkreadok, "param_local->p_child_circ");

        checkreadok = fscanf(param_file,"%lg", &(param_local->eff_circ_vmmc));
        check_if_cannot_read_param(checkreadok, "param_local->eff_circ_vmmc");

        checkreadok = fscanf(param_file,"%lg", &(param_local->eff_circ_tmc));
        check_if_cannot_read_param(checkreadok, "param_local->eff_circ_tmc");

        checkreadok = fscanf(param_file,"%lg", &(param_local->rr_circ_unhealed));
        check_if_cannot_read_param(checkreadok, "param_local->rr_circ_unhealed");

        checkreadok = fscanf(param_file,"%lg", &(param_local->t0_pmtct));
        check_if_cannot_read_param(checkreadok, "param_local->t0_pmtct");

        checkreadok = fscanf(param_file,"%lg", &(param_local->t50_pmtct));
        check_if_cannot_read_param(checkreadok, "param_local->t50_pmtct");

        checkreadok = fscanf(param_file,"%lg", &(param_local->average_log_viral_load));
        check_if_cannot_read_param(checkreadok, "param_local->average_log_viral_load");

        checkreadok = fscanf(param_file,"%lg", &(param_local->average_annual_hazard));
        check_if_cannot_read_param(checkreadok, "param_local->average_annual_hazard");

        checkreadok = fscanf(param_file,"%lg", &(param_local->RRacute_trans));
        check_if_cannot_read_param(checkreadok, "param_local->RRacute_trans");

        checkreadok = fscanf(param_file,"%lg", &(param_local->RRmale_to_female_trans));
        check_if_cannot_read_param(checkreadok, "param_local->RRmale_to_female_trans");

        checkreadok = fscanf(param_file, "%lg %lg %lg %lg", 
            &(param_local->RRCD4[0]), &(param_local->RRCD4[1]), 
            &(param_local->RRCD4[2]), &(param_local->RRCD4[3]));
        check_if_cannot_read_param(checkreadok, "param_local->RRCD4");

        checkreadok = fscanf(param_file,"%lg", &(param_local->SPVL_beta_k));
        check_if_cannot_read_param(checkreadok, "param_local->SPVL_beta_k");

        checkreadok = fscanf(param_file,"%lg", &(param_local->SPVL_beta_50));
        check_if_cannot_read_param(checkreadok, "param_local->SPVL_beta_50");

        // Calc RR in infectivity from effectiveness of initial ART, ART when VS and ART when not VS
        // Reduction in infectivity is 1 - (effectiveness)
        
        checkreadok = fscanf(param_file, "%lg", &temp);
        check_if_cannot_read_param(checkreadok, "param_local->RR_ART_INITIAL");
        param_local->RR_ART_INITIAL = 1.0 - temp;

        checkreadok = fscanf(param_file, "%lg", &temp);
        check_if_cannot_read_param(checkreadok, "param_local->RR_ART_VS");
        param_local->RR_ART_VS = 1.0 - temp;

        checkreadok = fscanf(param_file, "%lg", &temp);
        check_if_cannot_read_param(checkreadok, "param_local->RR_ART_VU");
        param_local->RR_ART_VU = 1.0 - temp;

        checkreadok = fscanf(param_file, "%lg %lg", 
            &(param_local->min_dur_acute), &(param_local->max_dur_acute));
        check_if_cannot_read_param(checkreadok,"param_local->dur_acute");

        for(spvl = 0; spvl < NSPVL; spvl++){
            checkreadok = fscanf(param_file, "%lg", &(param_local->p_initial_cd4_gt500[spvl]));
            check_if_cannot_read_param(checkreadok, "param_local->p_initial_cd4_gt500");

            checkreadok = fscanf(param_file,"%lg", &(param_local->p_initial_cd4_350_500[spvl]));
            check_if_cannot_read_param(checkreadok, "param_local->p_initial_cd4_350_500");

            checkreadok = fscanf(param_file,"%lg", &(param_local->p_initial_cd4_200_350[spvl]));
            check_if_cannot_read_param(checkreadok, "param_local->p_initial_cd4_200_350");

            checkreadok = fscanf(param_file,"%lg", &(param_local->p_initial_cd4_lt200[spvl]));
            check_if_cannot_read_param(checkreadok, "param_local->p_initial_cd4_lt200");

            // Ensure these are normalized to sum to 1
            normalise_four_quantities(
                &param_local->p_initial_cd4_gt500[spvl],
                &param_local->p_initial_cd4_350_500[spvl],
                &param_local->p_initial_cd4_200_350[spvl], 
                &param_local->p_initial_cd4_lt200[spvl]);
                
            // Now to save recomputing these, make cumulative sums
            cumulative_four_quantities(
                param_local->p_initial_cd4_gt500[spvl], 
                param_local->p_initial_cd4_350_500[spvl], 
                param_local->p_initial_cd4_200_350[spvl],
                param_local->p_initial_cd4_lt200[spvl],
                &param_local->cumulative_p_initial_cd4_gt500[spvl],
                &param_local->cumulative_p_initial_cd4_350_500[spvl],
                &param_local->cumulative_p_initial_cd4_200_350[spvl],
                &param_local->cumulative_p_initial_cd4_lt200[spvl]); 
                
            if(PRINT_DEBUG_INPUT == 1){
                printf("Cumulative values = %f %f %f %f\n",
                    param_local->cumulative_p_initial_cd4_gt500[spvl],
                    param_local->cumulative_p_initial_cd4_350_500[spvl],
                    param_local->cumulative_p_initial_cd4_200_350[spvl],
                    param_local->cumulative_p_initial_cd4_lt200[spvl]);
            }
        }

        checkreadok = fscanf(param_file, "%lg", &(param_local->initial_SPVL_mu));
        check_if_cannot_read_param(checkreadok, "param_local->initial_SPVL_mu");
        checkreadok = fscanf(param_file, "%lg", &(param_local->initial_SPVL_sigma));
        check_if_cannot_read_param(checkreadok, "param_local->initial_SPVL_sigma");
        checkreadok = fscanf(param_file, "%lg", &(param_local->SPVL_sigma_M));
        check_if_cannot_read_param(checkreadok, "param_local->SPVL_sigma_M");
        checkreadok = fscanf(param_file, "%lg", &(param_local->SPVL_sigma_E));
        check_if_cannot_read_param(checkreadok, "param_local->SPVL_sigma_E");

        // For a given true CD4, p_misclassify_cd4[j] stores the probability that the 
        // measured cd4 cat is j.
        // Note that icd4 is the true cd4, and the other index is the measured cd4.
        for(icd4 = 0; icd4 < NCD4; icd4++){
            checkreadok = fscanf(param_file, "%lg %lg %lg %lg", 
                &(p_misclassify_cd4[0]), &(p_misclassify_cd4[1]), 
                &(p_misclassify_cd4[2]), &(p_misclassify_cd4[3]));
            
            check_if_cannot_read_param(checkreadok, "param_local->p_misclassify_cd4");
            
            normalise_four_quantities(
                &(p_misclassify_cd4[0]), &(p_misclassify_cd4[1]), 
                &(p_misclassify_cd4[2]), &(p_misclassify_cd4[3]));
            
            cumulative_four_quantities(
                p_misclassify_cd4[0], p_misclassify_cd4[1], 
                p_misclassify_cd4[2], p_misclassify_cd4[3],
                &param_local->cumulative_p_misclassify_cd4[icd4][0],
                &param_local->cumulative_p_misclassify_cd4[icd4][1],
                &param_local->cumulative_p_misclassify_cd4[icd4][2],
                &param_local->cumulative_p_misclassify_cd4[icd4][3]);
                
            if(PRINT_DEBUG_INPUT == 1){
                printf("Cumulative CD4 misclassification values: %lg %lg %lg %lg\n", 
                    param_local->cumulative_p_misclassify_cd4[icd4][0], 
                    param_local->cumulative_p_misclassify_cd4[icd4][1], 
                    param_local->cumulative_p_misclassify_cd4[icd4][2], 
                    param_local->cumulative_p_misclassify_cd4[icd4][3]);
            }
        }   

        for(spvl = 0; spvl < NSPVL; spvl++){
            for(icd4 = 0; icd4 < NCD4; icd4++){
                checkreadok = fscanf(param_file, "%lg ", &param_local->time_hiv_event[icd4][spvl]);
                check_if_cannot_read_param(checkreadok, "param_local->time_hiv_event");
            }
        }

        for(icd4 = 0; icd4 < NCD4; icd4++){
            checkreadok = fscanf(param_file, "%lg ",
                &param_local->CoxPH_SPVL_RR_CD4progression[icd4]);
            
            check_if_cannot_read_param(checkreadok, "param_local->CoxPH_SPVL_RR_CD4progression");
        }
        for(icd4 = 0; icd4 < NCD4; icd4++){
            if(param_local->CoxPH_SPVL_RR_CD4progression[icd4] < 1){
                
                printf("param_local->CoxPH_SPVL_RR_CD4progression[icd4] = ");
                printf("%lf was below 1. Have forced it to be 1.\n",
                    param_local->CoxPH_SPVL_RR_CD4progression[icd4]);
                
                param_local->CoxPH_SPVL_RR_CD4progression[icd4] = 1.0;
            }
        }
        
        checkreadok = fscanf(param_file, "%lg", 
            &(param_local->factor_for_slower_progression_ART_VU));
        check_if_cannot_read_param(checkreadok, 
            "param_local->factor_for_slower_progression_ART_VU");
    }
    // Close parameter file
    fclose(param_file);
    return;
}


/**************************************************************************//**
 * @brief Read partnership related parameters
 * 
 * @details
 * The parameters p_age_per_gender are not directly read in the input files
 * what is read is p1, p2/(1-p1), p3/(1-p1-p2), etc...
 * this avoids having constraints on these probabilities summing to 1.
 * 
 * @param patch_tag pointer to a character array
 * @param allrunparameters pointer to an array of parameters structure 
 *     for storing all parameter values
 * @param n_runs Number of runs in the simulation
 *****************************************************************************/

void read_partnership_params(char *patch_tag, parameters *allrunparameters, int n_runs){
    FILE * param_file;
    int g, ag, r, bg, i_run;
    char param_file_name[LONGSTRINGLENGTH];
    
    // This is a temporary var so as not to keep writing allparameters+i_run
    // (or equivalently &allparameters[i_run]). 
    parameters *param_local;
    int checkreadok;
    double c_multiplier, breakup_scale_multiplier_overall, temp_int;
    double breakup_scale_multiplier_between_vs_within_patch;
    
    strncpy(param_file_name, patch_tag, LONGSTRINGLENGTH);
    strcat(param_file_name, "partnerships.csv");

    // Open parameter file
    if((param_file = fopen(param_file_name, "r")) == NULL){
        printf("Cannot open %s\n", param_file_name);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }else{
        if(VERBOSE_OUTPUT == 1){
            printf("Partnership parameters read from: %s:\n", param_file_name);
        }
    }

    // Throw away the first line (the header file)
    fscanf(param_file, "%*[^\n]\n");

    // Read parameters from each line (i_run), one for each simulation run
    for(i_run = 0; i_run < n_runs; i_run++){

        param_local = allrunparameters + i_run;
        
        checkreadok = fscanf(param_file,"%lg", &(param_local->assortativity));
        check_if_cannot_read_param(checkreadok, "param_local->assortativity");

        checkreadok = fscanf(param_file,"%lg", &(param_local->prop_compromise_from_males));
        check_if_cannot_read_param(checkreadok, "param_local->prop_compromise_from_males");

        for(ag = 0; ag < N_AGE ; ag++){
            checkreadok = fscanf(param_file, "%lg",
                &(param_local->c_per_gender_within_patch[FEMALE][ag]));
            check_if_cannot_read_param(checkreadok, "param_local->c_per_gender_within_patch");
        }

        for(ag = 0; ag < N_AGE; ag++){
            checkreadok = fscanf(param_file, "%lg", 
                &(param_local->c_per_gender_within_patch[MALE][ag]));
            check_if_cannot_read_param(checkreadok, "param_local->c_per_gender_within_patch");
        }

        checkreadok = fscanf(param_file, "%lg", &c_multiplier);
        check_if_cannot_read_param(checkreadok, "c_multiplier");
        
        for(ag = 0; ag < N_AGE; ag++){
            param_local->c_per_gender_within_patch[FEMALE][ag] *= c_multiplier;
            param_local->c_per_gender_within_patch[MALE][ag] *= c_multiplier;
        }

        checkreadok = fscanf(param_file, "%lg",
            &(param_local->rel_rate_partnership_formation_between_patches));
        check_if_cannot_read_param(checkreadok,
            "param_local->rel_rate_partnership_formation_between_patches");

        for(ag = 0; ag < N_AGE; ag++){
            param_local->c_per_gender_between_patches[FEMALE][ag] =
                param_local->rel_rate_partnership_formation_between_patches *
                param_local->c_per_gender_within_patch[FEMALE][ag] / (NPATCHES-1);
            
            param_local->c_per_gender_between_patches[MALE][ag] =
                param_local->rel_rate_partnership_formation_between_patches *
                param_local->c_per_gender_within_patch[MALE][ag] / (NPATCHES-1);
        }

        checkreadok = fscanf(param_file, "%lg", &(param_local->rr_hiv_between_vs_within_patch));
        check_if_cannot_read_param(checkreadok, "param_local->rr_hiv_between_vs_within_patch");

        for(r = 0; r < N_RISK ; r++){
            checkreadok = fscanf(param_file, "%lg", 
                &(param_local->relative_number_partnerships_per_risk[r]));
            check_if_cannot_read_param(checkreadok, 
                "param_local->relative_number_partnerships_per_risk");
            
            // Because we're now defining everything relative to the low risk group we need this to
            // be forced to 1
            if(r == 0){
                param_local->relative_number_partnerships_per_risk[r] = 1.0;
            }
        }

        // The parameters p_age_per_gender are not directly read in the input files 
        // what is read is p1, p2/(1-p1), p3/(1-p1-p2), etc... 
        // this avoids having constraints on these probabilities summing to 1.  

        double p_age_per_gender_wrong_scale[N_GENDER][N_AGE][N_AGE], fact, check_sum;

        for(g = 0; g < N_GENDER; g++){
            for(ag = 0; ag < N_AGE; ag++){
                //printf("--- INDEX PERSON IN AGE GROUP %d; probability of choosing partner in each age group is:\n",ag);
                fact = 1;
                
                check_sum = 0;
                for(bg = 0; bg < N_AGE; bg++){
                    
                    checkreadok = fscanf(param_file, "%lg", 
                        &(p_age_per_gender_wrong_scale[g][ag][bg]));
                    check_if_cannot_read_param(checkreadok, "param_local->p_age_per_gender");
                    
                    param_local->p_age_per_gender[g][ag][bg] =
                        p_age_per_gender_wrong_scale[g][ag][bg] * fact;
                    
                    fact -= param_local->p_age_per_gender[g][ag][bg];
                    
                    //printf("%lg\t",param_local->p_age_per_gender[g][ag][bg]);
                    check_sum += param_local->p_age_per_gender[g][ag][bg];
                }
                
                if(fabs(check_sum - 1.0) > 1e-12){
                    printf("In input.c the sum of param_local->p_age_per_gender[g][ag][bg] ");
                    printf("over bg is 1-%e for ag=%d. Exiting", check_sum - 1.0, ag);
                    printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                    fflush(stdout);
                    exit(1);
                }
            }
        }

        for(r = 0; r < N_RISK; r++){
            checkreadok = fscanf(param_file, "%lf", &(temp_int));
            param_local->max_n_part_noage[r] = (int) floor(temp_int);
            check_if_cannot_read_param(checkreadok, "param_local->max_n_part_noage");
        }

        for(r = 0; r < N_RISK; r++){
            checkreadok = fscanf(param_file, "%lg", 
                &(param_local->breakup_scale_lambda_within_patch[r]));
            check_if_cannot_read_param(checkreadok,
                "param_local->breakup_scale_lambda_within_patch");
        }
        for(r = 0; r < N_RISK; r++){
            checkreadok = fscanf(param_file, "%lg", &(param_local->breakup_shape_k[r]));
            check_if_cannot_read_param(checkreadok, "param_local->breakup_shape_k");
        }
        
        checkreadok = fscanf(param_file, "%lg", &breakup_scale_multiplier_overall);
        check_if_cannot_read_param(checkreadok, "breakup_scale_multiplier_overall");
        
        for(r = 0; r < N_RISK; r++){
            param_local->breakup_scale_lambda_within_patch[r] *= breakup_scale_multiplier_overall;
        }

        checkreadok = fscanf(param_file, "%lg", &breakup_scale_multiplier_between_vs_within_patch);
        check_if_cannot_read_param(checkreadok, "breakup_scale_multiplier_between_vs_within_patch");
        
        for(r = 0; r < N_RISK; r++){
            param_local->breakup_scale_lambda_between_patch[r] =
                param_local->breakup_scale_lambda_within_patch[r] *
                breakup_scale_multiplier_between_vs_within_patch;
        }
    } // end for i_run
    // Close parameter file
    fclose(param_file);
    return;
}


/**************************************************************************//**
 * @brief Read time related parameters
 * 
 * @details Nothing returned; the argument allrunparameters structure is 
 * populated (a @ref parameters structure)
 * 
 * @param patch_tag Pointer to a character array of the patch tag
 * @param allrunparameters Pointer to an array of @ref parameters structure for storing all parameter values
 * @param n_runs Number of simulation runs
 * @param p Patch number
 *****************************************************************************/

void read_time_params(char *patch_tag, parameters *allrunparameters, int n_runs, int p){
    FILE * param_file;
    char param_file_name[LONGSTRINGLENGTH];
    int i_run, i;
    
    // For rounding the time HIV starts to the nearest timestep
    double time_start_hiv_undiscretised;
    double CHIPS_START, CHIPS_END;
    
    // What is wanted is not the time of the end of the last timestep in a CHiPS round but the time
    // at the beginning of that timestep. i.e. CHIPS_END_STARTOFTIMESTEP = CHIPS_END-TIME_STEP;
    double CHIPS_END_STARTOFTIMESTEP;
    
    // This is a temporary var so as not to keep writing allparameters+i_run
    // (or equivalently &allparameters[i_run]). 
    parameters *param_local;
    int checkreadok;

    double temp_int;

    strncpy(param_file_name, patch_tag, LONGSTRINGLENGTH);
    strcat(param_file_name, "times.csv");

    // Open parameter file
    if((param_file = fopen(param_file_name, "r")) == NULL){
        printf("Cannot open %s\n", param_file_name);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }else{
        if(VERBOSE_OUTPUT == 1){
            printf("Times parameters read from: %s:\n", param_file_name);
        }
    }

    // Throw away first line of the file (the header line)
    fscanf(param_file, "%*[^\n]\n");

    // Read parameters from each line (i_run) of the file as used for each run of the simulation
    for(i_run = 0; i_run < n_runs; i_run++){
        param_local = allrunparameters + i_run;

        checkreadok = fscanf(param_file, "%lg", &time_start_hiv_undiscretised);
        check_if_cannot_read_param(checkreadok, "param_local->start_time_hiv");
        
        // Here we want start_time_hiv to be equal to a timestep value (this avoids issues when
        // seeding infections, where new_infections looks for the person in age_index[] using the
        // index calculated using the infection time. If not a timestep then that may produce the
        // wrong index and the person is not found (so HIV status not updated etc).
        
        // The year when HIV is seeded
        param_local->start_time_hiv_discretised_year = (int) floor(time_start_hiv_undiscretised);
        
        // The timestep within that year when HIV is seeded
        param_local->start_time_hiv_discretised_timestep = (int) floor((time_start_hiv_undiscretised - param_local->start_time_hiv_discretised_year)*N_TIME_STEP_PER_YEAR);
        
        if(param_local->start_time_hiv_discretised_timestep >= N_TIME_STEP_PER_YEAR){
            printf("Error in calculating when HIV is seeded. Exiting\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        param_local->start_time_hiv = param_local->start_time_hiv_discretised_year +
            param_local->start_time_hiv_discretised_timestep / (1.0 * N_TIME_STEP_PER_YEAR);
        
        checkreadok = fscanf(param_file, "%lf", &temp_int);
        param_local->start_time_simul = (int) floor(temp_int);
        check_if_cannot_read_param(checkreadok, "param_local->start_time_simul");
        
        checkreadok = fscanf(param_file, "%lf", &temp_int);
        param_local->end_time_simul = (int) floor(temp_int);
        check_if_cannot_read_param(checkreadok, "param_local->end_time_simul");

        // If this message comes up make sure that age_list is correctly indexed
        if(
        ((param_local->end_time_simul-param_local->start_time_simul) >= 2*(MAX_AGE - AGE_ADULT)) &&
        i_run == 0 && 
        p == 0
        ){
            printf("Potential issue: param->end_time_simul-param->start_time_simul is bigger than");
            printf(" 2*(MAX_AGE-AGE_ADULT)=%i.\n", 2 * (MAX_AGE - AGE_ADULT));
            printf("There may be issues with the age index for the age_list[] array.\n");
        }
        
        checkreadok = fscanf(param_file, "%lf", &temp_int);
        param_local->COUNTRY_HIV_TEST_START = (int) floor(temp_int);
        check_if_cannot_read_param(checkreadok, "param_local->COUNTRY_HIV_TEST_START");

        checkreadok = fscanf(param_file, "%lg", &(param_local->COUNTRY_ART_START));
        check_if_cannot_read_param(checkreadok, "param_local->COUNTRY_ART_START");

        checkreadok = fscanf(param_file, "%lg", &(param_local->COUNTRY_CD4_350_START));
        check_if_cannot_read_param(checkreadok, "param_local->COUNTRY_CD4_350_START");

        checkreadok = fscanf(param_file, "%lg", &(param_local->COUNTRY_CD4_500_START));
        check_if_cannot_read_param(checkreadok, "param_local->COUNTRY_CD4_500_START");
        
        checkreadok = fscanf(param_file, "%lg", &(param_local->COUNTRY_IMMEDIATE_ART_START));
        check_if_cannot_read_param(checkreadok, "param_local->COUNTRY_IMMEDIATE_ART_START");

        checkreadok = fscanf(param_file, "%lg", &(param_local->COUNTRY_VMMC_START));
        check_if_cannot_read_param(checkreadok, "param_local->COUNTRY_VMMC_START");

        for(i = 0; i < NCHIPSROUNDS; i++){
            checkreadok = fscanf(param_file, "%lg", &CHIPS_START);
            check_if_cannot_read_param(checkreadok, "CHIPS_START");
            
            param_local->CHIPS_START_YEAR[i] = (int) CHIPS_START;
            param_local->CHIPS_START_TIMESTEP[i] = (int) round((CHIPS_START - 
                param_local->CHIPS_START_YEAR[i])*N_TIME_STEP_PER_YEAR);

            checkreadok = fscanf(param_file, "%lg", &CHIPS_END);
            check_if_cannot_read_param(checkreadok, "CHIPS_END");
            
            // What is wanted is not the time of the end of the last timestep in a CHiPS round but
            // the time at the beginning of that timestep. i.e. CHIPS_END_STARTOFTIMESTEP =
            // CHIPS_END-TIME_STEP;
            
            CHIPS_END_STARTOFTIMESTEP = CHIPS_END - TIME_STEP;
            param_local->CHIPS_END_YEAR[i] = (int) CHIPS_END_STARTOFTIMESTEP;
            param_local->CHIPS_END_TIMESTEP[i] = (int) round((CHIPS_END_STARTOFTIMESTEP -
                param_local->CHIPS_END_YEAR[i]) * N_TIME_STEP_PER_YEAR);
            
            /* This is derived from the above times. */
            //printf("STart = %i %i End = %i %i\n",param_local->CHIPS_START_YEAR[i],param_local->CHIPS_START_TIMESTEP[i],param_local->CHIPS_END_YEAR[i],param_local->CHIPS_END_TIMESTEP[i]);
            param_local->chips_params->n_timesteps_per_round[i] =
                (param_local->CHIPS_END_YEAR[i]-param_local->CHIPS_START_YEAR[i]) *
                N_TIME_STEP_PER_YEAR +
                (param_local->CHIPS_END_TIMESTEP[i] - param_local->CHIPS_START_TIMESTEP[i]) + 1;
        }

        param_local->CHIPS_START_TIMESTEP_POSTTRIAL = 
            param_local->CHIPS_END_TIMESTEP[NCHIPSROUNDS - 1] + 1;
        
        if(param_local->CHIPS_START_TIMESTEP_POSTTRIAL >= N_TIME_STEP_PER_YEAR){
            param_local->CHIPS_START_TIMESTEP_POSTTRIAL =
                param_local->CHIPS_START_TIMESTEP_POSTTRIAL - N_TIME_STEP_PER_YEAR;
        }

        checkreadok = fscanf(param_file, "%lf", &temp_int);
        param_local->chips_params->n_timesteps_per_round_posttrial = (int) floor(temp_int);
        check_if_cannot_read_param(checkreadok, "param_local->n_timesteps_per_round_posttrial");

        if( (int) (param_local->start_time_simul) != (param_local->start_time_simul) || 
            (int) (param_local->end_time_simul) != param_local->end_time_simul){
            
            printf("start_time_simul and end_time_simul defined in '%s' ", param_file_name);
            printf("should be integers. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        if(param_local->end_time_simul - param_local->start_time_simul > MAX_N_YEARS){
            printf("MAX_N_YEARS < param_local->end_time_simul - param_local->start_time_simul.");
            printf(" Need to increase the value of MAX_N_YEARS. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        checkreadok = fscanf(param_file, "%lg", &temp_int);
        check_if_cannot_read_param(checkreadok, "param_local->DHS_params->NDHSROUNDS");
        param_local->DHS_params->NDHSROUNDS = (int) temp_int;

        for(i = 0; i < param_local->DHS_params->NDHSROUNDS; i++){
            checkreadok = fscanf(param_file, "%lg", &temp_int);
            check_if_cannot_read_param(checkreadok, "param_local->DHS_params->DHS_YEAR[round]");
            param_local->DHS_params->DHS_YEAR[i] = (int) temp_int;
        }

        checkreadok = fscanf(param_file, "%lf", &temp_int);
        param_local->write_annual_population_start = (int) floor(temp_int);
        check_if_cannot_read_param(checkreadok, "param_local->write_annual_population_start");

        checkreadok = fscanf(param_file, "%lf", &temp_int);
        param_local->write_annual_population_end = (int) floor(temp_int);
        check_if_cannot_read_param(checkreadok, "param_local->write_annual_population_end");
    }
    // Close parameter file
    fclose(param_file);
    return;
}


/**************************************************************************//**
 * @brief Read cascade parameters, except CHiPs visit rates
 * 
 * @param patch_tag Pointer to a character array of the patch tag
 * @param allrunparameters Pointer to an array of @ref parameters structure for storing all parameter values
 * @param n_runs Number of simulation runs
 *****************************************************************************/

void read_cascade_params(char *patch_tag, parameters *allrunparameters, int n_runs){
    FILE *param_file;
    int icd4, iround;
    char param_file_name[LONGSTRINGLENGTH];
    int i_run,i;
    double max_temp;  /* Local variable which allows us to get the range from the max value. */
    /* This is a local temp variable we use so we don't have to keep writing allparameters+i_run (or equivalently &allparameters[i_run]). */
    parameters *param_local;
    int checkreadok;

    strncpy(param_file_name,patch_tag,LONGSTRINGLENGTH);
    strcat(param_file_name, "cascade.csv");

    // Open parameter file
    if ((param_file=fopen(param_file_name,"r"))==NULL){
        printf("Cannot open %s\n",param_file_name);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }else{
        if(VERBOSE_OUTPUT==1)
            printf("Cascade parameters read from: %s:\n",param_file_name);
    }

    /* The first line of this file is a header - the text below reads in just that line and throws it away.
     * The * instructs fscanf (all of the scanf family, in fact) to parse the data out as presented in the format string,
     * but NOT to store it at any target address (which is good, because there is none provided in an argument list).
     * the [^\n] means take anything except a newline, so ALL data will be consumed up to (but not including) the newline.
     *  Finally, the final \n means "and consume (and ignore) the newline" (which we just stopped at when fulfilling the
     *  prior format spec). */
    fscanf(param_file, "%*[^\n]\n");

    // Read parameters from each line i_run
    for (i_run = 0; i_run<n_runs; i_run++){
        param_local = allrunparameters + i_run;

        checkreadok = fscanf(param_file,"%lg",&(param_local->time_to_background_HIVtestNOW));
        //printf("param_local->time_to_background_HIVtestNOW = %lf\n",param_local->time_to_background_HIVtestNOW);
        check_if_cannot_read_param(checkreadok,"param_local->time_to_background_HIVtestNOW");

        checkreadok = fscanf(param_file,"%lg",&(param_local->time_to_background_HIVtest_maxval));
        //printf("param_local->time_to_background_HIVtest_maxval = %lf\n",param_local->time_to_background_HIVtest_maxval);
        check_if_cannot_read_param(checkreadok,"param_local->time_to_background_HIVtest_maxval");

        checkreadok = fscanf(param_file,"%lg",&(param_local->time_to_background_HIVtest_exponent));
        //printf("param_local->time_to_background_HIVtest_exponent = %lf\n",param_local->time_to_background_HIVtest_exponent);
        check_if_cannot_read_param(checkreadok,"param_local->time_to_background_HIVtest_exponent");

        checkreadok = fscanf(param_file,"%lg",&(param_local->time_to_background_HIVtest_midpoint));
        //printf("param_local->time_to_background_HIVtest_midpoint = %lf\n",param_local->time_to_background_HIVtest_midpoint);
        check_if_cannot_read_param(checkreadok,"param_local->time_to_background_HIVtest_midpoint");

        /* Input probabilities for the cascade events: */
        checkreadok = fscanf(param_file,"%lg",&(param_local->p_HIV_background_testing_female_pre2006));
        check_if_cannot_read_param(checkreadok,"param_local->p_HIV_background_testing_female_pre2006");
        
        checkreadok = fscanf(param_file,"%lg",&(param_local->p_HIV_background_testing_female_current));
        check_if_cannot_read_param(checkreadok,"param_local->p_HIV_background_testing_female_current");

        checkreadok = fscanf(param_file,"%lg",&(param_local->RR_HIV_background_testing_male));
        check_if_cannot_read_param(checkreadok,"param_local->RR_HIV_background_testing_male");

        checkreadok = fscanf(param_file,"%lg",&(param_local->HIV_rapid_test_sensitivity_CHIPS));
        check_if_cannot_read_param(checkreadok,"param_local->HIV_rapid_test_sensitivity_CHIPS");

        /* Input probabilities for the cascade events: */
        checkreadok = fscanf(param_file,"%lg",&(param_local->p_collect_hiv_test_results_cd4_over200));
        check_if_cannot_read_param(checkreadok,"param_local->p_collect_hiv_test_results_cd4_over200");

        checkreadok = fscanf(param_file,"%lg",&(param_local->p_collect_hiv_test_results_cd4_under200));
        check_if_cannot_read_param(checkreadok,"param_local->p_collect_hiv_test_results_cd4_under200");
        //printf("param_local->p_collect_hiv_test_results_cd4_over200 = %lf\n",param_local->p_collect_hiv_test_results_cd4_over200);

        checkreadok = fscanf(param_file,"%lg",&(param_local->p_collect_cd4_test_results_cd4_nonpopart));
        check_if_cannot_read_param(checkreadok,"param_local->p_collect_cd4_test_results_cd4_nonpopart");

        checkreadok = fscanf(param_file,"%lg",&(param_local->p_collect_cd4_test_results_cd4_popartYEAR1));
        check_if_cannot_read_param(checkreadok,"param_local->p_collect_cd4_test_results_cd4_popartYEAR1");

        checkreadok = fscanf(param_file,"%lg",&(param_local->p_collect_cd4_test_results_cd4_popartYEAR2onwards));
        check_if_cannot_read_param(checkreadok,"param_local->p_collect_cd4_test_results_cd4_popartYEAR2onwards");
        //printf("param_local->p_collect_cd4_test_results_cd4_popartYEAR2onwards = %lf\n",param_local->p_collect_cd4_test_results_cd4_popartYEAR2onwards);

        for (icd4=0; icd4<NCD4; icd4++){
            checkreadok = fscanf(param_file,"%lg",&(param_local->p_dies_earlyart_cd4[icd4]));
            check_if_cannot_read_param(checkreadok,"param_local->p_dies_earlyart_cd4");
        }

        checkreadok = fscanf(param_file,"%lg",&(param_local->p_leaves_earlyart_cd4_over200_if_not_die_early));
        check_if_cannot_read_param(checkreadok,"param_local->p_leaves_earlyart_cd4_over200_if_not_die_early");

        checkreadok = fscanf(param_file,"%lg",&(param_local->p_leaves_earlyart_cd4_under200_if_not_die_early));
        check_if_cannot_read_param(checkreadok,"param_local->p_leaves_earlyart_cd4_under200_if_not_die_early");

        checkreadok = fscanf(param_file,"%lg",&(param_local->p_becomes_vs_after_earlyart_if_not_die_early_or_leave));
        check_if_cannot_read_param(checkreadok,"param_local->p_becomes_vs_after_earlyart_if_not_die_early_or_leave");

        checkreadok = fscanf(param_file,"%lg",&(param_local->p_stays_virally_suppressed));
        check_if_cannot_read_param(checkreadok,"param_local->p_stays_virally_suppressed");

        checkreadok = fscanf(param_file,"%lg",&(param_local->p_stays_virally_suppressed_male));
        check_if_cannot_read_param(checkreadok,"param_local->p_stays_virally_suppressed_male");

        checkreadok = fscanf(param_file,"%lg",&(param_local->p_stops_virally_suppressed));
        check_if_cannot_read_param(checkreadok,"param_local->p_stops_virally_suppressed");

        checkreadok = fscanf(param_file,"%lg",&(param_local->p_vu_becomes_virally_suppressed));
        check_if_cannot_read_param(checkreadok,"param_local->p_vu_becomes_virally_suppressed");
        /* Input times for the cascade events: */
        checkreadok = fscanf(param_file,"%lg",&(param_local->t_earlyart_dropout_min[NOTPOPART]));
        check_if_cannot_read_param(checkreadok,"param_local->t_earlyart_dropout_min");

        checkreadok = fscanf(param_file,"%lg",&max_temp);
        check_if_cannot_read_param(checkreadok,"param_local->t_earlyart_dropout_range");
        param_local->t_earlyart_dropout_range[NOTPOPART] = max_temp - param_local->t_earlyart_dropout_min[NOTPOPART];

        checkreadok = fscanf(param_file,"%lg",&(param_local->t_earlyart_dropout_min[POPART]));
        check_if_cannot_read_param(checkreadok,"param_local->t_earlyart_dropout_min");

        checkreadok = fscanf(param_file,"%lg",&max_temp);
        check_if_cannot_read_param(checkreadok,"param_local->t_earlyart_dropout_range");
        param_local->t_earlyart_dropout_range[POPART]    = max_temp - param_local->t_earlyart_dropout_min[POPART];

        checkreadok = fscanf(param_file,"%lg",&(param_local->t_dies_earlyart_min[NOTPOPART]));
        check_if_cannot_read_param(checkreadok,"param_local->t_dies_earlyart_min");

        checkreadok = fscanf(param_file,"%lg",&max_temp);
        check_if_cannot_read_param(checkreadok,"param_local->t_dies_earlyart_range");
        param_local->t_dies_earlyart_range[NOTPOPART] = max_temp - param_local->t_dies_earlyart_min[NOTPOPART];

        checkreadok = fscanf(param_file,"%lg",&(param_local->t_dies_earlyart_min[POPART]));
        check_if_cannot_read_param(checkreadok,"param_local->t_dies_earlyart_min");

        checkreadok = fscanf(param_file,"%lg",&max_temp);
        check_if_cannot_read_param(checkreadok,"param_local->t_dies_earlyart_range");
        param_local->t_dies_earlyart_range[POPART]    = max_temp - param_local->t_dies_earlyart_min[POPART];

        checkreadok = fscanf(param_file,"%lg",&(param_local->t_end_early_art));
        check_if_cannot_read_param(checkreadok,"param_local->t_end_early_art");

        checkreadok = fscanf(param_file,"%lg",&(param_local->t_cd4_retest_min[NOTPOPART]));
        check_if_cannot_read_param(checkreadok,"param_local->t_cd4_retest_min");

        checkreadok = fscanf(param_file,"%lg",&max_temp);
        check_if_cannot_read_param(checkreadok,"param_local->t_cd4_retest_range");
        param_local->t_cd4_retest_range[NOTPOPART] = max_temp - param_local->t_cd4_retest_min[NOTPOPART];

        checkreadok = fscanf(param_file,"%lg",&(param_local->t_cd4_retest_min[POPART]));
        check_if_cannot_read_param(checkreadok,"param_local->t_cd4_retest_min");

        checkreadok = fscanf(param_file,"%lg",&max_temp);
        check_if_cannot_read_param(checkreadok,"param_local->t_cd4_retest_range");
        param_local->t_cd4_retest_range[POPART]    = max_temp - param_local->t_cd4_retest_min[POPART]; 

        checkreadok = fscanf(param_file,"%lg",&(param_local->t_cd4_whenartfirstavail_min));
        check_if_cannot_read_param(checkreadok,"param_local->t_cd4_whenartfirstavail_min");

        checkreadok = fscanf(param_file,"%lg",&max_temp);
        check_if_cannot_read_param(checkreadok,"param_local->t_cd4_whenartfirstavail_range");
        param_local->t_cd4_whenartfirstavail_range = max_temp - param_local->t_cd4_whenartfirstavail_min;

        checkreadok = fscanf(param_file,"%lg",&(param_local->t_delay_hivtest_to_cd4test_min[NOTPOPART]));
        check_if_cannot_read_param(checkreadok,"param_local->t_delay_hivtest_to_cd4test_min");

        checkreadok = fscanf(param_file,"%lg",&max_temp);
        check_if_cannot_read_param(checkreadok,"param_local->t_delay_hivtest_to_cd4test_range");
        param_local->t_delay_hivtest_to_cd4test_range[NOTPOPART] = max_temp - param_local->t_delay_hivtest_to_cd4test_min[NOTPOPART];

        checkreadok = fscanf(param_file,"%lg",&(param_local->t_delay_hivtest_to_cd4test_min[POPART]));
        check_if_cannot_read_param(checkreadok,"param_local->t_delay_hivtest_to_cd4test_min");

        checkreadok = fscanf(param_file,"%lg",&max_temp);
        check_if_cannot_read_param(checkreadok,"param_local->t_delay_hivtest_to_cd4test_range");
        param_local->t_delay_hivtest_to_cd4test_range[POPART]    = max_temp - param_local->t_delay_hivtest_to_cd4test_min[POPART];

        // Uniform version
        //checkreadok = fscanf(param_file,"%lg",&(param_local->t_start_art_min[NOTPOPART]));
        //check_if_cannot_read_param(checkreadok,"param_local->t_start_art_min");
        //checkreadok = fscanf(param_file,"%lg",&max_temp);
        //check_if_cannot_read_param(checkreadok,"param_local->t_start_art_range");
        //param_local->t_start_art_range[NOTPOPART] = max_temp - param_local->t_start_art_min[NOTPOPART];
        //checkreadok = fscanf(param_file,"%lg",&(param_local->t_start_art_min[POPART]));
        //check_if_cannot_read_param(checkreadok,"param_local->t_start_art_min");
        //checkreadok = fscanf(param_file,"%lg",&max_temp);
        //check_if_cannot_read_param(checkreadok,"param_local->t_start_art_range");
        //param_local->t_start_art_range[POPART] = max_temp - param_local->t_start_art_min[POPART];

        // Exponential version
        checkreadok = fscanf(param_file,"%lg",&(param_local->t_start_art_mean_non_popart));
        check_if_cannot_read_param(checkreadok,"param_local->t_start_art_mean_non_popart");

        double temp;

        for(iround = 0 ; iround<NCHIPSROUNDS ; iround++)
        {
            /*print_here(iround);*/
            checkreadok = fscanf(param_file,"%lg",&temp);
            param_local->n_time_periods_art_popart_per_round[iround] = (int) temp;
            check_if_cannot_read_param(checkreadok,"param_local->n_time_periods_art_popart_per_round");
            /*printf("param_local->n_time_periods_art_popart_per_round[%d] is %d\n",iround,param_local->n_time_periods_art_popart_per_round[iround]);
            fflush(stdout);*/
        }
        
        for(iround = 0 ; iround < NCHIPSROUNDS; iround++){
            for(i = 0; i < param_local->n_time_periods_art_popart_per_round[iround]; i++){
                checkreadok = fscanf(param_file, "%lg", 
                    &(param_local->t_start_art_mean_fast_popart[iround][i]));
                
                check_if_cannot_read_param(checkreadok,
                    "param_local->t_start_art_mean_fast_popart_round1");
            }
            
            for(i = 0; i < param_local->n_time_periods_art_popart_per_round[iround] ; i++){
                checkreadok = fscanf(param_file, "%lg", 
                    &(param_local->t_start_art_mean_slow_popart[iround][i]));
                check_if_cannot_read_param(checkreadok,
                    "param_local->t_start_art_mean_slow_popart_round1");
            }
            
            for(i = 0; i < param_local->n_time_periods_art_popart_per_round[iround]; i++){
                checkreadok = fscanf(param_file, "%lg",
                    &(param_local->p_start_art_mean_fast_popart[iround][i]));
                check_if_cannot_read_param(checkreadok, 
                    "param_local->p_start_art_mean_fast_popart_round1");
            }
        }
        
        checkreadok = fscanf(param_file,"%lg",&(param_local->t_end_vs_becomevu_min[NOTPOPART]));
        check_if_cannot_read_param(checkreadok,"param_local->t_end_vs_becomevu_min");

        checkreadok = fscanf(param_file,"%lg",&max_temp);
        check_if_cannot_read_param(checkreadok,"param_local->t_end_vs_becomevu_range");
        param_local->t_end_vs_becomevu_range[NOTPOPART] = max_temp - param_local->t_end_vs_becomevu_min[NOTPOPART];

        checkreadok = fscanf(param_file,"%lg",&(param_local->t_end_vs_becomevu_min[POPART]));
        check_if_cannot_read_param(checkreadok,"param_local->t_end_vs_becomevu_min");

        checkreadok = fscanf(param_file,"%lg",&max_temp);
        check_if_cannot_read_param(checkreadok,"param_local->t_end_vs_becomevu_range");
        param_local->t_end_vs_becomevu_range[POPART] = max_temp - param_local->t_end_vs_becomevu_min[POPART];

        checkreadok = fscanf(param_file,"%lg",&(param_local->t_end_vs_dropout_min[NOTPOPART]));
        check_if_cannot_read_param(checkreadok,"param_local->t_end_vs_dropout_min");

        checkreadok = fscanf(param_file,"%lg",&max_temp);
        check_if_cannot_read_param(checkreadok,"param_local->t_end_vs_dropout_range");
        param_local->t_end_vs_dropout_range[NOTPOPART] = max_temp - param_local->t_end_vs_dropout_min[NOTPOPART];

        checkreadok = fscanf(param_file,"%lg",&(param_local->t_end_vs_dropout_min[POPART]));
        check_if_cannot_read_param(checkreadok,"param_local->t_end_vs_dropout_min");

        checkreadok = fscanf(param_file,"%lg",&max_temp);
        check_if_cannot_read_param(checkreadok,"param_local->t_end_vs_dropout_range");
        param_local->t_end_vs_dropout_range[POPART] = max_temp -  param_local->t_end_vs_dropout_min[POPART];

        checkreadok = fscanf(param_file,"%lg",&(param_local->t_end_vu_becomevs_min[NOTPOPART]));
        check_if_cannot_read_param(checkreadok,"param_local->t_end_vu_becomevs_min");

        checkreadok = fscanf(param_file,"%lg",&max_temp);
        check_if_cannot_read_param(checkreadok,"param_local->t_end_vu_becomevs_range");
        param_local->t_end_vu_becomevs_range[NOTPOPART] = max_temp - param_local->t_end_vu_becomevs_min[NOTPOPART];

        checkreadok = fscanf(param_file,"%lg",&(param_local->t_end_vu_becomevs_min[POPART]));
        check_if_cannot_read_param(checkreadok,"param_local->t_end_vu_becomevs_min");

        checkreadok = fscanf(param_file,"%lg",&max_temp);
        check_if_cannot_read_param(checkreadok,"param_local->t_end_vu_becomevs_range");
        param_local->t_end_vu_becomevs_range[POPART] = max_temp - param_local->t_end_vu_becomevs_min[POPART];

        checkreadok = fscanf(param_file,"%lg",&(param_local->t_end_vu_dropout_min[NOTPOPART]));
        check_if_cannot_read_param(checkreadok,"param_local->t_end_vu_dropout_min");

        checkreadok = fscanf(param_file,"%lg",&max_temp);
        check_if_cannot_read_param(checkreadok,"param_local->t_end_vu_dropout_range");
        param_local->t_end_vu_dropout_range[NOTPOPART] = max_temp - param_local->t_end_vu_dropout_min[NOTPOPART];  

        checkreadok = fscanf(param_file,"%lg",&(param_local->t_end_vu_dropout_min[POPART]));
        check_if_cannot_read_param(checkreadok,"param_local->t_end_vu_dropout_min");

        checkreadok = fscanf(param_file,"%lg",&max_temp);
        check_if_cannot_read_param(checkreadok,"param_local->t_end_vu_dropout_range");
        param_local->t_end_vu_dropout_range[POPART] = max_temp - param_local->t_end_vu_dropout_min[POPART];  

        /* This is the probability of getting back into the cascade during PopART for someone who 
         * dropped out before PopART. */
        //for (i=0; i<NCHIPSROUNDS; i++){
        for (i=0; i<1; i++){ /* Change after DSMB. */
            checkreadok = fscanf(param_file,"%lg",&(param_local->p_popart_to_cascade[i]));
            check_if_cannot_read_param(checkreadok,"param_local->p_popart_to_cascade[i]");
        }

        checkreadok = fscanf(param_file,"%lg",&(param_local->p_circ_nonpopart));
        check_if_cannot_read_param(checkreadok,"param_local->p_circ_nonpopart");

        //for (i=0; i<NCHIPSROUNDS; i++){
        for (i=0; i<1; i++){ /* Change after DSMB. */
            checkreadok = fscanf(param_file,"%lg",&(param_local->p_circ_popart[i]));
            check_if_cannot_read_param(checkreadok,"param_local->p_circ_popart[i]");
        }

        checkreadok = fscanf(param_file,"%lg",&(param_local->t_get_vmmc_min[NOTPOPART]));
        check_if_cannot_read_param(checkreadok,"param_local->t_get_vmmc_min");

        checkreadok = fscanf(param_file,"%lg",&max_temp);
        check_if_cannot_read_param(checkreadok,"param_local->t_get_vmmc_range");
        param_local->t_get_vmmc_range[NOTPOPART] = max_temp - param_local->t_get_vmmc_min[NOTPOPART];

        checkreadok = fscanf(param_file,"%lg",&(param_local->t_get_vmmc_min[POPART]));
        check_if_cannot_read_param(checkreadok,"param_local->t_get_vmmc_min");

        checkreadok = fscanf(param_file,"%lg",&max_temp);
        check_if_cannot_read_param(checkreadok,"param_local->t_get_vmmc_range");
        param_local->t_get_vmmc_range[POPART] = max_temp - param_local->t_get_vmmc_min[POPART];

        checkreadok = fscanf(param_file,"%lg",&(param_local->t_vmmc_healing));
        check_if_cannot_read_param(checkreadok,"param_local->t_vmmc_healing");

    }
    // Close parameter file
    fclose(param_file);
    return;
}


/**************************************************************************//**
 * @brief This function loops through CHiPs rounds and reads CHiPs uptake data
 * 
 * @details
 * This function populates the following arrays in the @ref parameters structure:\n
 *     `parameters.chips_params.prop_tested_by_chips_in_round` \n
 *     `parameters->chips_params->n_timesteps_per_round[NCHIPSROUNDS]`
 * 
 * And similarly the post-trial round equivalents of these arrays:\n
 *     `parameters->chips_params->prop_tested_by_chips_in_round_posttrial[][][]`\n
 *     `parameters->chips_params->n_timesteps_per_round_posttrial[NCHIPSROUNDS]`
 * 
 * The values saved within the `n_timesteps_per_round_posttrial` arrays are imported in the 
 * @ref read_time_params() function.
 * 
 * @param patch_tag Pointer to a char of the patch tag
 * @param allrunparameters pointer to an array of @ref parameters structs
 *****************************************************************************/

void read_chips_uptake_params(char *patch_tag, parameters *allrunparameters){
    int i_run, n_steps;
    double temp_time_continuous;
    int temp_time_discrete;
    int g, ac, i;  /* Indices for gender, age group and the timestep (within a chips round) */
    int chips_round; /* Index for the chips round. */
    
    /* This is a local temp variable we use so we don't have to keep writing 
    allparameters + i_run (or equivalently &allparameters[i_run]). */
    parameters *param_local;
    FILE *param_file;
    int checkreadok;
    // Add path before file name
    char param_file_name[LONGSTRINGLENGTH];
    char temp_tag[100];

    // Read in CHiPs uptake data
    for(chips_round = 0; chips_round <= NCHIPSROUNDS; chips_round++){
        
        // Add path before file name.
        strncpy(param_file_name,patch_tag,LONGSTRINGLENGTH);
        
        if(chips_round < NCHIPSROUNDS){
            // Chips round files are labelled 1,2,3,4
            sprintf(temp_tag, "chipsuptake_round%i.csv", chips_round + 1);
        }else{
            // Post-trial uptake.
            sprintf(temp_tag,"chipsuptake_roundposttrial.csv");
        }
        strcat(param_file_name, temp_tag);
        
        // Open parameter file
        if((param_file = fopen(param_file_name, "r")) == NULL){
            printf("Cannot open %s",param_file_name);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }else{
            if(VERBOSE_OUTPUT == 1)
                printf("CHiPs uptake parameters read from: %s:\n", param_file_name);
        }
        // Read parameters

        /* First line offile is a header, read that line and throw it away.*/
        fscanf(param_file, "%*[^\n]\n");

        /* Store outputs in the i_run = 0 parameter struct, and later copy it to every other run. */
        i_run = 0;
        /* No need for the i_run to be added - but just keep here for clarity. */
        param_local = allrunparameters + i_run;
        
        /* Set the cumulative counters to be zero. */
        for(g = 0; g < N_GENDER; g++){
            for(ac = 0; ac < (MAX_AGE - AGE_CHIPS + 1); ac++){
                if(chips_round < NCHIPSROUNDS){
                    param_local->chips_params->prop_tested_by_chips_in_round[g][ac][chips_round]=0;
                    n_steps = param_local->chips_params->n_timesteps_per_round[chips_round];
                }else{
                    param_local->chips_params->prop_tested_by_chips_in_round_posttrial[g][ac] = 0;
                    n_steps = param_local->chips_params->n_timesteps_per_round_posttrial;
                }
            }
        }
        
        /* Read in the data from the CHiPs uptake file*/
        for(i = 0; i < n_steps; i++){
            checkreadok = fscanf(param_file,"%lg", &temp_time_continuous);
            check_if_cannot_read_param(checkreadok, "atemp_time_continuous");
            
            checkreadok = fscanf(param_file,"%i", &temp_time_discrete);
            check_if_cannot_read_param(checkreadok, "temp_time_discrete");
            
            checkreadok = fscanf(param_file,"%i", &temp_time_discrete);
            check_if_cannot_read_param(checkreadok, "temp_time_discrete");
            
            for(g = 0; g < N_GENDER; g++){
                for(ac = 0; ac < (MAX_AGE - AGE_CHIPS + 1); ac++){
                    if(chips_round < NCHIPSROUNDS){
                        checkreadok = fscanf(param_file, "%lg", 
                            &(param_local->chips_params->prop_tested_by_chips_per_timestep[g][ac][i][chips_round]));
                        check_if_cannot_read_param(checkreadok,
                            "params->chips_params->prop_tested_...step[g][ac][i][chips_round]");

                        /* Update cumulative counters: */
                        param_local->chips_params->prop_tested_by_chips_in_round[g][ac][chips_round] += param_local->chips_params->prop_tested_by_chips_per_timestep[g][ac][i][chips_round];
                    }else{
                        checkreadok = fscanf(param_file, "%lg", &(param_local->chips_params->prop_tested_by_chips_per_timestep_posttrial[g][ac][i]));
                        check_if_cannot_read_param(checkreadok,
                            "params->chips_params->prop_tested_..._posttrial[g][ac][i]");

                        /* Update cumulative counters: */
                        param_local->chips_params->prop_tested_by_chips_in_round_posttrial[g][ac] += param_local->chips_params->prop_tested_by_chips_per_timestep_posttrial[g][ac][i];
                    }
                }
            }
        }
        /* Close chips uptake parameter file */
        fclose(param_file);
    }
}


/**************************************************************************//**
 * @brief Copy CHiPs parameters from one run to the next.
 * 
 * @param allrunparameters Pointer to an array of @ref parameters structure for storing all parameter values
 * @param n_runs Number of simulation runs
 *****************************************************************************/

void copy_chips_params( parameters **allrunparameters, int n_runs){
    int i, g, ac, i_run, chips_round, p;
    /* We assume that the CHiPs parameters are fixed from one run to the next: */
    for (i_run = 1; i_run<n_runs; i_run++){
        for (p=0; p<NPATCHES; p++){
            for (g=0; g<N_GENDER;g++){
                for (ac = 0; ac<(MAX_AGE-AGE_CHIPS+1); ac++){
                    for (chips_round=0;chips_round<NCHIPSROUNDS;chips_round++){
                        allrunparameters[p][i_run].chips_params->prop_tested_by_chips_in_round[g][ac][chips_round] = allrunparameters[p][0].chips_params->prop_tested_by_chips_in_round[g][ac][chips_round];
                        for(i=0; i<allrunparameters[p][0].chips_params->n_timesteps_per_round[chips_round];i++){
                            allrunparameters[p][i_run].chips_params->prop_tested_by_chips_per_timestep[g][ac][i][chips_round] = allrunparameters[p][0].chips_params->prop_tested_by_chips_per_timestep[g][ac][i][chips_round];
                        }
                    }
                    /* Now copy posttrial inputs. */
                    allrunparameters[p][i_run].chips_params->prop_tested_by_chips_in_round_posttrial[g][ac] = allrunparameters[p][0].chips_params->prop_tested_by_chips_in_round_posttrial[g][ac];
                    for(i=0; i<allrunparameters[p][0].chips_params->n_timesteps_per_round_posttrial;i++){
                        allrunparameters[p][i_run].chips_params->prop_tested_by_chips_per_timestep_posttrial[g][ac][i] = allrunparameters[p][0].chips_params->prop_tested_by_chips_per_timestep_posttrial[g][ac][i];
                    }
                }
            }
        }
    }
}


/**************************************************************************//**
 * @brief Read in the PC enrolment parameters
 * 
 * @param patch_tag Pointer to a char of the patch tag
 * @param community_id Community number
 * @param allrunparameters Pointer to an array of @ref parameters structure 
 *  for storing all parameter values
 * @param n_runs Number of simulation runs
 * @param p Patch number
 *****************************************************************************/

void read_pc0_enrolment_params(char *patch_tag, int community_id, parameters *allrunparameters, int n_runs, int p){
    int i_run;
    double temp_time_continuous;
    int temp_timestep_discrete, temp_year_discrete;
    int g, ap, i_pc_category, i;    /* Indices for gender, age group, HIV subgroup (HIV-, HIV+ aware, HIV+ unaware) and the timestep within a given PC round. */
    int pc_round; /* Index for the PC round. In this file we input PC0 numbers only *This will be true until the end of the trial*. */
    /* This is a local temp variable we use so we don't have to keep writing allparameters+i_run (or equivalently &allparameters[i_run]). */
    parameters *param_local;
    FILE *param_file;
    int checkreadok;
    // Add path before file name
    char param_file_name[LONGSTRINGLENGTH];
    char temp_tag[100];
    char buffer[100];

    int pc_enrolment_round = 0; /* We are reading in PC0 enrollee data only, not PC12N/24N. */

    // Read in PC enrolment data
    i_run=0; /* We store the outputs in the i_run=0 parameter struct, and later copy it to every other run. */
    param_local = allrunparameters + i_run; /* No need for the i_run to be added - but just keep here for clarity. */

    for (pc_round=0;pc_round<NPC_ROUNDS;pc_round++){
        strncpy(param_file_name,patch_tag,LONGSTRINGLENGTH);   // Adding path before file name.
        sprintf(temp_tag,"PC%i_community%i.csv",pc_round,community_id);   // NOTE - we use pc_round+1 so the round files are labelled 1,2,3,4.
        strcat(param_file_name, temp_tag);

        // Open parameter file
        if ((param_file=fopen(param_file_name,"r"))==NULL){    // Opening parameter file - print error if file not found.
            printf("Cannot open %s",param_file_name);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }else{
            if(VERBOSE_OUTPUT==1)
                printf("PC enrolment parameters read from: %s:\n",param_file_name);
        }
        // Read parameters

        /* The first line of this file is a header - the text below reads in just that line and throws it away.
         * The * instructs fscanf (all of the scanf family, in fact) to parse the data out as presented in the format string,
         * but NOT to store it at any target address (which is good, because there is none provided in an argument list).
         * the [^\n] means take anything except a newline, so ALL data will be consumed up to (but not including) the newline.
         *  Finally, the final \n means "and consume (and ignore) the newline" (which we just stopped at when fulfilling the
         *  prior format spec). Btw, the NULL isn't required in the argument list. */
        //fscanf(param_file, "%*[^\n]\n", NULL);
        fscanf(param_file, "%*[^\n]\n");

        /* Now read in the data from the given file: */
        i = 0;
        while(!feof(param_file)){

            /* We don't need these numbers - they are to make the csv files more readable. */
            checkreadok = fscanf(param_file,"%lg",&temp_time_continuous);
            if(feof(param_file)){
                /* Previous timestep was the last one. */
                param_local->PC_params->PC_END_TIMESTEP[pc_round] = temp_timestep_discrete;
                param_local->PC_params->PC_END_YEAR[pc_round] =  temp_year_discrete;
                break;
            }
            check_if_cannot_read_param(checkreadok,"btemp_time_continuous");
            checkreadok = fscanf(param_file,"%d",&temp_year_discrete);
            check_if_cannot_read_param(checkreadok,"temp_year_discrete");
            checkreadok = fscanf(param_file,"%d",&temp_timestep_discrete);
            check_if_cannot_read_param(checkreadok,"temp_timestep_discrete");

            /* Now start reading data.
             * First row after the header row gives the start time of the first round (note this is PC0 hence the '0' index). */
            if(i==0){
                param_local->PC_params->PC_START_TIMESTEP[pc_round] = temp_timestep_discrete;
                param_local->PC_params->PC_START_YEAR[pc_round] =  temp_year_discrete;
            }

            for (i_pc_category=0; i_pc_category<N_PC_HIV_STRATA; i_pc_category++){
                for (g=0; g<N_GENDER;g++){
                    for (ap = 0; ap<(AGE_PC_MAX-AGE_PC_MIN+1); ap++){
                        checkreadok = fscanf(param_file,"%d",&(param_local->PC_params->number_seen_by_PC_per_timestep[g][ap][i_pc_category][i][pc_round]));
                        if (checkreadok<1){
                            printf("ERROR: Format of PC enrolment file %s is not correct. Exiting\n",param_file_name);
                            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                            fflush(stdout);
                            exit(1);
                        }
                    }
                }
            }
            i = i+1;
        }

        param_local->PC_params->n_timesteps_per_round[pc_round] = (param_local->PC_params->PC_END_YEAR[pc_round]-param_local->PC_params->PC_START_YEAR[pc_round])*N_TIME_STEP_PER_YEAR + param_local->PC_params->PC_END_TIMESTEP[pc_round] - param_local->PC_params->PC_START_TIMESTEP[pc_round]+1;
        
        // Calculate PC midpoints from input data
        int midpoint_timestep = (param_local->PC_params->PC_START_TIMESTEP[pc_round] +
            param_local->PC_params->n_timesteps_per_round[pc_round]/2); 
        
        param_local->PC_params->PC_MIDPOINT_TIMESTEP[pc_round] = 
            midpoint_timestep % N_TIME_STEP_PER_YEAR;
        
        param_local->PC_params->PC_MIDPOINT_YEAR[pc_round] =
            param_local->PC_params->PC_START_YEAR[pc_round] + 
            midpoint_timestep / N_TIME_STEP_PER_YEAR;
        
        /* For debugging - check that total number of people seen agrees with the csv files, and check right duration of round. */
        if (param_local->PC_params->n_timesteps_per_round[pc_round]!=i){
            printf("ERROR: Mismatch in number of timesteps in PC round %i = %i %i. Exiting\n",pc_round,param_local->PC_params->n_timesteps_per_round[pc_round],i);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        /* Now check that there was nothing else in the file - ie there's no mismatch between the number of timesteps in the PC round we calculate and the number of lines of data we have. */
        int firsterror = 1;
        while (!feof(param_file)){
            /* Only print error statement once. */
            if (firsterror==1){
                printf("Error -extra data in %s\n",param_file_name);
                firsterror=0;
            }
            fgets(buffer, sizeof buffer, param_file);
            fputs(buffer, stdout);
        }
        if (firsterror==0){
            printf("Exiting\n");
            fclose(param_file);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        /* If no problems, then continue. */
        fclose(param_file);     // Closing PC uptake parameter file
        //printf("(pc_round=%i): Number of timesteps = %i . Start = %i %i. End = %i %i\n",pc_round,param_local->PC_params->n_timesteps_per_round[0], param_local->PC_params->PC_START_YEAR[0], param_local->PC_params->PC_START_TIMESTEP[0], param_local->PC_params->PC_END_YEAR[0], param_local->PC_params->PC_END_TIMESTEP[0]);
    }

    /* Set counters to zero: */
    for (g=0; g<N_GENDER;g++)
        for (ap = 0; ap<(AGE_PC_MAX-AGE_PC_MIN+1); ap++)
            for (i_pc_category=0; i_pc_category<N_PC_HIV_STRATA; i_pc_category++)
                param_local->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round] = 0;

    /* Cohort size is the same throughout (we deal with dropouts by ignoring some people in the cohort). */
    param_local->PC_params->cohort_size = 0;

    /* Now use the baseline PC0 to add up totals who should be seen in this cohort throughout the trial (if retention is 100%): */
    pc_round = 0;
    for (i_pc_category=0; i_pc_category<N_PC_HIV_STRATA; i_pc_category++){
        for (g=0; g<N_GENDER;g++){
            for (ap = 0; ap<(AGE_PC_MAX-AGE_PC_MIN+1); ap++){
                for (i=0; i<param_local->PC_params->n_timesteps_per_round[pc_round]; i++)
                    param_local->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round] += param_local->PC_params->number_seen_by_PC_per_timestep[g][ap][i_pc_category][i][pc_round];
                param_local->PC_params->cohort_size += param_local->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round];
            }
        }
    }
    /* Only need to worry about PC in patch 0: */
    if (p==0)
        printf("Cohort size = %i\n",param_local->PC_params->cohort_size);
    //printf("Number of timesteps = %i . Start = %i %i. End = %i %i\n",param_local->PC_params->n_timesteps_per_round[0], param_local->PC_params->PC_START_YEAR[0], param_local->PC_params->PC_START_TIMESTEP[0], param_local->PC_params->PC_END_YEAR[0], param_local->PC_params->PC_END_TIMESTEP[0]);

}


/**************************************************************************//**
 * @brief Copy PC parameters from one run to the next.
 * 
 * @param allrunparameters Pointer to an array of @ref parameters structure
 *  for storing all parameter values
 * @param n_runs Number of simulations
 ****************************************************************************/

void copy_pc_params( parameters **allrunparameters, int n_runs){
    int i, g, ap, i_run, i_pc_category, pc_round, p, pc_enrolment_round;
    /* We assume that the PC parameters are fixed from one run to the next: */
    for (i_run = 1; i_run<n_runs; i_run++){
        for (p=0; p<NPATCHES; p++){
            allrunparameters[p][i_run].PC_params->cohort_size = allrunparameters[p][0].PC_params->cohort_size;
            for (g=0; g<N_GENDER;g++){
                for (ap = 0; ap<(AGE_PC_MAX-AGE_PC_MIN+1); ap++){
                    for (i_pc_category=0; i_pc_category<N_PC_HIV_STRATA; i_pc_category++){
                        for (pc_enrolment_round=0;pc_enrolment_round<NPC_ENROLMENTS;pc_enrolment_round++){
                            allrunparameters[p][i_run].PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round] = allrunparameters[p][0].PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round];
                        }
                    }
                }
            }
            for (pc_round=0;pc_round<NPC_ROUNDS;pc_round++){
                allrunparameters[p][i_run].PC_params->PC_retention[pc_round] = allrunparameters[p][0].PC_params->PC_retention[pc_round];
                for (i_pc_category=0; i_pc_category<N_PC_HIV_STRATA; i_pc_category++){
                    for (g=0; g<N_GENDER;g++){
                        for (ap = 0; ap<(AGE_PC_MAX-AGE_PC_MIN+1); ap++){
                            for(i=0; i<allrunparameters[p][0].PC_params->n_timesteps_per_round[pc_round];i++){
                                allrunparameters[p][i_run].PC_params->number_seen_by_PC_per_timestep[g][ap][i_pc_category][i][pc_round] = allrunparameters[p][0].PC_params->number_seen_by_PC_per_timestep[g][ap][i_pc_category][i][pc_round];
                            }
                        }
                    }
                }
                allrunparameters[p][i_run].PC_params->n_timesteps_per_round[pc_round] = allrunparameters[p][0].PC_params->n_timesteps_per_round[pc_round];
                allrunparameters[p][i_run].PC_params->PC_START_TIMESTEP[pc_round] = allrunparameters[p][0].PC_params->PC_START_TIMESTEP[pc_round];
                allrunparameters[p][i_run].PC_params->PC_START_YEAR[pc_round] =  allrunparameters[p][0].PC_params->PC_START_YEAR[pc_round];
                allrunparameters[p][i_run].PC_params->PC_END_TIMESTEP[pc_round] = allrunparameters[p][0].PC_params->PC_END_TIMESTEP[pc_round];
                allrunparameters[p][i_run].PC_params->PC_END_YEAR[pc_round] =  allrunparameters[p][0].PC_params->PC_END_YEAR[pc_round];
                
                allrunparameters[p][i_run].PC_params->PC_MIDPOINT_TIMESTEP[pc_round] = allrunparameters[p][0].PC_params->PC_MIDPOINT_TIMESTEP[pc_round];
                allrunparameters[p][i_run].PC_params->PC_MIDPOINT_YEAR[pc_round] =  allrunparameters[p][0].PC_params->PC_MIDPOINT_YEAR[pc_round];
            }
        }
    }
}


/**************************************************************************//**
 * @brief Reads in dates of PC12, PC24, PC36 and expected retention rate
 * 
 * @param patch_tag Pointer to a char of the patch tag
 * @param allrunparameters Pointer to an array of @ref parameters structure for storing all parameter values
 * @param n_runs Number of simulation runs
 *****************************************************************************/

void read_pc_future_params(char *patch_tag, parameters *allrunparameters, int n_runs){
    FILE *param_file;
    int pc_round;   /* Index for PC round number. */
    //double PC_TIME; /* Temporary store of the (decimal) time when PC round starts/ends - used when we convert to discrete time. */

    int i_run;
    /* This is a local temp variable we use so we don't have to keep writing allparameters+i_run (or equivalently &allparameters[i_run]). */
    parameters *param_local;
    int checkreadok;
    // Add path before file name
    char param_file_name[LONGSTRINGLENGTH];

    strncpy(param_file_name,patch_tag,LONGSTRINGLENGTH);
    strcat(param_file_name, "PC.csv");


    // Open parameter file
    if ((param_file=fopen(param_file_name,"r"))==NULL){
        printf("Cannot open %s\n",param_file_name);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    else
        if(VERBOSE_OUTPUT==1)
            printf("PC parameters read from: %s:\n",param_file_name);

    /* The first line of this file is a header - the text below reads in just that line and throws it away.
     * The * instructs fscanf (all of the scanf family, in fact) to parse the data out as presented in the format string,
     * but NOT to store it at any target address (which is good, because there is none provided in an argument list).
     * the [^\n] means take anything except a newline, so ALL data will be consumed up to (but not including) the newline.
     *  Finally, the final \n means "and consume (and ignore) the newline" (which we just stopped at when fulfilling the
     *  prior format spec). */
    fscanf(param_file, "%*[^\n]\n");

    // Read parameters from each line i_run
    for (i_run = 0; i_run<n_runs; i_run++){
        param_local = allrunparameters + i_run;

        for(pc_round=0;pc_round<NPC_ROUNDS; pc_round++){
            checkreadok = fscanf(param_file,"%lg",&(param_local->PC_params->PC_retention[pc_round]));
            check_if_cannot_read_param(checkreadok,"param_local->PC_retention");
        }
    }
    fclose(param_file);
    return;
}


/**************************************************************************//**
 * @brief Read initial condition parameter values
 * 
 * @details Reads parameters related to initial conditions of the simulation
 * which are stored in the files named `param_processed_patch$p_init.csv`
 * @param patch_tag Pointer to name of the dir where the parameter file is stored
 * @param allrunparameters Pointer to an array of @ref parameters structure for storing all parameter values
 * @param n_runs Number of runs of the simulation
 *****************************************************************************/

void read_initial_params(char *patch_tag, parameters *allrunparameters, int n_runs){
    
    FILE *param_file;

    /* Indices for age group (0..N_AGE-1), risk and gender. */
    int ag, r, g;
    
    double checksum;
    double temp_int; /* Used when file input is a decimal but want to convert to an integer. */
    double temp_initial_propbyrisk[N_GENDER][N_RISK - 1];
    double temp_sum[N_GENDER];

    // A scaling factor (on log scale) which multiplies the size of the initial HIV seed 
    // (ie % of people who we seed to be HIV+).  
    double log_seed_multiplier; 
    double seed_multiplier; /* 10^log_seed_multiplier. */

    int i_run;
    
    // This is a local temp variable we use so we don't have to keep 
    // writing allparameters+i_run (or equivalently &allparameters[i_run]).
    parameters *param_local;
    int checkreadok;
    // Add path before file name
    char param_file_name[LONGSTRINGLENGTH];

    strncpy(param_file_name, patch_tag, LONGSTRINGLENGTH);
    strcat(param_file_name, "init.csv");

    // Open parameter file
    if((param_file = fopen(param_file_name, "r")) == NULL){
        printf("Cannot open %s\n", param_file_name);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }else{
        if(VERBOSE_OUTPUT == 1){
            printf("Init parameters read from: %s:\n", param_file_name);
        }
    }

    // Throw away the first line of the parameter file (the header)
    fscanf(param_file, "%*[^\n]\n");

    // Read parameters from each line i_run
    for (i_run = 0; i_run < n_runs; i_run++){
        param_local = allrunparameters + i_run;
        
        // Read the random seed
        checkreadok = fscanf(param_file,"%lf",&temp_int);
        param_local->rng_seed = (long) floor(temp_int);
        //checkreadok = fscanf(param_file,"%ld",&(param_local->rng_seed));
        //printf("Reading rng_seed = %li\n",param_local->rng_seed);
        check_if_cannot_read_param(checkreadok,"param_local->rng_seed");

        // Read the initial population size
        checkreadok = fscanf(param_file,"%lg",&(param_local->initial_population_size));
        check_if_cannot_read_param(checkreadok,"param_local->initial_population_size");

        // Read initial proportion of individuals in each age group
        for(ag = 0; ag < N_AGE; ag++){
            checkreadok = fscanf(param_file, "%lg", &(param_local->initial_prop_age[ag]));
            check_if_cannot_read_param(checkreadok, "param_local->initial_prop_age");
        }
        
        /* Load proportion of total popn who are of gender k and risk group r. */
        for(g = 0; g < N_GENDER; g++){
            temp_sum[g] = 0.0;
        }

        for(r = 0; r < N_RISK - 1; r++){
            for(g = 0; g < N_GENDER; g++){
                checkreadok = fscanf(param_file, "%lg", &(temp_initial_propbyrisk[g][r]));
                check_if_cannot_read_param(checkreadok, "param_local->initial_prop_gender_risk");
                
                if(r == 0){
                    param_local->initial_prop_gender_risk[g][r] = temp_initial_propbyrisk[g][r];
                }else{
                    param_local->initial_prop_gender_risk[g][r] =
                         temp_initial_propbyrisk[g][r] * 
                             (1.0 - param_local->initial_prop_gender_risk[g][r - 1]);
                }
                temp_sum[g] += param_local->initial_prop_gender_risk[g][r];
            }
        }
        if(N_RISK > 2){
            for(g = 0; g < N_GENDER; g++){
                param_local->initial_prop_gender_risk[g][N_RISK - 1] = 1.0 - temp_sum[g];
            }
        }

        /* Code if we want to fix male/female proportions to be the same: */
//      if (i_run==0) /* Only print out warning once. */
//          printf("WARNING: Setting female risk proportion to be same as male proportion by hand for now. May want to change this later. \n");
//      for (r=0; r<N_RISK; r++)
//          param_local->initial_prop_gender_risk[FEMALE][r] = param_local->initial_prop_gender_risk[MALE][r];


        //for (r=0; r<N_RISK; r++){
        //  for (g=0; g<N_GENDER; g++){
        //      printf("param_local->initial_prop_gender_risk[g=%i][r=%i] = %lf\n",g,r,param_local->initial_prop_gender_risk[g][r]);
        //  }
        //}

        // Load proportion of gender-risk specific population who are 
        // HIV+ at time of HIV introduction
        for(r = 0; r < N_RISK; r++){
            for(g = 0; g < N_GENDER; g++){
                checkreadok = fscanf(param_file, "%lg", 
                    &(param_local->initial_prop_infected_gender_risk[g][r]));
                    
                check_if_cannot_read_param(checkreadok,
                    "param_local->initial_prop_infected_gender_risk");
            }
        }

        checkreadok = fscanf(param_file, "%lg", &(log_seed_multiplier));
        check_if_cannot_read_param(checkreadok, "log_seed_multiplier");
        seed_multiplier = pow(10, log_seed_multiplier);
        
        /* Now scale up each % seeded: */
        for(r = 0; r < N_RISK; r++){
            for(g = 0; g < N_GENDER; g++){
                param_local->initial_prop_infected_gender_risk[g][r] =
                    param_local->initial_prop_infected_gender_risk[g][r] * seed_multiplier;
            }
        }
        
        checkreadok = fscanf(param_file, "%lf", &temp_int);
        param_local->n_years_HIV_seeding = (int) floor(temp_int);
        check_if_cannot_read_param(checkreadok, "param_local->n_years_HIV_seeding");
        
        ////// This is for debugging: 
        ////// Make sure proportions add to 1.
        ////// Shall we do it this way or have only 2 of the 3 read in a file and automatically calculating the third one? Could then print the third one to warn the user

        for(g = 0; g < N_GENDER; g++){
            checksum = 0.0;
            for(r = 0; r < N_RISK; r++){
                checksum += param_local->initial_prop_gender_risk[g][r];
            }
            if(fabs(checksum - 1.0) > 1e-12){
                printf("Sum of proportions of initial population by risk is not 1 %f.\n", checksum);
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
        }
        checksum = 0;
        for(ag = 0; ag < N_AGE; ag++){
            checksum += param_local->initial_prop_age[ag];
        }
        if(checksum != 1.0){
            printf("Sum of proportions of initial population by age is not 1 = %f.\n", checksum);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        /* End of debugging. */
    }
    fclose(param_file);
    return;
}


/**************************************************************************//**
 * @brief Read in the seed used by python to generate a Latin Hypercube Sample
 * 
 * @param file_directory Name of directory where parameter file 
 * `python_seed.txt` is stored.
 * 
 * @return The seed within the file `python_seed.txt`
 *****************************************************************************/

long get_python_seed(char *file_directory){
    long seed;
    FILE *PYTHON_SEED_FILE;
    PYTHON_SEED_FILE = fopen("", "r");
    char python_seed_filename[LONGSTRINGLENGTH];
    strncpy(python_seed_filename, file_directory, LONGSTRINGLENGTH);
    
    /* Adds a / or \ as needed if working in directory other than current local dir. */
    add_slash(python_seed_filename);
    strcat(python_seed_filename, "python_seed.txt");
    
    // Open parameter file
    if((PYTHON_SEED_FILE = fopen(python_seed_filename, "r")) == NULL){
        printf("Cannot open python_seed.txt");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    
    // Read the seed
    int checkreadok = fscanf(PYTHON_SEED_FILE, "%ld", &seed);
    
    // Close the parameter file
    fclose(PYTHON_SEED_FILE);
    
    if(checkreadok < 1){
        printf("Could not read python seed from python_seed.txt");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    return(seed);
}
