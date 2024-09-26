/**************************************************************************//**
 * @file pc.c
 * @brief Functions for simulating the PC sampling and PC visits
 * @details This code is not currently used in the model.  
 * PC = Population Cohort.
*****************************************************************************/

/* Contains functions relevant to PC sampling and visits.
 * 1) Create PC list, and PC reserve list (of 5 people  of each age) for PC0 (note that PC12N and PC24N can draw from reserve list apart from new age 18).
 * 2) Add PC variables to indiv structure -  round of entry to PC (-1 not in PC, 0 = PC0, 1=PC12N, 2=PC24N), date of visits by PC (0, 12, 24, 36, set to -1 if not visited in that round), date of exit from PC.
 * 3) Function to update PC lists - count number of deaths and then draw n_dropout people to drop out (by age and gender) (n_dropout = number lost in data - number died in model). If one of those people died, then draw another person.
 * 3b) Update reserve list for PC12N and PC24N - maybe easier to create a new reserve list?
 * 4) write input functions for number of people visited at each timestep in PC0
 * 5) Pull out data from PC0 on number of people visited at each timestep.
 */
#include "pc.h"
#include "structures.h"
#include "constants.h"


/**************************************************************************//**
 * @brief Determine how PC participants are stratified by HIV status
 * 
 * @param pc_participant Individual PC participant to be stratified
 
 * @returns Returns an integer indicating the stratification of the PC participant:\n
 * 0 if the individual is HIV-\n
 * 1 if the individual is HIV+ and aware\n
 * 2 if the individual is HIV+ but unaware
 *****************************************************************************/

int get_PC_HIV_stratum(individual *pc_participant){
    if (pc_participant->HIV_status==UNINFECTED)
        return 0;
    else if (
        (pc_participant->HIV_status>UNINFECTED) && 
        (pc_participant->ART_status>ARTNEG) && 
        (pc_participant->ART_status<ARTDEATH)
        )
        return 1;  /* If HIV+ and aware: */
    else
        return 2;   /* If HIV+ but unaware: */
}


/**************************************************************************//**
 * @brief Take the population of currently alive people (using age_list) 
 * and firstly sub-divides them into a 'chips_sampling_frame', 
 * e.g. dividing up men and women, as CHiPs tends to visit more women than men.
 * 
 * @param patch Pointer to the patch where the PC sampling is taking place
 * @param p Index of the patch where the PC sampling is taking place
 * @param pc_round Round of PC sampling.  
 *  The variable `pc_round` tells us what we are sampling: 
 * pc_round=0: PC0 entrant; pc_round=1: PC12N entrant; pc_round=2: PC24N entrant.
 * @param g Gender of group being visited
 * @param ap Age group of group being visited
 * @param i_pc_category HIV status of group being visited
 * @param original_size Original size of the group being visited
 * @param number_to_remove Number of people to remove from the group being visited
 *****************************************************************************/

void remove_extras_from_timestep_recruitment(patch_struct *patch, int p, int pc_enrolment_round, 
        int g, int ap, int i_pc_category, int original_size, int number_to_remove){
    int *TEMP_SAMPLE_FRAME_TO_REMOVE;
    int *TEMP_LIST_TO_REMOVE;

    TEMP_SAMPLE_FRAME_TO_REMOVE = malloc(original_size*sizeof(int));
    TEMP_LIST_TO_REMOVE = malloc(number_to_remove*sizeof(int));
    if (TEMP_SAMPLE_FRAME_TO_REMOVE==NULL || TEMP_LIST_TO_REMOVE==NULL){
        printf("Unable to allocate TEMP_SAMPLE_FRAME_TO_REMOVE or TEMP_LIST_TO_REMOVE in remove_extras_from_timestep_recruitment(). Execution aborted.");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    int n = 0;
    int i_dt,j;

    for (i_dt=0;i_dt<patch[p].param->PC_params->n_timesteps_per_round[pc_enrolment_round]; i_dt++){
        for (j=0; j<patch[p].param->PC_params->number_seen_by_PC_per_timestep[g][ap][i_pc_category][i_dt][pc_enrolment_round]; j++){
            TEMP_SAMPLE_FRAME_TO_REMOVE[n] = i_dt;
            n++;
        }
    }
    if (n!=original_size){
        printf("Error - can't match sample size in remove_extras_from_timestep_recruitment(). Exiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    gsl_ran_choose(rng, TEMP_LIST_TO_REMOVE, number_to_remove, TEMP_SAMPLE_FRAME_TO_REMOVE, original_size, sizeof (int));

    for (n=0; n<number_to_remove; n++){
        /* We remove one person who was supposed to be seen from timestep TEMP_LIST_TO_REMOVE[n]. */
        i_dt = TEMP_LIST_TO_REMOVE[n];
        patch[p].param->PC_params->number_seen_by_PC_per_timestep[g][ap][i_pc_category][i_dt][pc_enrolment_round]--;
    }

    for (i_dt=0;i_dt<patch[p].param->PC_params->n_timesteps_per_round[pc_enrolment_round]; i_dt++){
        if (patch[p].param->PC_params->number_seen_by_PC_per_timestep[g][ap][i_pc_category][i_dt][pc_enrolment_round]<0){
            printf("Error - have negative number of people to see per timestep in remove_extras_from_timestep_recruitment(). Exiting\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
    }
    free(TEMP_SAMPLE_FRAME_TO_REMOVE);
    free(TEMP_LIST_TO_REMOVE);
}


/**************************************************************************//**
 * @brief Draw up a list of people to be enrolled by PC (including reserves)
 * 
 * @param patch Pointer to the patch where the PC sampling is taking place
 * @param age_list Pointer to the age list
 * @param PC_sample Pointer to the PC sample object storing information
 * on the Population Cohort sample
 * @param param Pointer to the parameters object
 * @param pc_enrolment_round Round of PC sampling.
 * @param p Index of the patch where the PC sampling is taking place
 *****************************************************************************/

void create_popart_pc_sample(patch_struct *patch, age_list_struct *age_list, 
        PC_sample_struct *PC_sample, parameters *param, int pc_enrolment_round, int p){
    int g;
    int aa, ai,i, ap;
    int i_pc_category; /* Index splitting up the population by HIV status etc. */

    /* We are drawing our PC sample at the beginning of the round. Lots of things can happen in the mean time
     * (death, HIV infection, learning status) that make that person move outside the subpopulation - for example if someone was
     * uninfected when sampled but got HIV in the mean time they should no longer be in the 'uninfected' subpopulation.
     * So we need 'reserves' - people who are in the same population who can replace a person if needed. */
    int n_reserves;    /* Number of reserves to add for a given subpopulation. */
    double prop_reserves = 0.25; /* We want a minimum of 10% more people to account for deaths, movement to other subpopulations (e.g. HIV- gets infected). */
    int min_reserves = 15; /* Arbitrary minimum number of reserves. */

    /* For use with FOLLOW_INDIVIDUAL - we store these the first time we find them so we can find them easily next time: */
    int g_persontofollow = -1; /* Default value indicates that the FOLLOW_INDIVIDUAL did not turn up when going through - this is because they are too young to be visited by CHiPs. */
    int ap_persontofollow = -1;
    int i_pc_category_persontofollow = -1;

    individual *enrollee; /* Temporary pointer for the person we are enrolling at tthe time to make code more readable. As it points at an existing person, no need to malloc it. */
    int original_cohort_size = patch[p].param->PC_params->cohort_size;

    if (pc_enrolment_round<0 || pc_enrolment_round>2){
        printf("ERROR: pc_enrolment_round=%d takes range 0, 1 or 2 (PC0 entrant, PC12N entrant, PC24N entrant. Exiting\n",pc_enrolment_round);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    int max_number = 0;
    for (g=0;g<N_GENDER;g++){
        for (aa=(AGE_PC_MIN-AGE_ADULT); aa<=(AGE_PC_MAX-AGE_ADULT); aa++){
            ai = age_list->age_list_by_gender[g]->youngest_age_group_index + aa; /* ai is the index of the array age_list->number_per_age_group of the age group of people you want to be dead */
            while (ai>(MAX_AGE-AGE_ADULT-1))
                ai = ai - (MAX_AGE-AGE_ADULT);
            if (max_number<age_list->age_list_by_gender[g]->number_per_age_group[ai])
                max_number = age_list->age_list_by_gender[g]->number_per_age_group[ai];
        }
    }
    //printf("Maximum = %i\n",max_number);
    /* This acts as a temporary store for the PC sampling frame for a given subpopulation. We then draw the number of people
     * we want in the PC sample (+reserves) from this using gsl_ran_sample().
     * We assume that there are no more than 1/50th of the total population in a given age/gender/HIV status group.
     * Because of the way the sampling frame is constructed we need an array to store them. */
    //long TEMP_SAMPLE_FRAME[N_PC_HIV_STRATA][MAX_NUMBER_PC_PARTICIPANTS_PER_GROUP];
    long *TEMP_SAMPLE_FRAME[N_PC_HIV_STRATA];
    for (i_pc_category=0; i_pc_category<N_PC_HIV_STRATA; i_pc_category++){
        TEMP_SAMPLE_FRAME[i_pc_category] = malloc(max_number*sizeof(long));
        if (TEMP_SAMPLE_FRAME[i_pc_category]==NULL){
            printf("Unable to allocate TEMP_SAMPLE_FRAME in create_popart_pc_sample(). Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
    }
    int n_sample[N_PC_HIV_STRATA];
    int number_not_recruited; /* Variable used if we can't make the full cohort size. */

    for (g=0;g<N_GENDER;g++)
        for (ap=0; ap<(AGE_PC_MAX-AGE_PC_MIN+1); ap++)
            for (i_pc_category=0; i_pc_category<N_PC_HIV_STRATA; i_pc_category++)
                PC_sample->next_person_to_see[g][ap][i_pc_category] = 0;

    for (g=0;g<N_GENDER;g++){
        /* aa is age index for age_list by gender (not adjusting for youngest age group).
         * Here it is chosen to correspond to ages from AGE_PC_MIN (=18) to AGE_PC_MAX (=44) inclusive. */
        for (aa=(AGE_PC_MIN-AGE_ADULT); aa<=(AGE_PC_MAX-AGE_ADULT); aa++){
            ai = age_list->age_list_by_gender[g]->youngest_age_group_index + aa; /* ai is the index of the array age_list->number_per_age_group of the age group of people you want to be dead */
            while (ai>(MAX_AGE-AGE_ADULT-1))
                ai = ai - (MAX_AGE-AGE_ADULT);

            ap = aa-(AGE_PC_MIN-AGE_ADULT); /* This is the relationship between the index aa and ap (the index for the PC_sample_struct structures). */

            /* TEMP_SAMPLE_FRAME[][] is a sampling frame for a given age group and gender.
             * So set counter to zero, and list to 'blank' (-1) here. */
            for (i_pc_category=0; i_pc_category<N_PC_HIV_STRATA; i_pc_category++){
                n_sample[i_pc_category] = 0;
                for (i=0;i<MAX_NUMBER_PC_PARTICIPANTS_PER_GROUP;i++)
                    TEMP_SAMPLE_FRAME[i_pc_category][i] = -1;
            }


            for (i=0;i<age_list->age_list_by_gender[g]->number_per_age_group[ai];i++){
                enrollee = age_list->age_list_by_gender[g]->age_group[ai][i];
                if (enrollee->PC_cohort_index==-1){  /* Only consider people who aren't already in PC. */
                    /* For this person work out what PC HIV stratum they belong to: */
                    i_pc_category = get_PC_HIV_stratum(enrollee);

                    //if (n_sample[i_pc_category]<MAX_NUMBER_PC_PARTICIPANTS_PER_GROUP){
                    /* Add the person to the sampling frame: */
                    TEMP_SAMPLE_FRAME[i_pc_category][n_sample[i_pc_category]] = enrollee->id;
                    n_sample[i_pc_category]++;

                    /* For debugging. */
                    if(enrollee->id==FOLLOW_INDIVIDUAL && enrollee->patch_no==FOLLOW_PATCH){
                        printf("Possible PC enrolment %ld %d %d in round %d \n",enrollee->id,ai,i,pc_enrolment_round);
                        fflush(stdout);
                        /* Now store their characteristics so it's easier to find them: */
                        g_persontofollow = g;
                        ap_persontofollow = ap;
                        i_pc_category_persontofollow = i_pc_category;
                    }
                    /* Also for debugging. */
                    if(enrollee->cd4==DUMMYVALUE || enrollee->cd4==DEAD){
                        printf("Error -trying to schedule PC enrolment for dead/non-existent person %ld\n",enrollee->id);
                        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                        fflush(stdout);
                        exit(1);
                    }
                }
            }
            /* Now we have got the people of this age/gender to sample from, so we want to draw the people to visit. */
            for (i_pc_category=0; i_pc_category<N_PC_HIV_STRATA; i_pc_category++){
                /* First calculate the number of reserves needed: */
                if(param->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round]>0){
                    n_reserves = (int) round(prop_reserves*param->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round]);
                    if (n_reserves<min_reserves)
                        n_reserves = min_reserves;
                }else
                    n_reserves = 0; /* No need for reserves if we're not seeing anyone. */
                //printf("param->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round] = %i n_sample=%i\n",param->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round],n_sample[i_pc_category]);

                if (g==1 && ap==15 && i_pc_category==2){
                    printf("n_reserves = %i\n",n_reserves);
                    printf("param->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round] = %i n_sample=%i\n",param->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round],n_sample[i_pc_category]);
                }

                /* Make sure no error in counting number of people: */
                if (n_sample[i_pc_category]<param->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round]){
                    /* Right now  param->PC_params->number_enrolled_in_PC_round[] is not actually the number of people who are recruited as we cannot exceed the sample size.
                     So we will reduce param->PC_params->number_enrolled_in_PC_round[]: */

                    number_not_recruited = param->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round] - n_sample[i_pc_category];

                    remove_extras_from_timestep_recruitment(patch, p, pc_enrolment_round, g, ap, i_pc_category, param->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round], number_not_recruited);

                    param->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round] = n_sample[i_pc_category];

                    /* Also update the cohort size to reflect non-recruitment: */
                    patch[p].param->PC_params->cohort_size = patch[p].param->PC_params->cohort_size - number_not_recruited;
                }

                /* Make sure that the total sample size including reserves doesn't exceed the number of people available. */
                PC_sample->number_in_sample_including_reserves[g][ap][i_pc_category] = fmin(param->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round]+n_reserves,n_sample[i_pc_category]);
                if (g==1 && ap==15 && i_pc_category==2){
                    printf("PC_sample->number_in_sample_including_reserves[g][ap][i_pc_category] = %li\n",PC_sample->number_in_sample_including_reserves[g][ap][i_pc_category]);
                }

                /* Choose the people who will be in PC in this subpopulation: */
                if (PC_sample->number_in_sample_including_reserves[g][ap][i_pc_category]>0){
                    gsl_ran_choose(rng, PC_sample->list_ids_potential_enrollees[g][ap][i_pc_category], PC_sample->number_in_sample_including_reserves[g][ap][i_pc_category], TEMP_SAMPLE_FRAME[i_pc_category], n_sample[i_pc_category], sizeof (long));
                    /* Randomise the order (as gsl_ran_choose maintains the order of the original list). */
                    gsl_ran_shuffle(rng, PC_sample->list_ids_potential_enrollees[g][ap][i_pc_category], PC_sample->number_in_sample_including_reserves[g][ap][i_pc_category], sizeof (long));
                }
            }
        }
    }
    if (original_cohort_size>patch[p].param->PC_params->cohort_size){
        printf("Warning: IBM population in certain strata was smaller than PC sample needed in pc_enrolment_round=%i. Cohort size reduced from %i to %i.\n",pc_enrolment_round,original_cohort_size,patch[p].param->PC_params->cohort_size);
        fflush(stdout);
    }

    /* Check to see if this person was visited. Note that if g_persontofollow==-1 then they are too young to be in chips_sampling_frame->list_ids_to_visit. */
    if (p==FOLLOW_PATCH && g_persontofollow>-1){
        for(i=0;i<param->PC_params->number_enrolled_in_PC_round[g_persontofollow][ap_persontofollow][i_pc_category_persontofollow][pc_enrolment_round];i++){
            if(PC_sample->list_ids_potential_enrollees[g_persontofollow][ap_persontofollow][i_pc_category_persontofollow][i]==FOLLOW_INDIVIDUAL){
                printf("PC enrolment for adult %ld, gender %d PC cat=%d from patch %d now scheduled\n",PC_sample->list_ids_potential_enrollees[g_persontofollow][ap_persontofollow][i_pc_category_persontofollow][i],g_persontofollow,i_pc_category_persontofollow,p);
                fflush(stdout);
            }
        }
    }
    for (i_pc_category=0; i_pc_category<N_PC_HIV_STRATA; i_pc_category++)
        free(TEMP_SAMPLE_FRAME[i_pc_category]);
}


/**************************************************************************//**
 * @brief Carry out the PC enrolment visits for a given timestep at time t
 * 
 * @param t0 Current year
 * @param t_step Current timestep
 * @param patch Pointer to the patch where the PC enrolment is taking place
 * @param p Index of the patch where the PC enrolment is taking place
 * @param pc_enrolment_round Round of PC enrolment
*****************************************************************************/

void carry_out_PC_enrolment_per_timestep(int t0, int t_step, patch_struct *patch, int p, int pc_enrolment_round){
    int g;
    int ap;
    int j;
    int i_pc_category; /* Index splitting up the population by HIV status etc. */
    int i_pc_thisperson;

    //int number_deleted = 0; /* Count how many people we didn't get to put in PC.  For debugging. */
    //printf("Carrying out enrolment for round %i at time %i + %i\n",pc_enrolment_round,t0,t_step);
    /* i_dt tells us the index of the timestep within PC_sample corresponding to the current time t so we can work out how many people to visit.*/
    int i_dt;

    int n_dropped_from_cohort = 0; /* Keep track of the number of people we couldn't enroll. */

    individual *enrollee; /* Temporary pointer for the person we are enrolling at tthe time to make code more readable. As it points at an existing person, no need to malloc it. */

    //printf("Number of timesteps = %i in round %i. Start = %i %i. End = %i %i\n",patch[p].param->PC_params->n_timesteps_per_round[pc_enrolment_round],pc_enrolment_round, patch[p].param->PC_params->PC_START_YEAR[pc_enrolment_round], patch[p].param->PC_params->PC_START_TIMESTEP[pc_enrolment_round], patch[p].param->PC_params->PC_END_YEAR[pc_enrolment_round], patch[p].param->PC_params->PC_END_TIMESTEP[pc_enrolment_round]);

    i_dt = (t0 - patch[p].param->PC_params->PC_START_YEAR[pc_enrolment_round])*N_TIME_STEP_PER_YEAR + (t_step-patch[p].param->PC_params->PC_START_TIMESTEP[pc_enrolment_round]);

    /* For debugging: */
    if (i_dt>=MAX_N_TIMESTEPS_PER_PC_ROUND){
        printf("Problem - MAX_N_TIMESTEPS_PER_PC_ROUND is too small. We are %i timesteps into PC round %i. Exiting\n",i_dt,pc_enrolment_round);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if(i_dt<0||i_dt>=patch[p].param->PC_params->n_timesteps_per_round[pc_enrolment_round])
    {
        printf("ERROR i_dt=%i is outside range [0,%i] in PC round %i at time %i+%i*TIME_STEP. Exiting\n",i_dt,patch[p].param->PC_params->n_timesteps_per_round[pc_enrolment_round],pc_enrolment_round,t0,t_step);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    /* Now PC visits each sub-population (currently men/women by year age group ap and HIV status stratum) in turn:
     * Go through the list of id's of people to visit, and visit them: */

    for (g=0;g<N_GENDER;g++){
        /* Run from 18 to 44. */
        for (ap=0; ap<(AGE_PC_MAX-AGE_PC_MIN+1); ap++){
            for (i_pc_category=0; i_pc_category<N_PC_HIV_STRATA; i_pc_category++){
                j = 0; /* This keeps track of how many people we have successfully enrolled of this subpopulation.*/

                while (j<(patch[p].param->PC_params->number_seen_by_PC_per_timestep[g][ap][i_pc_category][i_dt][pc_enrolment_round])){

                    /* Check if there are enough people in the PC sample (including reserves) - otherwise reduce the number of people to be seen and print a warning. */
                    if (patch[p].PC_sample->next_person_to_see[g][ap][i_pc_category]<patch[p].PC_sample->number_in_sample_including_reserves[g][ap][i_pc_category]){

                        enrollee = &(patch[p].individual_population[patch[p].PC_sample->list_ids_potential_enrollees[g][ap][i_pc_category][patch[p].PC_sample->next_person_to_see[g][ap][i_pc_category]]]);
                        i_pc_thisperson = get_PC_HIV_stratum(enrollee);

                        /* We only keep this person if they are still in the right HIV stratum
                         * (so we get the correct number in the cohort who are HIV-, HIV+ unaware and HIV+ aware).
                         * Also make sure they are NOT already enrolled. */
                        //if (i_pc_thisperson==i_pc_category && enrollee->PC_cohort_index==-1){
                        if (i_pc_thisperson==i_pc_category){
                            //printf("Visiting %li\n",enrollee.id);
                            PC_enroll_person(enrollee, patch, t0+t_step*TIME_STEP, p, pc_enrolment_round, g, ap, i_pc_category);

                            patch[p].PC_sample->next_person_to_see[g][ap][i_pc_category]++;  /* Move on to the next person on the list. */
                            j++;  /* We have successfully enrolled someone so increase j. */
                        }
                        else
                            patch[p].PC_sample->next_person_to_see[g][ap][i_pc_category]++;  /* The current person was not eligible so move on to next person (will use reserves). Don't increment j as this was not successful. */
                    }
                    else{
                        patch[p].param->PC_params->number_seen_by_PC_per_timestep[g][ap][i_pc_category][i_dt][pc_enrolment_round]--;
                        /* Also update cohort size as it will be smaller. */
                        patch[p].param->PC_params->cohort_size--;
                        /* Update how many people have been dropped out. */
                        n_dropped_from_cohort++;
                        //printf("Cohort_size reduced by one to = %i\n",patch[p].param->PC_params->cohort_size);
                        patch[p].param->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round]--;
                    }
                }
            }
        }
    }

    //printf("At enrolment timestep %i : Number in [1][15][2] = %i \n",i_dt,patch[p].param->PC_params->number_enrolled_in_PC_round[1][15][2][0]);
    if (n_dropped_from_cohort>0)
        printf("Warning: had to reduce PC sample by %i.\n",n_dropped_from_cohort);
    //printf("number not enrolled = %i\n",number_deleted);
}


/**************************************************************************//**
 * @brief Enroll individual in the Population Cohort
 * 
 * @param indiv Individual to enroll in the PC
 * @param patch Pointer to the patch where the individual is located
 * @param t Time at which the individual is enrolled
 * @param p Index of the patch where the individual is located
 * @param pc_enrolment_round Round of PC enrolment
 * @param g Gender of the individual
 * @param ap Age group of the individual
 * @param i_pc_category HIV status of the individual
*****************************************************************************/

void PC_enroll_person(individual *indiv, patch_struct *patch, double t, int p, int pc_enrolment_round, int g, int ap, int i_pc_category){
    /* Because of the way we draw CHiPs visits at the beginning of the year, it is possible some people
     * die before they are visited. If this is the case then do nothing more.
     * They are deleted from age_list so won't be in next year's sample. */

    if (indiv->cd4==DUMMYVALUE){
        printf("Trying to PC-enrol a non-existent person %d %ld !!! Exiting\n",p,indiv->id);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    if(indiv->id==FOLLOW_INDIVIDUAL && p==FOLLOW_PATCH){
        printf("Enrolling %li at time %6.4lf g=%i ap=%i i_pc_category=%i\n",indiv->id,t,g,ap,i_pc_category);
        fflush(stdout);
    }

    /* Add this person to the cohort. */
    patch[p].PC_cohort->list_ids_in_cohort[g][ap][i_pc_category][patch[p].PC_cohort->number_in_cohort[g][ap][i_pc_category]] = indiv->id;

    /* Increment the counter of number of people. */
    patch[p].PC_cohort->number_in_cohort[g][ap][i_pc_category]++;
    indiv->PC_cohort_index = patch[p].PC_cohort_data->PC_cohort_counter[pc_enrolment_round] + 100000*pc_enrolment_round; /* Person is in PC cohort. */

    /* PC_cohort_counter acts as the index for people in PCX_cohort_data.
     * Set up a temporary local variable (i.e. with shorter name) to make code more readable. */
    int pc_counter = patch[p].PC_cohort_data->PC_cohort_counter[pc_enrolment_round];

    /* Add this person's data to PC_cohort_data: */

    if (pc_enrolment_round==0){
        patch[p].PC_cohort_data->PC0_cohort_data[pc_counter].gender = g;
        patch[p].PC_cohort_data->PC0_cohort_data[pc_counter].ap = ap;
        //printf("pc_counter = %i\n",pc_counter);
        patch[p].PC_cohort_data->PC0_cohort_data[pc_counter].RETAINED_IN_COHORT[0] = 1;
        patch[p].PC_cohort_data->PC0_cohort_data[pc_counter].PC_visit_dates[0] = t;
        patch[p].PC_cohort_data->PC0_cohort_data[pc_counter].id =  indiv->id;
        patch[p].PC_cohort_data->PC0_cohort_data[pc_counter].HIV_status[0] =  indiv->HIV_status;
        patch[p].PC_cohort_data->PC0_cohort_data[pc_counter].ART_status[0] =  indiv->ART_status;
        if (indiv->HIV_status==UNINFECTED)
            patch[p].PC_cohort_data->PC0_cohort_data[pc_counter].serodiscordant_status[0] = indiv->n_HIVpos_partners;
        else
            patch[p].PC_cohort_data->PC0_cohort_data[pc_counter].serodiscordant_status[0] = 0;
        //printf("PC0 Enrolled  %li with pc_id=%i\n",patch[p].PC_cohort_data->PC0_cohort_data[pc_counter].id,pc_counter);
    }
    else if (pc_enrolment_round==1){
        patch[p].PC_cohort_data->PC12N_cohort_data[pc_counter].gender = g;
        patch[p].PC_cohort_data->PC12N_cohort_data[pc_counter].ap = ap;
        patch[p].PC_cohort_data->PC12N_cohort_data[pc_counter].RETAINED_IN_COHORT[1] = 1;
        patch[p].PC_cohort_data->PC12N_cohort_data[pc_counter].PC_visit_dates[1] = t;
        patch[p].PC_cohort_data->PC12N_cohort_data[pc_counter].id =  indiv->id;
        patch[p].PC_cohort_data->PC12N_cohort_data[pc_counter].HIV_status[1] =  indiv->HIV_status;
        patch[p].PC_cohort_data->PC12N_cohort_data[pc_counter].ART_status[1] =  indiv->ART_status;
        //printf("*PC12 Enrolled  %li with pc_id=%i\n",patch[p].PC_cohort_data->PC12N_cohort_data[pc_counter].id,pc_counter);
    }
    else if (pc_enrolment_round==2){
        patch[p].PC_cohort_data->PC24N_cohort_data[pc_counter].gender = g;
        patch[p].PC_cohort_data->PC24N_cohort_data[pc_counter].ap = ap;
        patch[p].PC_cohort_data->PC24N_cohort_data[pc_counter].RETAINED_IN_COHORT[2] = 1;
        patch[p].PC_cohort_data->PC24N_cohort_data[pc_counter].PC_visit_dates[2] = t;
        patch[p].PC_cohort_data->PC24N_cohort_data[pc_counter].id =  indiv->id;
        patch[p].PC_cohort_data->PC24N_cohort_data[pc_counter].HIV_status[2] =  indiv->HIV_status;
        patch[p].PC_cohort_data->PC24N_cohort_data[pc_counter].ART_status[2] =  indiv->ART_status;
        //printf("*PC24 Enrolled  %li with pc_id=%i\n",patch[p].PC_cohort_data->PC24N_cohort_data[pc_counter].id,pc_counter);
    }
    /* Note this must be patch[p].PC_cohort_data->PC_cohort_counter, not the local variable pc_counter. */
    patch[p].PC_cohort_data->PC_cohort_counter[pc_enrolment_round]++;
}


/**************************************************************************//**
 * @brief Carry out the PC visits (AFTER enrolment round) for a given 
 * timestep at time t
 * 
 * @param t0 Current year
 * @param t_step Current timestep
 * @param patch Pointer to a @ref patch_struct
 * @param p Index of the patch
 * @param pc_round Index of the PC round
 * @param pc_enrolment_round Index of the PC enrolment round
 *****************************************************************************/

void carry_out_PC_visits_per_timestep(int t0, int t_step, patch_struct *patch, int p, int pc_round, int pc_enrolment_round){
    int g;
    int ap;
    int i_pc_category; /* Index splitting up the population by HIV status etc. */

    int i_n; /* Temporary store to reduce verbosity of line. */

    /* i_dt tells us the index of the timestep within PC_sample corresponding to the current time t.*/
    int i_dt;

    i_dt = (t0 - patch[p].param->PC_params->PC_START_YEAR[pc_round])*N_TIME_STEP_PER_YEAR + (t_step-patch[p].param->PC_params->PC_START_TIMESTEP[pc_round]);

    /* For debugging: */
    if (i_dt>=MAX_N_TIMESTEPS_PER_PC_ROUND){
        printf("Problem - MAX_N_TIMESTEPS_PER_PC_ROUND is too small. We are %i timesteps into PC round %i. Exiting\n",i_dt,pc_round);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if(i_dt<0||i_dt>=patch[p].param->PC_params->n_timesteps_per_round[pc_round]){
        printf("ERROR i_dt=%i is outside range [0,%i] in PC round %i at time %6.4lf. Exiting\n",i_dt,patch[p].param->PC_params->n_timesteps_per_round[pc_round],pc_round,t0+t_step*TIME_STEP);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    /* Now PC visits each sub-population (currently men/women by year age group ap and HIV status stratum) in turn:
     * Go through the list of id's of people to visit, and visit them: */
    for (g=0;g<N_GENDER;g++){
        /* Run from 18 to 44. */
        for (ap=0; ap<(AGE_PC_MAX-AGE_PC_MIN+1); ap++){
            for (i_pc_category=0; i_pc_category<N_PC_HIV_STRATA; i_pc_category++){
                while ((patch[p].PC_cohort->next_person_to_visit[g][ap][i_pc_category]<patch[p].PC_cohort->next_person_to_visit[g][ap][i_pc_category] + patch[p].param->PC_params->number_seen_by_PC_per_timestep[g][ap][i_pc_category][i_dt][pc_round])
                        && (patch[p].PC_cohort->next_person_to_visit[g][ap][i_pc_category]<patch[p].param->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round])){

                    i_n = patch[p].PC_cohort->next_person_to_visit[g][ap][i_pc_category];
                    /* Visit this person. */
                    if (patch[p].individual_population[patch[p].PC_cohort->list_ids_in_cohort[g][ap][i_pc_category][i_n]].id<=0)
                        printf("Zero ids: id = %li, g=%i ap=%i i_pc = %i, i_dt=%i\n",patch[p].PC_cohort->list_ids_in_cohort[g][ap][i_pc_category][i_n],g,ap,i_pc_category,i_dt);

                    PC_visit_person(&(patch[p].individual_population[patch[p].PC_cohort->list_ids_in_cohort[g][ap][i_pc_category][i_n]]),
                            patch, t0+t_step*TIME_STEP, p, pc_round, pc_enrolment_round); /* Send the address (ie pointer) to this person. */

                    patch[p].PC_cohort->next_person_to_visit[g][ap][i_pc_category]++;
                }

            }
        }
    }
}


/**************************************************************************//**
 * @brief Store new data on person at PC round pc_round representing 
 * a PC visit.
 * 
 * @details pc_enrolment_round reflects that people may have been recruited
 * at different PC rounds (PC0, PC12N, PC24N)
 * 
 * @param indiv Individual being visited as part of the PC
 * @ patch Pointer to a @ref patch_struct
 * @param t Current time
 * @param p Index of the patch
 * @param pc_round Index of the PC round
 * @param pc_enrolment_round Index of the PC enrolment round
 *****************************************************************************/

void PC_visit_person(individual *indiv, patch_struct *patch, double t, int p, 
        int pc_round, int pc_enrolment_round){
    /* Because of the way we draw CHiPs visits at the beginning of the year, it is possible some people
     * die before they are visited. If this is the case then do nothing more.
     * They are deleted from age_list so won't be in next year's sample. */
    if (indiv->cd4==DUMMYVALUE){
        printf("Trying to PC visit a non-existent person %d %ld !!! Exiting\n",p,indiv->id);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    int n = indiv->PC_cohort_index - 100000*pc_enrolment_round;

    if (n<0 || n>=patch[p].param->PC_params->cohort_size){
        printf("ERROR cohort index in PC_visit_person %i %i not defined when visiting %li - exiting\n",n,indiv->PC_cohort_index,indiv->id);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    if (p==0 && pc_enrolment_round==0 && pc_round==1){
        if (indiv->HIV_status>UNINFECTED){
            if (patch[p].PC_cohort_data->PC0_cohort_data[n].HIV_status[0]==UNINFECTED && VERBOSE_OUTPUT==1)
                printf("New PC12 infection id=%li\n",indiv->id);
        }
    }

    if (pc_enrolment_round==0){
        patch[p].PC_cohort_data->PC0_cohort_data[n].RETAINED_IN_COHORT[pc_round] = 1;
        patch[p].PC_cohort_data->PC0_cohort_data[n].PC_visit_dates[pc_round] = t;
        patch[p].PC_cohort_data->PC0_cohort_data[n].HIV_status[pc_round] =  indiv->HIV_status;
        patch[p].PC_cohort_data->PC0_cohort_data[n].ART_status[pc_round] =  indiv->ART_status;
        if (indiv->HIV_status==UNINFECTED)
            patch[p].PC_cohort_data->PC0_cohort_data[n].serodiscordant_status[pc_round] = indiv->n_HIVpos_partners;
        else
            patch[p].PC_cohort_data->PC0_cohort_data[n].serodiscordant_status[pc_round] = 0;
    }else if (pc_enrolment_round==1){
        patch[p].PC_cohort_data->PC12N_cohort_data[n].RETAINED_IN_COHORT[pc_round] = 1;
        patch[p].PC_cohort_data->PC12N_cohort_data[n].PC_visit_dates[pc_round] = t;
        patch[p].PC_cohort_data->PC12N_cohort_data[n].HIV_status[pc_round] =  indiv->HIV_status;
        patch[p].PC_cohort_data->PC12N_cohort_data[n].ART_status[pc_round] =  indiv->ART_status;
    }else if (pc_enrolment_round==2){
        patch[p].PC_cohort_data->PC24N_cohort_data[n].RETAINED_IN_COHORT[pc_round] = 1;
        patch[p].PC_cohort_data->PC24N_cohort_data[n].PC_visit_dates[pc_round] = t;
        patch[p].PC_cohort_data->PC24N_cohort_data[n].HIV_status[pc_round] =  indiv->HIV_status;
        patch[p].PC_cohort_data->PC24N_cohort_data[n].ART_status[pc_round] =  indiv->ART_status;
    }
}


/**************************************************************************//**
 * @brief Return PC round from time point
 * 
 * @param t0 Current year
 * @param t_step Current time step
 * @param patch Patch structure containing info on PC round times.
 * Assumed to have a `patch[0].param->PC_params->PC_END_TIMESTEP[NPC_ROUNDS]`
 * 
 * @returns Index of the current PC round:\n
 * 0 = PC0, 1 = PC12 etc given the current (discrete) time.  A value of
 * -1 is returned if the input time is not a PC round.
 *****************************************************************************/

int get_pc_round(int t0, int t_step, patch_struct *patch, int p){
    
    /* PC is only recruited in patch p=0. */
    if(p > 0){
        return -1;
    }

    /* If before/after PC years return -1. 
    Note that this is slightly more computationally efficient as we deal with most cases first. */
    if(t0 < patch[p].param->PC_params->PC_START_YEAR[0] ||
        t0 > patch[p].param->PC_params->PC_END_YEAR[NPC_ROUNDS - 1]){
        return -1;
    }
    int pc_round = 0;
    while(pc_round < NPC_ROUNDS){
        if((t0 == patch[p].param->PC_params->PC_START_YEAR[pc_round] &&
            t_step >= patch[p].param->PC_params->PC_START_TIMESTEP[pc_round]) ||
            (t0 > patch[p].param->PC_params->PC_START_YEAR[pc_round] &&
            t0 < patch[p].param->PC_params->PC_END_YEAR[pc_round]) ||
            (t0 == patch[p].param->PC_params->PC_END_YEAR[pc_round] &&
            t_step <= patch[p].param->PC_params->PC_END_TIMESTEP[pc_round])){
            return pc_round;
        }
        pc_round = pc_round +1;
    }
    return -1;
}


/**************************************************************************//**
 * @brief Reset PC visit counter for individuals
 * 
 * @param patch pointer to a patch_struct
 * @param p Patch index
 *****************************************************************************/

void reset_visit_counter(patch_struct *patch, int p){
    int g, ap, i_pc_category; /* Indices splitting up the cohort by gender, age, HIV status etc. */
    /* Reset the counter for visits: */
    for (g=0;g<N_GENDER;g++){
        /* Run from 18 to 44 (inclusive). */
        for (ap=0; ap<(AGE_PC_MAX-AGE_PC_MIN+1); ap++){
            for (i_pc_category=0; i_pc_category<N_PC_HIV_STRATA; i_pc_category++){
                /* This is the counter for who we visit next during visits for each round. */
                patch[p].PC_cohort->next_person_to_visit[g][ap][i_pc_category] = 0;
            }
        }
    }
}


/**************************************************************************//**
 * @brief Function used to debug PC sample
 * 
 * @details This function will check the number in each group and compare
 * it against the input.
 * 
 * @param patch pointer to a patch_struct
 * @param p Patch index
 *****************************************************************************/

void check_pc_after_enrollment(patch_struct *patch, int p){
    int Ncohort[N_GENDER][AGE_PC_MAX-AGE_PC_MIN+1][N_PC_HIV_STRATA];
    individual *person;
    int n_id;

    int g, ap, i_pc_category; /* Indices splitting up the cohort by gender, age, HIV status etc. */
    for (g=0; g<N_GENDER; g++){
        for (ap=0; ap<(AGE_PC_MAX-AGE_PC_MIN+1); ap++){
            for (i_pc_category=0; i_pc_category<N_PC_HIV_STRATA; i_pc_category++){
                Ncohort[g][ap][i_pc_category] = 0;
            }
        }
    }
    for (n_id=0; n_id<patch[p].id_counter; n_id++){
        person = &(patch[p].individual_population[n_id]);
        if (person->PC_cohort_index>=0){
            g = patch[p].PC_cohort_data->PC0_cohort_data[person->PC_cohort_index].gender;
            ap = patch[p].PC_cohort_data->PC0_cohort_data[person->PC_cohort_index].ap;
            /* Check status is OK: */
            if (patch[p].PC_cohort_data->PC0_cohort_data[person->PC_cohort_index].HIV_status[0]<UNINFECTED|| patch[p].PC_cohort_data->PC0_cohort_data[person->PC_cohort_index].HIV_status[0]>CHRONIC){
                printf("PANIC! Exiting check_pc_after_enrollment()\n");
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
            /* Look up that person's status in PC0: */
            if (patch[p].PC_cohort_data->PC0_cohort_data[person->PC_cohort_index].HIV_status[0]==UNINFECTED)
                i_pc_category = 0;
            else if (patch[p].PC_cohort_data->PC0_cohort_data[person->PC_cohort_index].HIV_status[0]>UNINFECTED && patch[p].PC_cohort_data->PC0_cohort_data[person->PC_cohort_index].ART_status[0]>ARTNEG && patch[p].PC_cohort_data->PC0_cohort_data[person->PC_cohort_index].ART_status[0]<ARTDEATH)
                i_pc_category = 1;
            else
                i_pc_category = 2;
            Ncohort[g][ap][i_pc_category]++;
        }
    }
}
