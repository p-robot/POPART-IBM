/**************************************************************************//**
 * @file debug.c
 * @brief Functions used for testing the model
*****************************************************************************/

#include "constants.h"
#include "output.h"
#include "debug.h"


/**************************************************************************//**
 * @brief Find individual in age list
 * 
 * @param t Current time
 * @param person_to_find Pointer to an @ref individual struct
 * @param age_list Pointer to an @ref age_list_struct struct
 * @param param Pointer to a @ref parameters struct
 ****************************************************************************/

void find_in_age_list(double t, individual* person_to_find, age_list_struct *age_list, 
    parameters *param){

    int ai_age_calc; /* This is the calculated ai index. */
    /* Index for manual search through age_list looking for indvidiuals. */
    int ai_age_manual_search;
    int found_person = 0;
    int g = person_to_find->gender;
    int aa;
    aa = (int) floor(floor(t) - person_to_find->DoB) - AGE_ADULT;
    
    ai_age_calc = age_list->age_list_by_gender[g]->youngest_age_group_index + aa;
    
    while (ai_age_calc>(MAX_AGE-AGE_ADULT-1)){
        ai_age_calc = ai_age_calc - (MAX_AGE-AGE_ADULT);
    }
    
    int ai_age_calcv1 = get_age_index(person_to_find->DoB, param->start_time_simul);
    int ai_age_calcv2 = get_age_indexv2(person_to_find->DoB, t, 
        age_list->age_list_by_gender[g]->youngest_age_group_index);
    
    /* Looking for the person_to_find in the age_list --> 
    presumably only for debugging, could get rid of this in final code to speed up. 
    */
    
    int n;
    for(ai_age_manual_search = 0; 
        ai_age_manual_search < (MAX_AGE - AGE_ADULT); 
        ai_age_manual_search++
    ){
        n = 0;
        //printf("ai_age_manual_search=%i %li\n",ai_age_manual_search,age_list->number_per_age_group[ai_age_manual_search]);
        while(
            (n < age_list->age_list_by_gender[g]->number_per_age_group[ai_age_manual_search]) && 
            (age_list->age_list_by_gender[g]->age_group[ai_age_manual_search][n]->id != person_to_find->id)
        ){
            n++;
        }
        if((n < age_list->age_list_by_gender[g]->number_per_age_group[ai_age_manual_search]) && 
            (age_list->age_list_by_gender[g]->age_group[ai_age_manual_search][n]->id == person_to_find->id)
        ){
            printf("*********DEBUG DATA at time t=%f : ", t);
            printf("Individual %li Gender %i DoB = %f found with ai = %i. ", 
                    person_to_find->id, g, person_to_find->DoB, ai_age_manual_search);
            
            printf("Calc ai=%i, Calcv1 =%i Calcv2 =%i aa=%i, youngest_index=%i\n",
                ai_age_calc, ai_age_calcv1, ai_age_calcv2, aa, 
                age_list->age_list_by_gender[g]->youngest_age_group_index);
            
            found_person = 1;
            break;
        }
    }
    
    if(found_person == 0){
        n = 0;
        while((n < age_list->age_list_by_gender[g]->number_oldest_age_group) &&
            (age_list->age_list_by_gender[g]->oldest_age_group[n]->id != person_to_find->id)){
            n++;
        }
        if((n < age_list->age_list_by_gender[g]->number_oldest_age_group) && 
            (age_list->age_list_by_gender[g]->oldest_age_group[n]->id == person_to_find->id)){
            
            printf("*********DEBUG DATA: Individual %li DoB = %f gender %d at time t=%f found in oldest age go. Calc ai=%i, aa=%i, youngest_index=%i\n",person_to_find->id,person_to_find->DoB,g,t,ai_age_calc,aa,age_list->age_list_by_gender[g]->youngest_age_group_index);
            found_person = 1;
        }
    }
    if(found_person == 0){
        printf("**********PROBLEM2: person %li from patch %d gender %d not found when searching by hand in age_list %li\n",person_to_find->id,person_to_find->patch_no,g,age_list->age_list_by_gender[g]->number_per_age_group[2]);
        fflush(stdout);
        //printf("age_list->age_group[2][53]->id = %li %li\n",age_list->age_group[2][53]->id,person_to_find->id);

    }
    
    if(found_person == 1){
        printf("Found person %li from patch %d gender %d when searching by hand in age_list\n",
            person_to_find->id, g, person_to_find->patch_no);
    }
}


/**************************************************************************//**
 * @brief For each year, print the ID of each individual in that age
 * @details This function is not currently used.
 * 
 * @param age_list pointer to an @ref age_list_struct struct
 ****************************************************************************/

void print_age_list(age_list_struct *age_list){
    
    int ai;
    int n, g;
    
    // Loop through genders
    for(g = 0; g < N_GENDER; g++){
        printf("Gender = %i\n", g);
        
        // Loop through each year of age
        for(ai = 0; ai < (MAX_AGE - AGE_ADULT); ai++){
            n = 0;
            printf("ai=%i :", ai);
            while((n < age_list->age_list_by_gender[g]->number_per_age_group[ai])){
                printf("%li ", age_list->age_list_by_gender[g]->age_group[ai][n]->id);
                n++;
            }
        }
    }
}


/**************************************************************************//**
 * @brief Count number of alive men and women in a patch, update counters
 * @details The counts of number of males and females are returned in the 
 * objects `n_m_indiv` and `n_f_indiv` respectively.  This function 
 * loops through individuals in the `individual` array of the 
 * @ref patch_struct object.  Used in @ref count_population_size_three_ways().
 * 
 * @param patch : patch_struct Patch structure in the simulation
 * @param n_m_indiv Long in which to store the number of alive males
 * @param n_f_indiv Long in which to store the number of alive females
 ****************************************************************************/

void count_population_by_going_through_indiv(patch_struct *patch, 
    long *n_m_indiv, long *n_f_indiv){
    long i;

    // Reset counters
    *n_m_indiv = 0;
    *n_f_indiv = 0;
    
    // Loop through each individual in the patch in question
    for(i = 0; i < patch->id_counter; i++){
        // Check the individual is alive
        if(patch->individual_population[i].cd4 != DEAD){
            // Check the age of the individual and update the counters
            if(patch->individual_population[i].gender == MALE){
                *n_m_indiv = *n_m_indiv + 1;
            }else{
                *n_f_indiv = *n_f_indiv + 1;
            }
        }
    }
}


/**************************************************************************//**
 * @brief Count the number of individuals in a patch, stratified by gender
 * 
 * @details In contrast to @ref count_population_by_going_through_indiv(), 
 * this function counts the number of individuals in each gender by traversing
 *  the arrays in @ref patch_struct objects:\n
 * `patch_struct->age_list->age_list_by_gender[GENDER]->number_per_age_group[AGEYEAR]`\n
 * `patch_struct->age_list->age_list_by_gender[GENDER]->number_oldest_age_group`.\n
 * Used in @ref count_population_size_three_ways().
 * 
 * @param patch Patch structure in the simulation
 * @param n_m_agelist Long in which to store the number of alive males
 * @param n_f_agelist Long in which to store the number of alive females
 ****************************************************************************/

void count_population_by_going_through_age_list(patch_struct *patch, 
    long *n_m_agelist, long *n_f_agelist){
    int aa;
    // Reset counters
    *n_m_agelist = 0;
    *n_f_agelist = 0;
    for(aa = 0; aa < MAX_AGE - AGE_ADULT; aa++){
        *n_m_agelist = *n_m_agelist + 
            patch->age_list->age_list_by_gender[MALE]->number_per_age_group[aa];
        
        *n_f_agelist = *n_f_agelist + 
            patch->age_list->age_list_by_gender[FEMALE]->number_per_age_group[aa];
    }
    
    *n_m_agelist = *n_m_agelist + 
        patch->age_list->age_list_by_gender[MALE]->number_oldest_age_group;
    
    *n_f_agelist = *n_f_agelist + 
        patch->age_list->age_list_by_gender[FEMALE]->number_oldest_age_group;
}


/**************************************************************************//**
 * @brief Count population of men and women using `pop_size_per_gender_age_risk`
 * 
 * @details Counts the popluation of men and women that is stored in the array
 * `pop_size_per_gender_age_risk[GENDER][NAGE][NRISK]` in a `population_size`
 * struct.  Used to check calculations are correct across different data 
 * structures in the model.  No return value, changes objects in place.
 * Used in @ref count_population_size_three_ways().
 * 
 * @param pop Structure for different stratifications of population size
 * @param n_m_pop Number of men in the population
 * @param n_f_pop Number of women in the population
 ****************************************************************************/

void count_population_using_n_population(population_size *pop, long *n_m_pop, long *n_f_pop){
    int ag, r;
    // Reset counters
    *n_m_pop = 0;
    *n_f_pop = 0;
    for(ag=0 ; ag<N_AGE ; ag++){
        for(r=0 ; r<N_RISK ; r++){
            *n_m_pop = *n_m_pop + pop->pop_size_per_gender_age_risk[MALE][ag][r];
            *n_f_pop = *n_f_pop + pop->pop_size_per_gender_age_risk[FEMALE][ag][r];
        }
    }
}


/**************************************************************************//**
 * @brief Count the patch population using three different methods and compare
 * 
 * @details This function counts the population of men and women who are alive
 * in a patch using the functions @ref count_population_by_going_through_indiv(), 
 * @ref count_population_by_going_through_age_list(), and 
 * @ref count_population_using_n_population().  If the population sizes differ
 * from one of these methods then an error is thrown and the program exist.
 * No return value.
 * 
 * @param patch Patch structure holding the population data
 * @param p Patch number in which to count the population
 * @param t Time step at which to count the population
 ****************************************************************************/

void count_population_size_three_ways(patch_struct *patch, int p, double t){
    long n_m_indiv, n_f_indiv;
    count_population_by_going_through_indiv(&patch[p], &n_m_indiv, &n_f_indiv);

    long n_m_agelist, n_f_agelist;
    count_population_by_going_through_age_list(&patch[p], &n_m_agelist, &n_f_agelist);

    long n_m_pop, n_f_pop;
    count_population_using_n_population(patch[p].n_population, &n_m_pop, &n_f_pop);

    if ((n_m_indiv!=n_m_agelist) || (n_m_indiv!=n_m_pop) || (n_f_indiv!=n_f_agelist) || (n_f_indiv!=n_f_pop)){
        printf("Error - population sizes not matching in patch %i t=%6.4lf male: %ld %ld %ld female: %ld %ld %ld\n",
            p,t,n_m_indiv,n_m_agelist,n_m_pop,n_f_indiv,n_f_agelist,n_f_pop);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
}


/**************************************************************************//**
 * @brief Find individual in patch by ID
 * 
 * @param id ID of individual to find
 * @param patch Patch structure in the simulation
 * @param p Patch number in which to search for the individual
 * 
 * @return int Flag for whether the individual was found (1) in the 
 * patch or not (0)
 ****************************************************************************/

int is_in_patch_individual_population(long id, patch_struct *patch, int p){
    long k = 0;
    int res = 0;

    while((patch[p].individual_population[k].id!=id) && (k<(patch[p].id_counter-1))){
        k++;
    }
    if(k<patch[p].id_counter-1){
        res = 1;
        printf("Individual %ld was found in patch %d at position %ld in list\n",id,p,k);
        fflush(stdout);
    }else{
        printf("Individual %ld was not found in patch %d\n",id,p);
        fflush(stdout);
    }
    return(res);
}


/**************************************************************************//**
 * @brief Function generates the text that goes in the header to the 
 * file `Age_distribution_check_CL....csv`
 * 
 * @details The files `Age_distribution_check_CL....csv` are used for 
 * validation of the demographics of the model.  This function generates 
 * the header for the file and therefore has no return value.
 * 
 * @param age_group_string Character array to store the header text
 * @param size_age_group_string Length of the character array
 ****************************************************************************/

void generate_demographics_byage_gender_file_header(char *age_group_string, int size_age_group_string){
    int g,ag;
    char temp_string[100];
    sprintf(age_group_string, "Time,");
    for (g=0;g<N_GENDER;g++){
        for (ag=0; ag<N_AGE_UNPD; ag++){
            if(g==MALE)
                sprintf(temp_string,"M:%i-%i,",AGE_GROUPS_UNPD[ag],AGE_GROUPS_UNPD[ag+1]-1);
            else
                sprintf(temp_string,"F:%i-%i,",AGE_GROUPS_UNPD[ag],AGE_GROUPS_UNPD[ag+1]-1);
            join_strings_with_check(age_group_string, temp_string, size_age_group_string, 
                "temp_string and age_group_string in generate_demographics_byage_gender_file_header");
        }
        if(g==MALE)
            sprintf(temp_string,"M%i+,",AGE_GROUPS_UNPD[N_AGE_UNPD]);
        else
            sprintf(temp_string,"F:%i+,",AGE_GROUPS_UNPD[N_AGE_UNPD]);
        join_strings_with_check(age_group_string, temp_string, size_age_group_string,
            "temp_string and age_group_string in generate_demographics_byage_gender_file_header");
    }
    for (g=0;g<N_GENDER;g++){
        for (ag=0; ag<N_AGE_UNPD; ag++){
            if(g==MALE)
                sprintf(temp_string,"DeadM:%i-%i,",AGE_GROUPS_UNPD[ag],AGE_GROUPS_UNPD[ag+1]-1);
            else
                sprintf(temp_string,"DeadF:%i-%i,",AGE_GROUPS_UNPD[ag],AGE_GROUPS_UNPD[ag+1]-1);
            join_strings_with_check(age_group_string, temp_string, size_age_group_string, 
                "temp_string and age_group_string in generate_demographics_byage_gender_file_header");
        }
        if(g==MALE)
            sprintf(temp_string,"DeadM%i+,",AGE_GROUPS_UNPD[N_AGE_UNPD]);
        else
            sprintf(temp_string,"DeadF:%i+,",AGE_GROUPS_UNPD[N_AGE_UNPD]);
        join_strings_with_check(age_group_string, temp_string, size_age_group_string, 
            "temp_string and age_group_string in generate_demographics_byage_gender_file_header");
    }
    sprintf(temp_string,"CumulativeDead\n");
    join_strings_with_check(age_group_string, temp_string, size_age_group_string, 
        "temp_string and age_group_string in generate_demographics_byage_gender_file_header");
}


/**************************************************************************//**
 * @brief Writes to file the number of people by age group and gender
 * 
 * @details Function goes through `individual_population` attribute within 
 * a @ref patch_struct object in a given patch at a given time 
 * (currently called yearly in @ref main.c).
 * Generates a csv file (for each run and patch) called 
 * `Age_distribution_check... .csv` which gives the number of adults in each 
 * 5 year UNPD age group by gender to compare against UNPD age distribution 
 * estimates for model validation.  Model also outputs number of dead people 
 * in each age group (their current age, not age ever).  Writes output data
 * so no return value.
 * 
 * @param patch Patch structure in the simulation
 * @param p Patch number in which to count the population
 * @param t Time step at which to count the population
 * @param file_data_store File structure for storing output data
 ****************************************************************************/

void write_demographics_byage_gender(patch_struct *patch, int p, double t, file_struct *file_data_store){
    long i;
    int ag;
    file_data_store->AGEDISTRIBUTIONFILE[p] = fopen(file_data_store->filename_debug_agedistribution[p],"a");
    if (file_data_store->AGEDISTRIBUTIONFILE[p]==NULL){
        printf("Cannot open output_file in write_demographics_byage_gender().\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    // Number of men/women in each age group (age groups correspond to UNPD 5 
    // year age groups apart from 13-14 year olds. Note that we DO NOT want to 
    // count children as childhood mortality is taken to occur at birth so the 
    // number of children in the model is the number of children who will 
    // survive to age 13, not the number of children alive at any given time.
    long  *num_age_people_m;
    long  *num_age_people_f;
    /* Number of deaths in the past year in men/women by age group : */
    long  *num_age_newdeaths_m;
    long  *num_age_newdeaths_f;
    /* This is a check that if we sum over 
    num_age_newdeaths_m and num_age_newdeaths_f we get the new deaths per 
    year. num_deaths_total is a cumulative measure so the difference 
    between one year and the next is incident deaths. */
    long num_deaths_total = 0;

    /* Use calloc() so that these are initialized to zero. */
    num_age_people_m = (long*)calloc(N_AGE_UNPD+1, sizeof(long));
    num_age_people_f = (long*)calloc(N_AGE_UNPD+1, sizeof(long));
    num_age_newdeaths_m = (long*)calloc(N_AGE_UNPD+1, sizeof(long));
    num_age_newdeaths_f = (long*)calloc(N_AGE_UNPD+1, sizeof(long));

    for (i=0; i<patch[p].id_counter; i++){
        if (patch[p].individual_population[i].cd4>DEAD){
            if (patch[p].individual_population[i].gender==MALE)
                num_age_people_m[get_age_group_unpd(patch[p].individual_population[i].DoB,t)]++;
            else if (patch[p].individual_population[i].gender==FEMALE)
                num_age_people_f[get_age_group_unpd(patch[p].individual_population[i].DoB,t)]++;
            else{
                printf("ERROR: Unknown gender. Exiting\n");
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
        }else{
            num_deaths_total += 1;
            if (patch[p].individual_population[i].DoD<1800){
                printf("Error: someone has died before 1800. Exiting.\n");
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
            /* write_demographics_byage_gender() is called at the 
            end of the year (so t is the start of the year), 
            this counts anyone who died in the last year. */
            else if(patch[p].individual_population[i].DoD>=t){
                if (patch[p].individual_population[i].gender==MALE)
                    num_age_newdeaths_m[get_age_group_unpd(patch[p].individual_population[i].DoB,t)]++;
                else if (patch[p].individual_population[i].gender==FEMALE)
                    num_age_newdeaths_f[get_age_group_unpd(patch[p].individual_population[i].DoB,t)]++;
                else{
                    printf("ERROR: Unknown gender in dead person. Exiting\n");
                    printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                    fflush(stdout);
                    exit(1);
                }
            }
        }
    }

    fprintf(file_data_store->AGEDISTRIBUTIONFILE[p],"%6.4f,",t+1);
    for (ag=0; ag<=N_AGE_UNPD; ag++)
        fprintf(file_data_store->AGEDISTRIBUTIONFILE[p],"%li,",num_age_people_m[ag]);
    for (ag=0; ag<=N_AGE_UNPD; ag++)
        fprintf(file_data_store->AGEDISTRIBUTIONFILE[p],"%li,",num_age_people_f[ag]);

    for (ag=0; ag<=N_AGE_UNPD; ag++)
        fprintf(file_data_store->AGEDISTRIBUTIONFILE[p], "%li,",num_age_newdeaths_m[ag]);
    for (ag=0; ag<=N_AGE_UNPD; ag++)
        fprintf(file_data_store->AGEDISTRIBUTIONFILE[p], "%li,",num_age_newdeaths_f[ag]);

    fprintf(file_data_store->AGEDISTRIBUTIONFILE[p], "%li\n",num_deaths_total);
    fclose(file_data_store->AGEDISTRIBUTIONFILE[p]);

    free(num_age_people_m);
    free(num_age_people_f);
    free(num_age_newdeaths_m);
    free(num_age_newdeaths_f);
}


/***************************************************************************//**
 * @brief Write blank file and header for the file called `OneYearAgeGp_*.csv`
 * 
 * @details For now we only use data from patch 0
 * 
 * @param file_data_store pointer to @ref file_struct structure that stores 
 * the file names of the files to write.  The file in
 * question has its file name stored in the attribute called 
 * `file_data_store->ONEYEARAGEDISTRIBUTIONFILE`
 ****************************************************************************/

void blank_one_year_age_groups_including_kids(file_struct *file_data_store){
    // Open a connection to the file
    file_data_store->ONEYEARAGEDISTRIBUTIONFILE[0] =
        fopen(file_data_store->filename_debug_one_yearage_dist_includingkids[0], "w");
    
    // Throw an error if we can't open the file.  
    if(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[0] == NULL){
        printf("Cannot open file_data_store->ONEYEARAGEDISTRIBUTIONFILE[0]\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    // Print header for year (time)
    fprintf(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[0], "t,");
    
    int i;
    // Write headers for children (they don't have a gender in the model)
    for(i = 0; i < AGE_ADULT + 1; i++){
        fprintf(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[0],"Nage%i-%i,", i, i + 1);
    }
    
    // Write headers for males
    for(i = AGE_ADULT + 1; i < MAX_AGE; i++){
        fprintf(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[0],"NMage%i-%i,", i, i + 1);
    }
    fprintf(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[0], "NMage80+,");
    
    // Write headers for females
    for(i = AGE_ADULT + 1; i < MAX_AGE; i++){
        fprintf(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[0],"NFage%i-%i,", i, i + 1);
    }
    fprintf(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[0], "NFage80+\n");
    
    fclose(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[0]);
    return;
}


/***************************************************************************//**
 * @brief Write population size per yearly age groups to file for a single year 
 * (append to file)
 * 
 * @details This function writes the total population size for children from 
 * ages 0-1 until 13-14 at yearly age gaps (children aren't split by gender
 * in the simulation), population sizes are then recorded from 14-15 until 80+
 * at yearly age groups split for males and females.  These values are taken
 * from the structures `child_population` and `age_list->age_list_by_gender` 
 * within the patch structure respectively.  The file in which these values
 * are appended starts with the prefix 'OneYearAgeGp_*.csv'.  This function
 * is called from @ref main.c and is only triggered if the macro (defined in 
 * @ref constants.h)` WRITE_DEBUG_DEMOGRAPHICS_AGE_DISTRIBUTION_ONEYEARINCKIDS`
 * is 1.  The function @ref blank_one_year_age_groups_including_kids() sets
 * up the file and header line that this function writes to.
 * 
 * The `child_population` array that's part of the @ref patch_struct structure
 * is indexed by HIV status of the children, 0 index for HIV- and 1 index for 
 * HIV-positive individuals.  However, at this stage (Feb 2018) all children
 * are HIV- so the second entry in the array is all zeros.
 * 
 * @param file_data_store pointer to a @ref file_struct structure that stores 
 * the file names of the files to write
 * @param patch pointer to a @ref patch_struct structure
 * @param p Patch number of interest
 * @param t Time step at which to write the data
 ****************************************************************************/

void write_one_year_age_groups_including_kids(file_struct *file_data_store, patch_struct *patch,
    int p, double t){
    int i_child, ai_m, ai_f, aa, hivstatus;
    
    // Open a connection to the file in question and append to it ("a")
    file_data_store->ONEYEARAGEDISTRIBUTIONFILE[p] =
        fopen(file_data_store->filename_debug_one_yearage_dist_includingkids[p], "a");
    
    // Throw an error if the file can't be opened
    if(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[p] == NULL){
        printf("Cannot open file_data_store->ONEYEARAGEDISTRIBUTIONFILE[p]\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    
    // Append the year
    fprintf(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[p], "%6.4lf,", t);

    // The 2 is because we have HIV- and HIV+ children
    int child_start_index[2];
    
    // If the oldest child group is at the end, then the youngest age group is [0]
    for(hivstatus = 0; hivstatus < 2; hivstatus++){
        
        if(
        patch[p].child_population[hivstatus].debug_tai == (AGE_ADULT + 1) * N_TIME_STEP_PER_YEAR - 1
        ){
            child_start_index[hivstatus] = 0;
        }else{
            // Otherwise the youngest age group is the one after the oldest age group
            child_start_index[hivstatus] = patch[p].child_population[hivstatus].debug_tai + 1;
        }
    }

    // Write to disk the number of kids by 1 year age groups
    long nkids_by_year = 0;
    
    for(i_child = 0; i_child < ((AGE_ADULT + 1) * N_TIME_STEP_PER_YEAR); i_child++){
        
        // Sum up number of kids by yearly age group
        for(hivstatus = 0; hivstatus < 2; hivstatus++){
            nkids_by_year +=
                patch[p].child_population[hivstatus].n_child[child_start_index[hivstatus]];
        }
        
        /* Print out number of kids and then reset counter for next yearly age group: */
        if((i_child % N_TIME_STEP_PER_YEAR) == (N_TIME_STEP_PER_YEAR - 1)){
            fprintf(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[p], "%li,", nkids_by_year);
            nkids_by_year = 0;
        }
        
        /* Now update pointer to loop over all kids: */
        for(hivstatus = 0; hivstatus < 2; hivstatus++){
            if(child_start_index[hivstatus] < (AGE_ADULT + 1) * N_TIME_STEP_PER_YEAR - 1){
                child_start_index[hivstatus]++;
            }else{
                child_start_index[hivstatus] = 0;
            }
        }
    }
    
    /* Write number of adults to disk.
    Note we don't print the age 13 group here as we call this routine at the start of the year -
    so no adults aged 13 yet.  We count the age 13-14 group from the kids above.  Hence the range
    of aa is 0 to (MAX_AGE-AGE_ADULT-1) instead of (MAX_AGE-AGE_ADULT), and we add 1 to ai_m and
    ai_f
    */
    // Output count of individuals per age for males
    for(aa = 0; aa < (MAX_AGE - AGE_ADULT - 1); aa++){
        ai_m = aa + patch[p].age_list->age_list_by_gender[MALE]->youngest_age_group_index + 1;
        
        while(ai_m > (MAX_AGE - AGE_ADULT - 1)){
            ai_m = ai_m - (MAX_AGE - AGE_ADULT);
        }
        fprintf(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[p], "%li,", 
            patch[p].age_list->age_list_by_gender[MALE]->number_per_age_group[ai_m]);
    }
    fprintf(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[p], "%li,", 
        patch[p].age_list->age_list_by_gender[MALE]->number_oldest_age_group);
    
    // Output count of individuals per age for females
    for(aa = 0; aa < (MAX_AGE - AGE_ADULT - 1); aa++){
        ai_f = aa + patch[p].age_list->age_list_by_gender[FEMALE]->youngest_age_group_index + 1;
        
        while(ai_f > (MAX_AGE - AGE_ADULT - 1)){
            ai_f = ai_f - (MAX_AGE - AGE_ADULT);
        }
        fprintf(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[p], "%li,", 
            patch[p].age_list->age_list_by_gender[FEMALE]->number_per_age_group[ai_f]);
    }
    fprintf(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[p], "%li\n", 
        patch[p].age_list->age_list_by_gender[FEMALE]->number_oldest_age_group);
    
    // Close the connection to the file
    fclose(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[p]);
}


/***************************************************************************//**
 * @brief Write the number of births, deaths, and new adults to file
 * 
 * @details Writes the number of births, deaths and new adults, kept in the 
 * attributes `DEBUG_NBIRTHS`, `DEBUG_NDEATHS`, and `DEBUG_NNEWADULTS` of a
 * @ref patch_struct structure, to a file.
 * @param file_data_store pointer to a @ref file_struct structure that stores
 * names and location of files
 * @param patch pointer to a @ref patch_struct structure
 * @param year Year at which to write the data
 ****************************************************************************/

void write_nbirths_nnewadults_ndeaths(file_struct *file_data_store, patch_struct *patch, int year){
    file_data_store->NBIRTHS_NNEWADULTS_NDEATHS_FILE = fopen(file_data_store->filename_debug_nnewadults_ndeaths_file ,"a");
    if(file_data_store->NBIRTHS_NNEWADULTS_NDEATHS_FILE==NULL){
        printf("Cannot open file_data_store->NBIRTHS_NNEWADULTS_NDEATHS_FILE\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    fprintf(file_data_store->NBIRTHS_NNEWADULTS_NDEATHS_FILE, "%d,", year);
    int p;
    for(p=0; p<NPATCHES; p++){
        fprintf(file_data_store->NBIRTHS_NNEWADULTS_NDEATHS_FILE,"%ld,%ld,%ld,",
            patch[p].DEBUG_NBIRTHS,patch[p].DEBUG_NNEWADULTS,patch[p].DEBUG_NDEATHS);
    }
    fprintf(file_data_store->NBIRTHS_NNEWADULTS_NDEATHS_FILE,"\n");
    fclose(file_data_store->NBIRTHS_NNEWADULTS_NDEATHS_FILE);
}


/***************************************************************************//**
 * @brief Write the duration of life of people in 5 year age cohorts to file.
 * 
 * @details Function needs to be run for a long simulation time (e.g. >2100).
 * Function not currently used, but is intended to write the number of
 * individuals in each age group to a file.  Used for debugging/testing.
 * 
 * @param output_file_directory Directory in which to write the file
 * @param patch Pointer to a @ref patch_struct structure
 * @param p Patch number of interest
 * @param i_run Number of the current run
 ****************************************************************************/

void output_life_expectancy(char *output_file_directory, patch_struct *patch, int p, int i_run){
    long i_id;
    int i_a;   /* Indexes which birth cohort you are in. For now make 5 year age cohorts (ie born 1900-1905, 1905-1910,...). */
    double age_at_death;
    /* If simulation run for too short a time this won't work: */
    if(patch[p].param->end_time_simul<2100){
        printf("Warning: need to make end_time_simul at least 2100 to run output_life_expectancy(). This function has not been run even though WRITE_DEBUG_DEMOGRAPHICS_LIFE_EXPECTANCY==1\n");
        return;
    }

    FILE *life_expectancy_file;
    char life_expectancy_filename[LONGSTRINGLENGTH];
    char templabel[LONGSTRINGLENGTH];     /* Length is arbitrary - but make sure we never have stupidly long filenames! */

    /* Make sure the arrays are blank: */
    memset(life_expectancy_filename, '\0', sizeof(life_expectancy_filename));
    memset(templabel, '\0', sizeof(templabel));

    /* Assume that country setting is same in all patches so use patch 0. */
    if (patch[0].country_setting==ZAMBIA)
        sprintf(templabel,"_Za.csv");
    else
        sprintf(templabel,"_SA.csv");

    strncpy(life_expectancy_filename,output_file_directory,LONGSTRINGLENGTH);
    add_slash(life_expectancy_filename); /* Adds a / or \ as needed if working in directory other than current local dir. */
    join_strings_with_check(life_expectancy_filename, "LifeExpectancy", LONGSTRINGLENGTH, 
        "'LifeExpectancy' and life_expectancy_filename in output_life_expectancy()");
    join_strings_with_check(life_expectancy_filename, templabel, LONGSTRINGLENGTH, 
        "templabel and life_expectancy_filename in output_life_expectancy()");
    printf("LE filename = %s\n",life_expectancy_filename);

    /* For simplicity we store all the data from every run in a single file. This means that we need to make the first run
     * (i_run=1) blank.
     */
    if (i_run==1){
        life_expectancy_file = fopen(life_expectancy_filename,"w");
        if (life_expectancy_file==NULL){
            printf("Cannot open life_expectancy_file\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        /* Print a header for the file: */
        fprintf(life_expectancy_file,"i_run,");
        for (i_a=0;i_a<20;i_a++){
            fprintf(life_expectancy_file,"%i-%i,",1900+i_a*5,1904+i_a*5);
        }
        fprintf(life_expectancy_file,"\n");
    }
    else{ /* if not the first run then just append to the existing file: */
        life_expectancy_file = fopen(life_expectancy_filename,"a");
        if (life_expectancy_file==NULL){
            printf("Cannot open life_expectancy_file\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
    }

    long n_by_age_cohort[20];            /* Store number of individual in each age cohort from 1900-2000 - so 20 age groups. */
    double cumulative_life_years[20];    /* Store number of life-years lived of people in cohort. */
    double adjusted_life_expectancy[20]; /* Life expectancy of each cohort adjusting for childhood mortality. */
    for (i_a=0;i_a<20;i_a++){
        n_by_age_cohort[i_a] = 0;
        cumulative_life_years[i_a] = 0;
    }

    /* Use this to estimate remaining life expectancy of someone still alive in 2000: */
    double annual_mortality_over80_in2000 = 0.5*(exp(patch[p].param->mortality_rate_by_gender_age_intercept[MALE][N_AGE_UNPD_MORTALITY-1] + patch[p].param->mortality_rate_by_gender_age_slope[MALE][N_AGE_UNPD_MORTALITY-1]*2000) + exp(patch[p].param->mortality_rate_by_gender_age_intercept[FEMALE][N_AGE_UNPD_MORTALITY-1] + patch[p].param->mortality_rate_by_gender_age_slope[FEMALE][N_AGE_UNPD_MORTALITY-1]*2000));
    //printf("annual_mortality_over80_in2000 = %lf\n",annual_mortality_over80_in2000); /* Checked in Zambia this was 0.16, so 6 years expected future life. */

    for (i_id=0;i_id<patch[p].id_counter;i_id++){
        /* Ignore anyone born before 1900 - ie the adults at initialisation - as we have incomplete mortality data for them.
         * Similarly assume that we may not have a good sample of those born after ~2000 as some may not die by the end of the simulation. */
        if ((patch[p].individual_population[i_id].DoB>1900) && (patch[p].individual_population[i_id].DoB<2000)){
            /* Is the person dead by the end of the simulation? */
            if (patch[p].individual_population[i_id].DoD>-1){
                age_at_death = patch[p].individual_population[i_id].DoD - patch[p].individual_population[i_id].DoB;
                i_a = (int) floor((patch[p].individual_population[i_id].DoB-1900)/5);
                n_by_age_cohort[i_a]++;                       /* Add one to counter in the age cohort i_a. */
                cumulative_life_years[i_a] += age_at_death;   /* Add the life years of that individual to age cohort i_a. */
            }else{
                i_a = (int) floor((patch[p].individual_population[i_id].DoB-1900)/5);
                //printf("Individual %li in age gp %i not dead yet DoB=%6.4f\n",patch[p].individual_population[i_id].id,i_a,patch[p].individual_population[i_id].DoB);
                n_by_age_cohort[i_a]++;                       /* Add one to counter in the age cohort i_a. */
                cumulative_life_years[i_a] += 100 + 1.0/annual_mortality_over80_in2000;   /* Add remaining expected number of life years of that individual to age cohort i_a. */
            }
        }
    }
    fprintf(life_expectancy_file,"Run%iUnadjusted,",i_run);
    for (i_a=0;i_a<20;i_a++){
        fprintf(life_expectancy_file,"%8.6lf,",cumulative_life_years[i_a]/(n_by_age_cohort[i_a]*1.0));
    }
    fprintf(life_expectancy_file,"\n");

    /* Now adjust estimate for childhood mortality: */
    fprintf(life_expectancy_file,"Run%iAdjusted,",i_run);
    double annual_mortality_rate_under5, mortality_rate_under5;
    double annual_mortality_rate_5to10,  mortality_rate_5to10;
    int g;
    double t;
    for (i_a=0;i_a<20;i_a++){
        t = 1902.5+i_a*5; /* Take mid-point of age cohort birth date and get corresponding mortality: */
        annual_mortality_rate_under5 = 0.0;
        annual_mortality_rate_5to10 = 0.0;
        for (g=0;g<N_GENDER;g++){
            /* We by default assume that mortality in under 5 is mostly perinatal mortality (so occurs at time t) and that mortality in 5-10 year olds occurs uniformly over that age so on average at time t+7.5.
             * However we need to adjust these times for the fact that we only have data from 1950-2100. */
            if (t<1950){
                /* Average mortality over genders. The [0] and [1] indices refer to age groups 0-4 and 5-9: */
                annual_mortality_rate_under5 += exp(patch[p].param->mortality_rate_by_gender_age_intercept[g][0] + patch[p].param->mortality_rate_by_gender_age_slope[g][0]*1950);
                annual_mortality_rate_5to10 += exp(patch[p].param->mortality_rate_by_gender_age_intercept[g][1] + patch[p].param->mortality_rate_by_gender_age_slope[g][1]*1957.5);
            }else if (t>=2100){ /* Note - should never need t>2000. This is copied from a function in demographics, and keep for completeness. */
                annual_mortality_rate_under5 += exp(patch[p].param->mortality_rate_by_gender_age_intercept[g][0] + patch[p].param->mortality_rate_by_gender_age_slope[g][0]*2100);
                annual_mortality_rate_5to10 += exp(patch[p].param->mortality_rate_by_gender_age_intercept[g][1] + patch[p].param->mortality_rate_by_gender_age_slope[g][1]*2100);
            }else if (t>=2092.5){
            /* for times 2092.5-2100 */
                annual_mortality_rate_under5 += exp(patch[p].param->mortality_rate_by_gender_age_intercept[g][0] + patch[p].param->mortality_rate_by_gender_age_slope[g][0]*t);
                annual_mortality_rate_5to10 += exp(patch[p].param->mortality_rate_by_gender_age_intercept[g][1] + patch[p].param->mortality_rate_by_gender_age_slope[g][1]*2100);
            }else{
            /* for times 1950-2092.5: */
                annual_mortality_rate_under5 += exp(patch[p].param->mortality_rate_by_gender_age_intercept[g][0] + patch[p].param->mortality_rate_by_gender_age_slope[g][0]*t);
                annual_mortality_rate_5to10 += exp(patch[p].param->mortality_rate_by_gender_age_intercept[g][1] + patch[p].param->mortality_rate_by_gender_age_slope[g][1]*(t+7.5));
            }
        }
        /* We take the average of the annual_mortality_rates over gender: */
        mortality_rate_under5 = 1- pow(1-annual_mortality_rate_under5/(N_GENDER*1.0),5);
        mortality_rate_5to10 = 1- pow(1-annual_mortality_rate_5to10/(N_GENDER*1.0),5);
        adjusted_life_expectancy[i_a] = cumulative_life_years[i_a]/(n_by_age_cohort[i_a]*1.0)*(1-mortality_rate_under5)*(1-mortality_rate_5to10) + 7.5*(1-mortality_rate_under5)*mortality_rate_5to10 + 0*mortality_rate_under5;
        fprintf(life_expectancy_file,"%8.6lf,",adjusted_life_expectancy[i_a]);
    }
    fprintf(life_expectancy_file,"\n");
    fclose(life_expectancy_file);
}


/***************************************************************************//**
 * @brief Check if an individual is alive, HIV-negative and has HIV-positive
 * partners, and if so whether they are in the right place in the list of 
 * susceptibles in a serodiscordant partnership
 * 
 * @param temp_ind Pointer to an individual
 * @param overall_partnerships Pointer to an @ref all_partnerships structure
 ****************************************************************************/

void check_if_individual_should_be_in_list_susceptibles_in_serodiscordant_partnership(
        individual *temp_ind, all_partnerships * overall_partnerships){
    int i;
    int isInSerodiscordantCouple = 0;
    if(temp_ind->cd4 > DEAD && temp_ind->HIV_status==0){
        if(temp_ind->n_partners>0){
            for(i=0 ; i<temp_ind->n_partners ; i++){ // loop over all partners in case there is something wrong with HIV+ partners tracking
                // is partner HIV+?
                if(temp_ind->partner_pairs[i]->ptr[1-temp_ind->gender]->HIV_status>0){
                    if(temp_ind->id==FOLLOW_INDIVIDUAL  && temp_ind->patch_no==FOLLOW_PATCH){
                        printf("seropositive partner: ");
                        print_individual(temp_ind->partner_pairs[i]->ptr[1-temp_ind->gender]);
                        fflush(stdout);
                    }
                    isInSerodiscordantCouple = 1;
                }
            }
        }
        if(isInSerodiscordantCouple>0){
            // check this person is in the right place in the list of susceptibles in a serodiscordant partnership
            if(temp_ind->idx_serodiscordant<0 || temp_ind->idx_serodiscordant>=overall_partnerships->n_susceptible_in_serodiscordant_partnership[0]){
                printf("PROBLEM: individual %ld from patch %d is not at all in the list of susceptible individuals in a serodiscordant couple\n",temp_ind->id,temp_ind->patch_no);
                print_individual(temp_ind);
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }else if(overall_partnerships->susceptible_in_serodiscordant_partnership[temp_ind->idx_serodiscordant]->id != temp_ind->id ||  overall_partnerships->susceptible_in_serodiscordant_partnership[temp_ind->idx_serodiscordant]->patch_no != temp_ind->patch_no){
                printf("PROBLEM: individual %ld from patch %d is not found where should be in the list of susceptible individuals in a serodiscordant couple\n",temp_ind->id,temp_ind->patch_no);
                print_individual(temp_ind);
                printf("BUT the person pointed to in the list is at idx %li:\n",temp_ind->idx_serodiscordant);
                print_individual(overall_partnerships->susceptible_in_serodiscordant_partnership[temp_ind->idx_serodiscordant]);
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
        }
    }
}


/***************************************************************************//**
 * @brief Check if individual is alive and has available partnerships, 
 * whether he is in the right place in the list of available partners
 * 
 * @param temp_ind Pointer to an individual
 * @param overall_partnerships Pointer to an @ref all_partnerships structure
 * @param t0 Start time of the simulation
 * @param t_step Time step of the simulation
 ****************************************************************************/

void check_if_individual_should_be_in_list_available_partners(
        individual *temp_ind, all_partnerships * overall_partnerships, int t0, int t_step){
    int i;
    long temp_id;
    int temp_patch;
    int ag;
    int n_times_found_in_list_available_partners;

    if(temp_ind->cd4 > DEAD && temp_ind->n_partners<temp_ind->max_n_partners){
        n_times_found_in_list_available_partners = 0;
        // find ag the age group of temp_ind
        ag = get_age_group(temp_ind->DoB,t0+t_step*TIME_STEP, AGE_GROUPS, N_AGE);

        for(i=0 ; i<overall_partnerships->n_pop_available_partners->pop_per_patch[temp_ind->patch_no].pop_size_per_gender_age_risk[temp_ind->gender][ag][temp_ind->sex_risk] ; i++){
            temp_id = overall_partnerships->pop_available_partners->pop_per_patch_gender_age_risk[temp_ind->patch_no][temp_ind->gender][ag][temp_ind->sex_risk][i]->id;
            temp_patch = overall_partnerships->pop_available_partners->pop_per_patch_gender_age_risk[temp_ind->patch_no][temp_ind->gender][ag][temp_ind->sex_risk][i]->patch_no;

            if(temp_id == temp_ind->id && temp_patch == temp_ind->patch_no){
                n_times_found_in_list_available_partners++;
            }
        }
        if(n_times_found_in_list_available_partners != temp_ind->max_n_partners - temp_ind->n_partners)
        {
            printf("PROBLEM: individual %ld from patch %d is not found as many time as expected in the list of available partners\n",temp_ind->id,temp_ind->patch_no);
            print_individual(temp_ind);
            printf("Individual was found %d times in the list but has %d available partnerships\n",n_times_found_in_list_available_partners,temp_ind->max_n_partners - temp_ind->n_partners);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
    }
}


/***************************************************************************//**
 * @brief Run two checks on all individuals in the simulation
 * 
 * @details First, check if the individual is alive, HIV- and has HIV+ partners, 
 * and if so whether they are in the right place in the list of susceptibles 
 * in a serodiscordant partnership.\n
 * Secondly, check if the individual is alive and has available partnerships, 
 * whether they are in the right place in the list of available partners.
 * @param patch Pointer to a @ref patch_struct structure
 * @param overall_partnerships Pointer to an @ref all_partnerships structure
 * @param t0 Start time of the simulation
 * @param t_step Time step of the simulation
 ****************************************************************************/

void sweep_through_all_and_check_lists_serodiscordant_and_available_partners(
        patch_struct *patch, all_partnerships * overall_partnerships, int t0, int t_step){
    int p, k;
    individual temp_ind;

    for(p=0 ; p<NPATCHES; p++){
        for(k=0 ; k<patch[p].id_counter; k++){
            temp_ind = patch[p].individual_population[k];
            if(temp_ind.id==FOLLOW_INDIVIDUAL  && temp_ind.patch_no==FOLLOW_PATCH){
                print_individual(&temp_ind);
            }
            // CHECK 1: check if individual is alive, HIV- and has HIV+ partners, 
            // and if so whether he is in the right place in the list of susceptibles 
            // in a serodiscordant partnership
            check_if_individual_should_be_in_list_susceptibles_in_serodiscordant_partnership(&temp_ind, overall_partnerships);
            // CHECK 2: check if individual is alive and has available partnerships, 
            // whether he is in the right place in the list of available partners
            check_if_individual_should_be_in_list_available_partners(&temp_ind, overall_partnerships, t0, t_step);
        }

    }
}


/***************************************************************************//**
 * @brief Check partners in inside and outside patches.
 * 
 * @param patch Pointer to a @ref patch_struct structure
 * @param overall_partnerships Pointer to an @ref all_partnerships structure
 * @param t0 Start time of the simulation
 * @param t_step Time step of the simulation.
 ****************************************************************************/

void sweep_through_all_and_check_n_partners_outside_n_HIVpos_partners_and_n_HIVpos_partners_outside(
        patch_struct *patch, all_partnerships * overall_partnerships, int t0, int t_step){
    int p, k, i;
    individual temp_ind;
    int temp_patch;
    int check_n_partners_outside, check_n_HIVpos_partners, check_n_HIVpos_partners_outside;

    for(p=0 ; p<NPATCHES; p++){
        for(k=0 ; k<patch[p].id_counter; k++){
            temp_ind = patch[p].individual_population[k];
            temp_patch = temp_ind.patch_no;
            if(temp_ind.cd4 > DEAD){ // only do if person still alive
                if(temp_ind.id==FOLLOW_INDIVIDUAL  && temp_patch==FOLLOW_PATCH){
                    print_individual(&temp_ind);
                }

                // is n_partners_outside what it should be?
                check_n_HIVpos_partners = 0;
                check_n_partners_outside = 0;
                check_n_HIVpos_partners_outside = 0;
                if(temp_ind.n_partners>0){
                    for(i=0 ; i<temp_ind.n_partners ; i++){
                        if(temp_ind.partner_pairs[i]->ptr[1-temp_ind.gender]->HIV_status>0){
                            check_n_HIVpos_partners++;
                        }
                        if(temp_ind.partner_pairs[i]->ptr[1-temp_ind.gender]->patch_no != temp_patch){
                            check_n_partners_outside++;
                            if(temp_ind.partner_pairs[i]->ptr[1-temp_ind.gender]->HIV_status>0){
                                check_n_HIVpos_partners_outside++;
                            }
                        }
                    }
                }
                if(check_n_partners_outside != temp_ind.n_partners_outside){
                    printf("PROBLEM: individual %ld from patch %d does not have the correct number of partners outside the patch\n",temp_ind.id,temp_patch);
                    print_individual(&temp_ind);
                    if(temp_ind.n_partners>0){
                        printf("--- PARTNERS OF THE PROBLEMATIC INDIVIDUAL ---\n");
                        fflush(stdout);
                        for(i=0 ; i<temp_ind.n_partners ; i++){
                            print_individual(temp_ind.partner_pairs[i]->ptr[1-temp_ind.gender]);
                        }
                        printf("------\n");
                        fflush(stdout);
                    }
                    printf("Individual has %d partners outside his/her patch but n_partners_outside is %d \n",check_n_partners_outside,temp_ind.n_partners_outside);
                    printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                    fflush(stdout);
                    exit(1);
                }
                if(temp_ind.HIV_status==0){ // we only keep track of HIV positive partners for the HIV negative individuals
                    if(check_n_HIVpos_partners != temp_ind.n_HIVpos_partners){
                        printf("PROBLEM: individual %ld from patch %d does not have the correct number of HIV+ partners\n",temp_ind.id,temp_patch);
                        print_individual(&temp_ind);
                        if(temp_ind.n_partners>0){
                            printf("--- PARTNERS OF THE PROBLEMATIC INDIVIDUAL ---\n");
                            fflush(stdout);
                            for(i=0 ; i<temp_ind.n_partners ; i++){
                                print_individual(temp_ind.partner_pairs[i]->ptr[1-temp_ind.gender]);
                            }
                            printf("------\n");
                            fflush(stdout);
                        }
                        printf("Individual has %d HIV+ partners but n_HIVpos_partners is %d \n",check_n_HIVpos_partners,temp_ind.n_HIVpos_partners);
                        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                        fflush(stdout);
                        exit(1);
                    }
                    if(check_n_HIVpos_partners_outside != temp_ind.n_HIVpos_partners_outside){
                        printf("PROBLEM: individual %ld from patch %d does not have the correct number of HIV+ partners outside the patch\n",temp_ind.id,temp_patch);
                        print_individual(&temp_ind);
                        if(temp_ind.n_partners>0){
                            printf("--- PARTNERS OF THE PROBLEMATIC INDIVIDUAL ---\n");
                            fflush(stdout);
                            for(i=0 ; i<temp_ind.n_partners ; i++){
                                print_individual(temp_ind.partner_pairs[i]->ptr[1-temp_ind.gender]);
                            }
                            printf("------\n");
                            fflush(stdout);
                        }
                        printf("Individual has %d HIV+ partners outside his/her patch but n_HIVpos_partners_outside is %d \n",check_n_HIVpos_partners,temp_ind.n_HIVpos_partners_outside);
                        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                        fflush(stdout);
                        exit(1);
                    }
                }
            }
        }
    }
}


/***************************************************************************//**
 * @brief Check age and risk groups of partnerships
 * 
 * @param patch Pointer to a @ref patch_struct structure
 * @param overall_partnerships Pointer to an @ref all_partnerships structure
 * @param t0 Start time of the simulation
 * @param t_step Time step of the simulation
 * @param debug Pointer to a @ref debug_struct structure
 ****************************************************************************/

void sweep_through_all_and_check_age_and_risk_of_partners(patch_struct *patch, 
    all_partnerships * overall_partnerships, int t0, int t_step, debug_struct *debug){

    int p, k, i;
    individual *ind1, *ind2;
    int age_f, age_m, risk_f, risk_m;

    for(age_f=0 ; age_f<N_AGE ; age_f++){
        for(age_m=0 ; age_m<N_AGE ; age_m++){
            debug->age_of_partners_cross_sectional[t0 - patch[0].param->start_time_simul][age_f][age_m] = 0.0;
        }
    }
    for(risk_f=0 ; risk_f<N_RISK ; risk_f++){
        for(risk_m=0 ; risk_m<N_RISK ; risk_m++){
            debug->risk_of_partners_cross_sectional[t0 - patch[0].param->start_time_simul][risk_f][risk_m] = 0.0;
        }
    }
    for(p=0 ; p<NPATCHES; p++){
        for(k=0 ; k<patch[p].id_counter; k++){
            ind1 = &patch[p].individual_population[k];
            if(ind1->id == FOLLOW_INDIVIDUAL && ind1->patch_no == FOLLOW_PATCH){
                print_individual(ind1);
                fflush(stdout);
            }
            if(ind1->cd4 > DEAD && ind1->n_partners>0){
                for(i=0 ; i<ind1->n_partners ; i++){
                    ind2 = ind1->partner_pairs[i]->ptr[1-ind1->gender];
                    if(ind1->gender == FEMALE){
                        age_f = get_age_group(ind1->DoB,t0+t_step*TIME_STEP, AGE_GROUPS, N_AGE);
                        age_m = get_age_group(ind2->DoB,t0+t_step*TIME_STEP, AGE_GROUPS, N_AGE);
                        risk_f = ind1->sex_risk;
                        risk_m = ind2->sex_risk;
                    }else{
                        age_f = get_age_group(ind2->DoB,t0+t_step*TIME_STEP, AGE_GROUPS, N_AGE);
                        age_m = get_age_group(ind1->DoB,t0+t_step*TIME_STEP, AGE_GROUPS, N_AGE);
                        risk_f = ind2->sex_risk;
                        risk_m = ind1->sex_risk;
                    }
                    debug->age_of_partners_cross_sectional[t0 - patch[0].param->start_time_simul][age_f][age_m] += 0.5; // will count each partnership twice otherwise
                    debug->risk_of_partners_cross_sectional[t0 - patch[0].param->start_time_simul][risk_f][risk_m] += 0.5; // will count each partnership twice otherwise
                }
            }
        }
    }
}


/***************************************************************************//**
 * @brief Create blank debug files, including headers, at beginning of simulation
 * 
 * @details Several checks are created in these debug files, including: number 
 * of HIV-related deaths, duration of HIV without ART by SPVL, duration of HIV
 * with ART, estimation of R0 and/or doubling time. \n
 * 
 * This function blanks the debugging files at the beginning of the run. 
 * The reason we need to do this is because debugging file are _inefficient_
 * as we write to them continuously using the "a" (ie append) write option so
 * that we don't have to keep a lot of potentially large arrays that we do not
 * need for non-debugging runs.  Essentially i/o (reading/writing from/to disk)
 * is inefficient so we want to avoid it as much as possible, but for debugging
 * files we make an exception.  For files which save output and only write 
 * occasionally (e.g. `Annual_output*.csv` files) we don't need to blank the files
 * at the beginning of a run.\n
 * From http://www.statssa.gov.za/?p=2973 (accessed 24 June 2016):\n
 * ....the number of AIDS related deaths is estimated to have decreased from 
 * 363 910 deaths in 2005 (51% of all deaths) to 171 733 deaths in 2014 
 * (31% of all deaths). This can be associated with the increased rollout of
 * antiretroviral therapy (ART).
 * 
 * @param file_data_store Pointer to a @ref file_struct structure
 ****************************************************************************/

void blank_debugging_files(file_struct *file_data_store){
    char age_group_string[1000];
    int p;
    for(p=0;p<NPATCHES;p++){
        if(WRITE_DEBUG_HIV_DURATION==1){
            file_data_store->HIVDURATIONFILE[p] = fopen(file_data_store->filename_debughivduration[p],"w");
            fprintf(file_data_store->HIVDURATIONFILE[p],"time,patch,id,Time_HIV+,Time_sc,ART_Status,SPVL_cat,CD4\n");
            fclose(file_data_store->HIVDURATIONFILE[p]);
        }
        if (WRITE_DEBUG_HIV_DURATION_KM==1){
            file_data_store->HIVDURATIONFILE_KM[p] = fopen(file_data_store->filename_debughivduration_km[p],"w");
            fprintf(file_data_store->HIVDURATIONFILE_KM[p],"time,time_pos,reason,gender,CD4,SPVLnum,SPVLcat\n");
            fclose(file_data_store->HIVDURATIONFILE_KM[p]);
        }
        if (WRITE_DEBUG_CD4_AFTER_SEROCONVERSION==1){
            file_data_store->HIVCD4_AFTER_SEROCONVERSION[p] = fopen(file_data_store->filename_debughivcd4_after_seroconversion[p],"w");
            fprintf(file_data_store->HIVCD4_AFTER_SEROCONVERSION[p],"id,CD4,SPVLcat\n");
            fclose(file_data_store->HIVCD4_AFTER_SEROCONVERSION[p]);
        }
        if (WRITE_DEBUG_INITIAL_SPVL_DISTRIBUTION==1){
            file_data_store->HIV_INITIAL_SPVL_DISTRIBUTION[p] = fopen(file_data_store->filename_debuginitial_spvl_distribution[p],"w");
            fprintf(file_data_store->HIV_INITIAL_SPVL_DISTRIBUTION[p],"SPVL_cat,SPVL_E,SPVL_G\n");
            fclose(file_data_store->HIV_INITIAL_SPVL_DISTRIBUTION[p]);
        }
        if (WRITE_DEBUG_HIV_STATES==1){
            file_data_store->HIVSTATEPOPULATIONFILE[p] = fopen(file_data_store->filename_debug_hivpopulation[p],"w");
            fprintf(file_data_store->HIVSTATEPOPULATIONFILE[p],"time,cd4_status,spvl,art_status,cumulative_t_earlyART,cumulative_t_ARTVS,cumulative_t_ARTVU,npartners\n");
            fclose(file_data_store->HIVSTATEPOPULATIONFILE[p]);
        }
        /* If WRITE_DEBUG_DEMOGRAPHICS_AGE_DISTRIBUTION_BY_GENDER==1 then we print age distribution at some specified times, so make sure that the file is initially blank. */
        if (WRITE_DEBUG_DEMOGRAPHICS_AGE_DISTRIBUTION_BY_GENDER==1){
            /* The "1000" refers to the size of the array age_group_string declared above. */
            generate_demographics_byage_gender_file_header(age_group_string, 1000);
            file_data_store->AGEDISTRIBUTIONFILE[p] = fopen(file_data_store->filename_debug_agedistribution[p],"w");
            fprintf(file_data_store->AGEDISTRIBUTIONFILE[p],"%s",age_group_string);
            fclose(file_data_store->AGEDISTRIBUTIONFILE[p]);
        }
        if(WRITE_DEBUG_ART_STATE==1){
            file_data_store->ARTPOPULATIONFILE[p] = fopen(file_data_store->filename_debug_artpopulation[p],"w");
            fprintf(file_data_store->ARTPOPULATIONFILE[p],"time,n_hivneg,n_hivpos_dontknowstatus,n_hivpos_knowposneverart,n_hivpos_earlyart,n_hivpos_artvs,n_hivpos_artvu,n_hivpos_dropout,n_hivpos_cascadedropout,n_artdeath,");
            fprintf(file_data_store->ARTPOPULATIONFILE[p],"cumulative_n_start_emergency_art_fromuntested,cumulative_n_start_emergency_art_fromartnaive,cumulative_n_start_emergency_art_fromartdroupout,cumulative_n_start_emergency_art_fromcascadedropout,");
            fprintf(file_data_store->ARTPOPULATIONFILE[p],"cumulative_n_learnhivpos_fromuntested,cumulative_n_startART_fromuntested,cumulative_n_startART_fromartnaive,cumulative_n_startART_fromartdropout,cumulative_n_startART_fromcascadedropout,cumulative_n_becomeVS_fromearlyart,cumulative_n_becomeVS_fromartvu,cumulative_n_becomeVU_fromearlyart,cumulative_n_becomeVU_fromartvs,cumulative_n_ARTdropout_fromearlyart,cumulative_n_ARTdropout_fromartvs,cumulative_n_ARTdropout_fromartvu,cumulative_n_cascadedropout_fromARTnaive,n_cascadedropout_fromARTneg,cumulative_n_aidsdeaths_fromuntested,cumulative_n_aidsdeaths_fromartnaive,cumulative_n_aidsdeaths_fromearlyart,cumulative_n_aidsdeaths_fromartvs,cumulative_n_aidsdeaths_fromartvu,cumulative_n_aidsdeaths_fromartdropout,cumulative_n_aidsdeaths_fromcascadedropout\n");
            fclose(file_data_store->ARTPOPULATIONFILE[p]);
        }
    }
    if (WRITE_DEBUG_DEMOGRAPHICS_NBIRTHS_NEWADULTS_DEATHS==1){
        file_data_store->NBIRTHS_NNEWADULTS_NDEATHS_FILE = fopen(file_data_store->filename_debug_nnewadults_ndeaths_file,"w");
        fprintf(file_data_store->NBIRTHS_NNEWADULTS_NDEATHS_FILE,"t,");
        for(p=0;p<NPATCHES;p++)
            fprintf(file_data_store->NBIRTHS_NNEWADULTS_NDEATHS_FILE,"NBirthsPatch%i,NNewAdultsPatch%i,NDeaths%i,",p,p,p);
        fprintf(file_data_store->NBIRTHS_NNEWADULTS_NDEATHS_FILE,"\n");
        fclose(file_data_store->NBIRTHS_NNEWADULTS_NDEATHS_FILE);
    }

    /* This header is quite complicated so put in a separate function. */
    if (WRITE_DEBUG_DEMOGRAPHICS_AGE_DISTRIBUTION_ONEYEARINCKIDS==1){
        blank_one_year_age_groups_including_kids(file_data_store);
    }

    if(WRITE_PHYLOGENETICS_OUTPUT==1)
        blank_phylo_transmission_data_file(file_data_store);

    if(DEBUG_PARTNERSHIP_DURATION ==1){
        file_data_store->DUR_BETWEEN_HIGHHIGH = fopen(file_data_store->filename_DUR_BETWEEN_HIGHHIGH,"w");
        if (file_data_store->DUR_BETWEEN_HIGHHIGH==NULL){
            printf("Cannot open DUR_BETWEEN_HIGHHIGH\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        fclose(file_data_store->DUR_BETWEEN_HIGHHIGH);

        file_data_store->DUR_BETWEEN_MEDMED = fopen(file_data_store->filename_DUR_BETWEEN_MEDMED,"w");
        if (file_data_store->DUR_BETWEEN_MEDMED==NULL){
            printf("Cannot open DUR_BETWEEN_MEDMED\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        fclose(file_data_store->DUR_BETWEEN_MEDMED);

        file_data_store->DUR_BETWEEN_LOWLOW = fopen(file_data_store->filename_DUR_BETWEEN_LOWLOW,"w");
        if (file_data_store->DUR_BETWEEN_LOWLOW==NULL){
            printf("Cannot open DUR_BETWEEN_LOWLOW\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        fclose(file_data_store->DUR_BETWEEN_LOWLOW);

        file_data_store->DUR_WITHIN_HIGHHIGH = fopen(file_data_store->filename_DUR_WITHIN_HIGHHIGH,"w");
        if (file_data_store->DUR_WITHIN_HIGHHIGH==NULL){
            printf("Cannot open DUR_WITHIN_HIGHHIGH\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        fclose(file_data_store->DUR_WITHIN_HIGHHIGH);

        file_data_store->DUR_WITHIN_MEDMED = fopen(file_data_store->filename_DUR_WITHIN_MEDMED,"w");
        if (file_data_store->DUR_WITHIN_MEDMED==NULL){
            printf("Cannot open DUR_WITHIN_MEDMED\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        fclose(file_data_store->DUR_WITHIN_MEDMED);

        file_data_store->DUR_WITHIN_LOWLOW = fopen(file_data_store->filename_DUR_WITHIN_LOWLOW,"w");
        if (file_data_store->DUR_WITHIN_LOWLOW==NULL){
            printf("Cannot open DUR_WITHIN_LOWLOW\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        fclose(file_data_store->DUR_WITHIN_LOWLOW);
    }

    if(CHECK_AGE_AND_RISK_ASSORTATIVITY ==1){
        file_data_store->age_assortativity = fopen(file_data_store->filename_age_assortativity,"w");
        if (file_data_store->age_assortativity==NULL){
            printf("Cannot open age_assortativity\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        fclose(file_data_store->age_assortativity);

        file_data_store->age_assortativity_cross_sectional = fopen(file_data_store->filename_age_assortativity_cross_sectional,"w");
        if (file_data_store->age_assortativity_cross_sectional==NULL){
            printf("Cannot open age_assortativity_cross_sectional\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        fclose(file_data_store->age_assortativity_cross_sectional);

        file_data_store->risk_assortativity = fopen(file_data_store->filename_risk_assortativity,"w");
        if (file_data_store->risk_assortativity==NULL){
            printf("Cannot open risk_assortativity\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        fclose(file_data_store->risk_assortativity);

        file_data_store->risk_assortativity_cross_sectional = fopen(file_data_store->filename_risk_assortativity_cross_sectional,"w");
        if (file_data_store->risk_assortativity_cross_sectional==NULL){
            printf("Cannot open risk_assortativity_cross_sectional\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        fclose(file_data_store->risk_assortativity_cross_sectional);
    }
    if (WRITE_HAZARDS==1)
        blank_hazard_file(file_data_store);
}


/***************************************************************************//**
 * @brief Print how long each person who is HIV-positive and dies of AIDS-related
 *  illness is alive for.
 * @details This function allows us to check that HIV duration has the right
 * distribution.  Can also be subset to ART-naive.
 * 
 * @param dead_person Pointer to an @ref individual structure
 * @param t Current time
 * @param file_data_store Pointer to a @ref file_struct structure of file names
 ****************************************************************************/

void write_hiv_duration(individual *dead_person, double t, file_struct *file_data_store){
    int p = dead_person->patch_no;
    file_data_store->HIVDURATIONFILE[p] = fopen(file_data_store->filename_debughivduration[p],"a");
    if (file_data_store->HIVDURATIONFILE[p]==NULL){
        printf("Cannot open output file in write_hiv_duration().\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    fprintf(file_data_store->HIVDURATIONFILE[p],
        "%6.4lf,%d,%ld,%8.6f,%6.4lf,%d,%d,%d\n",
        t,
        dead_person->patch_no,dead_person->id,
        dead_person->DEBUGTOTALTIMEHIVPOS,
        dead_person->t_sc,
        dead_person->ART_status,
        dead_person->SPVL_cat,
        dead_person->cd4);
    fclose(file_data_store->HIVDURATIONFILE[p]);
}


/***************************************************************************//**
 * @brief Write information related to HIV-positive individuals to file
 * 
 * @details The final argument is reason for being removed from survival cohort.
 * 1="AIDS death", 2="death from natural causes", 3="start ART".
 * Note we don't bother with the end of the simulation for now.
 * 
 * @param indiv Pointer to an @ref individual structure
 * @param t Current time
 * @param file_data_store Pointer to a @ref file_struct structure of file names
 * @param reason Reason for being removed from survival cohort
 ****************************************************************************/

void write_hiv_duration_km(individual *indiv, double t, file_struct *file_data_store, int reason){
    int p = indiv->patch_no;
    file_data_store->HIVDURATIONFILE_KM[p] = fopen(file_data_store->filename_debughivduration_km[p],"a");
    if (file_data_store->HIVDURATIONFILE_KM[p]==NULL){
        printf("Cannot open output file in write_hiv_duration_km().\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    fprintf(file_data_store->HIVDURATIONFILE_KM[p],
        "%6.4lf,%6.4lf,%d,%d,%d,%6.4lf,%d\n",
        t,
        t-indiv->t_sc,
        reason,
        indiv->gender,
        indiv->cd4,
        indiv->SPVL_num_E+indiv->SPVL_num_G,
        indiv->SPVL_cat);
    fclose(file_data_store->HIVDURATIONFILE_KM[p]);
}


/***************************************************************************//**
 * @brief Write information related to HIV-positive individuals to file
 * 
 * @param patch Pointer to a @ref patch_struct structure
 * @param t Current time
 * @param file_data_store Pointer to a @ref file_struct structure of file names
 ****************************************************************************/

/* The "4" is the "reason" (or censoring status) - 4 means reaching the end 
of the simulation before dying/starting ART. */
void write_hiv_duration_km_end_of_simulation(patch_struct *patch, double t, file_struct *file_data_store){
    int p;
    long n_id;
    for(p=0; p<NPATCHES; p++){
        printf("Opening file %s\n",file_data_store->filename_debughivduration_km[p]);
        file_data_store->HIVDURATIONFILE_KM[p] = fopen(file_data_store->filename_debughivduration_km[p],"a");
        if (file_data_store->HIVDURATIONFILE_KM[p]==NULL){
            printf("Cannot open output file in write_hiv_duration_km().\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        for(n_id=0; n_id<patch[p].id_counter; n_id++){
            if((patch[p].individual_population[n_id].cd4>DEAD) && (patch[p].individual_population[n_id].HIV_status>UNINFECTED) && (patch[p].individual_population[n_id].ART_status==ARTNAIVE || patch[p].individual_population[n_id].ART_status==ARTNEG)){
                fprintf(file_data_store->HIVDURATIONFILE_KM[p],
                    "%6.4lf,%6.4lf,4,%d,%d,%6.4lf,%d\n",
                    t,
                    t-patch[p].individual_population[n_id].t_sc,
                    patch[p].individual_population[n_id].gender,
                    patch[p].individual_population[n_id].cd4,
                    patch[p].individual_population[n_id].SPVL_num_E+patch[p].individual_population[n_id].SPVL_num_G,
                    patch[p].individual_population[n_id].SPVL_cat);
            }
        }
        fclose(file_data_store->HIVDURATIONFILE_KM[p]);
    }
}


/***************************************************************************//**
 * @brief Write CD4 count at seroconversion to file
 * 
 * @details Called in @ref hiv.c.  Nothing returned, just writes to file.
 * 
 * @param indiv Pointer to an @ref individual structure
 * @param file_data_store Pointer to a @ref file_struct structure of file
 ****************************************************************************/

void write_cd4_at_seroconversion(individual *indiv, file_struct *file_data_store){
    int p = indiv->patch_no;
    file_data_store->HIVCD4_AFTER_SEROCONVERSION[p] = fopen(file_data_store->filename_debughivcd4_after_seroconversion[p],"a");
    if (file_data_store->HIVCD4_AFTER_SEROCONVERSION[p]==NULL){
        printf("Cannot open output file in write_hiv_duration_km().\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    fprintf(file_data_store->HIVCD4_AFTER_SEROCONVERSION[p],"%li,%i,%i\n",indiv->id,indiv->cd4,indiv->SPVL_cat);
    fclose(file_data_store->HIVCD4_AFTER_SEROCONVERSION[p]);
}


/***************************************************************************//**
 * @brief Write initial SPVL distribution to file (i.e. at seeding of HIV)
 * 
 * @details Called in @ref hiv.c.
 * 
 * @param seeded_infection pointer to a @ref individual structure of the 
 * individuals who seed the HIV epidemic
 * @param file_data_store pointer to a @ref file_struct structure of file
 * names/locations
 ****************************************************************************/

void write_initial_spvl_distribution(individual *seeded_infection, file_struct *file_data_store){
    int p = seeded_infection->patch_no;
    file_data_store->HIV_INITIAL_SPVL_DISTRIBUTION[p] = fopen(file_data_store->filename_debuginitial_spvl_distribution[p],"a");
    if (file_data_store->HIV_INITIAL_SPVL_DISTRIBUTION[p]==NULL){
        printf("Cannot open output file in write_hiv_duration_km().\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    fprintf(file_data_store->HIV_INITIAL_SPVL_DISTRIBUTION[p],"%i,%6.4lf,%6.4lf\n",seeded_infection->SPVL_cat,seeded_infection->SPVL_num_E,seeded_infection->SPVL_num_G);
    fclose(file_data_store->HIV_INITIAL_SPVL_DISTRIBUTION[p]);
}


/***************************************************************************//**
 * @brief Write counts of CD4 and SPVL states to file
 * 
 * @param patch pointer to a @ref patch_struct structure
 * @param year Current year
 * @param file_data_store pointer to a @ref file_struct structure of file 
 * locations/names
 ****************************************************************************/

void write_cd4_spvl_states(patch_struct *patch, int year, file_struct *file_data_store){
    int p;
    long n_id;
    int art_status, cd4_status, npartners;
    double spvl, t_early_art, t_vs, t_vu;
    for(p=0; p<NPATCHES; p++){
        file_data_store->HIVSTATEPOPULATIONFILE[p] = fopen(file_data_store->filename_debug_hivpopulation[p],"a");
        for(n_id=0; n_id<patch[p].id_counter; n_id++){
            /* Only get data for HIV+ who are alive. */
            if((patch[p].individual_population[n_id].cd4>DEAD) && (patch[p].individual_population[n_id].HIV_status>UNINFECTED)){
                art_status = patch[p].individual_population[n_id].ART_status;
                cd4_status = patch[p].individual_population[n_id].cd4;
                spvl = patch[p].individual_population[n_id].SPVL_num_E+patch[p].individual_population[n_id].SPVL_num_G;
                npartners = patch[p].individual_population[n_id].n_partners;
                t_early_art = patch[p].individual_population[n_id].DEBUG_cumulative_time_on_ART_early;
                t_vs = patch[p].individual_population[n_id].DEBUG_cumulative_time_on_ART_VS;
                t_vu = patch[p].individual_population[n_id].DEBUG_cumulative_time_on_ART_VU;
                fprintf(file_data_store->HIVSTATEPOPULATIONFILE[p],"%d,%d,%6.4lf,%d,%6.4lf,%6.4lf,%6.4lf,%i\n",year,cd4_status,spvl,art_status,t_early_art,t_vs,t_vu,npartners);
            }
        }
        fclose(file_data_store->HIVSTATEPOPULATIONFILE[p]);
    }
}


/***************************************************************************//**
 * @brief Helper function for printing when something should be zero but isn't.
 * 
 * @param varname Name of variable that should be zero.
 ****************************************************************************/

void print_debugerror_shouldbezero_exit(char *varname){
    printf("ERROR: Variable %s should be zero. Exiting\n",varname);
    printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
    fflush(stdout);
    exit(1);
}


/***************************************************************************//**
 * @brief Write counts of individuals with different ART statuses to file
 * 
 * @param patch pointer to a @ref patch_struct structure
 * @param year Current year
 * @param debug pointer to a @ref debug_struct structure
 * @param file_data_store pointer to a @ref file_struct structure of file
 * locations.
 ****************************************************************************/

void write_art_states(patch_struct *patch, int year, debug_struct *debug, file_struct *file_data_store){
    int p;
    long n_id;
    int art_i;
    int art_status;

    long n_hivneg;
    long n_hivpos_dontknowstatus;
    long n_hivpos_knowposneverart;
    long n_hivpos_earlyart;
    long n_hivpos_artvs;
    long n_hivpos_artvu;
    long n_hivpos_dropout;
    long n_hivpos_cascadedropout;
    long n_artdeath;

    /****** Now look at transitions. ******/
    /* First AIDS deaths - should ONLY be possible if not ART VS. */
    long n_aidsdeaths_fromuntested,n_aidsdeaths_fromartnaive,n_aidsdeaths_fromearlyart,n_aidsdeaths_fromartvs,n_aidsdeaths_fromartvu,n_aidsdeaths_fromartdropout,n_aidsdeaths_fromcascadedropout;
    /* Next number of HIV+ who learn status (go from n_hivpos_dontknowstatus to n_hivpos_knowposneverart). */
    long n_learnhivpos_fromuntested;
    /* Number of people who get CD4 tested. */
    //long n_getcd4test_fromuntested,n_getcd4test_fromartnaive,n_getcd4test_fromartdropout,n_getcd4test_fromcascadedropout;
    /* Number of people who get ART tested: */
    long n_startART_fromuntested,n_startART_fromartnaive,n_startART_fromartdropout,n_startART_fromcascadedropout;
    /* Number of people becoming virally suppressed/virally unsuppressed: */
    long n_becomeVS_fromearlyart, n_becomeVS_fromartvu,n_becomeVU_fromearlyart, n_becomeVU_fromartvs;
    /* Number of people stopping ART: */
    long n_ARTdropout_fromearlyart,n_ARTdropout_fromartvs,n_ARTdropout_fromartvu;
    /* Number of people dropping out of cascade before ART: */
    long n_cascadedropout_fromARTnaive, n_cascadedropout_fromARTneg;

    /* ART_status runs from -1 to 6 (ARTDEATH) so need array to go from 0 to ARTDEATH+2. */
    long temp_state_counter[ARTDEATH+2];
    long temp_hivnegstate_counter;

    for (p=0;p<NPATCHES;p++){
        for (art_i=0;art_i<(ARTDEATH+2);art_i++)
            temp_state_counter[art_i] = 0;
        temp_hivnegstate_counter = 0;
        for (n_id=0;n_id<patch[p].id_counter;n_id++){
            if (patch[p].individual_population[n_id].cd4>DEAD){
                art_status = patch[p].individual_population[n_id].ART_status;
                if (art_status<ARTNEG || art_status>=ARTDEATH){
                    printf("ERROR - undefined ART status for individual %ld in patch %d. Exiting\n",patch[p].individual_population[n_id].id,p);
                    printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                    fflush(stdout);
                    exit(1);
                }

                temp_state_counter[art_status+1]++;
                if(patch[p].individual_population[n_id].HIV_status==UNINFECTED)
                    temp_hivnegstate_counter++;
            }
            else{ /* Check if dead person died while on ART: */
                if(patch[p].individual_population[n_id].ART_status==ARTDEATH)
                    temp_state_counter[ARTDEATH+1]++;
            }
        }
        /* Translate these into more readable variable names for clarity - this is not needed for coded. */
        n_hivneg = temp_hivnegstate_counter;
        n_hivpos_dontknowstatus = temp_state_counter[ARTNEG+1]-temp_hivnegstate_counter;
        n_hivpos_knowposneverart = temp_state_counter[ARTNAIVE+1];
        n_hivpos_earlyart = temp_state_counter[EARLYART+1];
        n_hivpos_artvs = temp_state_counter[LTART_VS+1];
        n_hivpos_artvu = temp_state_counter[LTART_VU+1];
        n_hivpos_dropout = temp_state_counter[ARTDROPOUT+1];
        n_hivpos_cascadedropout = temp_state_counter[CASCADEDROPOUT+1];
        n_artdeath = temp_state_counter[ARTDEATH+1];

        /* Now look at transitions between states. Again, translate first for clarity: */
        if (debug->art_vars[p].n_start_emergency_art != (debug->art_vars[p].n_start_emergency_art_fromuntested + debug->art_vars[p].n_start_emergency_art_fromartnaive+ debug->art_vars[p].n_start_emergency_art_fromartdroupout+ debug->art_vars[p].n_start_emergency_art_fromcascadedropout)){
            printf("Error - number starting emergency ART in a given timestep doesn't add up. Exiting\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        n_learnhivpos_fromuntested = debug->art_vars[p].cascade_transitions[ARTNEG+1][ARTNAIVE+1];
        if (debug->art_vars[p].cascade_transitions[ARTNAIVE+1][ARTNEG+1]>0) /* Impossible */
            print_debugerror_shouldbezero_exit("debug->art_vars[p].cascade_transitions[ARTNAIVE+1][ARTNEG+1]");
        for (art_i=EARLYART;art_i<=ARTDEATH;art_i++){
            if(debug->art_vars[p].cascade_transitions[art_i+1][ARTNEG+1]!=0){
                char errorstring[50];
                sprintf(errorstring,"art_vars[p].cascade_transitions[%i][ARTNEG+1] %li",art_i,debug->art_vars[p].cascade_transitions[art_i+1][ARTNEG+1]);
                print_debugerror_shouldbezero_exit(errorstring);
            }
        }

        n_startART_fromuntested = debug->art_vars[p].cascade_transitions[ARTNEG+1][EARLYART+1];
        n_startART_fromartnaive = debug->art_vars[p].cascade_transitions[ARTNAIVE+1][EARLYART+1];
        n_startART_fromartdropout = debug->art_vars[p].cascade_transitions[ARTDROPOUT+1][EARLYART+1];
        n_startART_fromcascadedropout = debug->art_vars[p].cascade_transitions[CASCADEDROPOUT+1][EARLYART+1];
        /* Check for impossible transitions: */
        for (art_i=EARLYART;art_i<ARTDROPOUT;art_i++){
            if(debug->art_vars[p].cascade_transitions[art_i+1][EARLYART+1]!=0){
                char errorstring[50];
                sprintf(errorstring,"art_vars[p].cascade_transitions[%i][EARLYART+1]",art_i);
                print_debugerror_shouldbezero_exit(errorstring);
            }
        }
        if(debug->art_vars[p].cascade_transitions[ARTDEATH+1][EARLYART+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTDEATH+1][EARLYART+1]");

        // Now becoming VS - should only be possible from early ART or VU:
        n_becomeVS_fromearlyart = debug->art_vars[p].cascade_transitions[EARLYART+1][LTART_VS+1];
        n_becomeVS_fromartvu = debug->art_vars[p].cascade_transitions[LTART_VU+1][LTART_VS+1];
        if(debug->art_vars[p].cascade_transitions[ARTNEG+1][LTART_VS+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTNEG+1][LTART_VS+1]");
        if(debug->art_vars[p].cascade_transitions[ARTNAIVE+1][LTART_VS+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTNAIVE+1][LTART_VS+1]");
        if(debug->art_vars[p].cascade_transitions[LTART_VS+1][LTART_VS+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[LTART_VS+1][LTART_VS+1]");
        if(debug->art_vars[p].cascade_transitions[ARTDROPOUT+1][LTART_VS+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTDROPOUT+1][LTART_VS+1]");
        if(debug->art_vars[p].cascade_transitions[CASCADEDROPOUT+1][LTART_VS+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[CASCADEDROPOUT+1][LTART_VS+1]");
        if(debug->art_vars[p].cascade_transitions[ARTDEATH+1][LTART_VS+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTDEATH+1][LTART_VS+1]");

        // Now becoming VU - should only be possible from early ART or VS:
        n_becomeVU_fromearlyart = debug->art_vars[p].cascade_transitions[EARLYART+1][LTART_VU+1];
        n_becomeVU_fromartvs = debug->art_vars[p].cascade_transitions[LTART_VS+1][LTART_VU+1];
        if(debug->art_vars[p].cascade_transitions[ARTNEG+1][LTART_VU+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTNEG+1][LTART_VU+1]");
        if(debug->art_vars[p].cascade_transitions[ARTNAIVE+1][LTART_VU+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTNAIVE+1][LTART_VU+1]");
        if(debug->art_vars[p].cascade_transitions[LTART_VU+1][LTART_VU+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[LTART_VU+1][LTART_VU+1]");
        if(debug->art_vars[p].cascade_transitions[ARTDROPOUT+1][LTART_VU+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTDROPOUT+1][LTART_VU+1]");
        if(debug->art_vars[p].cascade_transitions[CASCADEDROPOUT+1][LTART_VU+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[CASCADEDROPOUT+1][LTART_VU+1]");
        if(debug->art_vars[p].cascade_transitions[ARTDEATH+1][LTART_VU+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTDEATH+1][LTART_VU+1]");

        // Now dropping out from ART - should only be possible from early ART, VU or VS:
        n_ARTdropout_fromearlyart = debug->art_vars[p].cascade_transitions[EARLYART+1][ARTDROPOUT+1];
        n_ARTdropout_fromartvs = debug->art_vars[p].cascade_transitions[LTART_VS+1][ARTDROPOUT+1];
        n_ARTdropout_fromartvu = debug->art_vars[p].cascade_transitions[LTART_VU+1][ARTDROPOUT+1];
        if(debug->art_vars[p].cascade_transitions[ARTNEG+1][ARTDROPOUT+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTNEG+1][ARTDROPOUT+1]");
        if(debug->art_vars[p].cascade_transitions[ARTNAIVE+1][ARTDROPOUT+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTNAIVE+1][ARTDROPOUT+1]");
        if(debug->art_vars[p].cascade_transitions[ARTDROPOUT+1][ARTDROPOUT+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTDROPOUT+1][ARTDROPOUT+1]");
        if(debug->art_vars[p].cascade_transitions[CASCADEDROPOUT+1][ARTDROPOUT+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[CASCADEDROPOUT+1][ARTDROPOUT+1]");
        if(debug->art_vars[p].cascade_transitions[ARTDEATH+1][ARTDROPOUT+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTDEATH+1][ARTDROPOUT+1]");

        // Now dropping out from cascade before ART - should only be possible from ARTNAIVE or ARTNEG:
        n_cascadedropout_fromARTnaive = debug->art_vars[p].cascade_transitions[ARTNAIVE+1][CASCADEDROPOUT+1];
        n_cascadedropout_fromARTneg = debug->art_vars[p].cascade_transitions[ARTNEG+1][CASCADEDROPOUT+1];
        //if(debug->art_vars[p].cascade_transitions[ARTNEG+1][CASCADEDROPOUT+1]!=0)
        //  print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTNEG+1][CASCADEDROPOUT+1]");
        if(debug->art_vars[p].cascade_transitions[EARLYART+1][CASCADEDROPOUT+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[EARLYART+1][CASCADEDROPOUT+1]");
        if(debug->art_vars[p].cascade_transitions[LTART_VS+1][CASCADEDROPOUT+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[LTART_VS+1][CASCADEDROPOUT+1]");
        if(debug->art_vars[p].cascade_transitions[LTART_VU+1][CASCADEDROPOUT+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[LTART_VU+1][CASCADEDROPOUT+1]");
        if(debug->art_vars[p].cascade_transitions[ARTDROPOUT+1][CASCADEDROPOUT+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTDROPOUT+1][CASCADEDROPOUT+1]");
        if(debug->art_vars[p].cascade_transitions[CASCADEDROPOUT+1][CASCADEDROPOUT+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[CASCADEDROPOUT+1][CASCADEDROPOUT+1]");
        if(debug->art_vars[p].cascade_transitions[ARTDEATH+1][CASCADEDROPOUT+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTDEATH+1][CASCADEDROPOUT+1]");

        n_aidsdeaths_fromuntested = debug->art_vars[p].cascade_transitions[ARTNAIVE+1][ARTDEATH+1];
        n_aidsdeaths_fromartnaive = debug->art_vars[p].cascade_transitions[ARTNEG+1][ARTDEATH+1];
        n_aidsdeaths_fromearlyart = debug->art_vars[p].cascade_transitions[EARLYART+1][ARTDEATH+1];
        n_aidsdeaths_fromartvs = debug->art_vars[p].cascade_transitions[LTART_VS+1][ARTDEATH+1];
        n_aidsdeaths_fromartvu = debug->art_vars[p].cascade_transitions[LTART_VU+1][ARTDEATH+1];
        n_aidsdeaths_fromcascadedropout = debug->art_vars[p].cascade_transitions[CASCADEDROPOUT+1][ARTDEATH+1];
        n_aidsdeaths_fromartdropout = debug->art_vars[p].cascade_transitions[ARTDROPOUT+1][ARTDEATH+1];
        if(debug->art_vars[p].cascade_transitions[ARTDEATH+1][ARTDEATH+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTDEATH+1][ARTDEATH+1]");

        file_data_store->ARTPOPULATIONFILE[p] = fopen(file_data_store->filename_debug_artpopulation[p],"a");
        fprintf(file_data_store->ARTPOPULATIONFILE[p],"%d,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,",year,n_hivneg,n_hivpos_dontknowstatus,n_hivpos_knowposneverart,n_hivpos_earlyart,n_hivpos_artvs,n_hivpos_artvu,n_hivpos_dropout,n_hivpos_cascadedropout,n_artdeath);
        fprintf(file_data_store->ARTPOPULATIONFILE[p],"%ld,%ld,%ld,%ld,",debug->art_vars[p].n_start_emergency_art_fromuntested , debug->art_vars[p].n_start_emergency_art_fromartnaive, debug->art_vars[p].n_start_emergency_art_fromartdroupout, debug->art_vars[p].n_start_emergency_art_fromcascadedropout);
        //fprintf(file_data_store->ARTPOPULATIONFILE[p],"%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld\n",n_learnhivpos_fromuntested,n_getcd4test_fromuntested,n_getcd4test_fromartnaive,n_getcd4test_fromartdropout,n_getcd4test_fromcascadedropout,n_startART_fromuntested,n_startART_fromartnaive,n_startART_fromartdropout,n_startART_fromcascadedropout,n_becomeVS_fromearlyart, n_becomeVS_fromartvu,n_becomeVU_fromearlyart, n_becomeVU_fromartvs,n_ARTdropout_fromearlyart,n_ARTdropout_fromartvs,n_ARTdropout_fromartvu,n_cascadedropout_fromARTnaive,n_aidsdeaths_fromuntested,n_aidsdeaths_fromartnaive,n_aidsdeaths_fromearlyart,n_aidsdeaths_fromartvs,n_aidsdeaths_fromartvu,n_aidsdeaths_fromartdropout,n_aidsdeaths_fromcascadedropout);
        fprintf(file_data_store->ARTPOPULATIONFILE[p],"%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld\n",n_learnhivpos_fromuntested,n_startART_fromuntested,n_startART_fromartnaive,n_startART_fromartdropout,n_startART_fromcascadedropout,n_becomeVS_fromearlyart, n_becomeVS_fromartvu,n_becomeVU_fromearlyart, n_becomeVU_fromartvs,n_ARTdropout_fromearlyart,n_ARTdropout_fromartvs,n_ARTdropout_fromartvu,n_cascadedropout_fromARTnaive,n_cascadedropout_fromARTneg, n_aidsdeaths_fromuntested,n_aidsdeaths_fromartnaive,n_aidsdeaths_fromearlyart,n_aidsdeaths_fromartvs,n_aidsdeaths_fromartvu,n_aidsdeaths_fromartdropout,n_aidsdeaths_fromcascadedropout);
        fclose(file_data_store->ARTPOPULATIONFILE[p]);
    }
}


/***************************************************************************//**
 * @brief Reset counters on annual CHiPs visits
 * 
 * @details Used in @ref simul.c.  Nothing is returned, only counters reset.
 * 
 * @param age_list Pointer to @ref age_list_struct structure.
 ****************************************************************************/

void reset_annual_chips_visit_counter(age_list_struct *age_list){
    int g,ai,n;
    for(g=0; g<N_GENDER; g++){
        for(ai=0; ai<(MAX_AGE-AGE_ADULT); ai++)
            for(n=0;n<age_list->age_list_by_gender[g]->number_per_age_group[ai];n++)
                age_list->age_list_by_gender[g]->age_group[ai][n]->VISITED_BY_CHIPS_THISROUND = FALSE;
        for(n=0; n<age_list->age_list_by_gender[g]->number_oldest_age_group; n++)
            age_list->age_list_by_gender[g]->oldest_age_group[n]->VISITED_BY_CHIPS_THISROUND = FALSE;
    }
}


/***************************************************************************//**
 * @brief Print summaries associated with CHiPs visits.
 * 
 * @details Store the number of people who have been visited by CHiPs 
 * 0, 1, 2, 3, 4, 5+ times.The upper limit should be bigger than the 4 visits 
 * scheduled in PopART, but in principle this could be more than 5 visits - 
 * just change the size of MAXVISITSRECORDED.  This function records only people
 * who are alive at the end of the round.  Returns nothing, but prints to the
 * screen.
 * 
 * @param age_list Pointer to @ref age_list_struct structure.
 * @param t Current time.
 ****************************************************************************/

void print_chips_statistics_using_age_list(age_list_struct *age_list, double t){
    int g,aa,ai,i, v;
    long n_visited_by_chips_this_round = 0;
    int MAXVISITSRECORDED = 5;
    long *nvisits_distribution[N_GENDER];
    long n_chips_start_art = 0;
    long denom_chips_visits[N_GENDER]; /* Number of M/F eligible to be visited by ChiPs this round. */
    individual *person;            /* Temporary pointer to person currently being used, to make code more readable. As pointing to existing memory no malloc used. */
    for (g=0;g<N_GENDER;g++){
        nvisits_distribution[g] = malloc((MAXVISITSRECORDED+1)*sizeof(long)); /* Note that we go from 0..MAXVISITSRECORDED. */
        denom_chips_visits[g] = 0;
    }
    for (g=0;g<N_GENDER;g++)
        /* Again, we go from 0..MAXVISITSRECORDED - hence "<=" rather than "<". */
        for (v=0;v<=MAXVISITSRECORDED;v++)
            nvisits_distribution[g][v] = 0;
    int NDIED = 0;
    for (g=0;g<N_GENDER;g++){
        for (aa=(AGE_CHIPS-AGE_ADULT); aa<(MAX_AGE-AGE_ADULT); aa++){
            ai = age_list->age_list_by_gender[g]->youngest_age_group_index + aa; /* ai is the index of the array age_list->number_per_age_group of the age group of people you want to be dead */
            while (ai>(MAX_AGE-AGE_ADULT-1))
                ai = ai - (MAX_AGE-AGE_ADULT);
            for(i=0;i<age_list->age_list_by_gender[g]->number_per_age_group[ai];i++){
                person = age_list->age_list_by_gender[g]->age_group[ai][i];
                /* Note that there is a slight issue that people may be visited successfully by chips, then die (or be scheduled to receive a visit but die beforehand). To get around this use a special value for VISITED_BY_CHIPS_THISROUND. */
                if (person->VISITED_BY_CHIPS_THISROUND>DIEDBEFORECHIPSVISIT){
                    denom_chips_visits[g]++;
                    /* Since VISITED_BY_CHIPS_THISROUND=0 if not visited, and 1 if visited, can sum over this: */
                    if (person->VISITED_BY_CHIPS_THISROUND!=0 && person->VISITED_BY_CHIPS_THISROUND!=1)
                        printf("Error: at t=%6.4lf person.VISITED_BY_CHIPS_THISROUND=%i person.id=%li person.NCHIPSVISITS = %i\n",t,person->VISITED_BY_CHIPS_THISROUND,person->id,person->NCHIPSVISITS);
                    n_visited_by_chips_this_round += person->VISITED_BY_CHIPS_THISROUND;
                    /* v is the lifetime number of visits this individual has had by CHiPs. */
                    v = person->NCHIPSVISITS;
                    if (v>=MAXVISITSRECORDED)
                        nvisits_distribution[g][MAXVISITSRECORDED]++;
                    else
                        nvisits_distribution[g][v]++;
                    if(person->VISITEDBYCHIPS_TO_INIT_ART)
                        n_chips_start_art++;
                }else{
                    NDIED++;
                }
            }
        }
    }
    printf("Number died = %i\n",NDIED);
    printf("%6.4lf n_thisround=%ld ",t,n_visited_by_chips_this_round);
    for (v=0;v<=MAXVISITSRECORDED;v++)
        for (g=0;g<N_GENDER;g++)
            printf("%ld ",nvisits_distribution[g][v]);
    printf("\n");

    /* Free memory. Note - it's really important to NOT free person
    (as this memory is in use and was not locally allocated). */
    for (g=0; g<N_GENDER; g++)
        free(nvisits_distribution[g]);
}


/***************************************************************************//**
 * @brief Print summaries associated with CHiPs visits in the model
 * 
 * @details Print the number of people who have been visited 
 * 0, 1, 2, 3, 4, 5+ times by CHiPs in the model.  Used for testing/debugging.
 * The function does not return a value, but prints to the screen.  Not 
 * called in the model in typical runs.
 * 
 * @param patch Pointer to @ref patch_struct structure.
 * @param p Index of patch.
 * @param t Current time.
 ****************************************************************************/

void print_chips_statistics_using_chipsonly(patch_struct *patch, int p, double t){
    int g,ac,i, v;
    long n_visited_by_chips_this_round = 0;

    /* Store the number of people who have been visited 0, 1, 2, 3, 4, 5+ times.
    The upper limit should be bigger than the 4 visits scheduled in PopART, 
    but in principle this could be more than 5 visits - you just need
    to change the size of MAXVISITSRECORDED below. */
    int MAXVISITSRECORDED = 5;
    long *nvisits_distribution[N_GENDER];
    long n_chips_start_art = 0;
    /* Number of M/F eligible to be visited by ChiPs this round. */
    long denom_chips_visits[N_GENDER];
    for(g=0;g<N_GENDER;g++){
        /* Note: this runs from 0..MAXVISITSRECORDED. */
        nvisits_distribution[g] = malloc((MAXVISITSRECORDED+1)*sizeof(long));
        denom_chips_visits[g] = 0;
    }
    for(g=0; g<N_GENDER; g++)
        /* Again, we go from 0..MAXVISITSRECORDED - hence "<=" rather than "<". */
        for(v=0; v<=MAXVISITSRECORDED; v++)
            nvisits_distribution[g][v] = 0;

    for(g=0; g<N_GENDER; g++){
        /* Run from AGE_CHIPS to 80+. */
        for(ac=0; ac<(MAX_AGE-AGE_CHIPS+1); ac++){
            for(i=0; i<patch[p].chips_sample->number_to_visit[g][ac]; i++){
                denom_chips_visits[g]++;
                /* Since VISITED_BY_CHIPS_THISROUND=0 if not visited, and 1 if visited, can sum over this: */
                n_visited_by_chips_this_round += 
                    patch[p].individual_population[patch[p].chips_sample->list_ids_to_visit[g][ac][i]].VISITED_BY_CHIPS_THISROUND;
                /* v is the lifetime number of visits this individual has had by CHiPs. */
                v = patch[p].individual_population[patch[p].chips_sample->list_ids_to_visit[g][ac][i]].NCHIPSVISITS;
                if (v>=MAXVISITSRECORDED)
                    nvisits_distribution[g][MAXVISITSRECORDED]++;
                else
                    nvisits_distribution[g][v]++;
                if(patch[p].individual_population[patch[p].chips_sample->list_ids_to_visit[g][ac][i]].VISITEDBYCHIPS_TO_INIT_ART)
                    n_chips_start_art++;
            }
        }
    }
    printf("%6.4lf n_thisround=%ld ", t,n_visited_by_chips_this_round);
    for(v=0;v<=MAXVISITSRECORDED;v++)
        for(g=0; g<N_GENDER; g++)
            printf("%ld ", nvisits_distribution[g][v]);
    printf("\n");

    /* Free memory. Note - it's really important to NOT free person 
    (as this memory is in use and was not locally allocated). */
    for(g=0;g<N_GENDER;g++)
        free(nvisits_distribution[g]);
}


/***************************************************************************//**
 * @brief Output hazard data over time period
 * 
 * @details Function writes data to file so there is no return value.  Called 
 * within @ref hiv.c.
 * 
 * @param t Time at which hazard is calculated.
 * @param hazard Hazard value.
 * @param susceptible Pointer to susceptible individual.
 * @param pos_partner Pointer to positive partner.
 * @param file_data_store Pointer to @ref file_struct structure that holds file 
 * locations.
 * @param output Pointer to @ref output_struct structure that holds output data.
 ****************************************************************************/

void output_hazard_over_time_period(double t, double hazard, individual *susceptible, 
    individual *pos_partner, file_struct *file_data_store, output_struct *output){
    int infector_hiv_cd4_acute;
    char temp_string_hazard[100];

    /* Merge acute and CD4 to reduce output: */
    if (pos_partner->HIV_status==ACUTE)
        infector_hiv_cd4_acute = -1;
    else
        infector_hiv_cd4_acute = pos_partner->cd4;

    sprintf(temp_string_hazard,"%lf,%8.6lf,%i,%i,%i,%i,%i,%i,%i,%8.6lf\n",
        hazard,t,susceptible->circ,pos_partner->gender,infector_hiv_cd4_acute, 
        pos_partner->ART_status, pos_partner->patch_no, susceptible->sex_risk, 
        pos_partner->sex_risk,pos_partner->SPVL_num_E+pos_partner->SPVL_num_G);

    /* The -2 is because in C the last character in any string of length n is "\0" - 
    so we only have n-1 characters in the array we can write to. Make -2 instead 
    of -1 to be a bit more sure! */
    if((strlen(output->hazard_output_string)+strlen(temp_string_hazard))>(HAZARD_OUTPUT_STRING_LENGTH-2)){
        write_hazard_data(file_data_store, output->hazard_output_string);
        /* Empty output->phylogenetics_output_string. To do this we just need to set the first character to be '\0'. */
        (output->hazard_output_string)[0] = '\0';
    }
    /* Add to existing hazard output string. */
    strcat(output->hazard_output_string,temp_string_hazard);
}


/***************************************************************************//**
 * @brief Write hazard data to file
 * 
 * @param file_data_store pointer to @ref file_struct structure that holds
 * file names.
 * @param hazard_output_string Hazard data, in a string format, to be written to
 * file.
 ****************************************************************************/

void write_hazard_data(file_struct *file_data_store, char *hazard_output_string){
    file_data_store->HAZARD_FILE = fopen(file_data_store->filename_hazard_output, "a");
    fprintf(file_data_store->HAZARD_FILE,"%s", hazard_output_string);
    fclose(file_data_store->HAZARD_FILE);
}


/***************************************************************************//**
 * @brief Create a blank file for the hazard date and write header
 * 
 * @details Function writes data to file so there is no return value.
 * 
 * @param file_data_store pointer to @ref file_struct structure that 
 * holds file names.
 ****************************************************************************/

void blank_hazard_file(file_struct *file_data_store){
    file_data_store->HAZARD_FILE = fopen(file_data_store->filename_hazard_output,"w");
    fprintf(file_data_store->HAZARD_FILE,"Hazard,t,SusceptibleCircStatus,PartnerGender,PartnerHIVstatus,PartnerARTstatus,PartnerPatchNumber,SusceptibleRiskGp,PartnerRiskGp,SPVL\n");
    fclose(file_data_store->HAZARD_FILE);
}
