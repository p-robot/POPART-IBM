/**************************************************************************//**
 * @file constants.h
 * @brief Defines fixed constants and creates fixed-size global arrays
*****************************************************************************/

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include "compat.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

// For switching on and off different parts of the code

/** @brief Should the efficacy of VMMC be switched to zero (0) 
 * after the start of PopART or not (1)?
 * @details Default is 0. */
#define VMMC_EFF_ZERO_AT_POPART_START 0

/** @brief Run checks related to partnerships.
 * @details If set to 1, runs the functions @ref check_partnership_formation, 
 * @ref check_partnership_formation_and_HIV_acquisition, 
 * @ref check_partnership_dissolution, and 
 * @ref check_draw_number_partnership.  See call in @ref main.c.  Default is 0. */
#define SIMPLE_PARTNERSHIP_CHECK 0

/** @brief Run check on whole population annually to ensure different arrays are
 * aligned.
 * @details Checks concordance in the data structures holding the list of 
 * susceptibles, serodiscordant partnerships, and list of available partners.  
 * Note that this significantly slows the code (as all individuals and their 
 * partners are checked annually within the simulation).  This function does 
 * not produce any file outputs but will throw an error if someone is not in a list.
 * If it runs without throwing an error then the checks passed.  Default is 0.  */
#define SWEEP_THROUGH_TO_CHECK_LISTS 0

/** @brief Run various checks on number of partners outside the community
 * @details Run checks on number of partners outside the community, the number 
 * of HIV-positive partners and the number of HIV-positive partners outside
 * the community is recorded correctly for everyone this does not produce
 * any specific outputs, but stops with an error message if someone is 
 * not in a list it should belong to so if the code runs until the end with
 * no error message it means the numbers are updated as expected.  Default is 0. */
#define SWEEP_THROUGH_TO_CHECK_N_PARTNERS_OUTSIDE_AND_N_HIVPOS_PARTNERS_AND_N_HIVPOS_PARTNERS_OUTSIDE 0

/** @brief Run checks on age and risk assortativity
 * @details Runs checks at the time of partnership formation and once a year, 
 * of age and risk assortativity within partnerships.  Several output files
 * are generated:\n
 * `Age_assortativity_at_partnership_formation_runX.csv`\n
 * `Age_assortativity_cross_sectional_runX.csv`\n
 * `Risk_assortativity_at_partnership_formation_runX.csv`\n
 * `Risk_assortativity_cross_sectional_runX.csv`.\n  Default is 0. */
#define CHECK_AGE_AND_RISK_ASSORTATIVITY 0

/** @brief Run checks on partnership duration within the model
 * @details Checks are run at partnership formation and of partnership durations.
 * Several files are generated:\n
 * `duration_partnership_between_high_high.csv`\n
 * `duration_partnership_between_low_low.csv`\n
 * `duration_partnership_between_med_med.csv`\n
 * `duration_partnership_within_high_high.csv`\n
 * `duration_partnership_within_low_low.csv`\n
 * `duration_partnership_within_med_med.csv`.\n Default is 0. */
#define DEBUG_PARTNERSHIP_DURATION 0

/** @brief Write files of the form annual partnerships (`Annual_partnerships_outputs_*.csv`) 
 * @details Default is 0.*/
#define WRITE_ANNUAL_PARTNERSHIPS_OUTPUTS 0

/** @brief Generates the files on HIV-related survival (`HIVsurvival_individualdata.csv`)
 * @details The files generated contain DoB, DoD, gender, date first on ART etc for all 
 * HIV-positive individuals in simulation.  Default is 0.  */
#define WRITE_HIVSURVIVAL_OUTPUT 0

// Controlling the verbosity of printing
/** @brief Control the verbosity of printing to screen.  
 * @details Set to 1 if want to print normally to screen, 
 * set to 0 to minimise output to screen. */
#define VERBOSE_OUTPUT 0

/** @brief Should checks be run on parameters?  
 * @details Setting to 1 will run the function 
 * @ref check_if_parameters_plausible to check if parameters are in
 * a suitable range. Set to 0 if parameter values outside suitable ranges,
 * for example for debugging purposes. Default 1.*/
#define CHECKPARAMS 1

/** @brief Print extra demographic info to screen.  Set to 1 if want to print
 * extra demographic information to screen, 0 if not. */
#define PRINT_DEBUG_DEMOGRAPHICS 0

/** @brief Generates the files `NBirthsNNewAdultsNdeaths_Run*.csv` used 
 * for model validation
 * @details Outputs the number of births, new adults and deaths in a year. 
 * Compare with output of files produced by 
 * @ref write_one_year_age_groups_including_kids(). Default is 0. */
#define WRITE_DEBUG_DEMOGRAPHICS_NBIRTHS_NEWADULTS_DEATHS 0

/** @brief Generate the file `Age_distribution_check*.csv` for model validation.
 * @details Outputs the age distribution of the adult population 
 * using UNPD 5 year age-groups Used for for model validation.
 * Default is 0. */
#define WRITE_DEBUG_DEMOGRAPHICS_AGE_DISTRIBUTION_BY_GENDER 0

/** @brief Generates the files `OneYearAgeGp.csv` used for model validation
 * @details Outputs the age distribution of the population in one year age
 *  groups including kids. Note that this is not disaggregated by gender. 
 * Default is 0. */
#define WRITE_DEBUG_DEMOGRAPHICS_AGE_DISTRIBUTION_ONEYEARINCKIDS 0

/** @brief Generates the single file `LifeExpectancy_Za.csv` or `LifeExpectancy_SA.csv`
 * @details These files are only generated if the parameter `end_time_simul` 
 * is 2100 or later (so enough time for the simulated people to have
 * all deceased).  Default is 0.*/
#define WRITE_DEBUG_DEMOGRAPHICS_LIFE_EXPECTANCY 0

/** @brief Prints additional output to screen from within @ref input.c 
 * @details Default is 0. */
#define PRINT_DEBUG_INPUT 0

/** @brief Generates the files `DEBUG_HIV_initialSPVLdistribution*.csv` for 
 * model validation
 * @details Default is 0. */
#define WRITE_DEBUG_INITIAL_SPVL_DISTRIBUTION 0


/** @brief Generates the files `DEBUG_HIV_CD4_after_seroconversion*.csv`
 * @details Used for model validation. Default is 0. */
#define WRITE_DEBUG_CD4_AFTER_SEROCONVERSION 0

/** @brief Generates the files `DEBUG_HIVduration*.csv`
 * @details These files are used to validate the time each person 
 * spends as HIV-positive.  Default is 0. */
#define WRITE_DEBUG_HIV_DURATION 0

/** @brief Generates the files `DEBUG_HIVduration_KM*.csv`
 * @details This file is used to do a Kaplan-Meier survival analysis of 
 * time individuals spend as HIV-positive.  Default is 0. */
#define WRITE_DEBUG_HIV_DURATION_KM 0

/** @brief Generates the files `DEBUG_HIVstates_population*.csv`
 * @details Files generated are used in cross-validation of the model. 
 * Default is 0. */
#define WRITE_DEBUG_HIV_STATES 0

/** @brief Used in generated of `DEBUG_HIVstates_population*.csv`
 * @details Determines the maximum year for which the HIV state output is
 * when @ref WRITE_DEBUG_HIV_STATES is turned on (1).  Default is 2016. */
#define DEBUG_MAX_HIV_STATE_OUTPUT_TIME 2016

/** @brief Generates the files `DEBUG_ART_population_*.csv`
 * @details These files are used by `generate_art_distribution_files.py`
 * to make the files `ART_distribution.csv` and `ART_transition_dist.csv` files. 
 * Default is 0. */
#define WRITE_DEBUG_ART_STATE 0

/** @brief Write totals of number of individuals for each `ART_status`, stratified by sex and year of age.
 * @details Totals are written for each time step.  Default is 0. */
#define WRITE_ART_STATUS_BY_AGE_SEX 0

/** @brief Generate the files `CHIPS_outputs_annual*.csv` containing the 
 * annual data on people when they are visited by CHiPs. 
 * @details Default is 0. */
#define WRITE_DEBUG_CHIPS_STATES 0

/** @brief Write calibration output (`Calibration*.csv` files) to disk.
 * @details Set to 1 (default) if want to write calibration output to disk.*/
#define WRITE_CALIBRATION 1

/** @brief Print everything regardless of fitting
 * @details Default is 1. */
#define PRINT_ALL_RUNS 1

/** @brief Should an output file be generated for each run?
 * @details When calibrating the model this is set to 1.  Default is 1. */
#define PRINT_EACH_RUN_OUTPUT 1

/** @brief Generates the files `Timestep_outputs*.csv` of core indicators at each time step.
 * @details These files provide core indicators (population size, HIV-positive population,
 * cascade statistics, etc) at each time step of the simulation since the start of the 
 * introduction of HIV.  Default is 1. */
#define WRITE_EVERYTIMESTEP 1

/** @brief Generates the files `Timestep_age_outputs_*.csv` to disk.
 * @details Default is 0. */
#define TIMESTEP_AGE 0

/** @brief Patch identifier for which phylo information (transmission files) 
 * should be written for.
 * @details Used in combination with @ref WRITE_PHYLOGENETICS_OUTPUT.
 * Default is 0 (inside patch). */
#define PHYLO_PATCH 0

/** @brief Should phylogenetics output (transmission trees) be printed to file?
 * @details Only one patch is printed at a time (set by @ref PHYLO_PATCH), two
 * files are output - a file of HIV-positive individuals in the simulation 
 * and a file of transmission events.  Together these can be used to recreate
 * transmission trees.  For transmission trees including both patches (e.g. 
 * transmissions from outside the community) then this has to be run for both
 * @ref PHYLO_PATCH values of 0 and 1. */
#define WRITE_PHYLOGENETICS_OUTPUT 0

/** @brief Should the sexual network be written to file at fixed times to allow network plots?
 * @details Writes the files `Partnership_network_*.csv` to disk.  The years at which 
 * partnerships are output are hard-coded in @ref main.c. */
#define WRITE_PARTNERSHIP_NETWORK_SNAPSHOT 0

/** @brief Should the file `Partner_outside_inside_patch0.csv` be written to disk? */
#define WRITE_PARTNERS_OUTSIDE_COMMUNITY 0

/** @brief Generates the files `Hazards_*.csv` */
#define WRITE_HAZARDS 0

/** @brief Generates the files `Distr_n_lifetime_partners*.csv` and 
 * `Distr_n_partners_lastyear*.csv` files. 
 * @details These files are needed for generating 
 * `ReadAnnualOutputs-knitr.Rnw` */
#define WRITE_PARTNERSHIPS_AT_PC0 0

/** @brief Print additional outputs for a specific individual
 * @details Adjusting this variable to the ID of a specific individual will
 * print several specific events for this individual.  Setting this to -1 
 * will run the simulation normally.  Default is -1. */
#define FOLLOW_INDIVIDUAL -1

/** @brief Print additional outputs for a specific patch
 * Default is 0. */
#define FOLLOW_PATCH 0

/** @brief Generates a new file `cost_effectiveness_$.csv` used in the
 * cost-effectiveness analysis.
 * @details Default is 0. */
#define WRITE_COST_EFFECTIVENESS_OUTPUT 0

/** @brief Generates additional output for aligning the model 
 * with modelling work used in the TREATS clinical trial. 
 * @details Default is 0. */
#define WRITE_TREATS_OUTPUT 0

// Random number variables
/** @brief Used as part of GSL's RNG routines. */
const gsl_rng_type * TYPE_RNG;

/** @brief Used as part of GSL's RNG routines. */
gsl_rng * rng;

// General variables

/** @brief Macro for FALSE. */
#define FALSE 0

/** @brief Macro for TRUE. */
#define TRUE 1

/** @brief Macro for a dummy value, used in debugging. */
#define DUMMYVALUE -99

// Time variables
/** @brief Number of time steps within 1 year.
 * @details Default is 48, which is 1/4th of a month (approximately a week). */
#define N_TIME_STEP_PER_YEAR 48

/** @brief Timestep phrased as a fraction of a year.
 * @details Default is 1.0/@ref N_TIME_STEP_PER_YEAR. */
#define TIME_STEP 1.0/N_TIME_STEP_PER_YEAR

/** @brief Maximum number of years the simulation will run for */
#define MAX_N_YEARS 200

/** @brief When we want post-PopART CHiPs to roll out in the outside patches.
 * @details Default is 2020. */
#define T_ROLLOUT_CHIPS_EVERYWHERE 2020

/** @brief  */
#define ROLL_OUT_CHIPS_INSIDE_PATCH 1

/** @brief The year when to stop roll out of CHiPs to inside patch.
 * @details Default is 2100. */
#define T_STOP_ROLLOUT_CHIPS_INSIDE_PATCH 2100

/** @brief Should post-PopART rollout of CHiPs be allowed in counterfactual simulations?
 * @details Default is 0, so that post-trial CHiPs is not allowed in counterfactual
 * simulations */
#define ALLOW_COUNTERFACTUAL_ROLLOUT 0

// Settings

/** @brief Label for Zambia */
#define ZAMBIA 1

/** @brief Label for South Africa */
#define SOUTH_AFRICA 2

/** @brief Cut-off for Zambian community numbers
 * @details First 12 clusters in M&E reports are always Zambia, 
 * so cluster numbers 1-12 for Zambia and >12 for South Africa: */
#define IS_ZAMBIA 12

/** @brief The number of patches in the model */
#define NPATCHES 2

/** @brief Labels for the arms of the trial
 * @details Trial arm: 0=ARM C, 1=ARM A, 2=ARM B. */
#define ARM_C 0
#define ARM_A 1
#define ARM_B 2

/** @brief Time denoting the start of PC0 */
#define TIME_PC0 2014.5

// Demographics

/** @brief The maximum population size, used for memory allocation. 
 * @details Default is 500,000. */
#define MAX_POP_SIZE 500000

/** @brief Maximum number of people in each adult age year group (ie 13, 14, 15...)
 * @details The denominator is chosen to be really conservative (each age year will 
 * be <<5% of adult population)- as this is used for static memory allocation. 
 * It is also the maximum number of people who can die in each age group in a 
 * given timestep.  This variable may need to be increased if very long projections
 * are made.  Default is @ref MAX_POP_SIZE /6. */
#define MAX_N_PER_AGE_GROUP MAX_POP_SIZE/6

/** @brief Age at which individuals enter the simulation.
 * @details Default is 13. */
#define AGE_ADULT 13

/** @brief Number of age groups.
 * @details Default is 7. */
#define N_AGE 7

/** @brief Number of age groups when using UNPD 5-year agegroups
 * @details Default is 14. */
#define N_AGE_UNPD  14

/** @brief Age at which CHiPs visit can start to occur.
 * @brief Default is 18. */
#define AGE_CHIPS 18

/** @brief Maximum number of time steps per CHiPs round
 * @details Default is 96 (so 2 years) so no CHiPs round can
 * last greater than 2 years (without memory allocation issues). */
#define MAX_N_TIMESTEPS_PER_CHIPS_ROUND 96

/** @brief Number of CHiPs rounds used for calibration.
 * @details Default is 3. */
#define NCHIPSROUNDSFORFITTING 3

/** @brief Number of DHS rounds of data, used for fitting
 * @details This is used for allocating memory to the array storing the times of 
 * the DHS rounds.  There is also `param->DHS_params->NDHSROUNDS` that is read
 * from parameter files at input time.  Default is 4. */
#define NDHSROUNDS_MAX 4

/** @brief Minimum age in DHS data
 * @details Default is 15, used for determining size of arrays. */
#define AGE_DHS_MIN 15

/** @brief Maximum age in DHS data
 * @details Default is 59, used for determining size of arrays. */
#define AGE_DHS_MAX 59

/** @brief Maximum range of DHS ages (by single year)
 * @details Default is 45, (DHS typically runs from 15-59 
 * so that's 45 age groups) used for determining size of arrays. */
#define DHS_AGE_RANGE_MAX 45

// PC constants
/** @brief Flag to enroll and follow Population Cohort (PC) samples
 * @details Default is 0. */
#define RUN_PC 0

/** @brief Number of rounds in the Population Cohort (PC) (e.g. PC0, PC12, PC24, PC36).
 * @details Default is 4. */
#define NPC_ROUNDS 4

/** @brief Number of rounds of enrollment into the Population Cohort (PC)
 * @details Enrollment at PC0 only within the model. Default is 1. */
#define NPC_ENROLMENTS 1

/** @brief Minimum age at which enrollment into the Population Cohort can occur. 
 * @details Default is 18. */
#define AGE_PC_MIN 18

/** @brief Maximum age at which enrollment into the Population Cohort can occur. 
 * @details Default is 44. */
#define AGE_PC_MAX 44

/** @brief Age range (in years) of the Population Cohort age groups.
 * @details The Population Cohort (PC) enrolls individuals from 18-44 
 * so 27 age groups.  Default is 27. */
#define PC_AGE_RANGE_MAX 27

/** @brief Maximum number of steps in a Population Cohort round
 * @details No PC round can last >1.5 years (longest from data is 66 
 * timesteps in community 9).  Default is 72. */
#define MAX_N_TIMESTEPS_PER_PC_ROUND 72

/** @brief Number of HIV-related categories used for dividing up the 
 * Population Cohort (PC) sample.
 * @details To include the correct number of HIV-status categories
 * HIV-, HIV+ know status etc. Default is 3. */
#define N_PC_HIV_STRATA 3

/** @brief Maximum number of Population Cohort participants per sampled group
 * @details Fixes size of 
 * list_ids_in_cohort[g][ap][i_pc_category][MAX_NUMBER_PC_PARTICIPANTS_PER_GROUP].
 * Default is 200. */
#define MAX_NUMBER_PC_PARTICIPANTS_PER_GROUP 200

/** @brief Array for storing lower-bound cut-off ages for those considered for 
 * partnership formation, set in @ref constants.c */
extern const int AGE_GROUPS[N_AGE];

/** @brief Array for storing lower-bound cut-off ages for those considered for
 * partnership formation, with the oldest age group, set in @ref constants.c */
extern const int AGE_GROUPS_WITH_OLD[N_AGE+1];

/** @brief Array for storing lower-bound cut-off ages for those considered for
 * partnership formation, using UNPD age-groups, set in @ref constants.c */
extern const int AGE_GROUPS_UNPD[N_AGE_UNPD+1];

/** @brief Upper bound for age at start of simulation
 * @details Default is 80. */
#define MAX_AGE 80

/** @brief Array used to convert year-of-age to UNDP 5-year age group index.
 * @details For instance, FIND_AGE_GROUPS_UNPD[13] == 0 since a 
 * 13 yo individual would be in the first age group (with index
 * 0).  Defined in @ref constants.c. */
extern const int FIND_AGE_GROUPS_UNPD[MAX_AGE-AGE_ADULT+1];

/** @brief Array used to convert from (age - @ref AGE_ADULT) to the 
 * @ref AGE_GROUPS index
 * @details Each entry in the array is an @ref AGE_GROUPS index.  For instance, 
 * `FIND_AGE_GROUPS_UNPD[13] == 0` since a 13 yo individual would be in the 
 * first age group (with index 0).  Defined in @ref constants.c. */
extern const int FIND_AGE_GROUPS[MAX_AGE-AGE_ADULT+1];

/** @brief Number of genders. */
#define N_GENDER 2

/** @brief Indicator for male. */
#define MALE 0

/** @brief Indicator for female. */
#define FEMALE 1

/** @brief Indicator for deceased individual, used with CD4 attribute of 
 * @ref individual struct. */
#define DEAD -2

/** @brief Used to identify people who were scheduled to be visited by 
 * CHiPs in current round but died beforehand. */
#define DIEDBEFORECHIPSVISIT -2

/** @brief Gestation time, for birth rates, in timesteps
 * @details Default is 36 timesteps, which is 9 months. */
#define NSTEPS_GESTATION_TIME 36

/** @brief Number of time periods for which fertility data is given by the UNPD
 * @details Default is 30. */
#define N_UNPD_TIMEPOINTS 30

/** @brief Number of age groups in which fertility data is specified by the UNPD
 * @details Default is 7. */
#define N_AGE_UNPD_FERTILITY 7

/** @brief Number of age groups in which mortality data is specified by the UNPD
 * @details Default is 17. */
#define N_AGE_UNPD_MORTALITY 17

/** @brief UNPD estimates start in 1950-55.
 * @details Note that because of the way they are calculated the actual
 *  point time is 1952.5.*/
#define UNPD_START 1952.5

/** @brief As above end period is 2095-2100
 * @details Note we are currently using medium fertility variant projections
 *  for future demographic parameters. */
#define UNPD_END 2097.5

/** @brief UNPD fertility estimates only given for 15-49 year olds. */
#define UNPD_FERTILITY_YOUNGEST_AGE 15

/** @brief UNPD fertility estimates only given for 15-49 year olds. */
#define UNPD_FERTILITY_OLDEST_AGE 49

/** @brief Flag for whether the CHiPs sampling frame was established */
extern int POPART_SAMPLING_FRAME_ESTABLISHED;

// Partnership
/** @brief An individual can belong to up to 
 * `MAX_PARTNERSHIPS_PER_INDIVIDUAL` partnerships at any time point. */
#define MAX_PARTNERSHIPS_PER_INDIVIDUAL 15

/** @brief Maximum number of lifetime partners written in outputs. */
#define MAX_N_PARTNERS_IN_OUTPUTS 100

/** @brief Maximum number of breakups that can happen in a given time step.
 * @details Typically taken to be very conservative.  Default is 
 * @ref MAX_PARTNERSHIPS_PER_INDIVIDUAL * @ref MAX_POP_SIZE /2000 */
#define MAX_BREAKUPS_PER_TIME_STEP MAX_PARTNERSHIPS_PER_INDIVIDUAL*MAX_POP_SIZE/2000

/** @brief Number of groups of sexual risk taking behaviour */
#define N_RISK 3

/** @brief Codes for groups of riskiness of sexual behaviour of an individual
 * @details @ref LOW (0), @ref MEDIUM (1), @ref HIGH (2).  Note if these change 
 * then so does @ref RISK_GP_NAMES (defined in @ref constants.c). */
#define LOW 0
#define MEDIUM 1
#define HIGH 2

/** @brief Array for holding character names of the risk groups
 * @details Number 5 is the max length of each name +1 (ie length("High")+1) */
extern const char RISK_GP_NAMES[N_RISK][5];

// HIV-related variables
/** @brief Flag for determining the schedule of HIV testing.
 * @details Whether HIV test scheduling procedure happens for the whole population at fixed times (1),
 * or if each person has theirs scheduled sequentially (0).   Default is 1.*/
#define HIVTESTSCHEDULE 1

/** @brief How relative risk by risk group is handled in the model.
 * @details If 1, in @ref hiv_acquisition() the model boosts the high-high partnership transmission 
 * hazard by a factor of 2 and reduces the low-low partnership hazard by half.  This is based on the
 * Cori 2013 PLOS One paper.  However it may lead to there appearing to be distinct HIV epidemics in
 * low medium and high-risk.  Model sets the relative hazard by risk group to be 1 as runs with no 
 * difference in hazard seem to have a good spread, and we therefore use Occam's razor.  Default 
 * is 0*/
#define CHANGE_RR_BY_RISK_GROUP 0

/** @brief How background (non-CHiPs) HIV testing is carried out in the model
 * @details If 0, each person gets scheduled HIV tests sequentially.  If 1, we annually draw what
 * percentage of pop gets tested and then draw a time to next test (which is uniform).
 * Default is 1. */
#define DO_HIV_TESTING 1

/** @brief Flag for whether the PopART trial should be run.
 * @details Yes = 1, No = 0.  Default is 1. */
#define RUN_POPART 1

/** @brief Minimum age where HIV is introduced (during seeding).
 * @details Default is 18. */
#define YOUNGEST_AGE_SEED_WITH_HIV 18

/** @brief Maximum age where HIV is introduced (during seeding).
 * @details Default is 30. */
#define OLDEST_AGE_SEED_WITH_HIV 30

/** @brief Flag to allow people to start emergency ART
 * @details Yes = 1, No = 0.  Default is 1. */
#define ALLOW_EMERGENCY_ART 1

/** @brief Codes for HIV status 
 * @details @ref UNINFECTED (0), @ref ACUTE (1), @ref CHRONIC (2). */
#define UNINFECTED 0

/** @brief Codes for HIV status 
 * @details @ref UNINFECTED (0), @ref ACUTE (1), @ref CHRONIC (2). */
#define ACUTE 1

/** @brief Codes for HIV status 
 * @details @ref UNINFECTED (0), @ref ACUTE (1), @ref CHRONIC (2). */
#define CHRONIC 2

/** @brief Codes for ART_status of @ref individual
 * @details Codes run from -1..6 so 8 events (@ref NARTEVENTS).
 * The following codes are used:\n
 * @ref ARTNEG (-1): If never tested HIV positive (note that this is tested, 
 *  not serostatus).\n
 * @ref ARTNAIVE (0): Never been on ART\n
 * @ref EARLYART (1): First few weeks/months before achieve viral suppression. 
 *  Higher mortality and drop-out rate.\n
 * @ref LTART_VS (2): longer-term ART and Virally Suppressed (so low transmission, 
 *  no CD4 progression or drug resistance).\n
 * @ref LTART_VU (3): Longer-term ART and Virally Unsuppressed (so higher transmission, 
 *  could have (but not currently) CD4 progression and drug resistance).\n
 * @ref ARTDROPOUT (4): Has been on ART before but not currently.\n
 * @ref CASCADEDROPOUT (5): Dropped out of HIV care cascade prior to ever starting ART \n
 * @ref ARTDEATH (6): Signals that the person needs to be removed as they die while on ART.\n
 * Note that these are states and not processes. */
#define NARTEVENTS 8

/** @brief Codes for ART_status of @ref individual
 * @details Codes run from -1..6 so 8 events (@ref NARTEVENTS).
 * The following codes are used:\n
 * @ref ARTNEG (-1): If never tested HIV positive (note that this is tested, 
 *  not serostatus).\n
 * @ref ARTNAIVE (0): Never been on ART\n
 * @ref EARLYART (1): First few weeks/months before achieve viral suppression. 
 *  Higher mortality and drop-out rate.\n
 * @ref LTART_VS (2): longer-term ART and Virally Suppressed (so low transmission, 
 *  no CD4 progression or drug resistance).\n
 * @ref LTART_VU (3): Longer-term ART and Virally Unsuppressed (so higher transmission, 
 *  could have (but not currently) CD4 progression and drug resistance).\n
 * @ref ARTDROPOUT (4): Has been on ART before but not currently.\n
 * @ref CASCADEDROPOUT (5): Dropped out of HIV care cascade prior to ever starting ART \n
 * @ref ARTDEATH (6): Signals that the person needs to be removed as they die while on ART.\n
 * Note that these are states and not processes. */
#define ARTNEG  -1

/** @brief Codes for ART_status of @ref individual
 * @details Codes run from -1..6 so 8 events (@ref NARTEVENTS).
 * The following codes are used:\n
 * @ref ARTNEG (-1): If never tested HIV positive (note that this is tested, 
 *  not serostatus).\n
 * @ref ARTNAIVE (0): Never been on ART\n
 * @ref EARLYART (1): First few weeks/months before achieve viral suppression. 
 *  Higher mortality and drop-out rate.\n
 * @ref LTART_VS (2): longer-term ART and Virally Suppressed (so low transmission, 
 *  no CD4 progression or drug resistance).\n
 * @ref LTART_VU (3): Longer-term ART and Virally Unsuppressed (so higher transmission, 
 *  could have (but not currently) CD4 progression and drug resistance).\n
 * @ref ARTDROPOUT (4): Has been on ART before but not currently.\n
 * @ref CASCADEDROPOUT (5): Dropped out of HIV care cascade prior to ever starting ART \n
 * @ref ARTDEATH (6): Signals that the person needs to be removed as they die while on ART.\n
 * Note that these are states and not processes. */
#define ARTNAIVE 0

/** @brief Codes for ART_status of @ref individual
 * @details Codes run from -1..6 so 8 events (@ref NARTEVENTS).
 * The following codes are used:\n
 * @ref ARTNEG (-1): If never tested HIV positive (note that this is tested, 
 *  not serostatus).\n
 * @ref ARTNAIVE (0): Never been on ART\n
 * @ref EARLYART (1): First few weeks/months before achieve viral suppression. 
 *  Higher mortality and drop-out rate.\n
 * @ref LTART_VS (2): longer-term ART and Virally Suppressed (so low transmission, 
 *  no CD4 progression or drug resistance).\n
 * @ref LTART_VU (3): Longer-term ART and Virally Unsuppressed (so higher transmission, 
 *  could have (but not currently) CD4 progression and drug resistance).\n
 * @ref ARTDROPOUT (4): Has been on ART before but not currently.\n
 * @ref CASCADEDROPOUT (5): Dropped out of HIV care cascade prior to ever starting ART \n
 * @ref ARTDEATH (6): Signals that the person needs to be removed as they die while on ART.\n
 * Note that these are states and not processes. */
#define EARLYART 1

/** @brief Codes for ART_status of @ref individual
 * @details Codes run from -1..6 so 8 events (@ref NARTEVENTS).
 * The following codes are used:\n
 * @ref ARTNEG (-1): If never tested HIV positive (note that this is tested, 
 *  not serostatus).\n
 * @ref ARTNAIVE (0): Never been on ART\n
 * @ref EARLYART (1): First few weeks/months before achieve viral suppression. 
 *  Higher mortality and drop-out rate.\n
 * @ref LTART_VS (2): longer-term ART and Virally Suppressed (so low transmission, 
 *  no CD4 progression or drug resistance).\n
 * @ref LTART_VU (3): Longer-term ART and Virally Unsuppressed (so higher transmission, 
 *  could have (but not currently) CD4 progression and drug resistance).\n
 * @ref ARTDROPOUT (4): Has been on ART before but not currently.\n
 * @ref CASCADEDROPOUT (5): Dropped out of HIV care cascade prior to ever starting ART \n
 * @ref ARTDEATH (6): Signals that the person needs to be removed as they die while on ART.\n
 * Note that these are states and not processes. */
#define LTART_VS 2

/** @brief Codes for ART_status of @ref individual
 * @details Codes run from -1..6 so 8 events (@ref NARTEVENTS).
 * The following codes are used:\n
 * @ref ARTNEG (-1): If never tested HIV positive (note that this is tested, 
 *  not serostatus).\n
 * @ref ARTNAIVE (0): Never been on ART\n
 * @ref EARLYART (1): First few weeks/months before achieve viral suppression. 
 *  Higher mortality and drop-out rate.\n
 * @ref LTART_VS (2): longer-term ART and Virally Suppressed (so low transmission, 
 *  no CD4 progression or drug resistance).\n
 * @ref LTART_VU (3): Longer-term ART and Virally Unsuppressed (so higher transmission, 
 *  could have (but not currently) CD4 progression and drug resistance).\n
 * @ref ARTDROPOUT (4): Has been on ART before but not currently.\n
 * @ref CASCADEDROPOUT (5): Dropped out of HIV care cascade prior to ever starting ART \n
 * @ref ARTDEATH (6): Signals that the person needs to be removed as they die while on ART.\n
 * Note that these are states and not processes. */
#define LTART_VU 3

/** @brief Codes for ART_status of @ref individual
 * @details Codes run from -1..6 so 8 events (@ref NARTEVENTS).
 * The following codes are used:\n
 * @ref ARTNEG (-1): If never tested HIV positive (note that this is tested, 
 *  not serostatus).\n
 * @ref ARTNAIVE (0): Never been on ART\n
 * @ref EARLYART (1): First few weeks/months before achieve viral suppression. 
 *  Higher mortality and drop-out rate.\n
 * @ref LTART_VS (2): longer-term ART and Virally Suppressed (so low transmission, 
 *  no CD4 progression or drug resistance).\n
 * @ref LTART_VU (3): Longer-term ART and Virally Unsuppressed (so higher transmission, 
 *  could have (but not currently) CD4 progression and drug resistance).\n
 * @ref ARTDROPOUT (4): Has been on ART before but not currently.\n
 * @ref CASCADEDROPOUT (5): Dropped out of HIV care cascade prior to ever starting ART \n
 * @ref ARTDEATH (6): Signals that the person needs to be removed as they die while on ART.\n
 * Note that these are states and not processes. */
#define ARTDROPOUT 4

/** @brief Codes for ART_status of @ref individual
 * @details Codes run from -1..6 so 8 events (@ref NARTEVENTS).
 * The following codes are used:\n
 * @ref ARTNEG (-1): If never tested HIV positive (note that this is tested, 
 *  not serostatus).\n
 * @ref ARTNAIVE (0): Never been on ART\n
 * @ref EARLYART (1): First few weeks/months before achieve viral suppression. 
 *  Higher mortality and drop-out rate.\n
 * @ref LTART_VS (2): longer-term ART and Virally Suppressed (so low transmission, 
 *  no CD4 progression or drug resistance).\n
 * @ref LTART_VU (3): Longer-term ART and Virally Unsuppressed (so higher transmission, 
 *  could have (but not currently) CD4 progression and drug resistance).\n
 * @ref ARTDROPOUT (4): Has been on ART before but not currently.\n
 * @ref CASCADEDROPOUT (5): Dropped out of HIV care cascade prior to ever starting ART \n
 * @ref ARTDEATH (6): Signals that the person needs to be removed as they die while on ART.\n
 * Note that these are states and not processes. */
#define CASCADEDROPOUT 5

/** @brief Codes for ART_status of @ref individual
 * @details Codes run from -1..6 so 8 events (@ref NARTEVENTS).
 * The following codes are used:\n
 * @ref ARTNEG (-1): If never tested HIV positive (note that this is tested, 
 *  not serostatus).\n
 * @ref ARTNAIVE (0): Never been on ART\n
 * @ref EARLYART (1): First few weeks/months before achieve viral suppression. 
 *  Higher mortality and drop-out rate.\n
 * @ref LTART_VS (2): longer-term ART and Virally Suppressed (so low transmission, 
 *  no CD4 progression or drug resistance).\n
 * @ref LTART_VU (3): Longer-term ART and Virally Unsuppressed (so higher transmission, 
 *  could have (but not currently) CD4 progression and drug resistance).\n
 * @ref ARTDROPOUT (4): Has been on ART before but not currently.\n
 * @ref CASCADEDROPOUT (5): Dropped out of HIV care cascade prior to ever starting ART \n
 * @ref ARTDEATH (6): Signals that the person needs to be removed as they die while on ART.\n
 * Note that these are states and not processes. */
#define ARTDEATH 6

/** @brief CD4 value to identify that people are not infected with HIV.
 * @details Note: CD4==-2 means the person is dead, see @ref DEAD.
 * 0="CD4>500", 1="CD4 350-500", 2="CD4 200-350", 3="CD4 <200". */
#define CD4_UNINFECTED -1

/** @brief Number of CD3 categories (four)
 * @details Model includes four CD4 categories 
 * (0=">500", 1="350-500", 2= "200-350", 3="<200").
 * Default is 4. */
#define NCD4 4

// Recency assays
/** @brief Immune biomarker value if individual is uninfected with HIV
 * @details Dummy value used for simulation of recency assay.  Default is -1. */
#define BIOMARKER_UNINFECTED -1

/** @brief Should a proportion of the population be written to file annually?
 * @details Default is 0. */
#define WRITE_ANNUAL_SAMPLING 0

/** @brief Number of set-point viral load (SPVL) categories
 * @details 0="<4"; 1="4-4.5"; 2="4.5-5"; 3=">5" */
#define NSPVL 4

/** @brief Should set-point viral load be inherited from infector?
 * @details Yes (1), No (0).  Default is 0. */
#define SPVL_INHERITANCE 0

/** @brief HIV progression events for recording next cascade events of individuals
 * @details For coding the `next_HIV_event` or `next_cascade_event` on an 
 * @ref individual struct.
 * The following codes are used:\n
 * @ref HIVEVENT_ENDOFACUTE (0): Next event is the end of the acute phase \n
 * @ref HIVEVENT_CD4_PROGRESSION (1): Next event is CD4 progression\n
 * @ref HIVEVENT_STARTART_LOWCD4 (2): Next event is starting ART from low CD4 count\n
 * @ref HIVEVENT_AIDSDEATH (3): Next event is AIDS-related death\n
 * @ref NOEVENT (-1): There is no 'next' event for this individual.\n
 * @ref EVENTAFTERENDSIMUL (-2): The next event will take place after the end 
 *  of the simulation */
#define HIVEVENT_ENDOFACUTE 0

/** @brief HIV progression events for recording next cascade events of individuals
 * @details For coding the `next_HIV_event` or `next_cascade_event` on an 
 * @ref individual struct.
 * The following codes are used:\n
 * @ref HIVEVENT_ENDOFACUTE (0): Next event is the end of the acute phase \n
 * @ref HIVEVENT_CD4_PROGRESSION (1): Next event is CD4 progression\n
 * @ref HIVEVENT_STARTART_LOWCD4 (2): Next event is starting ART from low CD4 count\n
 * @ref HIVEVENT_AIDSDEATH (3): Next event is AIDS-related death\n
 * @ref NOEVENT (-1): There is no 'next' event for this individual.\n
 * @ref EVENTAFTERENDSIMUL (-2): The next event will take place after the end 
 *  of the simulation */
#define HIVEVENT_CD4_PROGRESSION 1

/** @brief HIV progression events for recording next cascade events of individuals
 * @details For coding the `next_HIV_event` or `next_cascade_event` on an 
 * @ref individual struct.
 * The following codes are used:\n
 * @ref HIVEVENT_ENDOFACUTE (0): Next event is the end of the acute phase \n
 * @ref HIVEVENT_CD4_PROGRESSION (1): Next event is CD4 progression\n
 * @ref HIVEVENT_STARTART_LOWCD4 (2): Next event is starting ART from low CD4 count\n
 * @ref HIVEVENT_AIDSDEATH (3): Next event is AIDS-related death\n
 * @ref NOEVENT (-1): There is no 'next' event for this individual.\n
 * @ref EVENTAFTERENDSIMUL (-2): The next event will take place after the end 
 *  of the simulation */
#define HIVEVENT_STARTART_LOWCD4 2

/** @brief HIV progression events for recording next cascade events of individuals
 * @details For coding the `next_HIV_event` or `next_cascade_event` on an 
 * @ref individual struct.
 * The following codes are used:\n
 * @ref HIVEVENT_ENDOFACUTE (0): Next event is the end of the acute phase \n
 * @ref HIVEVENT_CD4_PROGRESSION (1): Next event is CD4 progression\n
 * @ref HIVEVENT_STARTART_LOWCD4 (2): Next event is starting ART from low CD4 count\n
 * @ref HIVEVENT_AIDSDEATH (3): Next event is AIDS-related death\n
 * @ref NOEVENT (-1): There is no 'next' event for this individual.\n
 * @ref EVENTAFTERENDSIMUL (-2): The next event will take place after the end 
 *  of the simulation */
#define HIVEVENT_AIDSDEATH 3

/** @brief HIV progression events for recording next cascade events of individuals
 * @details For coding the `next_HIV_event`  or `next_cascade_event` on an
 * @ref individual struct.
 * The following codes are used:\n
 * @ref HIVEVENT_ENDOFACUTE (0): Next event is the end of the acute phase \n
 * @ref HIVEVENT_CD4_PROGRESSION (1): Next event is CD4 progression\n
 * @ref HIVEVENT_STARTART_LOWCD4 (2): Next event is starting ART from low CD4 count\n
 * @ref HIVEVENT_AIDSDEATH (3): Next event is AIDS-related death\n
 * @ref NOEVENT (-1): There is no 'next' event for this individual.\n
 * @ref EVENTAFTERENDSIMUL (-2): The next event will take place after the end 
 *  of the simulation */
#define NOEVENT -1

/** @brief HIV progression events for recording next cascade events of individuals
 * @details For coding the `next_HIV_event` or `next_cascade_event` on an 
 * @ref individual struct.
 * The following codes are used:\n
 * @ref HIVEVENT_ENDOFACUTE (0): Next event is the end of the acute phase \n
 * @ref HIVEVENT_CD4_PROGRESSION (1): Next event is CD4 progression\n
 * @ref HIVEVENT_STARTART_LOWCD4 (2): Next event is starting ART from low CD4 count\n
 * @ref HIVEVENT_AIDSDEATH (3): Next event is AIDS-related death\n
 * @ref NOEVENT (-1): There is no 'next' event for this individual.\n
 * @ref EVENTAFTERENDSIMUL (-2): The next event will take place after the end 
 *  of the simulation */
#define EVENTAFTERENDSIMUL -2


/** @brief Flag to decide if a cascade event (e.g. an HIV test) is from PopART or not.
 * @details If a process is part of PopART then processes may happen faster 
 * (e.g. time to CD4 test is quicker), and ART CD4 eligibility may be different. */
#define NOTPOPART 0

/** @brief Flag to decide if a cascade event (e.g. an HIV test) is from PopART or not.
 * @details If a process is part of PopART then processes may happen faster 
 * (e.g. time to CD4 test is quicker), and ART CD4 eligibility may be different. */
#define POPART 1

/** @brief Number of rounds of CHiPS visits to be simulated
 * @details Different to @ref NCHIPSROUNDSFORFITTING, which is the number of rounds
 * used for fitting (i.e. for which there is data for fitting).  Early in the trial
 * the number of rounds of data available for fitting would be different to the number
 * of rounds that the model needed to simulate.  Default is 3. */
#define NCHIPSROUNDS 3

/** @brief Denoting that the simulation is in post-trial. */
#define CHIPSROUNDPOSTTRIAL -1

/** @brief Number of time periods per round, used for ART initiation for which 
 * parameters are varied within a round over different time periods */
#define MAX_N_TIME_PERIODS_PER_ROUND 12

/** @brief Flag to indicate that CHiPs intervention is not currently running
 * @details This may occur between rounds of before the trial. */
#define CHIPSNOTRUNNING -2


/** @brief Number of cascade events that occur in the background cascade 
 * (or standard of care).
 * @details Note: a test is used (`EVENT<NCASCADEEVENTS_NONPOPART`) to 
 * determine if the event is due to popart (=1) or not (=0). Default is 7. */
#define NCASCADEEVENTS_NONPOPART 7

/** @brief Number of cascade events that occur in total
 * @details Default is 14. */
#define NCASCADEEVENTS 14

#define CASCADEEVENT_HIV_TEST_NONPOPART 0
#define CASCADEEVENT_CD4_TEST_NONPOPART 1
#define CASCADEEVENT_START_ART_NONPOPART 2
#define CASCADEEVENT_VS_NONPOPART 3
#define CASCADEEVENT_VU_NONPOPART 4
#define CASCADEEVENT_DROPOUT_NONPOPART 5
#define CASCADEEVENT_ARTDEATH_NONPOPART 6
#define CASCADEEVENT_HIV_TEST_POPART 7
#define CASCADEEVENT_CD4_TEST_POPART 8
#define CASCADEEVENT_START_ART_POPART 9
#define CASCADEEVENT_VS_POPART 10
#define CASCADEEVENT_VU_POPART 11
#define CASCADEEVENT_DROPOUT_POPART 12
#define CASCADEEVENT_ARTDEATH_POPART 13

/** @brief Flag to indicate individual has never tested for HIV
 * @details Used in the `time_last_hiv_test` attribute of @ref individual struct. */
#define NEVERHIVTESTED -1

/** @brief Sensitivity of HIV tests
 * @details Default is 1. */
#define HIVTESTSENSITIVITY 1

/** @brief Specificity of HIV tests
 * @details Default is 1. */
#define HIVTESTSPECIFICITY 1

/** @brief Determines how CD4 progression by SPVL is implemented
 * @details Using 4 SPVL categories (CD4PROGRESSIONMODEL=0) or using the 
 * RR from the Cox PH model (CD4PROGRESSIONMODEL=1). */
#define CD4PROGRESSIONMODEL 1

/** @brief Code for matching to `CD4PROGRESSIONMODEL` */
#define USEFOURSPVLCATEGORIES 0

/** @brief Code for matching to `CD4PROGRESSIONMODEL` */
#define USECOXPH 1

/** @brief Baseline SPVL in the Cox PH model */
#define COXMODELBASELINESPVL 4.0

/** @brief Indicator for whether this is NOT a counterfactual run
 * @details If a counterfactual run then set all arms to be arm C. 
 * Default is 0. */
#define NOT_COUNTERFACTUAL_RUN 0

/** @brief Indicator for whether this IS a counterfactual run
 * @details If a counterfactual run then set all arms to be arm C. 
 * Default is 1. */
#define IS_COUNTERFACTUAL_RUN 1

/** @brief Indicator for uncircumcised status
 * @details Used for comparison with `circ` arttribute of the 
 * @ref individual struct. */
#define UNCIRC 0

/** @brief Indicator for uncircumcised awaiting VMMC status
 * @details Used for comparison with `circ` arttribute of the 
 * @ref individual struct. */
#define UNCIRC_WAITING_VMMC 1

/** @brief Indicator for VMMC status
 * @details Used for comparison with `circ` arttribute of the 
 * @ref individual struct. */
#define VMMC 2

/** @brief Indicator for had VMMC and in healing phase
 * @details Used for comparison with `circ` arttribute of the 
 * @ref individual struct. */
#define VMMC_HEALING 3

/** @brief Indicator for having had a traditional circumcision (TMC)
 * @details Used for comparison with `circ` arttribute of the 
 * @ref individual struct. */
#define TRADITIONAL_MC 4

/** @brief Number of HIV+ people whose next progression event will 
 * happen in a given time step
 * @details Taken to be moderately conservative - in a 50% prevalence 
 * population assume average of 1 HIV event per year.  At some point 
 * in the future we will tune this to avoid too many `reallocs()`.
 * Default @ref MAX_POP_SIZE / 100. */
#define DEFAULT_N_HIV_PROGRESS_PER_TIME_STEP MAX_POP_SIZE/100

/** @brief Default number of people who will have an HIV cascade event (HIV test, start ART etc)
 * in a given time step.
 * @details Taken to be moderately conservative - in a 50% prevalence population assume average 
 * of 1 HIV event per year.  At some point in the future we will tune this to avoid too many 
 * `reallocs()`.  Default is @ref MAX_POP_SIZE / 50. */
#define DEFAULT_N_HIV_CASCADE_PER_TIME_STEP MAX_POP_SIZE/50

/** @brief Amount of extra memory to add if the model runs out of memory
 * @details Note that the value is arbitrary so can change to optimise. 
 * Default is 100. */
#define RESIZEMEM 100

/** @brief Amount of memory to add if depleted in the `planned_breakups` array.
 * Note that the value is arbitrary so can change to optimise.  Default is 100. */
#define RESIZEMEM_BREAKUP 100

/** @brief An array used by @ref hiv_acquisition() to store the per-partnership 
 * hazard from all HIV+ partners of an individual. */
double PER_PARTNERSHIP_HAZARD_TEMPSTORE[MAX_PARTNERSHIPS_PER_INDIVIDUAL];

/** @brief Defines amount that can be stored each year for a single patch to be 
 * stored in `AnnualOutputs` files.
 * @details If very long lines are to be written, this needs to increase. 
 * Default is 3000000. */
#define SIZEOF_annual_outputs_string 3000000

/** @brief Defines amount that can be stored each year for a single patch to be 
 * stored in `AnnualOutputs` PC only files.
 * @details If very long lines are to be written, this needs to increase. 
 * Default is 3000000. */
#define SIZEOF_annual_outputs_string_pconly 300000

/** @brief Size for an array storing a single indicator for a patch (e.g. prevalence, 
 * % know serostatus, % on ART). 
 * @details Default is 100000. */
#define SIZEOF_annual_outputs_tempstore 100000

/** @brief Defines amount that can be stored each year for a single patch to be 
 * stored in `Calibration` files.
 * @details Default is 3000000. */
#define SIZEOF_calibration_outputs 3000000

/** @brief Defines amount that can be stored each year for a single patch to be 
 * stored in `cost_effectiveness_*.csv` files.
 * @details Default is 3000000. */
#define SIZEOF_cost_effectiveness_outputs_string 3000000

/** @brief Flag used to avoid excessive writing to disk, only write out these 
 * data every @ref NRUNSPERWRITETOFILE runs.
 * @details *Note: Do not set this too large as it incurs excessive memory usage.*
 * Default is 100. */
#define NRUNSPERWRITETOFILE 100

/** @brief How many timesteps to store in `timestep_outputs_string`
 * @details Done to reduce size of file.
 * Default is 4. */
#define OUTPUTTIMESTEP 4

/** @brief This is the size of the arrays allocated for input/output filenames. 
 * @details Note that as we often specify full pathnames it should be quite big.
 * Default is 500. */
#define LONGSTRINGLENGTH 500

/** @brief This stores the output (all transmission info) for a run.
 * @details This stores the output we need for phylogenetics so that it is printed at the end of the simulation. 
 * It will contain the ID of infectee, infector, whether the infector was in acute infection etc.
 * Estimate that each line in phylogenetics_output_string contains roughly 30 characters, and that we may have up to 
 * 50,000 transmission events per simulation.  Default 400000.*/
#define PHYLO_OUTPUT_STRING_LENGTH 400000

/** @brief Macro used if we want to run several scenarios for PANGEA.
 * @details Default 0.*/
#define PHYLO_SCENARIO 0

/** @brief Length of array for storing hazard and associated factors.
 * @details Default 20000.*/
#define HAZARD_OUTPUT_STRING_LENGTH 20000

#endif
