/**************************************************************************//**
 * @file constants.c
 * @brief Defines fixed constants and creates fixed-size global arrays
*****************************************************************************/

#include "constants.h"

/** @brief Array for storing lower-bound cut-off ages for those considered for 
 * partnership formation.  */
const int AGE_GROUPS[N_AGE] = {13,18,23,30,40,50,60};

/** @brief Array for storing lower-bound cut-off ages for those considered for
 * partnership formation, with the oldest age group.
 * @details The 80 and over age group is only here for the ageing process. 
 * Model does not differentiate 61-79 and 80 and over for other processes.  */
const int AGE_GROUPS_WITH_OLD[N_AGE+1] = {13,18,23,30,40,50,60,80};

/** @brief Array for storing lower-bound cut-off ages for those considered for
 * partnership formation, using UNPD age-groups.  */
const int AGE_GROUPS_UNPD[N_AGE_UNPD+1] = {13,15,20,25,30,35,40,45,50,55,60,65,70,75,80};

/** @brief Array used to convert year-of-age to UNDP 5-year age group index.
 * @details For instance, FIND_AGE_GROUPS_UNPD[13] == 0 since a 
 * 13 yo individual would be in the first age group (with index
 * 0). */
const int FIND_AGE_GROUPS_UNPD[MAX_AGE - AGE_ADULT + 1] = {
    0,0, // 13, 14
    1,1,1,1,1, // 15, 16, 17, 18, 19
    2,2,2,2,2, // 20, 21, 22, 23, 24
    3,3,3,3,3, // 25, 26, 27, 28, 29
    4,4,4,4,4, // 30, 31, 32, 33, 34
    5,5,5,5,5, // 35, 36, 37, 38, 39
    6,6,6,6,6, // 40, 41, 42, 43, 44
    7,7,7,7,7, // 45, 46, 47, 48, 49
    8,8,8,8,8, // 50, 51, 52, 53, 54
    9,9,9,9,9, // 55, 56, 57, 58, 59
    10,10,10,10,10, // 60, 61, 62, 63, 64
    11,11,11,11,11, // 65, 66, 67, 68, 69
    12,12,12,12,12, // 70, 71, 72, 73, 74
    13,13,13,13,13, // 75, 76, 77, 78, 79
    14 // 80 +
};

/** @brief Array used to convert from (age - @ref AGE_ADULT) to the 
 * @ref AGE_GROUPS index
 * @details Each entry in the array is an @ref AGE_GROUPS index.  For instance, 
 * `FIND_AGE_GROUPS_UNPD[13] == 0` since a 13 yo individual would be in the 
 * first age group (with index 0).  */
const int FIND_AGE_GROUPS[MAX_AGE-AGE_ADULT+1] = {
    0,0,0,0,0,
    1,1,1,1,1,
    2,2,2,2,2,2,2,
    3,3,3,3,3,3,3,3,3,3,
    4,4,4,4,4,4,4,4,4,4,
    5,5,5,5,5,5,5,5,5,5,
    6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6
};

/** @brief Codes for groups of riskiness of sexual behaviour of an individual
 * @details @ref LOW (0), @ref MEDIUM (1), @ref HIGH (2). */
const char RISK_GP_NAMES[N_RISK][5] = {"Low","Med","High"};

/** @brief Flag for whether the CHiPs sampling frame was established */
int POPART_SAMPLING_FRAME_ESTABLISHED;
