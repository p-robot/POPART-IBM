/**************************************************************************//**
 * @file compat.h
 * @brief Handles issues related to compatibility between different compilers
*****************************************************************************/

#ifndef COMPAT_H
#define COMPAT_H

#ifdef _WIN32
   #define THISOS 0
   #define POPEN _popen
   #define PCLOSE _pclose
   #define _CRT_SECURE_NO_WARNINGS 1
   #define _CRT_SECURE_NO_DEPRECATION 1
#else
   #define THISOS 1
   #define POPEN popen
   #define PCLOSE pclose
#endif

#endif
