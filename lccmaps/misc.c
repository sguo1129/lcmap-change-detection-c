#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <getopt.h>
#include <dirent.h>
#include <fnmatch.h>
#include <string.h>
#include <limits.h>
#include <errno.h>

#include "utilities.h"
#include "lccmaps.h"

#define JULIAN_DATE_LAST_DAY_1972 720624 
#define LANDSAT_START_YEAR 1973
#define LEAP_YEAR_DAYS 366
#define NON_LEAP_YEAR_DAYS 365

/******************************************************************************
MODULE:  get_args

PURPOSE:  Gets the command-line arguments and validates that the required
arguments were specified.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
FAILURE         Error getting the command-line arguments or a command-line
                argument and associated value were not specified
SUCCESS         No errors encountered

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/5/2015    Song Guo         Original Development
20151203    Brian Davis      Added arguments for input and output
                             directories, and scene list file.

NOTES:
  1. Memory is allocated for the input and output files.  All of these should
     be character pointers set to NULL on input.  The caller is responsible
     for freeing the allocated memory upon successful return.
  2. chi2inv(T_cg, num_bands) = chi2inv(0.99, 5) = 15.0863 
  3. chi2inv(T_max_cg, num_bands) = chi2inv(1-1e-6, 5) = 35.8882 
******************************************************************************/
int get_args
(
    int argc,                 /* I: number of cmd-line args */
    char *argv[],             /* I: string of cmd-line args */
    int *start_year,          /* O: start year for generate maps */
    int *end_year,            /* O: end year for generate maps */
    int *month,               /* O: month of year for generate maps */
    int *day,                 /* O: day of year for generate maps */
    char *in_path,            /* O: directory locaiton for input data */
    char *out_path,           /* O: directory location for output files */
    char *scene_name,         /* O: one scene name to read in envi header info */
    bool *verbose             /* O: verbose flag */
)
{
    int c;                          /* current argument index */
    int option_index;               /* index for the command-line option */
    static int verbose_flag = 0;    /* verbose flag */
    char errmsg[MAX_STR_LEN];       /* error message */
    char FUNC_NAME[] = "get_args";  /* function name */
    static struct option long_options[] = {
        {"verbose", no_argument, &verbose_flag, 1},
        {"start_year", required_argument, 0, 's'},
        {"end_year", required_argument, 0, 'e'},
        {"month", required_argument, 0, 'm'},
        {"day", required_argument, 0, 'd'},
        {"in_path", required_argument, 0, 'i'},
        {"out_path", required_argument, 0, 'o'},
        {"scene_name", required_argument, 0, 'n'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    /* Loop through all the cmd-line options */
    opterr = 0; /* turn off getopt_long error msgs as we'll print our own */
    while (1)
    {
        /* optstring in call to getopt_long is empty since we will only
           support the long options */
        c = getopt_long (argc, argv, "", long_options, &option_index);
        if (c == -1)
        {                       /* Out of cmd-line options */
            break;
        }

        switch (c)
        {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
		{
                    break;
		}
		sprintf (errmsg, "option %s\n", long_options[option_index].name);
                if (optarg)
		{
		    sprintf (errmsg, "option %s with arg %s\n", 
                             long_options[option_index].name, optarg);
		}
                RETURN_ERROR (errmsg, FUNC_NAME, ERROR);
                break;

            case 'h':              /* help */
                usage ();
                exit(SUCCESS);
                break;

            case 's':             
                *start_year = atoi (optarg);
                break;

            case 'e':             
                *end_year = atoi (optarg);
                break;

            case 'm':             
                *month = atoi (optarg);
                break;

            case 'd':             
                *day = atoi (optarg);
                break;

            case 'i':
                strcpy (in_path, optarg);
                break;

            case 'o':
                strcpy (out_path, optarg);
                break;

            case 'n':
                strcpy (scene_name, optarg);
                break;

            case '?':
            default:
                sprintf (errmsg, "Unknown option %s", argv[optind - 1]);
                usage ();
                RETURN_ERROR (errmsg, FUNC_NAME, ERROR);
                break;
        }
    }

    /* If inputDir and out_path were not specified, assign local directory, */
    /* so that pre-pending a directory/path later on will not cause an    */
    /* error, but instead result in ./<filename> .                        */
    if (strlen(in_path) == 0)
        {
            strcpy (in_path, ".");
        }
    if (strlen(out_path) == 0)
        {
            strcpy (out_path, ".");
        }

    /* Check the verbose flag */
    if (verbose_flag)
        *verbose = true;
    else
        *verbose = false;

    if (*verbose)
    {
        printf ("start_year = %d\n", *start_year);
        printf ("end_year = %d\n", *end_year);
        printf ("month = %d\n", *month);
        printf ("day = %d\n", *day);
        printf ("verbose = %d\n", *verbose);
        printf ("in_path = %s\n", in_path);
        printf ("out_path = %s\n", out_path);
        printf ("scene_name = %s\n", scene_name);
    }

    return SUCCESS;
}

bool is_leap_year
(
    int year /*I: Year to test */
)
{
    if (((year % 4) != 0) || (((year % 100) == 0) && ((year % 400) != 0)))
        return false;
    else
        return true;
}


/* Calculate day of year given year, month, and day of month */
void convert_year_month_day_to_doy
(
    int year,  /* I: Year */
    int month, /* I: Month */
    int day,   /* I: Day of month */
    int *doy   /* O: Day of year */
)
{
    /* Days in month for non-leap years */
    static const int noleap[12] =
        {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    /* Days in month for leap years */
    static const int leap[12] =
        {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    int i;
    int doy_sum;

    /* Calculate day of year */
    doy_sum = 0;
    if (is_leap_year(year))
    {
        for (i = 0; i < month - 1; i++)
            doy_sum += leap[i];
    }
    else
    {
        for (i = 0; i < month - 1; i++)
            doy_sum += noleap[i];
    }
    doy_sum += day;

    *doy = doy_sum;
}

/******************************************************************************
MODULE:  convert_year_doy_to_jday_from_0000

PURPOSE:  convert day of year in a year to julian day counted from year 0000

RETURN VALUE: int
ERROR           Error for year less than 1973 as input
SUCCESS         No errors encountered

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/23/2015   Song Guo         Original Development

NOTES:
******************************************************************************/
int convert_year_mon_day_to_jday_from_0000
(
    int year,      /* I: year */
    int month,     /* I: month of the year */
    int day,       /* I: day of the year */
    int *jday      /* O: julian date since year 0000 */
)
{
    char FUNC_NAME[] = "convert_year_doy_to_jday_from_0000";
    int i;
    int doy;
    int status;

    convert_year_month_day_to_doy(year, month, day, &doy);

    if (year < 1973)
    {
        RETURN_ERROR ("Landsat data starts from 1973", FUNC_NAME, ERROR);
    }

    /* 12-31-1972 is 720624 in julian day since year 0000 */
    if (year != LANDSAT_START_YEAR)
    {
        *jday = JULIAN_DATE_LAST_DAY_1972;
        for (i = LANDSAT_START_YEAR; i < year; i++)
        {
            status = is_leap_year(i);
            if (status == true)
	    {
                *jday += LEAP_YEAR_DAYS;
	    }
            else
	    {
                *jday += NON_LEAP_YEAR_DAYS;
	    }
        }
    }    
    *jday += doy;

    return SUCCESS;
}

/******************************************************************************
MODULE:  convert_jday_from_0000_to_year_doy

PURPOSE:  convert julian day counted from year 0000 to year and doy

RETURN VALUE: int
ERROR           Error for year less than 1973 as input
SUCCESS         No errors encountered

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/23/2015   Song Guo         Original Development

NOTES:
******************************************************************************/
void convert_jday_from_0000_to_year_doy
(
    int jday,      /* I: julian date since year 0000 */
    int *year,     /* O: year */
    int *doy       /* O: month of the year */
)
{
    int yr = 1973;
    int status;

    /* 12-31-1972 is 720624 in julian day since year 0000 */
    jday -= JULIAN_DATE_LAST_DAY_1972;

    while ((is_leap_year(yr) && jday > 366) || (!(is_leap_year(yr)) && jday > 365))
    {
        status = is_leap_year(yr);
        if (status == true)
	{
            jday -= LEAP_YEAR_DAYS;
	}
        else
	{
            jday -= NON_LEAP_YEAR_DAYS;
	}
	yr++;
    }

    *year = yr;
    *doy = jday;
}

/******************************************************************************
MODULE:  matlab_norm

PURPOSE:  simulate matlab norm function 

RETURN VALUE:
Type = void
Value           Description
-----           -----------


HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
2/9/2015   Song Guo         Original Development

NOTES: 
******************************************************************************/
void matlab_norm
(
    float *array,        /* I: input array */
    int length,          /* I: number of input elements */
    float  *output_norm  /* O: output norm value */
)
{
    int i;
    float sum = 0.0;

    for (i = 0; i < length; i++)
    {
        sum += array[i] * array[i];
    }
    *output_norm = sqrt(sum);
}


