#include "lcmaps.h"
#include "utilities.h"

/******************************************************************************
MODULE: write_envi_header

PURPOSE: Reads envi header info into output envi header files
 
RETURN VALUE:
Type = None

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
5/26/2016    Song Guo         Original development

NOTES:
*****************************************************************************/
int write_envi_header
(
    char *filename,        /* I:header filename with full path*/
    int start_year,        /* I: start year of image */
    int end_year,          /* I: end year of image */
    int df,                /* I: envi data type value */ 
    Input_meta_t *meta     /* I: mata data info */
)
{
    FILE *fp; 
    int year;
    int nbands;
    char msgstr[MAX_STR_LEN];
    char FUNC_NAME[] = "write_envi_header"; /* function name */

    nbands = end_year - start_year + 1;
    /* Generate output map envi header files */
    fp = fopen(filename, "w");
    if (fp == NULL )
    {
        sprintf (msgstr, "Generate %s envi header \n", filename);
        RETURN_ERROR (msgstr, FUNC_NAME, FAILURE);
    }

    /* Write the header to the file */
    fprintf(fp,
            "ENVI\n"
            "description = {Landsat Scientific Data}\n"
            "samples = %d\n"
            "lines   = %d\n"
            "bands   = %d\n"
            "header offset = 0\n"
            "file type = ENVI Standard\n"
            "sensor type = Landsat\n"
            "data type = %d\n"
            "interleave = bil\n"
            "byte order = 0\n", meta->samples, meta->lines, nbands, df);

    fprintf(fp, "band names = {");
    for (year = start_year; year <= end_year; year++)
    {
        if (year != end_year)
            fprintf(fp, "%d,", year);
	else
            fprintf(fp, "%d", year);
    }
    fprintf(fp, "}\n");

    fprintf(fp,
            "map_info = %s\n"
            "projection_info = %s\n"
            "coroedinate system string = %s\n",meta->map_info,meta->projection_info,meta->coord_sys_str);


    /* Close all output map header files */
    fclose(fp);

    return SUCCESS;
}

