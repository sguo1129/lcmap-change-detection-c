#include "lcmaps.h"
#include "utilities.h"
#include <stdio.h>
#include <stdlib.h>
#include "string.h"
#include "ctype.h"

/******************************************************************************
MODULE: trimwhitespace

PURPOSE: Trim leading spaces of a sting
 
RETURN VALUE:
Type = string without trailing space

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
1/16/2015    Song Guo         Modified from online code

NOTES:
*****************************************************************************/
char *trimwhitespace(char *str)
{
  char *end;

  /* Trim leading space */
  while(isspace(*str)) str++;

  if(*str == 0)  
    return str;

  /* Trim trailing space */
  end = str + strlen(str) - 1;
  while(end > str && isspace(*end)) end--;

  /* Write new null terminator */
  *(end+1) = 0;

  return str;
}

/******************************************************************************
MODULE: read_envi_header

PURPOSE: Reads envi header info into input meta structure
 
RETURN VALUE:
Type = None

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
1/16/2015    Song Guo         Original development

NOTES:
*****************************************************************************/
int read_envi_header
(
    char *filename,        /* I:header filename with full path*/
    Input_meta_t *meta     /* O: saved header file info */
)
{
    char  buffer[MAX_STR_LEN] = "\0";
    char  *label = NULL;
    char  *tokenptr = NULL;
    char  *tokenptr2 = NULL;
    char  *seperator = "=";
    char  *seperator2 = ",";
    FILE *in;
    int ib;
    char map_info[10][MAX_STR_LEN]; 
    char FUNC_NAME[] = "read_envi_header"; /* function name */
#if 0
    if (landsat_number == 8)
        sprintf(filename, "%s_sr_band2.hdr", scene_name);
    else 
        sprintf(filename, "%s_sr_band1.hdr", scene_name);
    //    sprintf(filename, "%s.hdr", scene_name);
    printf("scene_name,filename=%s,%s\n",scene_name,filename);
#endif

    in=fopen(filename, "r");
    if (in == NULL)
    {
        RETURN_ERROR ("opening header file", FUNC_NAME, FAILURE);
    }

    /* process line by line */
    while(fgets(buffer, MAX_STR_LEN, in) != NULL) 
    {

        char *s;
        s = strchr(buffer, '=');
        if (s != NULL)
        {
            /* get string token */
            tokenptr = strtok(buffer, seperator);
            label=trimwhitespace(tokenptr);

            if (strcmp(label,"lines") == 0)
            {
                tokenptr = trimwhitespace(strtok(NULL, seperator));
                meta->lines = atoi(tokenptr);
            }

            if (strcmp(label,"data type") == 0)
            {
                tokenptr = trimwhitespace(strtok(NULL, seperator));
                meta->data_type = atoi(tokenptr);
            }

            if (strcmp(label,"byte order") == 0)
            {
                tokenptr = trimwhitespace(strtok(NULL, seperator));
                meta->byte_order = atoi(tokenptr);
            }

            if (strcmp(label,"samples") == 0)
            {
                tokenptr = trimwhitespace(strtok(NULL, seperator));
                meta->samples = atoi(tokenptr);
            }

            if (strcmp(label,"interleave") == 0)
            {
                tokenptr = trimwhitespace(strtok(NULL, seperator));
                strcpy(meta->interleave, tokenptr);
            }

            if (strcmp(label,"UPPER_LEFT_CORNER") == 0)
            {
                tokenptr = trimwhitespace(strtok(NULL, seperator));
            }

            if (strcmp(label,"map info") == 0)
            {
                tokenptr = trimwhitespace(strtok(NULL, seperator));
                strcpy(meta->map_info, tokenptr);
            }

            if (strcmp(label,"map info") == 0)
            {
                tokenptr2 = strtok(tokenptr, seperator2);
                ib = 0;
                while(tokenptr2 != NULL)
                {
                    strcpy(map_info[ib], tokenptr2);
                    if (ib == 3)
                        meta->upper_left_x = atoi(map_info[ib]);
                    if (ib == 4)
                        meta->upper_left_y = atoi(map_info[ib]);
                    if (ib == 5)
                        meta->pixel_size = atoi(map_info[ib]);
                    if(ib == 7)
                        meta->utm_zone = atoi(map_info[ib]);
                    tokenptr2 = strtok(NULL, seperator2);
                    ib++;
                }
            }

            if (strcmp(label,"projection info") == 0)
            {
                tokenptr = trimwhitespace(strtok(NULL, seperator));
                strcpy(meta->projection_info, tokenptr);
            }

            if (strcmp(label,"coordinate system string") == 0)
            {
                tokenptr = trimwhitespace(strtok(NULL, seperator));
                strcpy(meta->coord_sys_str, tokenptr);
            }
        }
    }
    fclose(in);

    return SUCCESS;
}
