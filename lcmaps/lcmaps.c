#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/timeb.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>

#include "utilities.h"
#include "matio.h"
#include "2d_array.h"
#include "lcmaps.h"

#define NUM_COEFFS 8
#define NUM_BANDS 7
#define NUM_CLASSES 12
#define NUM_LASSO_BANDS 5
#define TOTAL_IMAGE_BANDS 7
#define TOTAL_BANDS 8
#define MIN_NUM_C 4
#define MID_NUM_C 6
#define MAX_NUM_C 8
#define CONSE 6
#define MAX_YEARS 40
#define N_TIMES 3 /* number of clear observations/coefficients*/
#define NUM_YEARS 365.25 /* average number of days per year */
#define T_CONST 4.89 /* Threshold for cloud, shadow, and snow detection */
#define MIN_YEARS 1 /* minimum year for model intialization */
#define T_SN 0.75        /* no change detection for permanent snow pixels */ 
#define T_CLR 0.25       /* Fmask fails threshold */
#define T_CG 15.0863     /* chi-square inversed T_cg (0.99) for noise removal */
#define T_MAX_CG 35.8882 /* chi-square inversed T_max_cg (1e-6) for 
                            last step noise removal */
#define CFMASK_CLEAR 0
#define CFMASK_WATER 1
#define CFMASK_CLOUD 2
#define CFMASK_SNOW 3
#define CFMASK_SHADOW 4 
#define CFMASK_FILL  255 
#define IMAGE_FILL -9999

char *sub_string
(
    const char *source,
    size_t start,
    size_t length
) 
{
    int i;
    char *target;

    target = malloc(length*sizeof(char));

    for(i = 0; i != length; ++i) 
    {
        target[i] = source[start + i];
    }
    target[i] = 0;
    return target;
}

int stable_class[5] = {1, 2, 5, 8, 11}; 
/******************************************************************************
METHOD:  lcmaps

PURPOSE:  the main routine generating land cover map and its QA map in C

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           An error occurred during processing of the ccdc
SUCCESS         Processing was successful

PROJECT:  Land Change Monitoring, Assessment and Projection (LCMAP) Project

HISTORY:
Date        Programmer       Reason
--------    ---------------  ----------------------------------------------
4/15/2016   Song Guo         Original Development
 
NOTES: type ./lcmaps --help for information to run the code

Note: The commented out parts of code is the inputs using ESPA putputs,
      the current input part is only for reading inputs from Zhe's code
******************************************************************************/
int main (int argc, char *argv[])
{
    char FUNC_NAME[] = "main";    /* for printing error messages          */
    char msg_str[MAX_STR_LEN];    /* input data scene name                */
    Output_t *rec_cg = NULL;      /* output structure and metadata        */
    int i, j, k;                  /* loop counters                        */
    char in_path[MAX_STR_LEN];    /* directory location of input data/files     */
    char out_path[MAX_STR_LEN];   /* directory location of output data/files     */
    char input_mat[MAX_STR_LEN];  /* directory and file name for output.bin */
    char scene_name[MAX_STR_LEN]; /* input envi header filename */
    char filename[MAX_STR_LEN];   /* input filenames with path              */
    char tmpstr[MAX_STR_LEN];     /* char string for text manipulation */
    char *short_scene;            /* char string for text manipulation */
    int num_fc;
    float line_pt = 0.0;
    int count;
    int zero_count;
    int label;
    int non_unique;
    int *class_count;
    int maj_class;
    int len;
    int line;
    int verbose;
    Input_meta_t *meta;
    mat_t *matfp;
    matvar_t* mat_struct;
    char *structname = "rec_cg";
    void *fieldnames[12] = { "t_start","t_end","t_break","coefs","rmse","pos",
				   "change_prob","num_obs","category","magnitude","class","classQA"};
    matvar_t *field;
    uint8 **cover_map;
    uint8 **cover_qa_map;
    int ncols= 5000;
    int start_year;
    int end_year;
    int all_yrs;
    int month;
    int day;
    int year;
    int *jul_d;
    int col;
    int *n_band;
    int *non_zero_cover;
    int *zero_cover;
    FILE *fp_cover_map;
    FILE *fp_cover_qa_map;
    char fname_cover_map[MAX_STR_LEN];  
    char fname_cover_qa_map[MAX_STR_LEN];  
    int status;
    int sum_band;
    uint8 *year_class;
    int dummy_class;
    uint8 *year_class_qa;
    int non_stable_class[NUM_CLASSES];
    int start_class;
    int end_class;
    int t_id;
    int zero_start_year;

    time_t now;
    time (&now);
    snprintf (msg_str, sizeof(msg_str),
              "CCDC start_time=%s\n", ctime (&now));
    LOG_MESSAGE (msg_str, FUNC_NAME);
	
    /* because they are optional.                               */
    strcpy(in_path, "");
    strcpy(out_path, "");

    /* Read the command-line arguments */
    status = get_args (argc, argv, &start_year, &end_year, &month, &day, 
                       in_path, out_path, scene_name, &verbose);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("calling get_args", FUNC_NAME, EXIT_FAILURE);
    }
    all_yrs = end_year - start_year + 1;

    strcat(in_path, "/TSFitMap"); 
    strcat(out_path, "/CCDCMap"); 
    /* Create the Input metadata structure */
    meta = (Input_meta_t *)malloc(sizeof(Input_meta_t));
    if (meta == NULL) 
    {
        RETURN_ERROR("allocating Input data structure", FUNC_NAME, FAILURE);
    }

    len = strlen(scene_name);
    short_scene= strndup(scene_name, len-5);
    if (strncmp(short_scene, ".", 1) == 0)
    {
        strncpy(tmpstr, short_scene + 2, len - 2);
        sprintf(filename, "%s/%s_MTLstack.hdr", tmpstr, scene_name);
    }
    else
        sprintf(filename, "%s/%s_MTLstack.hdr", short_scene, scene_name);
    free(short_scene);

    printf("filename=%s\n",filename);
    status = read_envi_header(filename, meta);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Calling read_envi_header", 
                      FUNC_NAME, FAILURE);
    }

    if (verbose)
    {
        /* Print some info to show how the input metadata works */
        printf ("DEBUG: Number of input lines: %d\n", meta->lines);
        printf ("DEBUG: Number of input samples: %d\n", meta->samples);
        printf ("DEBUG: UL_MAP_CORNER: %d, %d\n", meta->upper_left_x,
                meta->upper_left_y);
        printf ("DEBUG: ENVI data type: %d\n", meta->data_type);
        printf ("DEBUG: ENVI byte order: %d\n", meta->byte_order);
        printf ("DEBUG: UTM zone number: %d\n", meta->utm_zone);
        printf ("DEBUG: Pixel size: %d\n", meta->pixel_size);
        printf ("DEBUG: Envi save format: %s\n", meta->interleave);
    }

    /* allocate memories */
    cover_map = (uint8 **) allocate_2d_array (all_yrs, ncols, sizeof (uint8));
    if (cover_map == NULL)
    {
        RETURN_ERROR ("Allocating memory: cover_map", FUNC_NAME, FAILURE);
    }

    cover_qa_map = (uint8 **) allocate_2d_array (all_yrs, ncols, sizeof (uint8));
    if (cover_qa_map == NULL)
    {
        RETURN_ERROR ("Allocating memory: int", FUNC_NAME, FAILURE);
    }

    jul_d = (int *)malloc(all_yrs * sizeof(int));
    if (jul_d == NULL)
    {
        RETURN_ERROR ("Allocating memory: jul_d", FUNC_NAME, FAILURE);
    }

    n_band = (int *)malloc(all_yrs * sizeof(int));
    if (n_band == NULL)
    {
        RETURN_ERROR ("Allocating memory: n_band", FUNC_NAME, FAILURE);
    }

    non_zero_cover = (int *)malloc(all_yrs * sizeof(int));
    if (non_zero_cover == NULL)
    {
        RETURN_ERROR ("Allocating memory: non_zero_cover", FUNC_NAME, FAILURE);
    }

    zero_cover = (int *)malloc(all_yrs * sizeof(int));
    if (zero_cover == NULL)
    {
        RETURN_ERROR ("Allocating memory: zero_cover", FUNC_NAME, FAILURE);
    }

    year_class = (uint8 *)malloc(all_yrs * sizeof(uint8));
    if (year_class == NULL)
    {
        RETURN_ERROR ("Allocating memory: year_class", FUNC_NAME, FAILURE);
    }

    year_class_qa = (uint8 *)malloc(all_yrs * sizeof(uint8));
    if (year_class_qa == NULL)
    {
        RETURN_ERROR ("Allocating memory: year_class_qa", FUNC_NAME, FAILURE);
    }

    class_count = (int *)malloc(NUM_CLASSES * sizeof(int));
    if (class_count == NULL)
    {
        RETURN_ERROR ("Allocating memory: class_count", FUNC_NAME, FAILURE);
    }

    for (i = 0; i < all_yrs; i++)
    {
        status = convert_year_mon_day_to_jday_from_0000(
                 i+start_year, month, day, &jul_d[i]);
	if (status != SUCCESS)
            RETURN_ERROR ("Calling convert_year_mon_day_to_jday_from_0000", 
                 FUNC_NAME, FAILURE);
    }

    printf("out_path=%s\n",out_path);
    /* Open the output map files */
    sprintf(fname_cover_map, "%s/CoverMap", out_path);
    fp_cover_map = fopen(fname_cover_map, "wb");
    if (fp_cover_map == NULL)
    {
        RETURN_ERROR ("Opening cover_map file\n", FUNC_NAME, FAILURE);
    }

    sprintf(fname_cover_qa_map, "%s/CoverQAMap", out_path);
    fp_cover_qa_map = fopen(fname_cover_qa_map, "wb");
    if (fp_cover_qa_map == NULL)
    {
        RETURN_ERROR ("Opening cover_qa_map file\n", FUNC_NAME, FAILURE);
    }

    /* Loop through all lines of the ARD data */
    for (line = 0; line < meta->lines; line++)
    {
        for (i = 0; i < all_yrs; i++)
        {
            for (j = 0; j < ncols; j++)
                cover_map[i][j] = 255;
                cover_qa_map[i][j] = 255;
        }

        if ((100.0*((float)line/(float)meta->lines) - line_pt) > 1.0)
	{
	    if (verbose)
                printf("Processing %f percent\n",ceil(100.0*((float)line/(float)meta->lines)));
            line_pt = 100.0*((float)line/(float)meta->lines);
	}

        sprintf(input_mat, "%s/record_change%d.mat", in_path, line + 1);
#if 0
        /* Check the existence of the TSFitMap directory, create one if not existed */
        if (!(stat(in_path, &sb) == 0 && S_ISDIR(sb.st_mode)))
        {
            status = mkdir(in_path, 0777);
            if (status != SUCCESS)
                RETURN_ERROR ("Creating TSFitMap directory", FUNC_NAME, FAILURE);
        }
    
        /* Check the existence of the input file */
        if (access(input_mat, F_OK) != 0) /* File does not exist */
            matfp = Mat_CreateVer(output_mat, NULL, MAT_FT_MAT5);
            //        matfp = Mat_CreateVer(output_mat, NULL, MAT_FT_MAT73);
        else
            matfp = Mat_Open(input_mat,MAT_ACC_RDWR);
        if (matfp == NULL) 
            RETURN_ERROR ("Opening MAT file", FUNC_NAME, FAILURE);
#endif
	//	printf("input_mat=%s\n",input_mat);
        matfp = Mat_Open(input_mat, MAT_ACC_RDONLY);
        if (matfp == NULL) 
            RETURN_ERROR ("Opening MAT file", FUNC_NAME, FAILURE);

        mat_struct = Mat_VarRead(matfp,(char*)structname);
        if (mat_struct == NULL)
        { 
            Mat_Close(matfp);
            RETURN_ERROR ("Reading MAT structure\n", FUNC_NAME, FAILURE);
        }
        num_fc = mat_struct->dims[1];
	//	printf("num_fc=%d\n",num_fc);

        /* Allocate memory for rec_cg */ 
        rec_cg = malloc(num_fc * sizeof(Output_t));
        if (rec_cg == NULL)
        {
            RETURN_ERROR("ERROR allocating rec_cg memory", FUNC_NAME, FAILURE);
        }

        /* Output rec_cg structure to the output file 
           note: can use fread to read out the structure from the output file */
        for (k = 0; k < num_fc; k++)
        { 	                    
            field = Mat_VarGetStructField(mat_struct,fieldnames[5],MAT_BY_NAME,k);
            if (field != NULL)
	    {
	        if (field->rank == 2)
	        {
                    size_t stride = Mat_SizeOf(field->data_type);
                    char *data = field->data;
                    for ( i = 0; i < field->dims[0] && i < 1; i++ ) 
                    {
                        for ( j = 0; j < field->dims[1] && j < 1; j++ ) 
                        {
                            size_t idx = field->dims[0]*j+i;
			    rec_cg[k].pos = (int)(*(double *)(data+idx*stride));
                        }
                    }
	        }
	    }

            field = Mat_VarGetStructField(mat_struct,fieldnames[10],MAT_BY_NAME,k);
            if (field != NULL)
            {
	      //	        Mat_VarPrint(field,1);
	      //	        printf("field->dims[0],field->dims[1]=%d,%d\n",field->dims[0],field->dims[1]);
	        if (field->rank == 2)
	        {
                    size_t stride = Mat_SizeOf(field->data_type);
                    char *data = field->data;
                    for ( i = 0; i < field->dims[0] && i < all_yrs; i++ ) 
                    {
                        for ( j = 0; j < field->dims[1] && j < 1; j++ ) 
                        {
                            size_t idx = field->dims[0]*j+i;
			    rec_cg[k].class[i] = (int)(*(double *)(data+idx*stride));
                        }
                    }
	        }
	    }

            field = Mat_VarGetStructField(mat_struct,fieldnames[11],MAT_BY_NAME,k);
            if (field != NULL)
            {
	        if (field->rank == 2)
	        {
                    size_t stride = Mat_SizeOf(field->data_type);
                    char *data = field->data;
                    for ( i = 0; i < field->dims[0] && i < all_yrs; i++ ) 
                    {
                        for ( j = 0; j < field->dims[1] && j < 1; j++ ) 
                        {
                            size_t idx = field->dims[0]*j+i;
			    rec_cg[k].class_qa[i] = (int)(*(double *)(data+idx*stride));
                        }
                    }
	        }
	    }
	}
        Mat_VarFree(mat_struct);
        Mat_Close(matfp);

        for (i = 0; i < num_fc; i++)
	{
            /* Put disturbed class (10) to grass/shrub (7) */
            for (year = 0; year < all_yrs; year++)
	    {
		  if (rec_cg[i].class[year] == 10)
		        rec_cg[i].class[year] = 7;
	    }

	    // col = (int)(((float)rec_cg[i].pos / (float)ncols - rec_cg[i].pos / ncols) * ncols);
	    col = rec_cg[i].pos - line * ncols - 1;
#if 0
            for (year = 0; year < all_yrs; year++)
		  printf("line,num_fc,col,year,rec_cg[i].class[year]=%d,%d,%d,%d,%d\n",
			 line,i,col,year,rec_cg[i].class[year]);
#endif
            /* initialize pixels have at least one model */
	    int cover_map_counter = 0;
            for (year = 0; year < all_yrs; year++)
	    {
         	if (cover_map[year][col] == 255)
		    cover_map_counter++;
	    }	
	    if (cover_map_counter == all_yrs)
	    {
        	for (year = 0; year < all_yrs; year++)
	        {
		    cover_map[year][col] = 0;
		    cover_qa_map[year][col] = 0;
		}
       	    }

            /* write only when there are more than one valid cover */
	    sum_band = 0;
            for (year = 0; year < all_yrs; year++)
	    {
		n_band[year] = 0;
	        if (rec_cg[i].class[year] > 0)
		{
		    n_band[year] = 1;
		    sum_band++;
		}
	    }

	    if (sum_band > 0)
	    {
                count = 0;
		for (year = 0; year < all_yrs; year++)
		{
    		    if (n_band[year] == 1)
		    {
         		year_class[count] = (uint8)rec_cg[i].class[year]; 
		        year_class_qa[count] = (uint8)rec_cg[i].class_qa[year]; 
			count++;
		    }
		}

		/* number of differenc class within a single curve */
		non_unique = 0;
		for (year = 0; year < count - 1; year++)
		{
		    for (k = year + 1; k < count; k++)
		    {
		        if (year_class[year] == year_class[k])
			{
		            non_unique = 1;
	      	            break;
			}
		    }
		    if (non_unique == 1)
			break;
		}

        	/* use majority class if it is water(1), developed(2), barren(5),
                   agriculture(8), or snow(11), as these cover will not change
                   without an abrupt change */
	        if (non_unique == 1)
	        {		  
	            for (k = 0; k < NUM_CLASSES; k++)
	            {
      		        class_count[k] = 0;
	            }

		    for (year = 0; year < count; year++)
		    {
             		if (year_class[year] != 0)
		            class_count[year_class[year]]++;
		    }

         	    int maj_class_count = 0;
	            for (k = 0; k < NUM_CLASSES; k++)
	            {
      		        if (class_count[k] > maj_class_count)
			{
		            maj_class_count = class_count[k];
                            maj_class = k;
			}
	            }
#if 0
		    printf("maj_class, maj_class_count=%d,%d\n",maj_class, maj_class_count);
#endif
		    label = 0;
        	    for ( k = 0; k < 5; k++)
		    {
         	        if (maj_class == stable_class[k])
		        {
			    label = 1;
			    break;
			}
		    }

		    if (label == 1)
		    {
		        for (year = 0; year < all_yrs; year++)
		        {
		            if (n_band[year] == 1)
			    {
  			        dummy_class = year_class[year];
			        year_class[year] = maj_class;
				/* give illogical QA (103) */
				if (year_class[year] != dummy_class)
				    year_class_qa[year] = 103.0;
			    }
		        }
		    }
		    else
		    {
		        count = 0;
		        for (year = 0; year < all_yrs; year++)
		        {
		            for ( k = 0; k < 5; k++)
		            {
			        if (year_class[year] != stable_class[k])
				{
    		                    if (n_band[year] == 1)
		                    {
				        non_stable_class[count] = year_class[year];
				        count++;
				        break;
				    }
				}
			    } 
			}
			start_class = non_stable_class[0];
			end_class = non_stable_class[count-1];

			if (start_class == end_class)
			{
			    /* no transition */
		            for (year = 0; year < all_yrs; year++)
			    {
    		                if (n_band[year] == 1)
		                {
  			            dummy_class = year_class[year];
				    year_class[year] = start_class;
				    /* give illogical QA (103) */
				    if (year_class[year] != dummy_class)
				        year_class_qa[year] = 103.0;
				}
			    }
			}
			else
		        {
			    /* add trasition: S/G<=>F; S/G <=>W; F<=>W 
                               transition id = first end class */
		            for (year = 0; year < all_yrs; year++)
		            {
    		                if (n_band[year] == 1)
		                {
			            if (year_class[year] != end_class)
				    {
  			                dummy_class = year_class[year];
				        year_class[year] = start_class;
				        /* give vegetation QA (104) */
				        if (year_class[year] != dummy_class)
				            year_class_qa[year] = 104.0;
				    }
  				    else
		                    {
				        t_id = year;
					break;
				    }
				}
			    }

			    for (year = t_id; year < all_yrs; year++)
			    {
    		                if (n_band[year] == 1)
		                {
  			            dummy_class = year_class[year];
			            year_class[year] = end_class;
				    /* give vegetation QA (104) */
				    if (year_class[year] != dummy_class)
				        year_class_qa[year] = 104.0;
				}
			    }
			}
		    }
		}
	    
	        /* write land cover to CoverMap */
	        for (year = 0; year < all_yrs; year++)
	        {
    		    if (n_band[year] == 1)
		    {
		        cover_map[year][col] = year_class[year];
		        cover_qa_map[year][col] = year_class_qa[year];
		    }
	        }

        	/* Fill in zeros with the next valid obs */
	        count = 0;
                for (year = 0; year < all_yrs; year++)
	        {
	            if (cover_map[year][col] > 0)
		    {
		        non_zero_cover[count] = year;
		        count++;
		    }
		}

		/* all ids of zeros before the last nozero value */
		zero_count = 0;
		for (year = 0; year <= non_zero_cover[count-1]; year++)
		{
		    if (cover_map[year][col] == 0)
		    {
		        zero_cover[zero_count] = year;
		        zero_count++;
		    }
	        }

        	/* there are zeros before the last nonzero value */
		if (zero_count > 0)
		{
                    for (year = 0; year <= zero_cover[zero_count-1]; year++)
	            {
	                if (cover_map[year][col] == 0)
		        {
  		            cover_map[year][col] = 
                                cover_map[zero_cover[zero_count-1]+1][col];
			    /* gap filled with next curve: 102 */
		            cover_qa_map[year][col] = 102;
			}
		    }
	        }

        	/* fill in zeros with the previous valid obs */
	        if (i != num_fc - 1)
	        {
                    if (rec_cg[i].pos != rec_cg[i+1].pos) 
		    {
		        int zero_count = 0;
                        for (year = 0; year < all_yrs; year++)
	                {
	                    if (cover_map[year][col] == 0)
			    {
         		        zero_start_year = year;
				zero_count++;
		                break;
			    }
		        }

       			if (zero_count > 0)
			{
                            for (year = 0; year < all_yrs; year++)
	                    {
	                        if (cover_map[year][col] == 0)
				{
  		                    cover_map[year][col] = cover_map[zero_start_year-1][col];
			            /* gap filled with previous curve: 101 */
		                    cover_qa_map[year][col] = 101;
				}
			    }
	                }
		    }
	        }
	        else
	        {
		    zero_count = 0;
                    for (year = 0; year < all_yrs; year++)
	            {
	                if (cover_map[year][col] == 0)
		        {
         		    zero_start_year = year;
		            zero_count++;
		        }
		    }

		    if (zero_count > 0)
		    {
                        for (year = 0; year < all_yrs; year++)
        	        {
	                    if (cover_map[year][col] == 0)
			    {
  	                        cover_map[year][col] = cover_map[non_zero_cover[0]-1][col];
			        /* gap filled with previous curve: 101 */
		                cover_qa_map[year][col] = 101;
			    }
			}
	            }
	        }
	    }	    
	}

        /* Write out one line of output map files */
        for (year = 0; year < all_yrs; year++)
        {
            fwrite(cover_map[year], sizeof(uint8), ncols, fp_cover_map);
            fwrite(cover_qa_map[year], sizeof(uint8), ncols, fp_cover_qa_map);
        }
	free(rec_cg);
    }

    /* Close all output map files */
    fclose(fp_cover_map);
    fclose(fp_cover_qa_map);

    /* Create and write out all output envi header files */
    sprintf(filename, "%s/CoverMap.hdr", out_path);
    status = write_envi_header(filename, start_year, end_year, 1, meta);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Calling read_envi_header for CoverMap.hdr", 
                      FUNC_NAME, FAILURE);
    }
    sprintf(filename, "%s/CoverQAMap.hdr", out_path);
    status = write_envi_header(filename, start_year, end_year, 1, meta);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Calling read_envi_header for CoverQAMap", 
                      FUNC_NAME, FAILURE);
    }

    /* Free memories outside line loop */
    free(meta);
    free(jul_d);
    free(n_band);
    free(non_zero_cover);
    free(zero_cover);
    free(year_class);
    free(year_class_qa);
    free(class_count);
    status = free_2d_array ((void **) cover_map);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: cover_map\n", FUNC_NAME,
                      FAILURE);
    }
    status = free_2d_array ((void **) cover_qa_map);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: cover_qa__map\n", FUNC_NAME,
                      FAILURE);
    }

    time (&now);
    snprintf (msg_str, sizeof(msg_str),
              "CCDC end_time=%s\n", ctime (&now));
    LOG_MESSAGE (msg_str, FUNC_NAME);

    return SUCCESS;
}

/******************************************************************************
MODULE:  usage

PURPOSE:  Prints the usage information for this application.

RETURN VALUE:
Type = None

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
2/15/2015   Song Guo         Original Development

******************************************************************************/
void
usage ()
{
    printf ("Continuous Change Detection and Classification\n");
    printf ("\n");
    printf ("usage:\n");
    printf ("lcmaps"
            " --start_year=<start of image year>"
            " --end_year=<end of image year>"
            " --month=<month for image in each year>"
            " --day=<day for image in each year>"
            " --in_path=<input directory>"
            " --out_path=<output directory>"
            " --scene_name=<one scene_name file>"
            " [--verbose]\n");

    printf ("\n");
    printf ("where the following parameters are required:\n");
    printf ("    --start_year=:start of image year\n");
    printf ("    --end_year=:end of image year\n");
    printf ("    --month=:month for image in each year\n");
    printf ("    --day=:day for image in each year\n");
    printf ("    --in_path=: input data directory location\n");
    printf ("    --out_path=: directory location for output files\n");
    printf ("    --scene_name=: one scene_name file\n");
    printf ("\n");
    printf ("and the following parameters are optional:\n");
    printf ("    -verbose: should intermediate messages be printed?"
            " (default is false)\n");
    printf ("\n");
    printf ("lcmaps --help will print the usage statement\n");
    printf ("\n");
    printf ("Example:\n");
    printf ("lcmaps"
            " --start_year=1985"
            " --month=7"
            " --day=1"
            " --end_year=2014"
            " --in_path=/data2/sguo/CCDC/ARDMatlab/grid07"
            " --out_path=/data2/sguo/CCDC/ARDMatlab/grid07"
            " --scene_name=LT50480272010286PAC02"
            " --verbose\n");
    printf ("Note: Default is ccdc must run from the directory"
            " where the input data are located.\n\n");
}
