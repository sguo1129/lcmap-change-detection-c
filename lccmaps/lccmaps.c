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
#include "lccmaps.h"

#define NUM_COEFFS 8
#define NUM_BANDS 7
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

/******************************************************************************
METHOD:  lccmaps

PURPOSE:  the main routine generating land cover change and associate maps in C

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
 
NOTES: type ./lccmaps --help for information to run the code

The current category in output structure: 
first digit:          
0: normal model (no change)     
1: change at the beginning of time series model                  
2: change at the end of time series model                    
3: disturbance change in the middle               
4: fmask fail scenario
5: permanent snow scenario
second digit:
1: model has only constant term
4: model has 3 coefs + 1 const 
6: model has 5 coefs + 1 const
8: model has 7 coefs + 1 const

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
    float **condition_map;
    uint8 **qa_map;
    uint16 **change_map;
    float **change_mag_map;
    uint16 **number_map;
    int ncols= 5000;
    int start_year;
    int end_year;
    int all_yrs;
    int month;
    int day;
    int year;
    int *jul_d;
    int col;
    float b_start;
    float b_end;
    float r_start;
    float r_end;
    float n_start;
    float n_end;
    float evi_start;
    float evi_end;
    float evi_slope;
    int n_year;
    int n_doy;
    FILE *fp_qa_map;
    FILE *fp_change_map;
    FILE *fp_change_mag_map;
    FILE *fp_condition_map;
    FILE *fp_number_map;
    char fname_qa_map[MAX_STR_LEN];  
    char fname_change_map[MAX_STR_LEN];  
    char fname_change_mag_map[MAX_STR_LEN];  
    char fname_condition_map[MAX_STR_LEN];  
    char fname_number_map[MAX_STR_LEN];  
    int status;

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
    condition_map = (float **) allocate_2d_array (all_yrs, ncols,
                                         sizeof (float));
    if (condition_map == NULL)
    {
        RETURN_ERROR ("Allocating memory: condition_map", FUNC_NAME, FAILURE);
    }

    qa_map = (uint8 **) allocate_2d_array (all_yrs, ncols,
                                         sizeof (uint8));
    if (qa_map == NULL)
    {
        RETURN_ERROR ("Allocating memory: qa_map", FUNC_NAME, FAILURE);
    }

    change_map = (uint16 **) allocate_2d_array (all_yrs, ncols,
                                         sizeof (uint16));
    if (change_map == NULL)
    {
        RETURN_ERROR ("Allocating memory: change_map", FUNC_NAME, FAILURE);
    }

    change_mag_map = (float **) allocate_2d_array (all_yrs, ncols,
                                         sizeof (float));
    if (change_mag_map == NULL)
    {
        RETURN_ERROR ("Allocating memory: change_mag_map", FUNC_NAME, FAILURE);
    }

    number_map = (uint16 **) allocate_2d_array (ncols, all_yrs,
                                         sizeof (uint16));
    if (number_map == NULL)
    {
        RETURN_ERROR ("Allocating memory: number_map", FUNC_NAME, FAILURE);
    }

    jul_d = (int *)malloc(all_yrs * sizeof(int));
    if (jul_d == NULL)
    {
        RETURN_ERROR ("Allocating memory: jul_d", FUNC_NAME, FAILURE);
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
    sprintf(fname_qa_map, "%s/QAMap", out_path);
    fp_qa_map = fopen(fname_qa_map, "wb");
    if (fp_qa_map == NULL)
    {
        RETURN_ERROR ("Opening qa_map file\n", FUNC_NAME, FAILURE);
    }

    sprintf(fname_change_map, "%s/ChangeMap", out_path);
    fp_change_map = fopen(fname_change_map, "wb");
    if (fp_change_map == NULL)
    {
        RETURN_ERROR ("Opening change_map file\n", FUNC_NAME, FAILURE);
    }

    sprintf(fname_change_mag_map, "%s/ChangeMagMap", out_path);
    fp_change_mag_map = fopen(fname_change_mag_map, "wb");
    if (fp_change_mag_map == NULL)
    {
        RETURN_ERROR ("Opening change_mag_map file\n", FUNC_NAME, FAILURE);
    }

    sprintf(fname_condition_map, "%s/ConditionMap", out_path);
    fp_condition_map = fopen(fname_condition_map, "wb");
    if (fp_condition_map == NULL)
    {
        RETURN_ERROR ("Opening condition_map file\n", FUNC_NAME, FAILURE);
    }

    sprintf(fname_number_map, "%s/NumberMap", out_path);
    fp_number_map = fopen(fname_number_map, "wb");
    if (fp_number_map == NULL)
    {
        RETURN_ERROR ("Opening number_map file\n", FUNC_NAME, FAILURE);
    }

    /* Loop through all lines of the ARD data */
    for (line = 0; line < meta->lines; line++)
    {
        for (i = 0; i < all_yrs; i++)
        {
            for (j = 0; j < ncols; j++)
	    {
	        condition_map[i][j] = 9999.0;
                qa_map[i][j] = 255;
	        change_map[i][j] = 9999;
                change_mag_map[i][j] = 9999.0;
		number_map[i][j] = 9999;
	    }
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
	//        printf("num_fc=%d\n",num_fc);

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
            field = Mat_VarGetStructField(mat_struct,fieldnames[0],MAT_BY_NAME,k);
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
			    rec_cg[k].t_start = (int)(*(double *)(data+idx*stride));
                        }
                    }
	        }
	    }

            field = Mat_VarGetStructField(mat_struct,fieldnames[1],MAT_BY_NAME,k);
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
			    rec_cg[k].t_end = (int)(*(double *)(data+idx*stride));
                        }
                    }
	        }
	    }

            field = Mat_VarGetStructField(mat_struct,fieldnames[2],MAT_BY_NAME,k);
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
			    rec_cg[k].t_break = (int)(*(double *)(data+idx*stride));
                        }
                    }
	        }
	    }

            field = Mat_VarGetStructField(mat_struct,fieldnames[3],MAT_BY_NAME,k);
            //            Mat_VarPrint(field,1);
            if (field != NULL)
	    {
	        if (field->rank == 2)
	        {
                    size_t stride = Mat_SizeOf(field->data_type);
                    char *data = field->data;
                    for ( i = 0; i < field->dims[0] && i < NUM_BANDS; i++ ) 
                    {
                        for ( j = 0; j < field->dims[1] && j < NUM_COEFFS; j++ ) 
                        {
                            size_t idx = field->dims[0]*j+i;
			    rec_cg[k].coefs[j][i] = (float)(*(double *)(data+idx*stride));
                        }
                    }
	        }
	    }

            field = Mat_VarGetStructField(mat_struct,fieldnames[4],MAT_BY_NAME,k);
            if (field != NULL)
            {
	        if (field->rank == 2)
	        {
                    size_t stride = Mat_SizeOf(field->data_type);
                    char *data = field->data;
                    for ( i = 0; i < field->dims[0] && i < NUM_BANDS; i++ ) 
                    {
                        for ( j = 0; j < field->dims[1] && j < 1; j++ ) 
                        {
                            size_t idx = field->dims[0]*j+i;
			    rec_cg[k].rmse[i] = (float)(*(double *)(data+idx*stride));
                        }
                    }
	        }
	    }

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

            field = Mat_VarGetStructField(mat_struct,fieldnames[6],MAT_BY_NAME,k);
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
			    rec_cg[k].change_prob = (float)(*(double *)(data+idx*stride));
                        }
                    }
	        }
	    }

            field = Mat_VarGetStructField(mat_struct,fieldnames[7],MAT_BY_NAME,k);
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
			    rec_cg[k].num_obs = (int)(*(double *)(data+idx*stride));
                        }
                    }
	        }
	    }

            field = Mat_VarGetStructField(mat_struct,fieldnames[8],MAT_BY_NAME,k);
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
			    rec_cg[k].category = (int)(*(double *)(data+idx*stride));
                        }
                    }
	        }
	    }

            field = Mat_VarGetStructField(mat_struct,fieldnames[9],MAT_BY_NAME,k);
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
			    rec_cg[k].magnitude[i] = (float)(*(double *)(data+idx*stride));
                        }
                    }
	        }
	    }
	}
        Mat_VarFree(mat_struct);
        Mat_Close(matfp);

	/* initialize pixels have at least one model */
        for (j = 0; j < ncols; j++)
	{
	    int qa_map_counter = 0;
            for (year = 0; year < all_yrs; year++)
	    {
         	if (qa_map[year][j] == 255)
		    qa_map_counter++;
	    }	
	    if (qa_map_counter == all_yrs)
	    {
        	for (year = 0; year < all_yrs; year++)
	        {
		    condition_map[year][j] = 0.0;
		    qa_map[year][j] = 0;
		    change_map[year][j] = 0;
		    change_mag_map[year][j] = 0.0;
		    number_map[year][j] = 0;
		}
       	    }
	}

	for (i = 0; i < num_fc; i++)
	{
	    // col1 = (int)(((float)rec_cg[i].pos / (float)ncols - rec_cg[i].pos / ncols) * ncols);
	    col = rec_cg[i].pos - line * ncols - 1; 
            /* Get start and end EVI values and EVI slope */ 
            b_start = rec_cg[i].coefs[0][0] + rec_cg[i].t_start * rec_cg[i].coefs[0][1];
	    r_start = rec_cg[i].coefs[2][0] + rec_cg[i].t_start * rec_cg[i].coefs[2][1];
	    n_start = rec_cg[i].coefs[3][0] + rec_cg[i].t_start * rec_cg[i].coefs[3][1];
	    b_end = rec_cg[i].coefs[0][0] + rec_cg[i].t_end * rec_cg[i].coefs[0][1];
	    r_end = rec_cg[i].coefs[2][0] + rec_cg[i].t_end * rec_cg[i].coefs[2][1];
	    n_end = rec_cg[i].coefs[3][0] + rec_cg[i].t_end * rec_cg[i].coefs[3][1];

            evi_start = 2.5*(n_start - r_start)/(n_start + 6.0*r_start - 7.5*b_start + 10000.0);
            evi_end = 2.5*(n_end - r_end)/(n_end + 6.0*r_end - 7.5*b_end + 10000.0);
	    evi_slope = 10000.0 * (evi_end - evi_start) / (rec_cg[i].t_end - rec_cg[i].t_start);

	    /* year (band) the curve belongs to */
	    for (year = 0; year < all_yrs; year++)
	    {
                if (jul_d[year] >= rec_cg[i].t_start && (jul_d[year] <= rec_cg[i].t_end 
                    || jul_d[year] < rec_cg[i].t_break))
		{
         	    /* write EVI slope to ConditionMap */
  		    condition_map[year][col] = evi_slope;
		    qa_map[year][col] = rec_cg[i].category;
		    number_map[year][col] = rec_cg[i].num_obs;
		}

	        if (rec_cg[i].change_prob == 1.0)
	        {
                    convert_jday_from_0000_to_year_doy(rec_cg[i].t_break, &n_year, &n_doy);
        	    if (year + start_year == n_year)
		    {
  		        change_map[year][col] = n_doy;
		        matlab_norm(rec_cg[i].magnitude, NUM_BANDS, &change_mag_map[year][col]);
		    }
		}
	    }

	}

        /* Write out one line of output map files */
        for (year = 0; year < all_yrs; year++)
        {
            fwrite(condition_map[year], sizeof(float), ncols, fp_condition_map);
            fwrite(number_map[year], sizeof(uint16), ncols, fp_number_map);
            fwrite(qa_map[year], sizeof(uint8), ncols, fp_qa_map);
            fwrite(change_map[year], sizeof(uint16), ncols, fp_change_map);
            fwrite(change_mag_map[year], sizeof(float), ncols, fp_change_mag_map);

#if 0
   	        for (k = 0; k < ncols; k++)
		    printf("year,k,condition_map[year][col],number_map[year][col],qa_map[year][col],"
                       "change_map[year][col],change_mag_map[year][col]=%d,%d,%f,%d,%d,%d,%f\n",
		        year,k,condition_map[year][col],number_map[year][col],qa_map[year][col],
                        change_map[year][col],change_mag_map[year][col]);
#endif
        }
	free(rec_cg);
    }

    /* Close all output map files */
    fclose(fp_qa_map);
    fclose(fp_change_map);
    fclose(fp_change_mag_map);
    fclose(fp_condition_map);
    fclose(fp_number_map);

    /* Create and write out all output envi header files */
    sprintf(filename, "%s/QAMap.hdr", out_path);
    status = write_envi_header(filename, start_year, end_year, 1, meta);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Calling read_envi_header for QAMap.hdr", 
                      FUNC_NAME, FAILURE);
    }
    sprintf(filename, "%s/ChangeMap.hdr", out_path);
    status = write_envi_header(filename, start_year, end_year, 12, meta);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Calling read_envi_header for ChangeMap.hdr", 
                      FUNC_NAME, FAILURE);
    }
    sprintf(filename, "%s/ChangeMagMap.hdr", out_path);
    status = write_envi_header(filename, start_year, end_year, 4, meta);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Calling read_envi_header for ChangeMagMap.hdr", 
                      FUNC_NAME, FAILURE);
    }
    sprintf(filename, "%s/NumberMap.hdr", out_path);
    status = write_envi_header(filename, start_year, end_year, 12, meta);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Calling read_envi_header for NumberMap.hdr", 
                      FUNC_NAME, FAILURE);
    }
    sprintf(filename, "%s/ConditionMap.hdr", out_path);
    status = write_envi_header(filename, start_year, end_year, 4, meta);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Calling read_envi_header for QualityMap.hdr", 
                      FUNC_NAME, FAILURE);
    }

    /* Free memories outside line loop */
    free(meta);
    free(jul_d);
    status = free_2d_array ((void **) condition_map);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: condition_map\n", FUNC_NAME,
                      FAILURE);
    }
    status = free_2d_array ((void **) qa_map);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: qa_map\n", FUNC_NAME,
                      FAILURE);
    }
    status = free_2d_array ((void **) change_map);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: change_map\n", FUNC_NAME,
                      FAILURE);
    }
    status = free_2d_array ((void **) change_mag_map);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: change_mag_map\n", FUNC_NAME,
                      FAILURE);
    }
    status = free_2d_array ((void **) number_map);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: number_map\n", FUNC_NAME,
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
    printf ("lccmaps"
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
    printf ("lccmaps --help will print the usage statement\n");
    printf ("\n");
    printf ("Example:\n");
    printf ("lccmaps"
            " --start_year=1985"
            " --month=7"
            " --day=1"
            " --end_year=2014"
            " --in_path=/data2/sguo/CCDC/ARDMatlab/grid08"
            " --out_path=/data2/sguo/CCDC/ARDMatlab/grid08"
            " --scene_name=LC80440272013106LGN01"
            " --verbose\n");
    printf ("Note: Default is ccdc must run from the directory"
            " where the input data are located.\n\n");
}
