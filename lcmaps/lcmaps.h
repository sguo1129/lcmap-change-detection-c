typedef int bool;
#define true 1
#define false 0

typedef unsigned short uint16;
typedef unsigned char uint8;

#ifndef SUCCESS
    #define SUCCESS  0
#endif

#ifndef ERROR
    #define ERROR -1
#endif


#ifndef FAILURE
    #define FAILURE 1
#endif

#ifndef TRUE
    #define TRUE 1
#endif

#ifndef FALSE
    #define FALSE 0
#endif

#define MAX_STR_LEN 512
#define MAX_NUM_YEARS 40
#define NUM_BANDS 7
#define NUM_COEFFS 8

/* Structure for the metadata */
typedef struct {
    int lines;            /* number of lines in a scene */ 
    int samples;          /* number of samples in a scene */
    int data_type;        /* envi data type */
    int byte_order;       /* envi byte order */
    int utm_zone;         /* UTM zone; use a negative number if this is a
                             southern zone */
    int pixel_size;       /* pixel size */
    char interleave[MAX_STR_LEN];  /* envi save format */ 
    int  upper_left_x;    /* upper left x coordinates */                         
    int  upper_left_y;    /* upper left y coordinates */ 
    char map_info[MAX_STR_LEN];         /* map info */ 
    char projection_info[MAX_STR_LEN];  /* projection info */ 
    char coord_sys_str[MAX_STR_LEN];    /* coordinating system string */ 
} Input_meta_t;

typedef struct
{
    int t_start;           /* time when series model gets started */
    int t_end;             /* time when series model gets ended */
    int t_break;           /* time when the first break (change) is observed */
    float coefs[NUM_BANDS][NUM_COEFFS];
                           /*  coefficients for each time series model for each 
                               spectral band*/    
    float rmse[NUM_BANDS];
                           /*  RMSE for each time series model for each 
                               spectral band*/    
    int pos;               /* the pixel location of each time series model */
    float change_prob;     /* the probability of a pixel that have undergone 
                              change (between 0 and 100) */
    int num_obs;           /* the number of "good" observations used for model 
                              estimation */
    int category;          /* the quality of the model estimation (what model 
                              is used, what process is used) 
                              1x: persistent snow    2x: persistent water 
                              3x: Fmask fails        4x: normal precedure
                              x1: mean value (1)     x4: simple fit (4)
                              x6: basic fit (6)      x8: full fit (8) */
    float magnitude[NUM_BANDS];/* the magnitude of change (difference between model 
                                  prediction and observation for each spectral band)*/
    int class[MAX_NUM_YEARS];     /* land cover classification type */
    int class_qa[MAX_NUM_YEARS];/* use unsupervised ensemble margin as QA */
} Output_t;

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
);

int read_envi_header
(
    char *filename,        /* I:header filename with full path*/
    Input_meta_t *meta     /* O: saved header file info */
);

int convert_year_mon_day_to_jday_from_0000
(
    int year,      /* I: year */
    int month,     /* I: month of the year */
    int day,       /* I: day of the year */
    int *jday      /* O: julian date since year 0000 */
);

void convert_jday_from_0000_to_year_doy
(
    int jday,      /* I: julian date since year 0000 */
    int *year,     /* O: year */
    int *doy       /* O: month of the year */
);

void matlab_norm
(
    float *array,        /* I: input array */
    int length,          /* I: number of input elements */
    float  *output_norm  /* O: output norm value */
);

int write_envi_header
(
    char *filename,        /* I:header filename with full path*/
    int start_year,        /* I: start year of image */
    int end_year,          /* I: end year of image */
    int df,                /* I: envi data type value */ 
    Input_meta_t *meta     /* I: mata data info */
);

void usage ();

