// 3d_visualization.h

#include <stdio.h>
#include <GL/glut.h>
#include <limits.h> 
#include <float.h>
#include <string.h>
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <math.h>

#define TRUE 1
#define FALSE 0

//#define SWITCH_PPM_IMAGE

// visualization parameters
#define WINDOW_WIDTH 500 //640
#define WINDOW_HEIGHT 500 //480
#define POINT_SIZE 10
#define SPHERE_SIZE 0.015
#define DEFAULT_PLAY_SPEED 0.04

// simulation setting
#define N_REGIONS 1
#define R_REGION 0
#define N_NEURON_TYPES_M1 18
#define N_NEURON_TYPES_M1_THALAMUS 4
#define N_NEURON_TYPES_CEREBELLUM 8
#define N_NEURON_TYPES_S1 14
#define N_SPIKES 4000000
#define N_CONNECTIONS_PER_CELL 10000
#define N_MAX_CELL 40000


unsigned int n_proc = 4;

unsigned int n_subregions[N_REGIONS] = { 5, 2, 7};
unsigned int n_neuron_types[N_REGIONS][7] =
  {
    {2, 3, 5, 5, 3},
    {3, 1},
    {1, 1, 1, 2, 1, 1, 1},
  };

// visualization mode 
// 0: none, 1: full cell body, 2: gray cell body
#define N_TYPE_SHOW_CELL_LOCATION 3
// 0: all, 1:all cell in one region, 2:one cell in one region
#define N_TYPE_SWITCH_SHOW_DEGREE 3

#define SHOW_CELL_LOCATION_FULL 1
#define SHOW_CELL_LOCATION_GRAY 2

// vertex and edge array for frame drawing
GLdouble vertex[8][3];

int edge[][2] = {
  { 0, 1 },
  { 1, 2 },
  { 2, 3 },
  { 3, 0 },
  { 4, 5 },
  { 5, 6 },
  { 6, 7 },
  { 7, 4 },
  { 0, 4 },
  { 1, 5 },
  { 2, 6 },
  { 3, 7 }
};

struct region
{
  struct subregion *sr;
};

struct subregion
{
  struct neuron_type *nt;
};

struct neuron_type
{
  float **neuron_positions;
  float *spike_times;
  unsigned int *Idx;
  unsigned int n_spikes;
  unsigned int n_neurons;
  unsigned int count_spikes;
};

// light information
float light0_position[4]={0., 0.0, 0.0, 1.0};
float light1_position[4]={0., 0.0, 0.0, 1.0};
float light0_diffuse[4]={1.0, 1.0, 1.0, 0.8};
float light0_ambient[4]={0.5, 0.5, 0.5, 0.5};

// RGB color parameters for each cell type
GLfloat neuron_type_colors[N_REGIONS][N_NEURON_TYPES_M1][4] = {
  {
    { 0.5,  0.5,  0.5,  1. }, //1
    { 1.,  1.,  1.,  1. }, //1
    { 1.,  0.,  0.,  1. }, //2/3
    { 0.,  1.,  0.,  1. }, //2/3
    { 0.,  0.,  1.,  1. }, //2/3
    { 1.,  0.5,  0.,  1. }, //5a
    { 1.,  0.5,  0.,  1. }, //5a
    { 1.,  0.5,  0.,  1. }, //5a
    { 0.,  0.,  1.,  1. }, //5a
    { 0.,  1.,  0.,  1. }, //5a
    { 1.,  1.,  0.,  1. }, //5b
    { 1.,  1.,  0.,  1. },  //5b
    { 1.,  1.,  0.,  1. }, //5b
    { 0.,  0.,  1.,  1. }, //5b
    { 0.,  1.,  0.,  1. },  //5b
    { 1.,  0.,  1.,  1. }, //6
    { 0. , 0. , 1.,  1. }, //6
    { 0. , 1. , 0.,  1. }, //6
  },
  {
    { 1.,  0.,  0.,  1. },
    { 0.,  1.,  0.,  1. },
    { 0.,  0.,  1.,  1. },
    { 1.,  0.,  1.,  1. },
    { 0.,  0.,  0.,  0. },
    { 0.,  0.,  0.,  0. },
    { 0.,  0.,  0.,  0. },
    { 0.,  0.,  0.,  0. },
    { 0.,  0.,  0.,  0. },
    { 0.,  0.,  0.,  0. },
    { 0.,  0.,  0.,  0. },
    { 0.,  0.,  0.,  0. },
  },
  {
    { 1.,  0.,  0.,  1. },
    { 0.,  1.,  0.,  1. },
    { 0.,  0.,  1.,  1. },
    { 1.,  0.,  1.,  1. },
    { 1.,  0.,  0.,  1. },
    { 0.,  1.,  0.,  1. },
    { 0.,  0.,  1.,  1. },
    { 1.,  0.,  1.,  1. },
  }
};
  
  // labels for reginons and cells
/*
const char regions[N_REGIONS][256] = {"M1", "M1_Th", "cerebellum"}; 
const char neuron_types[N_REGIONS][N_NEURON_TYPES_M1][512] = { 
  {"L1_SBC", "L1_ENGC", "L2_3_CC", "L2_3_FS", "L2_3_LTS", "L5A_CS", "L5A_CC", "L5A_CT", "L5A_FS", "L5A_LTS", "L5B_PT", "L5B_CS", "L5B_CC", "L5B_FS", "L5B_LTS", "L6_CT", "L6_FS", "L6_LTS"}, 
  {"TN_HT", "TN_TC", "TN_IN", "TRN_RE"},
  {"molecular_layer_upper_ST", "molecular_layer_deep_BA", "PurkinjeCell_layer_PC", "granular_layer_GR", "granular_layer_GO", "deep_cerebellar_nuclei_DCN_E", "IO_CF", "Pons_MF" }
};
*/
const char regions[N_REGIONS][256] = {"cerebellum", "M1", "M1_Th"}; 
const char neuron_types[N_REGIONS][N_NEURON_TYPES_M1][512] = { 
  {"molecular_layer_upper_ST", "molecular_layer_deep_BA", "PurkinjeCell_layer_PC", "granular_layer_GR", "granular_layer_GO", "deep_cerebellar_nuclei_DCN_E", "IO_CF", "Pons_MF" },
  {"L1_SBC", "L1_ENGC", "L2_3_CC", "L2_3_FS", "L2_3_LTS", "L5A_CS", "L5A_CC", "L5A_CT", "L5A_FS", "L5A_LTS", "L5B_PT", "L5B_CS", "L5B_CC", "L5B_FS", "L5B_LTS", "L6_CT", "L6_FS", "L6_LTS"}, 
  {"TN_HT", "TN_TC", "TN_IN", "TRN_RE"},
};

unsigned int n_neuron_types_per_region[N_REGIONS] = {N_NEURON_TYPES_M1, N_NEURON_TYPES_M1_THALAMUS, N_NEURON_TYPES_S1};

FILE *fp, *fp_output;

double play_speed = (double)DEFAULT_PLAY_SPEED;

char data_directory_name[255];
char job_index[255];

// scalebar length: 500 micron
const double scalebar_length = 0.5;

// variables for view (not used now)
double xrot = 0., zrot = 0., xtrans = 0., ztrans = 0.;
double drag_L_rotate = 0., drag_L_zoom =0.;

double previous_x_drag_L = 0.0, previous_y_drag_L = 0.0;
double previous_x_drag_R = 0.0, previous_y_drag_R = 0.0;

// camera, center, and upper direction
const double default_camera_position[3] = {-4., -4., -1.5};
double camera_position[3] = {-4., -4., -1.5};
double center_position[3] = {1.6/2., 1.6/2., 2.4/2.};
double upper_direction[3] = {0., 0., -1.,};

//mouse button
unsigned int left_button = 0, right_button = 0;

// simulation time
unsigned int t = 0; // time

// memory consumption
float total_memory = 0.;

// switches
unsigned int switch_show_spikes = TRUE;
unsigned int switch_show_connections = TRUE;
//0: all, 1:all cell in one region, 2:one cell in one region
unsigned int switch_show_degree = 0;
unsigned int switch_show_cell_location = FALSE;
unsigned int switch_show_frame = TRUE;
unsigned int switch_mouse_mode = TRUE;
unsigned int switch_save_ppm_image = FALSE;

unsigned int i_shown_neuron_type = 0;
unsigned int i_shown_region = 1;

struct data
{
  double length_on_a_side;
  unsigned int simulation_time;
  double subregion_thicknesses[5];
  double subregion_z_upper_limits[5];
  double M1_thalamus_subregion_thicknesses[2];
  double M1_thalamus_subregion_z_upper_limits[2];
  double S1_subregion_thicknesses[7];
  double S1_subregion_z_upper_limits[7];
  
  double model_center_position[3];
  double model_opposing_corner_position[3];

  //double positions[N_REGIONS][N_NEURON_TYPES_M1][N_MAX_CELL][4];
  unsigned int n_cell[N_REGIONS][N_NEURON_TYPES_M1];
  double spike_times[N_REGIONS][N_NEURON_TYPES_M1][N_SPIKES];
  unsigned int count_spikes[N_REGIONS][N_NEURON_TYPES_M1];
  unsigned int n_connections[N_REGIONS][N_NEURON_TYPES_M1][N_REGIONS][N_NEURON_TYPES_M1];
  unsigned int connections[N_REGIONS][N_NEURON_TYPES_M1][N_REGIONS][N_NEURON_TYPES_M1][N_CONNECTIONS_PER_CELL][2];
  unsigned int running_numbers[N_REGIONS][N_NEURON_TYPES_M1][4000000];

  struct region *rg;

};

struct data *sd;

//function declarization
void draw_subregion_frame();
void DrawString(const char *str,void *font,float x,float y,float z);
void get_value_with_key_from_SLI_file(char * line);
void draw_scalebar();
void SaveImage_PPM(char* fname);
void get_model_position(struct data *sd);
void prepare_data_structure(struct data **sd);
void get_neuron_positions(struct data *sd);
void get_spike_times(struct data *sd);
void get_connection_patterns(struct data *sd);
void init(void);
void DrawString(const char *str,void *font,float x,float y,float z);
void special_input(int key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void SaveImage_PPM(char* fname);
void get_arguments(int argc, char *argv[]);
void *malloc2d(size_t size, unsigned int n_row, unsigned int n_col);
