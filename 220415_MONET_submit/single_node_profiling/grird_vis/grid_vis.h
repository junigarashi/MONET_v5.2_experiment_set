#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <malloc.h>
#include <stdbool.h>
#include <GL/glut.h>
#include <math.h>
#include <stdbool.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>
#include <getopt.h>

#define readline_length 4194304
#define file_name_character_length 256

// widow parameters
#define WINDOW_WIDTH 640
#define WINDOW_HEIGHT 480
#define POINT_SIZE 10
#define DEFAULT_PLAY_SPEED 0.04
// 0: all, 1:all cell in one region, 2:one cell in one region
#define N_TYPE_SWITCH_SHOW_DEGREE 3
//0:sg_nt , 1:sg , 2:process
#define N_TYPE_SWITCH_SPATIAL_GRANULARITY 3
//0: false, 1:countrclockwise, 2:clockwise
#define N_CAMERA_MOVING 3


// switch
extern bool switch_mouse_mode;
extern bool switch_save_ppm_image;

// light information
extern float light0_position[4];
extern float light0_ambient[4];

//mouse button
extern unsigned int left_button;
extern double previous_x_drag_L;
extern double previous_x_drag_R;

// camera
//extern double camera_position[3];
//extern double center_position[3];
//extern double upper_direction[3];
//extern double model_center_position[3];

extern double play_speed;
extern unsigned int t;

extern char data_directory_name[255];


struct system_parameters
{
  float dt;
  unsigned int n_steps_per_dt;
  unsigned int n_calculation_steps;
  float start_time, end_time;
  unsigned int n_min_delay;
  unsigned int n_max_delay;
  unsigned int communication_interval;
  float conductance_buffer_time_window;
  unsigned int n_conductance_buffer_time_window;
  float memory_per_node;
  unsigned long n_temp_connection_buffer;
  unsigned long n_temp_connection_buffer_omp;
  unsigned int n_total_processes;
  unsigned int n_total_regions;
  unsigned int n_OpenMP_threads;
  unsigned int n_threads;
  unsigned int PRNG_seed;
  unsigned int n_steps_for_spike_counts_per_time_bin;
  //unsigned int n_recording_time_points_for_spike_counts;
  unsigned int n_recording_time_points;
  char ** region_GID_to_region_name;
  char *** neuron_type_names;
  unsigned int * region_GID_to_n_neuron_types;
  unsigned int * region_GID_to_n_subregions;

  // cache settings
  float L1_size;
  unsigned int n_L1_cache_blocks;

  unsigned int n_cache_blocks;
  float L1_size_KB;
  float L2_size_MB;
  float DRAM_size_GB;
  float GFLOPS_per_processor;
  float memory_bandwidth_per_processor_GBs;

  // for measurement
  unsigned int count_intra_regional_communication;
  bool sw_record_calculate_time;
  bool sw_record_build_time;
  char record_calculate_time_process_rank[file_name_character_length];
  char record_build_time_process_rank[file_name_character_length];
  bool sw_record_neuron_positions_per_process;
  bool sw_record_spike_times_per_process;
  bool sw_record_spike_times_per_neuron_type;
  bool sw_record_neuron_positions_per_neuron_type;
  bool sw_record_neuron_state_varibles_per_neuron_type;
  bool sw_record_synaptic_conductance_buffers_per_neuron_type;
  bool sw_record_spike_info;
  bool sw_record_intra_regional_connection_positions;
  bool sw_record_intra_regional_connection_samples;
  bool sw_record_inter_regional_connection_positions;
  bool sw_record_inter_regional_connection_samples;
  bool sw_record_settings;
  bool sw_record_spike_time_communications;
  bool sw_all_time;
  bool sw_time_change_all;
  bool sw_time_change_region;
  bool sw_time_change_neuron_type;
  bool sw_time_change_region_per_tile;
  bool sw_time_change_neuron_type_per_tile;
  bool sw_time_change_region_per_subgrid_in_tile;
  bool sw_time_change_neuron_type_per_subgrid_in_tile;
  bool sw_time_change_region_per_cluster_in_tile;
  bool sw_time_change_neuron_type_per_cluster_in_tile;
  bool sw_FR_subgrid_neuron_type_per_process;
  bool sw_cache_block_manual_setting;
  bool sw_terminal_output;
  //float record_firing_rate_time_bin;
  
  unsigned int n_indegree_connections_per_post_synaptic_neuron;
  unsigned int n_outdegree_connections_per_pre_synaptic_neuron;
  unsigned int n_record_target_neuron_variables;
  unsigned int n_record_connection_weights;
  
  struct global_shared_info *gsi;

  // position
  float max_position[3];
  float min_position[3];
  
  bool switch_FR_sg_nt;
  bool switch_FR_sg;
  bool switch_FR_process;
  unsigned int switch_camera_moving;
  unsigned int switch_loop;
  unsigned int switch_show_degree;
  unsigned int switch_time_stop;
  unsigned int switch_spatial_granularity;

  unsigned int i_shown_region;
  unsigned int i_shown_neuron_type;
  char current_path[file_name_character_length];

  float camera_position[3];
  float view_point_position[3];
  float upper_direction[3];
  float model_center_position[3];

};

struct FR_with_time
{
  unsigned int *time_steps;
  float *FRs;
};

struct process_info
{
  unsigned int n_regions;
  unsigned int * region_GIDs;
  
  float ** top_bottom_diagonal_2vertices;
  float *** FRs_sg_nt;
  float *** FRs_sg;
  float ** FRs_process;
  char ** region_names;
};

struct system_parameters sp;
struct process_info PI;


void read_global_shared_info(struct system_parameters *sp);
void lntrim(char *str);
void trim(char *buf);
void get_firing_rates_with_time(struct system_parameters *sp, struct process_info * PI);
void resize(int w, int h);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void keyboard(unsigned char key, int x, int y);
void special_input(int key, int x, int y);
void timer(int value);
void display(void);
void read_region_info(struct system_parameters *sp, struct process_info * PI);
void read_process_info(struct system_parameters *sp, struct process_info * PI);
void get_arguments(int argc, char *argv[]);
void SaveImage_PPM(char* fname);
void moving_scene(void);
