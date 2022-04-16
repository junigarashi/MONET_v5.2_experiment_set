#include "./grid_vis.h"

bool switch_mouse_mode = true;
bool switch_save_ppm_image = false;

float light0_position[4] = {0., 0.0, 0.0, 1.0};
float light0_ambient[4] = {0.5, 0.5, 0.5, 0.5};

//mouse button
unsigned int left_button = 0, right_button = 0;
double previous_x_drag_L = 0.0, previous_y_drag_L = 0.0;
double previous_x_drag_R = 0.0, previous_y_drag_R = 0.0;
double drag_L_rotate = 0., drag_L_zoom = 0.;

// camera
const double default_camera_position[3] = {-4., -4., -1.5};

double play_speed = (double)DEFAULT_PLAY_SPEED;
unsigned int t = 0;

unsigned int n_subgrid = 25;
unsigned int n_nt = 18;

char data_directory_name[255];

void get_firing_rates_with_time(struct system_parameters *sp, struct process_info *PI)
{
  FILE *fp;
  char file_path[file_name_character_length];
  unsigned int count_lines = 0;
  char *readline = (char*) malloc(readline_length * sizeof(char));
  char *tok;
  
  printf("# reading firing rate files\n");

  if (sp->switch_FR_sg_nt == true)
    {
      for (unsigned int i_proc = 0; i_proc < sp->n_total_processes; i_proc++)
	{
	  //set file path
	  sprintf(file_path, "../FR_P%d_sg_nt.dat", i_proc);
	  printf("%s\n", file_path);
	  
	  //file open again
	  if ((fp = fopen(file_path, "r")) != NULL) 
	    {
	      count_lines = 0;
	      while ( fgets(readline, readline_length, fp) != NULL ) 
		{
		  tok = strtok(readline, " " );
		  unsigned int process_GID = atoi(tok);
		  tok = strtok(NULL, " " );
		  unsigned int i_nt_sg = atoi(tok);
		  
		  for(unsigned int i_tp = 0; i_tp < sp->n_recording_time_points; i_tp++)
		    {
		      tok = strtok( NULL, " " );
		      tok = strtok( NULL, " " );
		      PI->FRs_sg_nt[i_proc][count_lines][i_tp] = atof(tok);
		      //printf("t:%d p:%d sgnt:%d FR:%f\n", i_tp, i_proc, count_lines, PI->FRs[i_proc][count_lines][i_tp]);
		    }
		  count_lines++;
		}

	      fclose(fp);
	    }
	  else
	    {
	      printf("file open error: spike times: %s\n", file_path);
	      exit(EXIT_FAILURE);
	    }
	}
    }

  printf("# get_spike_times: finish\n");
  printf("total_memory:%f GB\n",
	 2 * 32 * sp->n_total_processes * count_lines * sp->n_recording_time_points / pow(10,9));
  
  if(sp->switch_FR_sg == true)
    {
      for (unsigned int i_proc = 0; i_proc < sp->n_total_processes; i_proc++)
	{
	  //set file path
	  sprintf(file_path, "../FR_P%d_sg.dat", i_proc);
	  printf("%s\n", file_path);
	  
	  //file open again
	  if ((fp = fopen(file_path, "r")) != NULL) 
	    {
	      count_lines = 0;
	      
	      while ( fgets(readline, readline_length, fp) != NULL ) 
		{
		  tok = strtok(readline, " " );
		  unsigned int process_GID = atoi(tok);
		  
		  for(unsigned int i_tp = 0; i_tp < sp->n_recording_time_points; i_tp++)
		    {
		      tok = strtok( NULL, " " );
		      tok = strtok( NULL, " " );
		      PI->FRs_sg[i_proc][count_lines][i_tp] = atof(tok); /////////////////////////
		      //printf("t:%d p:%d FR:%s \n", i_tp, i_proc, tok); //, PI->FRs[i_proc][count_lines][i_tp]);
		    }
		  count_lines++;
		}
	      
	      fclose(fp);
	    }
	  else
	    {
	      printf("file open error: spike times: %s\n", file_path);
	      exit(EXIT_FAILURE);
	    }
	}       
    }
  
  if(sp->switch_FR_process == true)
    {
      for (unsigned int i_proc = 0; i_proc < sp->n_total_processes; i_proc++)
	{
	  //set file path
	  sprintf(file_path, "../FR_P%d_process.dat", i_proc);
	  printf("%s\n", file_path);
	  
	  //file open again
	  if ((fp = fopen(file_path, "r")) != NULL) 
	    {
	      while ( fgets(readline, readline_length, fp) != NULL ) 
		{
		  tok = strtok(readline, " " );
		  unsigned int process_GID = atoi(tok);
		  
		  for(unsigned int i_tp = 0; i_tp < sp->n_recording_time_points; i_tp++)
		    {
		      tok = strtok( NULL, " " );
		      tok = strtok( NULL, " " );
		      PI->FRs_process[i_proc][i_tp] = atof(tok); /////////////////////////
		      //printf("t:%d p:%d FR:%s \n", i_tp, i_proc, tok); //, PI->FRs[i_proc][count_lines][i_tp]);
		    }
		}
	      
	      fclose(fp);
	    }
	  else
	    {
	      printf("file open error: spike times: %s\n", file_path);
	      exit(EXIT_FAILURE);
	    }
	}       
    }
    
}


void read_global_shared_info(struct system_parameters *sp)
{
  FILE *fp;
  char filename[file_name_character_length];
  char *readline = (char*) malloc(readline_length * sizeof(char));
  char *tok;
    
  //read reagion info
  sprintf(filename, "../global_shared_info.dat");
  
  if((fp=fopen(filename, "r"))==NULL)
    {
      fprintf(stderr, "Error: %s: Please prepare global_shared_info.dat by python script.\n", __func__);
      exit(EXIT_FAILURE);
    }
  else{
    while ( fgets(readline, readline_length, fp) != NULL ) 
      {
	
	tok = strtok(readline, "," );
	lntrim(tok);
	trim(tok);

	if (!strcmp(tok, "simulation_time"))
	  {
	    tok = strtok( NULL, "," );
	    sp->end_time = atof(tok);
	  }
	if (!strcmp(tok, "dt"))
	  {
	    tok = strtok( NULL, "," );
	    sp->dt = atof(tok);
	    sp->n_steps_per_dt = (unsigned int) (1.00001/sp->dt);
	    printf("%f %d\n", sp->dt, sp->n_steps_per_dt );
	  }
	
	if (!strcmp(tok, "n_total_process"))
	  {
	    tok = strtok( NULL, "," );
	    sp->n_total_processes = atoi(tok);
	  }
	
 	if (!strcmp(tok, "gsi_n_regions"))
	  {
	    tok = strtok( NULL, "," );
	    sp->n_total_regions = atoi(tok);
	  }
	
 	if (!strcmp(tok, "record_FR_subgrid_neuron_type_per_process"))
	  {
	    tok = strtok( NULL, "," );
	    lntrim(tok);
	    trim(tok);

	    if (!strcmp(tok, "True"))
	      {
		sp->switch_FR_sg_nt = true;
	      }
	    else
	      {
		sp->switch_FR_sg_nt = false;
	      }
	  }

 	if (!strcmp(tok, "record_FR_subgrid_per_process"))
	  {
	    tok = strtok( NULL, "," );
	    lntrim(tok);
	    trim(tok);
	    
	    if (!strcmp(tok, "True"))
	      {
		sp->switch_FR_sg = true;
	      }
	    else
	      {
		sp->switch_FR_sg = false;
	      }
	  }
	
 	if (!strcmp(tok, "record_FR_process"))
	  {
	    tok = strtok( NULL, "," );
	    lntrim(tok);
	    trim(tok);
	    
	    if (!strcmp(tok, "True"))
	      {
		sp->switch_FR_process = true;
	      }
	    else
	      {
		sp->switch_FR_process = false;
	      }
	  }
	
     }

    sp->region_GID_to_region_name = (char **) malloc(sizeof(char*) * sp->n_total_regions);
    for(unsigned int i_region = 0; i_region < sp->n_total_regions; i_region++)
      {
	sp->region_GID_to_region_name[i_region] = (char *) malloc(sizeof(char) * file_name_character_length);
      }
    sp->region_GID_to_n_neuron_types = (unsigned int *) malloc(sizeof(unsigned int) * sp->n_total_regions);
    sp->region_GID_to_n_subregions = (unsigned int *) malloc(sizeof(unsigned int) * sp->n_total_regions);
    for(unsigned int i_region = 0; i_region < sp->n_total_regions; i_region++)
      {
	sp->region_GID_to_n_neuron_types[i_region] = 0;
      }
  
    if((fp=fopen(filename, "r"))==NULL)
      {
	fprintf(stderr, "Error: %s: Please prepare global_shared_info.dat by python script.\n", __func__);
	exit(EXIT_FAILURE);
      }
    else{
      while ( fgets(readline, readline_length, fp) != NULL ) 
	{
	  
	  tok = strtok(readline, "," );
	  lntrim(tok);
	  trim(tok);
	  
	  if (!strcmp(tok, "gsi_region_names"))
	    {
	      for(unsigned int i_region = 0; i_region < sp->n_total_regions; i_region++)
		{
		  tok = strtok( NULL, "," );
		  lntrim(tok);
		  trim(tok);
		  sprintf(sp->region_GID_to_region_name[i_region], "%s", tok);
		  printf("%d %s\n", i_region, sp->region_GID_to_region_name[i_region]);
		}
	    }
	}
    }
      
    sp->start_time = 0;  //ms
    sp->n_recording_time_points = (unsigned int) ((sp->end_time - sp->start_time) / 5 ); //step

    printf("end:%f dt:%f n:%d start:%f steps:%d tp:%d\n",
	   sp->end_time, sp->dt, sp->n_steps_per_dt, sp->start_time, sp->n_calculation_steps,
	   sp->n_recording_time_points);
    
    fclose(fp); 
  }

  // initialization
  sp->camera_position[0] = -4;
  sp->camera_position[1] = -4;
  sp->camera_position[2] = -1.5;
  sp->view_point_position[0] = 1.6/2.;
  sp->view_point_position[1] = 1.6/2.;
  sp->view_point_position[2] = 2.4/2.;
  sp->upper_direction[0] = 0.;
  sp->upper_direction[1] = 0.;
  sp->upper_direction[2] = -1.;
  sp->model_center_position[3] = 1.6/2.;
  sp->model_center_position[3] = 1.6/2.;
  sp->model_center_position[3] = 2.4/2.;
  
  // initialization switch
  sp->i_shown_region = 0;
  sp->i_shown_neuron_type = 0;
  sp->switch_show_degree = 0;
  sp->switch_loop = 1;
  sp->switch_time_stop = 0;
  sp->switch_camera_moving = 0;
}


void read_region_info(struct system_parameters *sp, struct process_info * PI)
{
  FILE *fp;
  char filename[file_name_character_length];
  char *readline = (char*) malloc(readline_length * sizeof(char));
  char *tok;
  
  //read reagion info
  for(unsigned int i_region = 0; i_region < sp->n_total_regions; i_region++)
    {
      sprintf(filename, "../region_%s.dat", sp->region_GID_to_region_name[i_region]);
      printf("%s\n", filename);
    

      if((fp = fopen(filename, "r")) == NULL)
	{
	  fprintf(stderr, "Error: %s: Please prepare region_%s.dat by python script.\n", __func__);
	  exit(EXIT_FAILURE);
	}
      else{
	while ( fgets(readline, readline_length, fp) != NULL ) 
	  {
	    tok = strtok(readline, "," );
	    lntrim(tok);
	    trim(tok);
	    if (!strcmp(tok, "n_subregions"))
	      {
		tok = strtok( NULL, "," );
		sp->region_GID_to_n_subregions[i_region] = atoi(tok);
	      }
	  }
	fclose(fp); 
      }
      printf("sp->region_GID_to_n_subregions[%d]:%d\n", i_region, sp->region_GID_to_n_subregions[i_region]);
      
      if((fp = fopen(filename, "r")) == NULL)
	{
	  fprintf(stderr, "Error: %s: Please prepare region_%s.dat by python script.\n", __func__);
	  exit(EXIT_FAILURE);
	}
      else{
	while ( fgets(readline, readline_length, fp) != NULL ) 
	  {
	    tok = strtok(readline, "," );
	    lntrim(tok);
	    trim(tok);
	    
	    if (!strcmp(tok, "n_neuron_types_per_subregion"))
	      {
		for (unsigned int i_sr = 0; i_sr < sp->region_GID_to_n_subregions[i_region]; i_sr++)
		  {
		    tok = strtok( NULL, "," );
		    sp->region_GID_to_n_neuron_types[i_region] += atoi(tok);
		  }
	      }
	  }
	fclose(fp); 
      }
      printf("sp->region_GID_to_n_neuron_types[%d]:%d\n", i_region, sp->region_GID_to_n_neuron_types[i_region]);
    }

  sp->neuron_type_names = (char ***) malloc(sizeof(char**) * sp->n_total_regions);
  for(unsigned int i_region = 0; i_region < sp->n_total_regions; i_region++)
    {
      sp->neuron_type_names[i_region] = (char **) malloc(sizeof(char*) * sp->region_GID_to_n_neuron_types[i_region]);
      for(unsigned int i_nt = 0; sp->region_GID_to_n_neuron_types[i_region] > i_nt; i_nt++)
	{
	  sp->neuron_type_names[i_region][i_nt] = (char *) malloc(sizeof(char) * file_name_character_length);
	}
    }
  
  for(unsigned int i_region = 0; i_region < sp->n_total_regions; i_region++)
    {
      if((fp=fopen(filename, "r"))==NULL)
	{
	  fprintf(stderr, "Error: %s: Please prepare region_%s.dat by python script.\n", __func__);
	  exit(EXIT_FAILURE);
	}
      else{
	while ( fgets(readline, readline_length, fp) != NULL ) 
	  {
	    
	    tok = strtok(readline, "," );
	    lntrim(tok);
	    trim(tok);
	    
	    if (!strcmp(tok, "neuron_type_names"))
	      {
		for(unsigned int i_nt = 0; i_nt < sp->region_GID_to_n_neuron_types[i_region]; i_nt++)
		  {
		    tok = strtok( NULL, "," );
		    lntrim(tok);
		    trim(tok);
		    sprintf(sp->neuron_type_names[i_region][i_nt], "%s", tok);
		    printf("%d %s\n", i_region, sp->neuron_type_names[i_region][i_nt]);
		  }
	      }
	  }
      }
    }
        
  if (sp->switch_FR_sg_nt == true)
    {
      PI->FRs_sg_nt = (float***) malloc(sizeof(float**) * sp->n_total_processes);
      for (unsigned int i_proc = 0; i_proc < sp->n_total_processes; i_proc++)
	{
	  printf("i_proc:%d region_GID:%d n_neuron_type:%d\n", i_proc, PI->region_GIDs[i_proc], sp->region_GID_to_n_neuron_types[ PI->region_GIDs[i_proc]] );
	  unsigned int n_neuron_types = sp->region_GID_to_n_neuron_types[ PI->region_GIDs[i_proc] ];
	  PI->FRs_sg_nt[i_proc] = (float**) malloc(sizeof(float*) * n_subgrid * n_neuron_types );
	  for (unsigned int i_gs_nt = 0; i_gs_nt < n_subgrid * n_neuron_types; i_gs_nt++)
	    {
	      PI->FRs_sg_nt[i_proc][i_gs_nt] = (float*) malloc(sizeof(float) * sp->n_recording_time_points);
	    }
	}
    }
  
  if (sp->switch_FR_sg == true)
    {
      PI->FRs_sg = (float***) malloc(sizeof(float**) * sp->n_total_processes);
      for (unsigned int i_proc = 0; i_proc < sp->n_total_processes; i_proc++)
	{
	  printf("i_proc:%d region_GID:%d \n", i_proc, PI->region_GIDs[i_proc]);
	  PI->FRs_sg[i_proc] = (float**) malloc(sizeof(float*) * n_subgrid);
	  for (unsigned int i_gs = 0; i_gs < n_subgrid; i_gs++)
	    {
	      PI->FRs_sg[i_proc][i_gs] = (float*) malloc(sizeof(float) * sp->n_recording_time_points);
	    }
	}
    }
  
  if (sp->switch_FR_process == true)
    {
      PI->FRs_process = (float**) malloc(sizeof(float*) * sp->n_total_processes);
      for (unsigned int i_proc = 0; i_proc < sp->n_total_processes; i_proc++)
	{
	  printf("i_proc:%d \n", i_proc);
	  PI->FRs_process[i_proc] = (float*) malloc(sizeof(float) * sp->n_recording_time_points);
	}
    }
}


void read_process_info(struct system_parameters *sp, struct process_info * PI)
{
  FILE *fp;
  char filename[file_name_character_length];
  char *readline = (char*) malloc(readline_length * sizeof(char));
  char *tok;
  
  // allocate memory to process_info
  PI->top_bottom_diagonal_2vertices = (float **) malloc(sizeof(float*) * sp->n_total_processes);
  for (unsigned int i_proc = 0; i_proc < sp->n_total_processes; i_proc++)
    {
      PI->top_bottom_diagonal_2vertices[i_proc] = (float *) malloc(sizeof(float) * 6 ); // x1, y1, z1 x2, y2, z2
    }

  PI->region_GIDs = (unsigned int *) malloc(sizeof(unsigned int) * sp->n_total_processes);

  PI->region_names = (char **) malloc(sizeof(char*) * sp->n_total_processes);
  for (unsigned int i_proc = 0; i_proc < sp->n_total_processes; i_proc++)
    {
      PI->region_names[i_proc] = (char *) malloc(sizeof(char) * file_name_character_length);
    }

  for(unsigned int i_element = 0; i_element < 3; i_element++)
    {
      sp->max_position[i_element] = 0. ;
      sp->min_position[i_element] = 0.;
    }
  
  //read process_info file
  sprintf(filename, "../process_info.dat");

  if((fp = fopen(filename, "r")) == NULL)
    {
      fprintf(stderr, "Error: %s: Please prepare process_info.dat by python script.\n", __func__);
      exit(EXIT_FAILURE);
    }
  else{
    while ( fgets(readline, readline_length, fp) != NULL ) 
      {
	// process GID
	tok = strtok(readline, "," );
	tok = strtok(NULL, "," );
	lntrim(tok);
	unsigned int process_GID = atoi(tok);

	// region GID
	tok = strtok(NULL, "," );
	tok = strtok(NULL, "," );
	lntrim(tok);
	PI->region_GIDs[process_GID] = atoi(tok);

	// region name
	tok = strtok(NULL, "," );
	tok = strtok(NULL, "," );
	lntrim(tok);
	trim(tok);
	
	sprintf(PI->region_names[process_GID], "%s", tok);
	//printf("%s\n", PI->region_names[process_GID]);
	
	// xy plane vertex
	tok = strtok(NULL, "," );
	for(unsigned int i_vertex = 0; i_vertex < 4; i_vertex++)
	  {
	    tok = strtok(NULL, "," );
	    lntrim(tok);
	    //PI->xy_plane_4vertices[i_vertex].x[process_GID] = atof(tok);
	    tok = strtok(NULL, "," );
	    lntrim(tok);
	    //PI->xy_plane_4vertices[i_vertex].y[process_GID] = atof(tok);
	    tok = strtok(NULL, "," );
	    lntrim(tok);
	    //PI->xy_plane_4vertices[i_vertex].z[process_GID] = atof(tok);
	  }
	
	// get spatial extent info of own tile 
	// xyz top bottom edge vertex
	tok = strtok(NULL, "," );
	for(unsigned int i_vertex = 0; i_vertex < 6; i_vertex++)
	  {
	    tok = strtok(NULL, "," );
	    lntrim(tok);
	    PI->top_bottom_diagonal_2vertices[process_GID][i_vertex] = atof(tok);
	    //printf("%d %d %f \n", process_GID, i_vertex, PI->top_bottom_diagonal_2vertices[process_GID][i_vertex]);
	  }

	// maximum position
	if (PI->top_bottom_diagonal_2vertices[process_GID][0] > sp->max_position[0])
	  {
	    sp->max_position[0] = PI->top_bottom_diagonal_2vertices[process_GID][0];
	  }
	if (PI->top_bottom_diagonal_2vertices[process_GID][3] > sp->max_position[0])
	  {
	    sp->max_position[0] = PI->top_bottom_diagonal_2vertices[process_GID][3];
	  }
	if (PI->top_bottom_diagonal_2vertices[process_GID][1] > sp->max_position[1])
	  {
	    sp->max_position[1] = PI->top_bottom_diagonal_2vertices[process_GID][1];
	  }
	if (PI->top_bottom_diagonal_2vertices[process_GID][4] > sp->max_position[1])
	  {
	    sp->max_position[1] = PI->top_bottom_diagonal_2vertices[process_GID][4];
	  }
	if (PI->top_bottom_diagonal_2vertices[process_GID][2] > sp->max_position[2])
	  {
	    sp->max_position[2] = PI->top_bottom_diagonal_2vertices[process_GID][2];
	  }
	if (PI->top_bottom_diagonal_2vertices[process_GID][5] > sp->max_position[2])
	  {
	    sp->max_position[2] = PI->top_bottom_diagonal_2vertices[process_GID][5];
	  }
	    
	// minimum position
	if (PI->top_bottom_diagonal_2vertices[process_GID][0] < sp->min_position[0])
	  {
	    sp->min_position[0] = PI->top_bottom_diagonal_2vertices[process_GID][0];
	  }
	if (PI->top_bottom_diagonal_2vertices[process_GID][3] < sp->min_position[0])
	  {
	    sp->min_position[0] = PI->top_bottom_diagonal_2vertices[process_GID][3];
	  }
	if (PI->top_bottom_diagonal_2vertices[process_GID][1] < sp->min_position[1])
	  {
	    sp->min_position[1] = PI->top_bottom_diagonal_2vertices[process_GID][1];
	  }
	if (PI->top_bottom_diagonal_2vertices[process_GID][4] < sp->min_position[1])
	  {
	    sp->min_position[1] = PI->top_bottom_diagonal_2vertices[process_GID][4];
	  }
	if (PI->top_bottom_diagonal_2vertices[process_GID][2] < sp->min_position[2])
	  {
	    sp->min_position[2] = PI->top_bottom_diagonal_2vertices[process_GID][2];
	  }
	if (PI->top_bottom_diagonal_2vertices[process_GID][5] < sp->min_position[2])
	  {
	    sp->min_position[2] = PI->top_bottom_diagonal_2vertices[process_GID][5];
	  }

	/*
	// subregion cube vertices
	tok = strtok(NULL, "," );
	tok = strtok(NULL, "," );
	//PI->n_my_subregions = atoi(tok);
	//PI->my_subregion_top_bottom_diagonal_2vertices = (numpre * ) malloc(sizeof(numpre) * PI->n_my_subregions * 6);
	
	tok = strtok(NULL, "," );
	for(unsigned int i_sr = 0; i_sr < PI->n_my_subregions; i_sr++)
	  {
	    //printf("%s PID:%d i_sr:%d ", PI->my_region_name, process_GID, i_sr);
	    for(unsigned int i_vertex = 0; i_vertex < 6; i_vertex++)
	      {
		tok = strtok(NULL, "," );
		lntrim(tok);
		//PI->my_subregion_top_bottom_diagonal_2vertices[ i_sr * 6 + i_vertex] = atof(tok);
		//printf("%f, ", PI->my_subregion_top_bottom_diagonal_2vertices[ i_sr * 6 + i_vertex]);
	      }
	    //printf("\n");
	  }
	*/
	
	
      }
    fclose(fp); 
  }
  free(readline);

  // model center position
  for(unsigned int i_element = 0; i_element < 3; i_element++)
    {
      sp->model_center_position[i_element] = (sp->max_position[i_element] + sp->min_position[i_element]) / 2. * 0.001;
      // set initial view pont to model center
      sp->view_point_position[i_element] = sp->model_center_position[i_element];
      printf("%f\n", sp->model_center_position[i_element]);
    }
  
} //read_process_info
 

void lntrim(char *str) {  
  char *p;  
  p = strchr(str, '\n');  
  if(p != NULL) {  
    *p = '\0';  
  }  
}

void trim(char *buf)
{
    char *start, *end;
    char *p = buf;
    
    while(isspace(*p))
    {
        p++;
    }
    start = p;
    
    end = NULL;
    while(*p != '\0')
    {
        if(!isspace(*p))
        {
            end = p;
        }
        p++;
    }
    
    if (end != NULL)
    {
        memmove(buf, start, end - start + 1);
        buf[end - start + 1] = '\0';
    }
    else
        buf[0] = '\0';
}

void resize(int w, int h)
{
  glViewport(0, 0, w, h);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(30.0, (double)w / (double)h, 1.0, 100.0);

  glMatrixMode(GL_MODELVIEW);
}


void mouse(int button, int state, int x, int y)
{
  if (switch_mouse_mode)
    {
      //left button
      if (button == GLUT_LEFT_BUTTON)
	{
	  if (state == GLUT_DOWN) 
	    {
	      left_button = 1;
	      // initialize history of x, y for left drag
	      previous_x_drag_L = x;
	      previous_y_drag_L = y;
	    }
	  else if(state == GLUT_UP) {
	    left_button = 0;
	  }
	}
      //right button
      else if (button == GLUT_RIGHT_BUTTON)
	{
	  if (state == GLUT_DOWN) 
	    {
	      right_button = 1;
	      // initialize history of x, y for right drag
	      previous_x_drag_R = x;
	      previous_y_drag_R = y;
	    }
	  else if(state == GLUT_UP) {
	    right_button = 0;
	  }
	}
    }
}

void motion(int x, int y)
{
  //left drag
  if(left_button == 1)
    {
      sp.camera_position[2] -= (previous_y_drag_L - y) * 0.01;
      previous_y_drag_L = y;
    }

  //no assignment for right drag now
  if(right_button == 1) 
    {

    }
}

void timer(int value) 
{
  GLfloat top = -0.9;
  static GLboolean isUp = GL_TRUE;
  
  if (top > 0.9F) isUp = GL_FALSE;
  else if (top <= -0.9F) isUp = GL_TRUE;
  top += (isUp == GL_TRUE ? 0.01 : -0.01);
  
  glutPostRedisplay();
  glutTimerFunc(1./(double)(play_speed) , timer , 0);
}

void keyboard(unsigned char key, int x, int y)
{

  unsigned int i_rg = 2, i_nt = 0;
  switch (key) {
    /*
  case 'c':
    switch_show_connections  = !switch_show_connections;
    switch(switch_show_connections)
      {
      case 0:
	printf("\rTurn off connection                 ");fflush(stdout);
	break;
      case 1:
	printf("\rTurn on connection                  ");fflush(stdout);
	break;
      }
     break;
    */
  case 'c':
    sp.switch_camera_moving  = ( sp.switch_camera_moving + 1 ) % N_CAMERA_MOVING;
    switch(sp.switch_camera_moving)
      {
      case 0:
	printf("\rTurn off campera moving             ");fflush(stdout);
	break;
      case 1:
	printf("\rCamera moving counterclockwisely    ");fflush(stdout);
	break;
      case 2:
	printf("\rCamera moving clockwisely           ");fflush(stdout);
	break;
      }
    break;

  case 'n':
    if (sp.switch_show_degree == 2) 
      {
	unsigned int n_neuron_types = sp.region_GID_to_n_neuron_types[ PI.region_GIDs[sp.i_shown_region]];
					 
	sp.i_shown_neuron_type = (sp.i_shown_neuron_type + 1) % n_neuron_types; 
	printf("\rCell type:%d                  ", sp.i_shown_neuron_type);
	printf("\rNeuon type: %s               ",  sp.neuron_type_names[sp.i_shown_region][sp.i_shown_neuron_type]);
	fflush(stdout);
      }
    break;
  case 'r':
    if (sp.switch_show_degree == 1 || sp.switch_show_degree == 2) 
      {
	sp.i_shown_region = (sp.i_shown_region + 1) % sp.n_total_regions;
	// reset index of cell type for avoidng different number of cell types between cortex and thalamus
	//printf("\rRegion: %s               ",  regions[i_shown_region]);
	//fflush(stdout);
	sp.i_shown_neuron_type = 0;
      }
    break;
  case 'a':
    //0: all, 1:all cell in one region, 2:one cell in one region
    sp.switch_show_degree = ( sp.switch_show_degree + 1 ) % N_TYPE_SWITCH_SHOW_DEGREE;
    
    switch(sp.switch_show_degree)
      {
      case 0:
	printf("\rShow all cells                   ");
	fflush(stdout);
	break;
      case 1:
	printf("\rShow one region                  ");
	fflush(stdout);
	break;
      case 2:
	printf("\rShow one cell types              ");
	fflush(stdout);
	break;
      }
    break;
    
  case 'g':
    //0:sg_nt , 1:sg , 2:process
    sp.switch_spatial_granularity = ( sp.switch_spatial_granularity + 1 ) % N_TYPE_SWITCH_SPATIAL_GRANULARITY;
    
    switch(sp.switch_spatial_granularity)
      {
      case 0:
	printf("\rShow average FR by neuron type and subgrid ");
	fflush(stdout);
	break;
      case 1:
	printf("\rShow average FR by subgrid                 ");
	fflush(stdout);
	break;
      case 2:
	printf("\rShow average FR by process                 ");
	fflush(stdout);
	break;
      }
    
    break;
    /*
      case 's':
    switch_show_spikes = !switch_show_spikes;
    switch(switch_show_spikes)
      {
      case 0:
	printf("\rTurn off spike                  ");fflush(stdout);
	break;
      case 1:
	printf("\rTurn on spike                  ");fflush(stdout);
	break;
      }

    fflush(stdout);
    break;
    */

    // 0: none, 1: full cell body, 2: gray cell body
  case 'l':
    sp.switch_loop = ( sp.switch_loop + 1 ) % 2;
    switch(sp.switch_loop)
      {
      case 0:
	printf("\rTurn off Loop mode               ");
	fflush(stdout);
	break;
      case 1:
	printf("\rTurn on Loop mode               ");
	fflush(stdout);
	break;
      }    
    break;
    
  case 't':
    sp.switch_time_stop = ( sp.switch_time_stop + 1 ) % 2;
    switch(sp.switch_time_stop)
      {
      case 0:
 	printf("\r Run visualization               ");
	fflush(stdout);
	break;
      case 1:
	printf("\r Stoped time at t = %d               ", t);
	fflush(stdout);
	break;
      }    
    break;

    /*
  // 0: none, 1: full cell body, 2: gray cell body
  case 'l':
    switch_show_cell_location = ( switch_show_cell_location + 1 ) % N_TYPE_SHOW_CELL_LOCATION;
    switch(switch_show_cell_location)
      {
      case 0:
	printf("\rShow no cell bodies               ");
	fflush(stdout);
	break;
      case 1:
	printf("\rShow full colored cell bodies               ");
	fflush(stdout);
	break;
      case 2:
	printf("\rShow small gray cell bodies               ");
	fflush(stdout);
	break;
      }    
    break;

  case 'f':
    switch_show_frame  = !switch_show_frame;
    printf("\rTurn on/off frame               ");
    fflush(stdout);
    break;
  case 'm':
    switch_mouse_mode = !switch_mouse_mode;
    switch(switch_mouse_mode)
      {
      case 0:
	printf("\rTurn off mouse mode                 ");fflush(stdout);
	break;
      case 1:
	printf("\rTurn on mouse mode                  ");fflush(stdout);
	break;
      }
    
    break;
    */
    
  case 'h':
    printf("\n### Help message of 3d_visualization ###\n");
    printf("a: Switch shown extent - 0: all, 1: region mode, 2: neuron type mode \n");
    printf("r: Switch shown region in region mode and neuron type mode \n");
    printf("n: Switch shown neuron types in neuron type mode \n");
    //printf("s: Turn on/off spikes \n");
    //printf("l: Turn on/off cell body positions - 0: no, 1: full color, 2: small gray \n");
    //printf("c: Turn on/off connections if there is preserved connection data \n");
    //printf("f: Turn on/off frame of regions \n");
    printf("b: Set t = 0 \n");
    //printf("Right: Step back by 500 ms \n");
    //printf("Left: Step forward by 500 ms \n");
    //printf("Upper: Increase movie speed x1.25 \n");
    //printf("Lower: Decrease movie speed x1.25 \n");
    printf("1: Default upper diagonal view point \n");
    printf("2: Top view point \n");
    printf("3: Bottom view point \n");
    printf("3: Horizontal view point \n");
    //printf("m: Turn on/off mouse mode\n");
    printf("q or ESC: Quit this application \n");
    printf("h: Help message \n");
    break;

  case 'q':
  case 'Q':
  case '\033': 
    exit(0);
    
  // go back to default view
  case '1':
    drag_L_rotate = 0., drag_L_zoom = 0.;
    sp.camera_position[0] = default_camera_position[0];
    sp.camera_position[1] = default_camera_position[1];
    sp.camera_position[2] = default_camera_position[2];
    sp.view_point_position[0] = sp.model_center_position[0];
    sp.view_point_position[1] = sp.model_center_position[1];
    sp.view_point_position[2] = sp.model_center_position[2];
    sp.upper_direction[0] = 0.;
    sp.upper_direction[1] = 0.;
    sp.upper_direction[2] = -1.;
    printf("\rDefalut view                ");
    fflush(stdout);
    break;

  // change to top view 
  case '2':
    drag_L_rotate = 0., drag_L_zoom = 0.;
    sp.camera_position[0] = sp.model_center_position[0];
    sp.camera_position[1] = sp.model_center_position[1];
    sp.camera_position[2] = -4.5;
    sp.view_point_position[0] = sp.model_center_position[0];
    sp.view_point_position[1] = sp.model_center_position[1];
    sp.view_point_position[2] = sp.model_center_position[2];
    sp.upper_direction[0] = 0.;
    sp.upper_direction[1] = 1.;
    sp.upper_direction[2] = 0.;
    printf("\rTop view                ");
    fflush(stdout);
    break;

  case '3':
    drag_L_rotate = 0., drag_L_zoom = 0.;
    sp.camera_position[0] = sp.model_center_position[0];
    sp.camera_position[1] = sp.model_center_position[1];
    sp.camera_position[2] = 6.;
    sp.view_point_position[0] = sp.model_center_position[0];
    sp.view_point_position[1] = sp.model_center_position[1];
    sp.view_point_position[2] = sp.model_center_position[2];
    sp.upper_direction[0] = 1.;
    sp.upper_direction[1] = 0.;
    sp.upper_direction[2] = 0.;
    printf("\rBottom view                "); 
    fflush(stdout);
    break;

    // lateral view
  case '4': 
    drag_L_rotate = 0., drag_L_zoom = 0.;
    sp.camera_position[0] = sp.model_center_position[0];
    sp.camera_position[1] = 12.;
    sp.camera_position[2] = sp.model_center_position[2] + 1.;
    sp.view_point_position[0] = sp.model_center_position[0];
    sp.view_point_position[1] = sp.model_center_position[1];
    sp.view_point_position[2] = 1.;
    sp.upper_direction[0] = 0.;
    sp.upper_direction[1] = 0.;
    sp.upper_direction[2] = -1.;
    printf("\rLateral view                ");
    fflush(stdout);
    break;

    /*
  case GLUT_KEY_DOWN:
    sp.camera_position[0] += 1;
    sp.camera_position[1] += 1;
    sp.camera_position[2] += 1;
    printf("\rUpper                "); 
    fflush(stdout);
    break;
    */
    
  case 'b':
    t = 0;
    printf("\rGo back to t = 0               "); 
    fflush(stdout);
    break;

  default:
    break;
  }
}

void special_input(int key, int x, int y)
{
  
  switch(key)
    {
      
      // no assignment yet for left right
    case GLUT_KEY_LEFT:
      /*
	if ((int)t - 500 > 0 )
	{
	for ( i_rg = 0; i_rg < N_REGIONS; i_rg++ )
	{
	for ( i_nt = 0; i_nt < n_neuron_types_per_region[i_rg]; i_nt++ )
	{
		sd->count_spikes[i_rg][i_nt] = 0;
	      }
	  }
	t -= 500;
	}
      */
      sp.camera_position[0] += 1;
      printf("\rt = %d               ", t);
      fflush(stdout);
      break;
      
    case GLUT_KEY_RIGHT:
      /*
	(t < sd->simulation_time - 500) ? ( t += 500 ): (t);
      */
      sp.camera_position[0] -= 1;
      printf("\rt = %d               ", t);
      fflush(stdout);
      break;
      
      // change play speed
    case GLUT_KEY_UP:    play_speed *= 1.25; 
      printf("\rPlay speed = %f               ", play_speed);
      fflush(stdout);
      break;
      
    case GLUT_KEY_DOWN:  play_speed /= 1.25;
      printf("\rPlay speed = %f               ", play_speed);
      fflush(stdout);
      break;
    default: return;
    }
  glutPostRedisplay();
}

 
void display(void)
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0); 
  glEnable(GL_COLOR_MATERIAL);
  glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
  glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
  glEnable(GL_DEPTH_TEST);

  // show time
  char str[256];
  sprintf(str, "%d", t);
  glColor4d(1., 1., 1., 0.8); 
  //DrawString(str, GLUT_BITMAP_TIMES_ROMAN_24, 1., 1., 3.);

  glLoadIdentity();
  //gluLookAt: camera position x, y, z, target position x, y, z, camera upper direction x, y, z
  gluLookAt(sp.camera_position[0], sp.camera_position[1], sp.camera_position[2], 
	    sp.view_point_position[0], sp.view_point_position[1], sp.view_point_position[2], 
	    sp.upper_direction[0], sp.upper_direction[1], sp.upper_direction[2]);
  
  // rotation
  //glRotated(drag_L_rotate, 0., 0., 1.);
  
  // draw frame of layers
  //draw_subregion_frame();
  
  // draw axis lines
  //draw_scalebar();

  /*
  unsigned int i_rg = 0, i_sr = 0, i_nt = 0, i_nr = 0, nt_SerLID=0;
  // show cell location
  if (switch_show_cell_location != FALSE)
    { 
      for ( i_rg = 0; i_rg < N_REGIONS; i_rg++ )
	{
	  nt_SerLID = 0;
	  for ( i_sr = 0; i_sr < n_subregions[i_rg]; i_sr++ )
	    {
	      //printf("n_subregions[%d]:%d\n", i_sr, n_subregions[i_rg]);
	      for ( i_nt = 0; i_nt < n_neuron_types[i_rg][i_sr]; i_nt++ )
		{
		  //printf("n_neuron_types[%d][%d]:%d\n", i_nt, n_neuron_types[i_rg][i_sr]);
		  //printf("switch_show_degree:%d, switch_show_degree:%d switch_show_cell_location:%d \n", switch_show_degree, switch_show_degree, switch_show_cell_location);
		  if (switch_show_degree == 0 || (switch_show_degree == 1 && i_rg == i_shown_region) || (i_rg == i_shown_region && nt_SerLID == i_shown_neuron_type ) )
		    {
		      double ps, ss;
		      //printf("switch_show_degree:%d, switch_show_degree:%d switch_show_cell_location:%d \n", switch_show_degree, switch_show_degree, switch_show_cell_location);

		      if (switch_show_cell_location == SHOW_CELL_LOCATION_FULL)
			{
			  ps = POINT_SIZE;
			  ss = SPHERE_SIZE;
			  glColor4d(neuron_type_colors[i_rg][nt_SerLID][0], 
				    neuron_type_colors[i_rg][nt_SerLID][1], 
				    neuron_type_colors[i_rg][nt_SerLID][2], 
				    neuron_type_colors[i_rg][nt_SerLID][3]);
			}
		      else if (switch_show_cell_location == SHOW_CELL_LOCATION_GRAY)
			{
			  ps = POINT_SIZE / 2.;
			  ss = SPHERE_SIZE / 2.;
			  glColor4d(neuron_type_colors[i_rg][nt_SerLID][0] * 0.1, 
				    neuron_type_colors[i_rg][nt_SerLID][1] * 0.1, 
				    neuron_type_colors[i_rg][nt_SerLID][2] * 0.1, 
				    neuron_type_colors[i_rg][nt_SerLID][3] * 0.01);
			}
		      glPointSize(ps);
		      
		      for ( i_nr = 0; i_nr < sd->rg[i_rg].sr[i_sr].nt[i_nt].n_neurons; i_nr++)
			{
			  glPushMatrix();
			  
			  double x = sd->rg[i_rg].sr[i_sr].nt[i_nt].neuron_positions[i_nr][0]*0.001;
			  double y = sd->rg[i_rg].sr[i_sr].nt[i_nt].neuron_positions[i_nr][1]*0.001;
			  double z = sd->rg[i_rg].sr[i_sr].nt[i_nt].neuron_positions[i_nr][2]*0.001;
			  //printf("rg:%d sr:%d nt:%d %f %f %f %d %d %d\n", i_rg, i_sr, i_nt, x, y, z, nt_SerLID, neuron_type_colors[i_rg][nt_SerLID][0], neuron_type_colors[i_rg][nt_SerLID][1], neuron_type_colors[i_rg][nt_SerLID][2]);
			  glTranslated(x, y, z);//平行移動値の設定
			  glutSolidSphere(ss, 5, 5);//引数：(半径, Z軸まわりの分割数, Z軸に沿った分割数)
			  glPopMatrix();
			}
		    }
		  
		  nt_SerLID += 1;
		}
	    }
	}
    }
  */
  /*
  // show spiking 
  if (switch_show_spikes)
    {
      // show spiking 
      for ( i_rg = 0; i_rg < N_REGIONS; i_rg++ )
	{
	  nt_SerLID = 0;
	  for ( i_sr = 0; i_sr < n_subregions[i_rg]; i_sr++ )
	    {
	      for ( i_nt = 0; i_nt < n_neuron_types[i_rg][i_sr]; i_nt++ )
		{
		  if (switch_show_degree == 0 || (switch_show_degree == 1 && i_rg == i_shown_region) || (i_rg == i_shown_region && nt_SerLID == i_shown_neuron_type ) )
		    {
		      double previous_spike_time = 0.;

		      if ( (sd->rg[i_rg].sr[i_sr].nt[i_nt].n_spikes > 0 ) && ( sd->rg[i_rg].sr[i_sr].nt[i_nt].n_spikes > sd->rg[i_rg].sr[i_sr].nt[i_nt].count_spikes ) )
			{
			  //printf("t:%d rg:%d sr:%d nt:%d ns:%d nc:%d \n", t, i_rg, i_sr, i_nt, sd->rg[i_rg].sr[i_sr].nt[i_nt].n_spikes, sd->rg[i_rg].sr[i_sr].nt[i_nt].count_spikes);
			  
			  // skip past spike information
			  while( (double) sd->rg[i_rg].sr[i_sr].nt[i_nt].spike_times[sd->rg[i_rg].sr[i_sr].nt[i_nt].count_spikes] < (double)t*0.025 ){
			    sd->rg[i_rg].sr[i_sr].nt[i_nt].count_spikes++;
			  }
			  //printf("t:%d rg:%d sr:%d nt:%d ns:%d nc:%d \n", t, i_rg, i_sr, i_nt, sd->rg[i_rg].sr[i_sr].nt[i_nt].n_spikes, sd->rg[i_rg].sr[i_sr].nt[i_nt].count_spikes);

			  while( (double)(t * 0.025 ) <= sd->rg[i_rg].sr[i_sr].nt[i_nt].spike_times[ sd->rg[i_rg].sr[i_sr].nt[i_nt].count_spikes] && sd->rg[i_rg].sr[i_sr].nt[i_nt].spike_times[sd->rg[i_rg].sr[i_sr].nt[i_nt].count_spikes] < (double)(t+40)*0.025 )
			    {
			      //printf("t:%d rg:%d sr:%d nt:%d ns:%d nc:%d \n", t, i_rg, i_sr, i_nt, sd->rg[i_rg].sr[i_sr].nt[i_nt].n_spikes, sd->rg[i_rg].sr[i_sr].nt[i_nt].count_spikes);
			      if ( previous_spike_time > sd->rg[i_rg].sr[i_sr].nt[i_nt].spike_times[sd->rg[i_rg].sr[i_sr].nt[i_nt].count_spikes] )
				{
				  printf("Spike times does not seeme to be sorted correctly: %f, -> %f (should became larger) \n", previous_spike_time, sd->rg[i_rg].sr[i_sr].nt[i_nt].spike_times[sd->rg[i_rg].sr[i_sr].nt[i_nt].count_spikes] );
				  exit(EXIT_FAILURE);
				}
			      
			      glPointSize(POINT_SIZE);
			      glColor4d(1., 0., 0., 0.8);
			      glColor4d(neuron_type_colors[i_rg][nt_SerLID][0], 
					neuron_type_colors[i_rg][nt_SerLID][1], 
					neuron_type_colors[i_rg][nt_SerLID][2], 
					neuron_type_colors[i_rg][nt_SerLID][3]);
			      glPushMatrix();
			      
			      unsigned int i_nr = sd->rg[i_rg].sr[i_sr].nt[i_nt].Idx[ sd->rg[i_rg].sr[i_sr].nt[i_nt].count_spikes ];
			      //printf("rg:%d sr:%d nt:%d c:%d n:%d idx:%d %f \n", i_rg, i_sr, i_nt, sd->rg[i_rg].sr[i_sr].nt[i_nt].count_spikes, sd->rg[i_rg].sr[i_sr].nt[i_nt].n_spikes,
			      //sd->rg[i_rg].sr[i_sr].nt[i_nt].Idx[ sd->rg[i_rg].sr[i_sr].nt[i_nt].count_spikes ],
			      //sd->rg[i_rg].sr[i_sr].nt[i_nt].spike_times[ sd->rg[i_rg].sr[i_sr].nt[i_nt].count_spikes]);
			      
			      //printf("t:%d rg:%d sr:%d nt:%d i_nr:%d SerLID:%d c:%d n:%d\n", t, i_rg, i_sr, i_nt, i_nr, nt_SerLID, sd->rg[i_rg].sr[i_sr].nt[i_nt].count_spikes, sd->rg[i_rg].sr[i_sr].nt[i_nt].n_spikes);
			      
			      double x = sd->rg[i_rg].sr[i_sr].nt[i_nt].neuron_positions[i_nr][0]*0.001;
			      double y = sd->rg[i_rg].sr[i_sr].nt[i_nt].neuron_positions[i_nr][1]*0.001;
			      double z = sd->rg[i_rg].sr[i_sr].nt[i_nt].neuron_positions[i_nr][2]*0.001;
			      
			      
			      glTranslated(x, y, z);//平行移動値の設定
			      glutSolidSphere(SPHERE_SIZE, 5, 5);//引数：(半径, Z軸まわりの分割数, Z軸に沿った分割数)
			      glPopMatrix();
			      
			      sd->rg[i_rg].sr[i_sr].nt[i_nt].count_spikes++;
			    }
			}
		      
		    }
		  nt_SerLID += 1;
		}
	    }
	}
    }
  */

  //unsigned int current_step = (unsigned int) (t / sp.n_steps_for_spike_counts_per_time_bin);
  if(sp.switch_spatial_granularity == 0)
    {
      if (sp.switch_show_degree == 0) 
	{
	  for (unsigned int i_proc = 0; i_proc < sp.n_total_processes; i_proc++)
	    {
	      //float proc_x = (i_proc / 2) * 0.6;
	      //float proc_y = (i_proc % 2) * 0.6;
	      // process interval
	      float x_process_origin = PI.top_bottom_diagonal_2vertices[i_proc][0] * 0.001;
	      float y_process_origin = PI.top_bottom_diagonal_2vertices[i_proc][1] * 0.001;
	      float z_process_origin = PI.top_bottom_diagonal_2vertices[i_proc][2] * 0.001 - 1;
	      float x_process_width = (PI.top_bottom_diagonal_2vertices[i_proc][3] - PI.top_bottom_diagonal_2vertices[i_proc][0]) * 0.001;
	      float y_process_width = (PI.top_bottom_diagonal_2vertices[i_proc][4] - PI.top_bottom_diagonal_2vertices[i_proc][1]) * 0.001;
	      float z_process_width = (PI.top_bottom_diagonal_2vertices[i_proc][6] - PI.top_bottom_diagonal_2vertices[i_proc][2]) * 0.001;
	      
	      unsigned int n_neuron_types = sp.region_GID_to_n_neuron_types[ PI.region_GIDs[i_proc] ];
	      
	      for (unsigned int i_sg_nt = 0; i_sg_nt < n_subgrid * n_neuron_types; i_sg_nt++)
		{
		  unsigned int i_nt = (i_sg_nt / 25);
		  unsigned int i_nr = (i_sg_nt % 25);
		  float r, g, b, o;
		  
		  float x_sg_nt = (float) x_process_width * (i_nr / 5) * 0.2;
		  float y_sg_nt = (float) y_process_width * (i_nr % 5) * 0.2;
		  float z_sg_nt = (float) i_nt * 0.2;
		  
		  float x = x_process_origin + x_sg_nt;
		  float y = y_process_origin + y_sg_nt;
		  float z = z_process_origin + z_sg_nt;
		  
		  if (i_nt == 0)
		    {
		      r = 0;
		      g = 0;
		      b = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      o = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      glColor4d(r, g, b, o);
		      glPushMatrix();
		      glTranslated(x, y, z);//平行移動値の設定
		      glutSolidSphere(0.03, 5, 5);//引数：(半径, Z軸まわりの分割数, Z軸に沿った分割数)
		      glPopMatrix();
		      
		    }
		  else if (i_nt == 1)
		    {
		      r = 0;
		      g = 0;
		      b = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      o = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      glColor4d(r, g, b, o);
		      glPushMatrix();
		      glTranslated( x, y, z);//平行移動値の設定
		      glutSolidSphere(0.03, 5, 5);//引数：(半径, Z軸まわりの分割数, Z軸に沿った分割数)
		      glPopMatrix();
		      
		    }
		  else if (i_nt == 2 )
		    {
		      r = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      g = 0;
		      b = 0;
		      o = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      glColor4d(r, g, b, o);
		      glPushMatrix();
		      glTranslated( x, y, z);//平行移動値の設定
		      glutSolidSphere(0.03, 5, 5);//引数：(半径, Z軸まわりの分割数, Z軸に沿った分割数)
		      glPopMatrix();
		      
		    }
		  else if (i_nt == 3 )
		    {
		      r = 0;
		      g = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      b = 0;
		      o = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      glColor4d(r, g, b, o);
		      glPushMatrix();
		      glTranslated( x, y, z);//平行移動値の設定
		      glutSolidSphere(0.03, 5, 5);//引数：(半径, Z軸まわりの分割数, Z軸に沿った分割数)
		      glPopMatrix();
		      
		    }
		  else if (i_nt == 4 || i_nt == 5)
		    {
		      r = 0;
		      g = 0;
		      b = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      o = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      glColor4d(r, g, b, o);
		      glPushMatrix();
		      glTranslated( x, y, z);//平行移動値の設定
		      glutSolidSphere(0.03, 5, 5);//引数：(半径, Z軸まわりの分割数, Z軸に沿った分割数)
		      glPopMatrix();
		      
		    }
		  
		  else if (i_nt == 6)
		    {
		      r = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      g = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      b = 0;
		      o = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      glColor4d(r, g, b, o);
		      glPushMatrix();
		      glTranslated( x, y, z);//平行移動値の設定
		      glutSolidSphere(0.03, 5, 5);//引数：(半径, Z軸まわりの分割数, Z軸に沿った分割数)
		      glPopMatrix();
		      
		    }
		  else if (i_nt == 7)
		    {
		      r = 0;
		      g = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      b = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      o = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      glColor4d(r, g, b, o);
		      glPushMatrix();
		      glTranslated( x, y, z);//平行移動値の設定
		      glutSolidSphere(0.03, 5, 5);//引数：(半径, Z軸まわりの分割数, Z軸に沿った分割数)
		      glPopMatrix();
		    }
		  /*
		    glPushMatrix();
		    glTranslated( x, y, z);//平行移動値の設定
		    glutSolidSphere(0.04, 5, 5);//引数：(半径, Z軸まわりの分割数, Z軸に沿った分割数)
		    glPopMatrix();
		  */
		}
	    }
	}
      else if (sp.switch_show_degree == 1) 
	{
	  
	}
      
      else if (sp.switch_show_degree == 2) 
	{
	  for (unsigned int i_proc = 0; i_proc < sp.n_total_processes; i_proc++)
	    {
	      //float proc_x = (i_proc / 2) * 0.6;
	      //float proc_y = (i_proc % 2) * 0.6;
	      // process interval
	      float x_process_origin = PI.top_bottom_diagonal_2vertices[i_proc][0] * 0.001;
	      float y_process_origin = PI.top_bottom_diagonal_2vertices[i_proc][1] * 0.001;
	      float z_process_origin = PI.top_bottom_diagonal_2vertices[i_proc][2] * 0.001 - 1;
	      float x_process_width = (PI.top_bottom_diagonal_2vertices[i_proc][3] - PI.top_bottom_diagonal_2vertices[i_proc][0]) * 0.001;
	      float y_process_width = (PI.top_bottom_diagonal_2vertices[i_proc][4] - PI.top_bottom_diagonal_2vertices[i_proc][1]) * 0.001;
	      float z_process_width = (PI.top_bottom_diagonal_2vertices[i_proc][6] - PI.top_bottom_diagonal_2vertices[i_proc][2]) * 0.001;
	      
	      //sp.i_shown_neuron_type
	      //unsigned int n_neuron_types = sp.region_GID_to_n_neuron_types[ PI.region_GIDs[i_proc] ];
	      
	      for (unsigned int i_sg = 0; i_sg < n_subgrid; i_sg++)
		{
		  unsigned int i_nt = sp.i_shown_neuron_type;
		  unsigned int i_sg_nt = n_subgrid * sp.i_shown_neuron_type + i_sg;
		  
		  float r, g, b, o;
		  
		  float x_sg = (float) x_process_width * (i_sg / 5) * 0.2;
		  float y_sg = (float) y_process_width * (i_sg % 5) * 0.2;
		  float z_sg = (float) i_nt * 0.2;
		  
		  float x = x_process_origin + x_sg;
		  float y = y_process_origin + y_sg;
		  float z = z_process_origin + z_sg;
		  
		  if (i_nt == 0)
		    {
		      r = 0;
		      g = 0;
		      b = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      o = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      glColor4d(r, g, b, o);
		      glPushMatrix();
		      glTranslated(x, y, z); //平行移動値の設定
		      glutSolidSphere(0.03, 5, 5);//引数：(半径, Z軸まわりの分割数, Z軸に沿った分割数)
		      glPopMatrix();
		    }
		  
		  else if (i_nt == 1)
		    {
		      r = 0;
		      g = 0;
		      b = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      o = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      glColor4d(r, g, b, o);
		      glPushMatrix();
		      glTranslated( x, y, z);//平行移動値の設定
		      glutSolidSphere(0.03, 5, 5);//引数：(半径, Z軸まわりの分割数, Z軸に沿った分割数)
		      glPopMatrix();
		    }
		  else if (i_nt == 2 )
		    {
		      r = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      g = 0;
		      b = 0;
		      o = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      glColor4d(r, g, b, o);
		      glPushMatrix();
		      glTranslated( x, y, z);//平行移動値の設定
		      glutSolidSphere(0.03, 5, 5);//引数：(半径, Z軸まわりの分割数, Z軸に沿った分割数)
		      glPopMatrix();
		    }
		  else if (i_nt == 3 )
		    {
		      r = 0;
		      g = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      b = 0;
		      o = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      glColor4d(r, g, b, o);
		      glPushMatrix();
		      glTranslated( x, y, z);//平行移動値の設定
		      glutSolidSphere(0.03, 5, 5);//引数：(半径, Z軸まわりの分割数, Z軸に沿った分割数)
		      glPopMatrix();
		    }
		  else if (i_nt == 4 || i_nt == 5)
		    {
		      r = 0;
		      g = 0;
		      b = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      o = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      glColor4d(r, g, b, o);
		      glPushMatrix();
		      glTranslated( x, y, z);//平行移動値の設定
		      glutSolidSphere(0.03, 5, 5);//引数：(半径, Z軸まわりの分割数, Z軸に沿った分割数)
		      glPopMatrix();
		    }
		  else if (i_nt == 6)
		    {
		      r = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      g = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      b = 0;
		      o = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      glColor4d(r, g, b, o);
		      glPushMatrix();
		      glTranslated( x, y, z);//平行移動値の設定
		      glutSolidSphere(0.03, 5, 5);//引数：(半径, Z軸まわりの分割数, Z軸に沿った分割数)
		      glPopMatrix();
		      
		    }
		  else if (i_nt == 7)
		    {
		      r = 0;
		      g = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      b = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      o = PI.FRs_sg_nt[i_proc][i_sg_nt][t] / 100.;
		      glColor4d(r, g, b, o);
		      glPushMatrix();
		      glTranslated( x, y, z);//平行移動値の設定
		      glutSolidSphere(0.03, 5, 5);//引数：(半径, Z軸まわりの分割数, Z軸に沿った分割数)
		      glPopMatrix();
		    }
		}
	    }
	}
    }
  else if(sp.switch_spatial_granularity == 1)
    {
      for (unsigned int i_proc = 0; i_proc < sp.n_total_processes; i_proc++)
	{
	  float x_process_origin = PI.top_bottom_diagonal_2vertices[i_proc][0] * 0.001;
	  float y_process_origin = PI.top_bottom_diagonal_2vertices[i_proc][1] * 0.001;
	  float z_process_origin = PI.top_bottom_diagonal_2vertices[i_proc][2] * 0.001 - 1;
	  float x_process_width = (PI.top_bottom_diagonal_2vertices[i_proc][3] - PI.top_bottom_diagonal_2vertices[i_proc][0]) * 0.001;
	  float y_process_width = (PI.top_bottom_diagonal_2vertices[i_proc][4] - PI.top_bottom_diagonal_2vertices[i_proc][1]) * 0.001;
	  float z_process_width = (PI.top_bottom_diagonal_2vertices[i_proc][6] - PI.top_bottom_diagonal_2vertices[i_proc][2]) * 0.001;

	  float r, g, b, o;
	  
	  for (unsigned int i_sg = 0; i_sg < n_subgrid; i_sg++)
	    {
	      r = PI.FRs_sg[i_proc][i_sg][t] / 100.;
	      g = 0;
	      b = 0;
	      o = PI.FRs_sg[i_proc][i_sg][t] / 100.;

	      float x_sg = (float) x_process_width * (i_sg / 5) * 0.2;
	      float y_sg = (float) y_process_width * (i_sg % 5) * 0.2;
	      
	      float x = x_process_origin + x_sg;
	      float y = y_process_origin + y_sg;
	      float z = z_process_origin;
	      
	      glColor4d(r, g, b, o);
	      glPushMatrix();
	      glTranslated(x, y, z); //平行移動値の設定
	      glutSolidSphere(0.03, 5, 5);//引数：(半径, Z軸まわりの分割数, Z軸に沿った分割数)
	      glPopMatrix();
	    }
	}
    }
  else if(sp.switch_spatial_granularity == 2)
    {
      for (unsigned int i_proc = 0; i_proc < sp.n_total_processes; i_proc++)
	{
	  float x = (PI.top_bottom_diagonal_2vertices[i_proc][3] + PI.top_bottom_diagonal_2vertices[i_proc][0]) / 2. * 0.001;
	  float y = (PI.top_bottom_diagonal_2vertices[i_proc][4] + PI.top_bottom_diagonal_2vertices[i_proc][1]) / 2. * 0.001;
	  float z = (PI.top_bottom_diagonal_2vertices[i_proc][6] + PI.top_bottom_diagonal_2vertices[i_proc][2]) / 2. * 0.001;
	  
	  float r, g, b, o;
	  
	  r = PI.FRs_process[i_proc][t] / 100.;
	  g = 0;
	  b = 0;
	  o = PI.FRs_process[i_proc][t] / 100.;

	  glColor4d(r, g, b, o);
	  glPushMatrix();
	  glTranslated(x, y, z); //平行移動値の設定
	  glutSolidSphere(0.1, 5, 5);//引数：(半径, Z軸まわりの分割数, Z軸に沿った分割数)
	  glPopMatrix(); 
	}
    }
  
  // moving scene
  moving_scene();
 
  // save PPM image
  if (switch_save_ppm_image)
    {
      char filename[256];
      sprintf(filename, "./ppm_data/%d.ppm", t);
      SaveImage_PPM(&filename);
    }
  
  glutSwapBuffers();

  if (sp.switch_time_stop)
    {
    }
  else
    {
      if (t  < sp.end_time / 5 )
	{
	  //printf("\r %d               ", t);
	  t += 1;
	}
      else if (sp.switch_loop)
	{
	  t = 0;
	}
      else
	{
	  glutPostRedisplay();
	}
    }
}

void get_arguments(int argc, char *argv[])
{
  int c;
  char filename[256];

  if( argc < 2 )
    {
      printf("Usage: %s data_directory_name\n", argv[0]);
      exit(EXIT_FAILURE);
    }
  else if (argc == 2)
    {
      sprintf(data_directory_name, "%s", argv[1]);
    }
  
  else if (argc > 2 )
    {
      sprintf(data_directory_name, "%s", argv[1]);
      
      while( ( c = getopt(argc, argv, "i")) != -1)
	{
	  switch (c) 
	    {
	    case 'i':
	      fprintf(stdout, "###Save image mode (-i option) ###\n");
	      sprintf(filename, "./ppm_data/", data_directory_name);

	      if (mkdir(filename, 0777) )
		{
		  fprintf(stdout, "?");
		}
	      else{
		fprintf(stdout, "?");
	      }
	      switch_save_ppm_image = true;

	      break;
	    }
	}
    }
  
}

void moving_scene(void)
{
  if(sp.switch_camera_moving == 1)
    {
      // camera
      sp.camera_position[0] = 5 * sinf( 2.0 * 3.14 * (t % 1000) / 1000.) ;
      sp.camera_position[1] = 5 * cosf( 2.0 * 3.14 * (t % 1000) / 1000.) ;
    }
  else if (sp.switch_camera_moving == 2)
    {
      sp.camera_position[0] = 5 * sinf( -2.0 * 3.14 * (t % 1000) / 1000.) ;
      sp.camera_position[1] = 5 * cosf( -2.0 * 3.14 * (t % 1000) / 1000.) ;
     }
}


void SaveImage_PPM(char* fname)
{ // save current screen buffe image to the ppm file
  FILE* fp_output;
  int viewport[4];
  glGetIntegerv( GL_VIEWPORT, viewport);   // get location and size of viewport
  const unsigned int width = viewport[2];   // viewport width
  const unsigned int height = viewport[3];  //  viewport height
  void* image = malloc(3 * width * height);  // allocate buffer (rgb)*width*height
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadPixels(0, 0, width,height, GL_RGB, GL_UNSIGNED_BYTE, image);   // get image
  
  // write to file
  if ((fp_output = fopen(fname, "w")) == NULL)
    {
      printf("File open error: SaveImage_PPM: %s\n", fname);
      exit(EXIT_FAILURE);
    }
  fprintf(fp_output, "P3\n");
  fprintf(fp_output, "%d %d\n", width, height);
  fprintf(fp_output, "255\n");
  char* img = (char*)image;
  for(unsigned int ih = 0; ih < height; ih++){
    for(unsigned int iw = 0; iw < width; iw++){
      unsigned int i = (height-1-ih)*width+iw;
      int r = (unsigned char)img[i * 3 + 0];
      int g = (unsigned char)img[i * 3 + 1];
      int b = (unsigned char)img[i * 3 + 2];
      fprintf(fp_output, "%d %d %d\n", r, g, b);
    }
  }
  fclose(fp_output);
}
