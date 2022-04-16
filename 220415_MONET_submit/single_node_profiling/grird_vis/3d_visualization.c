// 3d_visualization.c

#include "3d_visualization.h"

void idle(void)
{
  glutPostRedisplay();
}

void display(void)
{
  //printf("%d\n", t);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);//光源0を利用
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
  gluLookAt(camera_position[0], camera_position[1], camera_position[2], 
	    center_position[0], center_position[1], center_position[2], 
	    upper_direction[0], upper_direction[1], upper_direction[2]);
  
  // rotation
  //glRotated(drag_L_rotate, 0., 0., 1.);
  
  // draw frame of layers
  //draw_subregion_frame();
  
  // draw axis lines
  //draw_scalebar();

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
  
  // show connections
  unsigned int i_pre_rg = 0, i_post_rg = 0, i_pre_sr = 0, i_post_sr = 0, i_pre_nt = 0, i_post_nt = 0, i_pre_nr = 0, i_post_nr = 0, i_connection = 0;
  if (switch_show_connections)
    {
      if (switch_show_degree == 0)
	{
	  // loop presynaptic region
	  for ( i_pre_rg = 0; i_pre_rg < N_REGIONS; i_pre_rg++)
	    {
	      for ( i_pre_sr = 0; i_pre_sr < n_subregions[i_pre_rg]; i_pre_sr++)
		{
		  for ( i_pre_nt = 0; i_pre_nt < n_neuron_types[i_pre_rg][i_pre_sr]; i_pre_nt++)
		    {
		      // loop postsynaptic region
		      for ( i_post_rg = 0; i_post_rg < N_REGIONS; i_post_rg++)
			{ 
			  for ( i_post_sr = 0; i_post_sr < n_subregions[i_post_rg]; i_post_sr++)
			    {
			      for ( i_post_nt = 0; i_post_nt < n_neuron_types[i_post_rg][i_post_sr]; i_post_nt++)
				{
				  /* need change n_connections, connections
				     for (i_connection = 0; i_connection < sd->n_connections[i_pre_rg][i_pre_nt][i_post_rg][i_post_nt]; i_connection++)
				     {
				     int i_pre_runnning_number = sd->connections[i_pre_rg][i_pre_nt][i_post_rg][i_post_nt][i_connection][0];
				     double pre_x = sd->rg[i_pre_rg].sr[i_pre_sr].nt[i_pre_nt].neuron_positions[i_pre_nr][1]*0.001;
				      double pre_y = sd->rg[i_pre_rg].sr[i_pre_sr].nt[i_pre_nt].neuron_positions[i_pre_nr][2]*0.001;
				      double pre_z = sd->rg[i_pre_rg].sr[i_pre_sr].nt[i_pre_nt].neuron_positions[i_pre_nr][3]*0.001;

				      
				      int i_post_running_number = sd->connections[i_pre_rg][i_pre_nt][i_post_rg][i_post_nt][i_connection][1];
				      double post_x = sd->rg[i_post_rg].sr[i_post_sr].nt[i_post_nt].neuron_positions[i_post_nr][1]*0.001;
				      double post_y = sd->rg[i_post_rg].sr[i_post_sr].nt[i_post_nt].neuron_positions[i_post_nr][2]*0.001;
				      double post_z = sd->rg[i_post_rg].sr[i_post_sr].nt[i_post_nt].neuron_positions[i_post_nr][3]*0.001;

				      glLineWidth(1.);
				      glBegin( GL_LINES );
				      glColor4d(neuron_type_colors[i_pre_rg][i_pre_nt][0], 
						neuron_type_colors[i_pre_rg][i_pre_nt][1], 
						neuron_type_colors[i_pre_rg][i_pre_nt][2], 
						neuron_type_colors[i_pre_rg][i_pre_nt][3]);
				      glVertex3f( pre_x, pre_y, pre_z );
				      glVertex3f( post_x, post_y, post_z );
				      glEnd();
				    }
				  */
				}
			    }
			}
		    }
		}
	    }
	}
      else 
	{
	  for ( i_post_rg = 0; i_post_rg < N_REGIONS; i_post_rg++)
	    {
	      for ( i_post_sr = 0; i_post_sr < n_subregions[i_post_rg]; i_post_sr++)
		{
		  for ( i_post_nt = 0; i_post_nt < n_neuron_types[i_pre_rg][i_pre_sr]; i_post_nt++)
		    {
		      // change sd reladed data
		      /*
		      for (i_connection = 0; i_connection < sd->n_connections[i_pre_rg][i_pre_nt][i_post_rg][i_post_nt]; i_connection++)
			{
			  int i_pre_runnning_number = sd->connections[i_pre_rg][i_pre_nt][i_post_rg][i_post_nt][i_connection][0];
			  double pre_x = sd->rg[i_pre_rg].sr[i_pre_sr].nt[i_pre_nt].neuron_positions[i_pre_nr][1]*0.001;
			  double pre_y = sd->rg[i_pre_rg].sr[i_pre_sr].nt[i_pre_nt].neuron_positions[i_pre_nr][2]*0.001;
			  double pre_z = sd->rg[i_pre_rg].sr[i_pre_sr].nt[i_pre_nt].neuron_positions[i_pre_nr][3]*0.001;
			  
			  int i_post_running_number = sd->connections[i_pre_rg][i_pre_nt][i_post_rg][i_post_nt][i_connection][1];
			  double post_x = sd->rg[i_post_rg].sr[i_post_sr].nt[i_post_nt].neuron_positions[i_post_nr][1]*0.001;
			  double post_y = sd->rg[i_post_rg].sr[i_post_sr].nt[i_post_nt].neuron_positions[i_post_nr][2]*0.001;
			  double post_z = sd->rg[i_post_rg].sr[i_post_sr].nt[i_post_nt].neuron_positions[i_post_nr][3]*0.001;

			  glLineWidth(1.);
			  glBegin( GL_LINES );
			  glColor4d(neuron_type_colors[i_pre_rg][i_pre_nt][0],
				    neuron_type_colors[i_pre_rg][i_pre_nt][1],
				    neuron_type_colors[i_pre_rg][i_pre_nt][2],
				    neuron_type_colors[i_pre_rg][i_pre_nt][3]);
			  glVertex3f( pre_x, pre_y, pre_z );
			  glVertex3f( post_x, post_y, post_z );
			  glEnd();
			  }
		      */
		    }
		}
	    }
	}
    }

  // save PPM image 
  if (switch_save_ppm_image)
    {
      char filename[256];
      sprintf(filename, "../../ppm_data/%d.ppm", t/40);
      SaveImage_PPM(&filename);
    }
  
  glutSwapBuffers();
  
  if (t*0.025 < sd->simulation_time){
    t += 40;
  }
  else{
    glutPostRedisplay();
  }
}


void draw_subregion_frame()
{
  if (switch_show_frame)
    {
      unsigned int i_subregion, i;
      
      for (i_subregion = 0; i_subregion < 5; i_subregion++)
	{
	  double O_X = 0., O_Y = 0., upper_Z = sd->subregion_z_upper_limits[i_subregion];
	  double layer_depth = sd->subregion_thicknesses[i_subregion];
	  
	  vertex[0][0] = O_X;
	  vertex[0][1] = O_Y;
	  vertex[0][2] = upper_Z;
	  vertex[1][0] = O_X+sd->length_on_a_side;
	  vertex[1][1] = O_Y;
	  vertex[1][2] = upper_Z;
	  vertex[2][0] = O_X+sd->length_on_a_side;
	  vertex[2][1] = O_Y+sd->length_on_a_side;
	  vertex[2][2] = upper_Z;
	  vertex[3][0] = O_X;
	  vertex[3][1] = O_Y+sd->length_on_a_side;
	  vertex[3][2] = upper_Z;
	  vertex[4][0] = O_X;
	  vertex[4][1] = O_Y;
	  vertex[4][2] = upper_Z + layer_depth;
	  vertex[5][0] = O_X+sd->length_on_a_side;
	  vertex[5][1] = O_Y;
	  vertex[5][2] = upper_Z + layer_depth;
	  vertex[6][0] = O_X+sd->length_on_a_side;
	  vertex[6][1] = O_Y+sd->length_on_a_side;
	  vertex[6][2] = upper_Z + layer_depth;
	  vertex[7][0] = O_X;
	  vertex[7][1] = O_Y+sd->length_on_a_side;
	  vertex[7][2] = upper_Z + layer_depth;
	  
	  glLineWidth(1.);
	  glBegin(GL_LINES);
	  glColor4d(0., 0.0, 0., 0.5);
	  for (i = 0; i < 12; ++i) {
	    glVertex3dv(vertex[edge[i][0]]);
	    glVertex3dv(vertex[edge[i][1]]);
	  }
	  glEnd();
	  
	}      

      for (i_subregion = 0; i_subregion < 2; i_subregion++)
	{
	  double O_X = 0., O_Y = 0., upper_Z = sd->M1_thalamus_subregion_z_upper_limits[i_subregion];
	  double layer_depth = sd->M1_thalamus_subregion_thicknesses[i_subregion];

	  vertex[0][0] = O_X;
	  vertex[0][1] = O_Y;
	  vertex[0][2] = upper_Z;
	  vertex[1][0] = O_X+sd->length_on_a_side;
	  vertex[1][1] = O_Y;
	  vertex[1][2] = upper_Z;
	  vertex[2][0] = O_X+sd->length_on_a_side;
	  vertex[2][1] = O_Y+sd->length_on_a_side;
	  vertex[2][2] = upper_Z;
	  vertex[3][0] = O_X;
	  vertex[3][1] = O_Y+sd->length_on_a_side;
	  vertex[3][2] = upper_Z;
	  vertex[4][0] = O_X;
	  vertex[4][1] = O_Y;
	  vertex[4][2] = upper_Z + layer_depth;
	  vertex[5][0] = O_X+sd->length_on_a_side;
	  vertex[5][1] = O_Y;
	  vertex[5][2] = upper_Z + layer_depth;
	  vertex[6][0] = O_X+sd->length_on_a_side;
	  vertex[6][1] = O_Y+sd->length_on_a_side;
	  vertex[6][2] = upper_Z + layer_depth;
	  vertex[7][0] = O_X;
	  vertex[7][1] = O_Y+sd->length_on_a_side;
	  vertex[7][2] = upper_Z + layer_depth;
	  
	  glLineWidth(1.);
	  glBegin(GL_LINES);
	  glColor4d(0., 0.0, 0., 0.5);
	  for (i = 0; i < 12; ++i) {
	    glVertex3dv(vertex[edge[i][0]]);
	    glVertex3dv(vertex[edge[i][1]]);
	  }
	  glEnd();
	}
    }
}

void draw_scalebar()
{
  glLineWidth(5.);
  glBegin(GL_LINES);
  glColor3d(1., 0., 0.);
  glVertex3d(-scalebar_length - 0.1,  sd->model_opposing_corner_position[1] + scalebar_length, sd->model_opposing_corner_position[2] + scalebar_length);
  glVertex3d(-0.1, sd->model_opposing_corner_position[1] + scalebar_length, sd->model_opposing_corner_position[2] + scalebar_length);
  glEnd();
  
  glLineWidth(5.);
  glBegin(GL_LINES);
  glColor3d(0., 1., 0.);
  glVertex3d(-scalebar_length - 0.1, sd->model_opposing_corner_position[1] + scalebar_length, sd->model_opposing_corner_position[2] + scalebar_length);
  glVertex3d(-scalebar_length - 0.1, sd->model_opposing_corner_position[1], sd->model_opposing_corner_position[2] + scalebar_length);
  glEnd();
  
  glLineWidth(5.);
  glBegin(GL_LINES);
  glColor3d(0., 0., 1.);
  glVertex3d(-scalebar_length - 0.1, sd->model_opposing_corner_position[1] + scalebar_length, sd->model_opposing_corner_position[2] + scalebar_length);
  glVertex3d(-scalebar_length - 0.1, sd->model_opposing_corner_position[1] + scalebar_length, sd->model_opposing_corner_position[2]);
  glEnd();
}

void resize(int w, int h)
{
  glViewport(0, 0, w, h);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(30.0, (double)w / (double)h, 1.0, 100.0);

  glMatrixMode(GL_MODELVIEW);
}

void keyboard(unsigned char key, int x, int y)
{
  unsigned int i_rg = 2, i_nt = 0;

  switch (key) {
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
  case 'n':
    if (switch_show_degree == 2) 
      {
	i_shown_neuron_type = (i_shown_neuron_type + 1) % (n_neuron_types_per_region[i_shown_region]);
	printf("\rCell type: %s               ",  neuron_types[i_shown_region][i_shown_neuron_type]);
	fflush(stdout);
      }
    break;
  case 'r':
    if (switch_show_degree == 1 || switch_show_degree == 2) 
      {
	i_shown_region = (i_shown_region + 1) % N_REGIONS;
	// reset index of cell type for avoidng different number of cell types between cortex and thalamus
	printf("\rRegion: %s               ",  regions[i_shown_region]);
	fflush(stdout);
	i_shown_neuron_type = 0;
      }
    break;
  //0: all, 1:all cell in one region, 2:one cell in one region
  case 'a':
    switch_show_degree = ( switch_show_degree + 1 ) % N_TYPE_SWITCH_SHOW_DEGREE;
    switch(switch_show_degree)
      {
      case 0:
	printf("\rShow all cells                   ");fflush(stdout);
	break;
      case 1:
	printf("\rShow one region                  ");fflush(stdout);
	break;
      case 2:
	printf("\rShow one cell types              ");fflush(stdout);
	break;
      }
    break;
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
  case 'h':
    printf("\n### Help message of 3d_visualization ###\n");
    printf("a: Switch shown extent - 0: all, 1: region mode, 2: cell type mode \n");
    printf("r: Switch shown region in region mode and cell type modefor key a - motor cortex or thalamus \n");
    printf("n: Switch shown cell types in cell type mode \n");
    printf("s: Turn on/off spikes \n");
    printf("l: Turn on/off cell body positions - 0: no, 1: full color, 2: small gray \n");
    printf("c: Turn on/off connections if there is preserved connection data \n");
    printf("f: Turn on/off frame of regions \n");
    printf("b: Set t = 0 \n");
    printf("Right: Step back by 500 ms \n");
    printf("Left: Step forward by 500 ms \n");
    printf("Upper: Increase movie speed x1.25 \n");
    printf("Lower: Decrease movie speed x1.25 \n");
    printf("1: Default upper diagonal view point \n");
    printf("2: Top view point \n");
    printf("3: Bottom view point \n");
    printf("m: Turn on/off mouse mode\n");
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
    camera_position[0] = default_camera_position[0];
    camera_position[1] = default_camera_position[1];
    camera_position[2] = default_camera_position[2];
    center_position[0] = sd->model_center_position[0];
    center_position[1] = sd->model_center_position[1];
    center_position[2] = sd->model_center_position[2];
    upper_direction[0] = 0.;
    upper_direction[1] = 0.;
    upper_direction[2] = -1.;
    printf("\rDefalut view                ");
    fflush(stdout);
    break;

  // change to top view 
  case '2':
    drag_L_rotate = 0., drag_L_zoom = 0.;
    camera_position[0] = sd->model_center_position[0];
    camera_position[1] = sd->model_center_position[1];
    camera_position[2] = -4.5;
    center_position[0] = sd->model_center_position[0];
    center_position[1] = sd->model_center_position[1];
    center_position[2] = sd->model_center_position[2];
    upper_direction[0] = 1.;
    upper_direction[1] = 0.;
    upper_direction[2] = 0.;
    printf("\rTop view                ");
    fflush(stdout);
    break;

  case '3':
    drag_L_rotate = 0., drag_L_zoom = 0.;
    camera_position[0] = sd->model_center_position[0];
    camera_position[1] = sd->model_center_position[1];
    camera_position[2] = 6.;
    center_position[0] = sd->model_center_position[0];
    center_position[1] = sd->model_center_position[1];
    center_position[2] = sd->model_center_position[2];
    upper_direction[0] = 1.;
    upper_direction[1] = 0.;
    upper_direction[2] = 0.;
    printf("\rBottom view                "); 
    fflush(stdout);
    break;

    // lateral view
  case '4': 
    drag_L_rotate = 0., drag_L_zoom = 0.;
    camera_position[0] = sd->model_center_position[0];
    camera_position[1] = 6.;
    camera_position[2] = sd->model_center_position[2] +1.;
    center_position[0] = sd->model_center_position[0];
    center_position[1] = sd->model_center_position[1];
    center_position[2] = 1.;
    upper_direction[0] = 0.;
    upper_direction[1] = 0.;
    upper_direction[2] = -1.;
    printf("\rLateral view                ");
    fflush(stdout);
    break;

  case 'b':
    for ( i_rg = 0; i_rg < N_REGIONS; i_rg++ )
      {
	for ( i_nt = 0; i_nt < n_neuron_types_per_region[i_rg]; i_nt++ )
	  {
	    sd->count_spikes[i_rg][i_nt] = 0;
	  }
      }
    t = 0;
    printf("\rGo back to t = 0               "); 
    fflush(stdout);
    break;
    
  default:
    break;
  }
}

// prepare data structure
void prepare_data_structure(struct data **sd)
{
  int i, j, k;

  (*sd) = (struct data * ) malloc(sizeof(struct data) * N_REGIONS);
  (*sd)->rg = (struct region * ) malloc(sizeof(struct region) * N_REGIONS);
  
  for (i = 0; i < N_REGIONS; i++)
    {
      (*sd)->rg[i].sr = (struct subregion * ) malloc(sizeof(struct subregion) * n_subregions[i]);
      
      for (j = 0; j < n_subregions[i]; j++)
	{
	  (*sd)->rg[i].sr[j].nt = (struct neuron_type * ) malloc(sizeof(struct neuron_type ) * n_neuron_types[i][j]);

	}
    }
}


// get information of simulation time and positions of frame
void get_model_position(struct data *sd)
{
  char file_path[256], s[256], s1[256], s2[256], s3[256], s4[256], s5[256], s6[256], s7[256], s8[256], *tok; 
  char input_words[10][256];
  double d1, d2, d3, d4, d5;
  int i, ret, count_words=0;
  char readline[256] = {'\0'};

  //sd = (struct data *) malloc(sizeof(struct data));
  
  //display message
  printf("# reading model position\n");

  // path to setting file
  sprintf(file_path, "../inputs/%s/preprocessed_data/system_info.dat", data_directory_name);

  //file open
  /*
  if ((fp = fopen(file_path, "r")) == NULL) 
    {
      printf("file open error: get_GID_and_position: %s\n", file_path);
      exit(EXIT_FAILURE);
    }
  double temp_length_on_a_side;

  
  if((fp=fopen(file_path, "r"))!=NULL)
    {
      //just count 
      while ( fgets(readline, 256, fp) != NULL ) 
	{
	  tok=strtok(readline, "," );
	  lntrim(tok);
	  trim(tok);
	  
	  if (!strcmp(tok, "length_on_a_side"))
	    {
	      tok = strtok( NULL, "," );
	      printf("length_on_a_side:%f\n", atof(tok));	  
	      sd->length_on_a_side = atof(tok)/1000.;
	      sd->model_opposing_corner_position[0] = sd->length_on_a_side;
	      sd->model_opposing_corner_position[1] = sd->length_on_a_side;
	      sd->model_center_position[0] = sd->length_on_a_side/2.;
	      sd->model_center_position[1] = sd->length_on_a_side/2.;
	    }
	  else if (!strcmp(tok, "simulation_time"))
	    {
	      tok = strtok( NULL, "," );
	      sd->simulation_time = (int)(atof(tok) * 40);
	      printf("simulation_time:%d\n", sd->simulation_time);
	    }
	}
    }
  fclose(fp);
  */
  sd->length_on_a_side = 1300./1000.;
  sd->model_opposing_corner_position[0] = sd->length_on_a_side;
  sd->model_opposing_corner_position[1] = sd->length_on_a_side;
  sd->model_center_position[0] = sd->length_on_a_side/2.;
  sd->model_center_position[1] = sd->length_on_a_side/2.;
  sd->simulation_time = 1000;
  
  /*
  // path to M1 parameter sli file
  sprintf(file_path, "../../../outputs/simulations/%s/M1_parameters.sli", data_directory_name);

  //file open 
  if ((fp = fopen(file_path, "r")) == NULL) 
    {
      printf("file open error: get_GID_and_position: %s\n", file_path);
      exit(EXIT_FAILURE);
    }

  while( ( ret = fscanf( fp, "%s %s %lf %lf %lf %lf %lf %s", s1, s2, &d1, &d2, &d3, &d4, &d5, s3 ) ) != EOF )
    {
      if (! strcmp(s1, "/subregion_thicknesses"))
	{
	  sd->subregion_thicknesses[0] = d1/1000.;
	  sd->subregion_thicknesses[1] = d2/1000.;
	  sd->subregion_thicknesses[2] = d3/1000.;
	  sd->subregion_thicknesses[3] = d4/1000.;
	  sd->subregion_thicknesses[4] = d5/1000.;
	}
      if (! strcmp(s1, "/subregion_z_upper_limits"))
	{
	  sd->subregion_z_upper_limits[0] = d1/1000.;
	  sd->subregion_z_upper_limits[1] = d2/1000.;
	  sd->subregion_z_upper_limits[2] = d3/1000.;
	  sd->subregion_z_upper_limits[3] = d4/1000.;
	  sd->subregion_z_upper_limits[4] = d5/1000.;
	}
    }
  fclose(fp);

  // path to S1 parameter sli file
  sprintf(file_path, "../../../outputs/simulations/%s/S1_parameters.sli", data_directory_name);

  //file open 
  if ((fp = fopen(file_path, "r")) == NULL) 
    {
      printf("file open error: get_GID_and_position: %s\n", file_path);
      exit(EXIT_FAILURE);
    }

  while( ( ret = fscanf( fp, "%s %s %lf %lf %lf %lf %lf %s", s1, s2, &d1, &d2, &d3, &d4, &d5, s3 ) ) != EOF )
    {
      if (! strcmp(s1, "/subregion_thicknesses"))
	{
	  sd->S1_subregion_thicknesses[0] = d1/1000.;
	  sd->S1_subregion_thicknesses[1] = d2/1000.;
	  sd->S1_subregion_thicknesses[2] = d3/1000.;
	  sd->S1_subregion_thicknesses[3] = d4/1000.;
	  sd->S1_subregion_thicknesses[4] = d5/1000.;
	  sd->S1_subregion_thicknesses[5] = d5/1000.;
	  sd->S1_subregion_thicknesses[6] = d5/1000.;
	}
      if (! strcmp(s1, "/subregion_z_upper_limits"))
	{
	  sd->S1_subregion_z_upper_limits[0] = d1/1000.;
	  sd->S1_subregion_z_upper_limits[1] = d2/1000.;
	  sd->S1_subregion_z_upper_limits[2] = d3/1000.;
	  sd->S1_subregion_z_upper_limits[3] = d4/1000.;
	  sd->S1_subregion_z_upper_limits[4] = d5/1000.;
	  sd->S1_subregion_z_upper_limits[5] = d5/1000.;
	  sd->S1_subregion_z_upper_limits[6] = d5/1000.;
	}
    }
  fclose(fp);

  // path to M1 thalamus parameter sli file
  sprintf(file_path, "../../../outputs/simulations/%s/M1_thalamus_parameters.sli", data_directory_name);

  //file open 
  if ((fp = fopen(file_path, "r")) == NULL) 
    {
      printf("file open error: get_GID_and_position: %s\n", file_path);
      exit(EXIT_FAILURE);
    }

  while( fgets( s, 256, fp ) != NULL )
    {
      if(strstr(s,"/subregion_thicknesses")!=NULL)
	{
	  tok = strtok(s, " ");
	  for(i=0; tok[i] != '\0';i++)
	    {
	      input_words[0][i] = tok[i];
	    }

	  count_words = 1;
	  while( (tok = strtok( NULL, " " ) ) != NULL && count_words < 10)
	    {
	      //tok = strtok( NULL, " " );  
	      for(i=0; tok[i] != '\0';i++)
		{
		  input_words[count_words][i] = tok[i];
		}
	      count_words++;
	    }
	  sd->M1_thalamus_subregion_thicknesses[0] = atof(input_words[2])/1000.;
	  sd->M1_thalamus_subregion_thicknesses[1] = atof(input_words[3])/1000.;
	}
      
      if(strstr(s,"/subregion_z_upper_limits")!=NULL)
	{
	  tok = strtok(s, " ");
	  for(i=0; tok[i] != '\0';i++)
	    {
	      input_words[0][i] = tok[i];
	    }

	  count_words = 1;
	  while( (tok = strtok( NULL, " " ) ) != NULL && count_words < 10)
	    {
	      //tok = strtok( NULL, " " );  
	      for(i=0; tok[i] != '\0';i++)
		{
		  input_words[count_words][i] = tok[i];
		}
	      count_words++;
	    }
	  sd->M1_thalamus_subregion_z_upper_limits[0] = atof(input_words[2])/1000.;
	  sd->M1_thalamus_subregion_z_upper_limits[1] = atof(input_words[3])/1000.;
	}
    }
  fclose(fp);
  
  // set model position for z direction
  sd->model_opposing_corner_position[2] =  sd->M1_thalamus_subregion_z_upper_limits[1] + sd->M1_thalamus_subregion_thicknesses[1];
  sd->model_center_position[2] = sd->model_opposing_corner_position[2]/2.;
  */
  
  // set camera position for first default view
  center_position[0] = sd->model_center_position[0];
  center_position[1] = sd->model_center_position[1];
  center_position[2] = 2500/1000.;//sd->model_center_position[2];
}


void get_neuron_positions(struct data *sd)
{

  char file_path[256], s[256], buf[256];
  int i_nt = 0, i_rg = 0, i_sr, i_nr = 0, n = 0, n_lines = 0, count_lines = 0, i_proc = 0;

  //display message
  printf("# reading neuron positions\n");

  // only read for count lines to get numbers of neurons
  for ( i_rg = 0; i_rg < N_REGIONS; i_rg++) 
    {
      unsigned int count_nt = 0;
      for ( i_sr = 0; i_sr < n_subregions[i_rg]; i_sr++ )
	{
	  for ( i_nt = 0; i_nt < n_neuron_types[i_rg][i_sr]; i_nt++ )
	    {
	      sd->rg[i_rg].sr[i_sr].nt[i_nt].n_neurons = 0;
	      count_lines = 0;

	      for (i_proc = 0; i_proc < n_proc; i_proc++)
		{
		  //set file path
		  sprintf(file_path, "../../np_p%d_%s_%s.dat", 
			  i_proc, regions[R_REGION], neuron_types[R_REGION][count_nt]);
		  
		  //file open again
		  if ((fp = fopen(file_path, "r")) == NULL) 
		    {
		      printf("file open error: get_GID_and_position: %s\n", file_path);
		      exit(EXIT_FAILURE);
		    }
		  		  
		  // change sd reladed data
		  while ( fgets(buf, sizeof(buf), fp) != NULL)
		    {
		      count_lines++;
		    }
		  fclose(fp);
		}
	      
	      sd->rg[i_rg].sr[i_sr].nt[i_nt].n_neurons = count_lines;
	      total_memory += sizeof(float) * sd->rg[i_rg].sr[i_sr].nt[i_nt].n_neurons * 4;
	      sd->rg[i_rg].sr[i_sr].nt[i_nt].neuron_positions = (float**) malloc2d( sizeof(float), sd->rg[i_rg].sr[i_sr].nt[i_nt].n_neurons, 4);
	      count_nt++;
	    }
	}
    }
  
  // get neuron positions
  for ( i_rg = 0; i_rg < N_REGIONS; i_rg++) 
    {
      unsigned int count_nt = 0;
      for ( i_sr = 0; i_sr < n_subregions[i_rg]; i_sr++ )
	{
	  for ( i_nt = 0; i_nt < n_neuron_types[i_rg][i_sr]; i_nt++ )
	    {
	      count_lines = 0;

	      for (i_proc = 0; i_proc < n_proc; i_proc++)
		{
		  //set file path
		  sprintf(file_path, "../../np_p%d_%s_%s.dat", 
			  i_proc, regions[R_REGION], neuron_types[R_REGION][count_nt]);
		  //file open again
		  if ((fp = fopen(file_path, "r")) == NULL) 
		    {
		      printf("file open error: get_GID_and_position: %s\n", file_path);
		      exit(EXIT_FAILURE);
		    }
		  
		  //count_lines = 0;
		  
		  // change sd reladed data
		  while ( fgets(buf, sizeof(buf), fp) != NULL)
		    {
		      // load to memory
		      n = sscanf(buf, "%f, %f, %f, %f", 
				 &(sd->rg[i_rg].sr[i_sr].nt[i_nt].neuron_positions[count_lines][0]), 
				 &(sd->rg[i_rg].sr[i_sr].nt[i_nt].neuron_positions[count_lines][1]), 
				 &(sd->rg[i_rg].sr[i_sr].nt[i_nt].neuron_positions[count_lines][2]), 
				 &(sd->rg[i_rg].sr[i_sr].nt[i_nt].neuron_positions[count_lines][3]) 
				 );
		      
		    if(n != 4 && n != -1)
		      {
			printf("n = %d, error\n", n);
		      }
		    
		    count_lines++;
		    }
		  fclose(fp);
		}
	      count_nt++;
	    }
	}
    }
  
  printf("# get_neuron_positions: finish\n");

}

void get_spike_times(struct data *sd)
{
  char file_path[256], s[256], buf[256];
  int i_rg = 0, i_sr = 0, i_nt = 0, i_spike = 0, i_cl = 0, n = 0, n_lines = 0, count_lines = 0, i_proc;

  printf("# reading spike times\n");

  // only read for count lines to get numbers of neurons
  for ( i_rg = 0; i_rg < N_REGIONS; i_rg++) 
    {
      unsigned int count_nt = 0;
      for ( i_sr = 0; i_sr < n_subregions[i_rg]; i_sr++ )
	{
	  for ( i_nt = 0; i_nt < n_neuron_types[i_rg][i_sr]; i_nt++ )
	    {
	      sd->rg[i_rg].sr[i_sr].nt[i_nt].n_spikes = 0;
	      count_lines = 0;

	      for (i_proc = 0; i_proc < n_proc; i_proc++)
		{
		  //set file path
		  sprintf(file_path, "../../st_p%d_%s_%s.dat", 
			  i_proc, regions[R_REGION], neuron_types[R_REGION][count_nt]);
		  
		  //file open again
		  if ((fp = fopen(file_path, "r")) == NULL) 
		    {
		      printf("file open error: spike times: %s\n", file_path);
		      exit(EXIT_FAILURE);
		    }
		  		  
		  // change sd reladed data
		  while ( fgets(buf, sizeof(buf), fp) != NULL)
		    {
		      count_lines++;
		    }
		  fclose(fp);
		}
	      
	      sd->rg[i_rg].sr[i_sr].nt[i_nt].n_spikes = count_lines;
	      sd->rg[i_rg].sr[i_sr].nt[i_nt].count_spikes = 0;
	      total_memory += sizeof(float) * sd->rg[i_rg].sr[i_sr].nt[i_nt].n_spikes + sizeof(unsigned int) * sd->rg[i_rg].sr[i_sr].nt[i_nt].n_spikes;
	      
	      printf("a rg:%d sr:%d nt:%d ns:%d\n", i_rg, i_sr, i_nt, sd->rg[i_rg].sr[i_sr].nt[i_nt].n_spikes);
	      
	      if (sd->rg[i_rg].sr[i_sr].nt[i_nt].n_spikes > 0)
		{
		  sd->rg[i_rg].sr[i_sr].nt[i_nt].spike_times = (float*) malloc( sizeof(float) * sd->rg[i_rg].sr[i_sr].nt[i_nt].n_spikes);
		  sd->rg[i_rg].sr[i_sr].nt[i_nt].Idx = (unsigned int*) malloc( sizeof(unsigned int) * sd->rg[i_rg].sr[i_sr].nt[i_nt].n_spikes);
		  for( i_cl = 0; i_cl < sd->rg[i_rg].sr[i_sr].nt[i_nt].n_spikes; i_cl++)
		    {
		      sd->rg[i_rg].sr[i_sr].nt[i_nt].spike_times[i_cl] = FLT_MAX;
		      sd->rg[i_rg].sr[i_sr].nt[i_nt].Idx[i_cl] = INT_MAX;
		    }
		}
	      count_nt++;

	    }
	}
    }
  
  // get spike times
  for ( i_rg = 0; i_rg < N_REGIONS; i_rg++) 
    {
      unsigned int count_nt = 0;
      for ( i_sr = 0; i_sr < n_subregions[i_rg]; i_sr++ )
	{
	  for ( i_nt = 0; i_nt < n_neuron_types[i_rg][i_sr]; i_nt++ )
	    {
	      count_lines = 0;
	      
	      for (i_proc = 0; i_proc < n_proc; i_proc++)
		{
		  //set file path
		  sprintf(file_path, "../../st_p%d_%s_%s.dat", 
			  i_proc, regions[R_REGION], neuron_types[R_REGION][count_nt]);
		  printf("b: %s\n", file_path);

		  //file open again
		  if ((fp = fopen(file_path, "r")) == NULL) 
		    {
		      printf("file open error: spike times: %s\n", file_path);
		      exit(EXIT_FAILURE);
		    }
		  
		  // change sd reladed data
		  while ( fgets(buf, sizeof(buf), fp) != NULL)
		    {
		      n = sscanf(buf, "%f, %d",
				 &(sd->rg[i_rg].sr[i_sr].nt[i_nt].spike_times[count_lines]),
				 &(sd->rg[i_rg].sr[i_sr].nt[i_nt].Idx[count_lines])
				 );
		      //printf("b rg:%d sr:%d nt:%d %d %f, %d\n", i_rg, i_sr, i_nt, count_lines, sd->rg[i_rg].sr[i_sr].nt[i_nt].spike_times[count_lines], sd->rg[i_rg].sr[i_sr].nt[i_nt].Idx[count_lines]);
		      
		      if(n != 2)
			{
			  printf("Error:%s is broken \n", file_path);
			  exit(EXIT_FAILURE);
			}
		      
		      /*
			if( (count_lines >= 1 )&& (sd->rg[i_rg].sr[i_sr].nt[i_nt].spike_times[count_lines-1] > sd->rg[i_rg].sr[i_sr].nt[i_nt].spike_times[count_lines]))
			{
			printf("Error: spike times are not ordeder according to time: %s %f > %f \n", file_path, sd->spike_times[i_rg][i_nt][count_lines-1], sd->spike_times[i_rg][i_nt][count_lines]);
			exit(EXIT_FAILURE);
			}
		      */
		      count_lines++;
		    }
		  fclose(fp);
		}
	      printf("count_lines:%d\n", count_lines);
	      count_nt++;
	    }
	}
    }
  
  printf("# get_spike_times: finish\n");
  printf("total_memory:%f\n", total_memory / pow(10,9)); 
}


void get_connection_patterns(struct data *sd)
{
  char file_path[256], s[256], buf[256];
  int i_pre_rg = 0,i_post_rg = 0, i_pre_nt = 0, i_post_nt = 0, n = 0, n_lines = 0, count_lines = 0;
  
  printf("# reading connections\n");

  for ( i_pre_rg = 0; i_pre_rg < N_REGIONS; i_pre_rg++)
    {
      for ( i_pre_nt = 0; i_pre_nt < n_neuron_types_per_region[i_pre_rg]; i_pre_nt++)
	{
	  for ( i_post_rg = 0; i_post_rg < N_REGIONS; i_post_rg++)
	    {
	      for ( i_post_nt = 0; i_post_nt < n_neuron_types_per_region[i_post_rg]; i_post_nt++)
		{
		  //initialize number of connections
		  sd->n_connections[i_pre_rg][i_pre_nt][i_post_rg][i_post_nt] = 0;
		  
		  //set file path
		  sprintf(file_path, "../../../input/%s/connection_%s_%s_%s_%s.dat",
			  data_directory_name,
			  regions[i_pre_rg], neuron_types[i_pre_rg][i_pre_nt], 
			  regions[i_post_rg], neuron_types[i_post_rg][i_post_nt] );
		  
		  if ((fp = fopen(file_path, "r")) != NULL)
		    {
		      count_lines=0;
		      while ( fgets(buf, sizeof(buf), fp) != NULL)
		  	{
		  	  //load to memory
		  	  n = sscanf(buf, "%d, %d",
		  		     &sd->connections[i_pre_rg][i_pre_nt][i_post_rg][i_post_nt][count_lines][0],
		  		     &sd->connections[i_pre_rg][i_pre_nt][i_post_rg][i_post_nt][count_lines][1]
		  		     );
		  	  if(n!=2 && n!=-1)
		  	    {
		  	      printf("n = %d, error\n", n);
		  	    }
		  	  count_lines++;
		  	}
		      sd->n_connections[i_pre_rg][i_pre_nt][i_post_rg][i_post_nt] = count_lines;

		      fclose(fp);

		    }
		}
	    }
	}
    }
  printf("# get_connection_patterns: finish\n");
}

void init(void)
{
  glClearColor(1., 1., 1., 1.0);
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

void DrawString(const char *str,void *font,float x,float y,float z)
{
  glRasterPos3f(x,y,z);
  glColor4d(0.7, 0.7, 0.5, 0.9);
  while(*str){
    glutBitmapCharacter(font, *str);
    ++str;
  }
}

void special_input(int key, int x, int y)
{
  unsigned int i_rg = 0, i_nt = 0;
  switch(key){
    
    // no assignment yet for left right
  case GLUT_KEY_LEFT:  
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
    printf("\rt = %d               ", t);
    fflush(stdout);
    break;

  case GLUT_KEY_RIGHT: (t < sd->simulation_time - 500) ? ( t += 500 ): (t);
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
      camera_position[2] -= (previous_y_drag_L - y) * 0.01;
      previous_y_drag_L = y;
    }

  //no assignment for right drag now
  if(right_button == 1) 
    {

    }
}

void SaveImage_PPM(char* fname)
{ // save current screen buffe image to the ppm file
  unsigned int ih, iw;
  int viewport[4];
  glGetIntegerv( GL_VIEWPORT, viewport);   // get location and size of viewport
  const unsigned int width = viewport[2];   // viewport width
  const unsigned int height = viewport[3];  //  viewport height
  void* image = malloc(3*width* height);  // allocate buffer (rgb)*width*height
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
  for(ih=0; ih < height; ih++){
    for(iw=0; iw < width; iw++){
      unsigned int i = (height-1-ih)*width+iw;
      int r = (unsigned char)img[i*3+0];
      int g = (unsigned char)img[i*3+1];
      int b = (unsigned char)img[i*3+2];
      fprintf(fp_output, "%d %d %d\n", r, g, b);
    }
  }
  fclose(fp_output);
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
	      sprintf(filename, "../../ppm_data/", data_directory_name);

	      if (mkdir(filename, 0777) )
		{
		  fprintf(stdout, "?");
		}
	      else{
		fprintf(stdout, "?");
	      }
	      switch_save_ppm_image = TRUE;

	      break;
	    }
	}
    }
  
}

void *malloc2d(size_t size, unsigned int n_row, unsigned int n_col)
{
    char **a, *b;
    unsigned long t = size * n_col;
    unsigned long i;
 
    a = (char**)malloc((sizeof(*a) + t) * n_row);
     
    if (a) {
        b = (char*)(a + n_row);
 
        for (i = 0; i < n_row; i++) {
            a[i] = b;
            b += t; 
        }
         return a;
    }
     return NULL;
}

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
 
 
int main(int argc, char *argv[])
{
  get_arguments(argc, argv);
  prepare_data_structure(&sd);
  get_model_position(sd);
  get_neuron_positions(sd);
  get_spike_times(sd);
  //get_connection_patterns(sd);
  
  glutInit(&argc, argv);
  
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
  glutInitWindowSize(WINDOW_WIDTH , WINDOW_HEIGHT);
  glutInitWindowPosition(100, 100);
  glutCreateWindow(argv[0]);
  glutDisplayFunc(display);

  glutReshapeFunc(resize);
  glutSpecialFunc(special_input);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
  glutKeyboardFunc(keyboard);
  glutReshapeFunc(resize);
  glutTimerFunc(1./(double)(play_speed), timer, 0);
  init();
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);	
  glutMainLoop();
  
  return 0;
}
