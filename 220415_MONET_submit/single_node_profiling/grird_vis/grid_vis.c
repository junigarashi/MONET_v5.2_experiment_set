#include "./grid_vis.h"

int main(int argc, char **argv)
{
  get_arguments(argc, argv);
  read_global_shared_info(&sp);
  read_process_info(&sp, &PI);
  read_region_info(&sp, &PI);

  get_firing_rates_with_time(&sp, &PI);

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
  glClearColor(1., 1., 1., 1.0);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);	
  glutMainLoop();
 
}


