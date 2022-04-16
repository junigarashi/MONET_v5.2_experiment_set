CC = gcc
MPICC = mpicc
CFLAGS = -std=c11 -Wall -O2 -fopenmp -DMPI -Wno-unknown-pragmas -g3 -I/usr/X11R6/include
LDLIBS = -L/usr/X11R6/lib -lglut -lGLU -lGL -lXmu -lXi -lXext -lX11 -lm -lpthread 


grid_vis: grid_vis.o grid_vis_tool.o
	$(MPICC) $(CFLAGS) -o $@ $^ -lm $(LDLIBS)

grid_vis_tool.o: grid_vis_tool.c
	$(MPICC) $(CFLAGS)  -o $@ -c $^ $(LDLIBS)

clean:
	rm -f grid_vis	 *~ *.o *.lst
