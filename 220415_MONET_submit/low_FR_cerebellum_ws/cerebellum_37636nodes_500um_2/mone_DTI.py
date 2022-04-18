#!/bin/env python
import json
#import numpy as np
#import matplotlib.pylab as pl

class voxel_sections_on_flatmap:
    """ voxel sections on flatmap """
    def __init__(self):   
        self.position=[]
        self.terminals=[]

def read_flatmap_data():
    pass

def read_DTI_data(sd, rd_array, process):

    # tetntativly generate pseudo data
    length_on_a_side_of_voxel = 50.

    # visualization
    """
    fig = pl.figure()
    fig.add_subplot(111, aspect='equal')
    """
    
    # generate artificial voxel surface on plane per MPI process.
    for i_proc, proc in enumerate(process):
        n_voxel_on_a_side = int(rd_array[proc.region_GID]["structure_info"]["xy_length_per_tile"]/length_on_a_side_of_voxel)

        if i_proc==0:
            color="black"
            vis_gap=0.
        else:
            color="red"
            vis_gap=500.
            
        proc.voxel_sections = [ [] for x in xrange(n_voxel_on_a_side*n_voxel_on_a_side)]
        # generate coordinates of artificial voxel segments on flatmap by n_voxels x n_voxels
        i_voxel_per_process=0
        for xv in range(n_voxel_on_a_side):
            for yv in range(n_voxel_on_a_side):
                proc.voxel_sections[i_voxel_per_process] = voxel_sections_on_flatmap()
                
                x0_element = xv * length_on_a_side_of_voxel + proc.spatial_extent[0]
                y0_element = yv * length_on_a_side_of_voxel + proc.spatial_extent[1]
                x1_element = (xv+1) * length_on_a_side_of_voxel + proc.spatial_extent[0]
                y1_element = yv * length_on_a_side_of_voxel + proc.spatial_extent[1]
                x2_element = (xv+1) * length_on_a_side_of_voxel + proc.spatial_extent[0]
                y2_element = (yv+1) * length_on_a_side_of_voxel + proc.spatial_extent[1]
                x3_element = xv * length_on_a_side_of_voxel + proc.spatial_extent[0]
                y3_element = (yv+1) * length_on_a_side_of_voxel + proc.spatial_extent[1]

                proc.voxel_sections[i_voxel_per_process].position = [[x0_element, y0_element],
                                                             [x1_element, y1_element],
                                                             [x2_element, y2_element],
                                                             [x3_element, y3_element]]
                i_voxel_per_process+=1

                #visualization
                """
                pl.plot([x0_element+vis_gap, x1_element+vis_gap], [y0_element+vis_gap, y1_element+vis_gap], "-", color=color)
                pl.plot([x1_element+vis_gap, x2_element+vis_gap], [y1_element+vis_gap, y2_element+vis_gap], "-", color=color)
                pl.plot([x2_element+vis_gap, x3_element+vis_gap], [y2_element+vis_gap, y3_element+vis_gap], "-", color=color)
                pl.plot([x3_element+vis_gap, x0_element+vis_gap], [y3_element+vis_gap, y0_element+vis_gap], "-", color=color)
                """
                
    #pl.xlim([-100, 2600])
    #pl.ylim([-100, 1600])

    # Pseudo DTI connections
    # connection process 0 -> process 1
    for i_proc, proc in enumerate(process):
        if proc.region_GID == 1:
            for i_vs, voxel_sections in enumerate(proc.voxel_sections):
                proc.voxel_sections[i_vs].terminals.append(
                    [(process[1].voxel_sections[i_vs].position[0][0] + process[1].voxel_sections[i_vs].position[1][0]) / 2.,
                     (process[1].voxel_sections[i_vs].position[0][1] + process[1].voxel_sections[i_vs].position[2][1]) / 2.]
                )
                
                # voxel combination for searching 
                
                
                # visualization
                """
                pl.plot([proc.voxel_sections[i_vs].position[0][0], proc.voxel_sections[i_vs].terminals[0][0][0]+vis_gap],
                        [proc.voxel_sections[i_vs].position[0][1], proc.voxel_sections[i_vs].terminals[0][0][1]+vis_gap], "-", color="blue", lw=0.5)
                """
    # visualization
    #pl.show()
