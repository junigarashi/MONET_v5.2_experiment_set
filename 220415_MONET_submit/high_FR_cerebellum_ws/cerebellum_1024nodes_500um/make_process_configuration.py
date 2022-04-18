#!/bin/env python

#import matplotlib.pyplot as plt
import multiprocessing as mp
from multiprocessing import Process
import time
#import numpy as np
import math
import os
import argparse

# mone_in.py
import mone_in

# mone_conn.py
import mone_conn

# mone_io.py
import mone_out

# mone_DTI.py
import mone_DTI

# mone_vis.py
try:
    import vtk
    switch_vtk_availablity = True
    import mone_vis

except:
    switch_vtk_availablity = False

# Abbribiations and term
# GID: Global ID,  ID numbered across process
# LID: Local ID, local ID numbered within group within one process without numbering across groups
# SerLID: Serial local ID, serial ID numbered across groups within one process

# region: rg, typical brain regions as broadmann's 50 brain regions
# subregion: srg, sr, typical modules or sub-structures in brain region as layers or nucleus. Not for column in cortex
# neuron type: nt, neuron types, which are clearly classified, as pyramidal cell, palvalbmin interneuron

# bundle: collection of connections from one region to another region
# conenction type: classification for connection group defined by specific combination of presynaptic and postsynaptic neuron types

# sd: system data read from system.json
# rd: region data read from region_X.json

# MPI process
class process_info:
    """ MPI process info class"""
    def __init__(self):   
        self.spatial_extent = []
        self.center_position = []
        self.xy_plane_4vertices = []
        self.region_name = ""
        self.process_GID = []
        self.intra_regional_process_pairs = []
        self.intra_regional_pre_post_address = []
        self.intra_regional_pre_post_address2 = []
        self.intra_regional_pre_post_address3 = []
        self.inter_regional_process_pairs = []
        self.inter_regional_pre_post_address = []
        self.bundle = {}
        self.connection_type_SerLIDs_per_post = []
        self.i_connection_type_SerLIDs_per_post = 0
        self.ct = {}
        self.bdl = {}
        self.count_ct = 0
        self.intra_regional_connection_nt_pairs = []
        self.inter_regional_connection_nt_pairs = []
        self.intra_regional_connection_process = []
        self.inter_regional_connection_process = []
        self.inter_regional_connection_pre_process = []
        self.inter_regional_connection_post_process = []
        self.n_bundles = 0
        self.count_connection_type_LID_in_bundle = 0
        self.n_indegree_inter_regional_connection_types_per_process = 0
        self.n_indegree_intra_regional_connection_types_per_process = 0
        self.cluster_center = []
        self.cluster_upper_plane = []
        self.cluster_lower_plane = []
        
def set_process_info(sd, rd_array):
    process = []
    if sd["mode"] == "artificial_square" or sd["mode"] == "minimum_artificial_square":
        process = set_process_info_artificial_square_sheet_1d_array(sd, rd_array)
    elif sd["mode"] == "human_cortical_sheet":
        process = set_process_info_human_cortical_sheet(sd, rd_array)
    elif sd["mode"] == "corticothalamic_circuit_for_cerebellum":
        process = set_process_info_corticothalamic_circuit_for_cerebellum(sd, rd_array)
    elif sd["mode"] == "cerebellum":
        process = set_process_info_artificial_square_sheet_1d_array(sd, rd_array)
    elif sd["mode"] == "minimum_test" or sd["mode"] == "example_colocalization_channels" or sd["mode"] == "example_STDP" or sd["mode"] == "example_HTC":
        process = set_process_info_minimum_test(sd, rd_array)
    elif sd["mode"] == "example_minimum_thalamic_neuron":
        process = set_process_info_artificial_square_sheet_1d_array(sd, rd_array)
    elif sd["mode"] == "cortico_thalamo_cerebello_sheets":
        process = set_process_info_cortico_thalamo_cerebello_sheets(sd, rd_array)
    elif sd["mode"] == "DTI_test":
        process = set_process_info_artificial_square_sheet_1d_array(sd, rd_array)
    elif sd["mode"] == "cluster_test":
        process = set_process_info_artificial_square_sheet_1d_array(sd, rd_array)
    elif sd["mode"] == "ID_connect_example":
        print("set_process_info: ID_connect")
        process = set_process_info_artificial_square_sheet_1d_array(sd, rd_array)
    elif sd["mode"] == "PFCFP_test":
        process = set_process_info_artificial_square_sheet_1d_array(sd, rd_array)
    else:
        print(sd["mode"] + " is not supported" )
    return process

def set_process_info_artificial_square_sheet_1d_array(sd, rd_array):
    process = [ [] for x in range(sd["n_process"])]
    
    # set corresponding process IDs to each region data
    x_region_width = sd["x_points"] / sd["n_x_regions"]
    y_region_width = sd["y_points"] / sd["n_y_regions"]

    # loop regions
    for i_rn, rn in enumerate(sd["region_names"]):

        # x y rank of regions
        i_x_region = i_rn % sd["n_x_regions"]
        i_y_region = (i_rn / sd["n_y_regions"]) % sd["n_y_regions"]

        # rank of tile in x and y axis 
        min_x_rank_tile = i_x_region * x_region_width
        max_x_rank_tile = (i_x_region + 1) * x_region_width
        min_y_rank_tile = i_y_region * y_region_width
        max_y_rank_tile = (i_y_region + 1) * y_region_width
        
        rd_array[i_rn]["structure_info"]["n_tiles_on_a_x_side"] = x_region_width 
        rd_array[i_rn]["structure_info"]["n_tiles_on_a_y_side"] = y_region_width 
        
        # when tile is larger than 32x32 tile, tile group is prepared
        if max_x_rank_tile - min_x_rank_tile >= 32 and max_y_rank_tile - min_y_rank_tile >= 32:
            
            minimum_tile_group_length = int(math.ceil(rd_array[i_rn]["structure_info"]["tile_link_limit"] / float(rd_array[i_rn]["structure_info"]["xy_length_per_tile"])))
            #tile_group_length = rd_array[i_rn]["structure_info"]["tile_link_limit"] / float(rd_array[i_rn]["structure_info"]["xy_length_per_tile"])
            # take both direction of length of tile_link_limit and one tile length 
            #tile_group_length = 2 * int(math.ceil( rd_array[i_rn]["structure_info"]["xy_length_per_tile"] / float(rd_array[i_rn]["structure_info"]["tile_link_limit"]))) + 1
            #print "mtgl, tgl:", minimum_tile_group_length, tile_group_length
            
            #if tile_group_length < minimum_tile_group_length:
            #    print "Please set larger tile_group_length than ", minimum_tile_group_length, ", current value is", tile_group_length
            #    exit()
                
            # take both direction of length of tile_link_limit and one tile length.
            # The one tile length is for center one, multiply 2 is both sides from edges of the center one.
            # to take sufficient tiles covering the tile link limit,
            # round up tile link limit by length on a side of tiel by math.ceil
            rd_array[i_rn]["structure_info"]["n_points_on_a_side_of_tile_group"] = 2 * int(math.ceil( rd_array[i_rn]["structure_info"]["tile_link_limit"] / float(rd_array[i_rn]["structure_info"]["xy_length_per_tile"]) )) + 1

            # round up n_tiles on a x (y) side / n tiles on a side of tile group
            # to get full length and not full length (at edge) tile group  
            rd_array[i_rn]["structure_info"]["n_tile_groups"] = int(math.ceil( rd_array[i_rn]["structure_info"]["n_tiles_on_a_x_side"] / float(rd_array[i_rn]["structure_info"]["n_points_on_a_side_of_tile_group"]))
                                                                   * math.ceil( rd_array[i_rn]["structure_info"]["n_tiles_on_a_y_side"] / float(rd_array[i_rn]["structure_info"]["n_points_on_a_side_of_tile_group"] )) )

            rd_array[i_rn]["structure_info"]["n_tile_groups_on_a_x_side"] = int(math.ceil( rd_array[i_rn]["structure_info"]["n_tiles_on_a_x_side"] / float(rd_array[i_rn]["structure_info"]["n_points_on_a_side_of_tile_group"] )))
            rd_array[i_rn]["structure_info"]["n_tile_groups_on_a_y_side"] = int(math.ceil( rd_array[i_rn]["structure_info"]["n_tiles_on_a_y_side"] / float(rd_array[i_rn]["structure_info"]["n_points_on_a_side_of_tile_group"] )))
            #print "ntg, ntg_os", rd_array[i_rn]["structure_info"]["n_tile_groups"], rd_array[i_rn]["structure_info"]["n_tile_groups_on_a_x_side"]
            rd_array[i_rn]["structure_info"]["tile_group"] = [ [] for i in range(rd_array[i_rn]["structure_info"]["n_tile_groups"])]
            rd_array[i_rn]["structure_info"]["tile_group_pair"] = [ [] for i in range(rd_array[i_rn]["structure_info"]["n_tile_groups"])]
            rd_array[i_rn]["structure_info"]["tile_group_center_position"] = [ [] for i in range(rd_array[i_rn]["structure_info"]["n_tile_groups"])]
            rd_array[i_rn]["structure_info"]["xy_plane_4vertices"] = [ [ [] for j in range(4) ] for i in range(rd_array[i_rn]["structure_info"]["n_tile_groups"]) ]
            
            # positions of tile groups
            for i_tg in range(rd_array[i_rn]["structure_info"]["n_tile_groups"]):
                # center position of tile group
                rd_array[i_rn]["structure_info"]["tile_group_center_position"][i_tg] = [ (i_tg % rd_array[i_rn]["structure_info"]["n_tile_groups_on_a_x_side"])
                                                                                         * rd_array[i_rn]["structure_info"]["n_points_on_a_side_of_tile_group"]
                                                                                         + rd_array[i_rn]["structure_info"]["n_points_on_a_side_of_tile_group"]/2
                                                                                         + min_x_rank_tile * rd_array[i_rn]["structure_info"]["xy_length_per_tile"],
                                                                                         i_tg / rd_array[i_rn]["structure_info"]["n_tile_groups_on_a_x_side"]
                                                                                         * rd_array[i_rn]["structure_info"]["n_points_on_a_side_of_tile_group"]
                                                                                         + rd_array[i_rn]["structure_info"]["n_points_on_a_side_of_tile_group"]/2
                                                                                         + min_y_rank_tile * rd_array[i_rn]["structure_info"]["xy_length_per_tile"] ]
                
                # four edge vetices of tile group
                rd_array[i_rn]["structure_info"]["xy_plane_4vertices"][i_tg][0] = [ (( 0 + (i_tg % rd_array[i_rn]["structure_info"]["n_tile_groups_on_a_x_side"]))
                                                                                    * rd_array[i_rn]["structure_info"]["n_points_on_a_side_of_tile_group"]
                                                                                    + min_x_rank_tile) * rd_array[i_rn]["structure_info"]["xy_length_per_tile"],
                                                                                    (( 0 + (i_tg / rd_array[i_rn]["structure_info"]["n_tile_groups_on_a_x_side"]) )
                                                                                    * rd_array[i_rn]["structure_info"]["n_points_on_a_side_of_tile_group"]
                                                                                    + min_y_rank_tile) * rd_array[i_rn]["structure_info"]["xy_length_per_tile"] ]
                rd_array[i_rn]["structure_info"]["xy_plane_4vertices"][i_tg][1] = [ (( 1 + (i_tg % rd_array[i_rn]["structure_info"]["n_tile_groups_on_a_x_side"]))
                                                                                    * rd_array[i_rn]["structure_info"]["n_points_on_a_side_of_tile_group"]
                                                                                    + min_x_rank_tile) * rd_array[i_rn]["structure_info"]["xy_length_per_tile"],
                                                                                    (( 0 + (i_tg / rd_array[i_rn]["structure_info"]["n_tile_groups_on_a_x_side"]) )
                                                                                    * rd_array[i_rn]["structure_info"]["n_points_on_a_side_of_tile_group"]
                                                                                    + min_y_rank_tile) * rd_array[i_rn]["structure_info"]["xy_length_per_tile"] ]
                rd_array[i_rn]["structure_info"]["xy_plane_4vertices"][i_tg][2] = [ (( 1 + (i_tg % rd_array[i_rn]["structure_info"]["n_tile_groups_on_a_x_side"]))
                                                                                    * rd_array[i_rn]["structure_info"]["n_points_on_a_side_of_tile_group"]
                                                                                    + min_x_rank_tile) * rd_array[i_rn]["structure_info"]["xy_length_per_tile"],
                                                                                    (( 1 + (i_tg / rd_array[i_rn]["structure_info"]["n_tile_groups_on_a_x_side"]) )
                                                                                    * rd_array[i_rn]["structure_info"]["n_points_on_a_side_of_tile_group"]
                                                                                    + min_y_rank_tile) * rd_array[i_rn]["structure_info"]["xy_length_per_tile"] ]
                rd_array[i_rn]["structure_info"]["xy_plane_4vertices"][i_tg][3] = [ (( 0 + (i_tg % rd_array[i_rn]["structure_info"]["n_tile_groups_on_a_x_side"]))
                                                                                    * rd_array[i_rn]["structure_info"]["n_points_on_a_side_of_tile_group"]
                                                                                    + min_x_rank_tile) * rd_array[i_rn]["structure_info"]["xy_length_per_tile"],
                                                                                    (( 1 + (i_tg / rd_array[i_rn]["structure_info"]["n_tile_groups_on_a_x_side"]) )
                                                                                    * rd_array[i_rn]["structure_info"]["n_points_on_a_side_of_tile_group"]
                                                                                    + min_y_rank_tile) * rd_array[i_rn]["structure_info"]["xy_length_per_tile"] ]
                

        else:
            rd_array[i_rn]["structure_info"]["n_tile_groups"] = 0
        
        # match x y coorditates of tiles and spatial extent of region, and register corresponding process ID   
        rd_array[i_rn]["process_ID"]= []
        for x in range(sd["x_points"]):
            for y in range(sd["y_points"]):
                if (min_x_rank_tile <= x < max_x_rank_tile) and (min_y_rank_tile <= y < max_y_rank_tile):
                    process_ID = y * sd["x_points"] + x
                   
                    rd_array[i_rn]["process_ID"].append(process_ID)
        rd_array[i_rn]["process_ID"].sort()    
        
    # assign region to process and set positions of process in cortical sheet
    sd["region_y_points"] = sd["y_points"] / sd["n_y_regions"]
    sd["region_x_points"] = sd["x_points"] / sd["n_x_regions"]

    # set process spatial info 
    for y in range(sd["y_points"]):
        i_y_region = int(y / sd["region_y_points"])
        for x in range(sd["x_points"]):
            i_x_region = int(x / sd["region_x_points"])
            
            i_proc = y * sd["x_points"] + x

            
           # class instance
            process[i_proc] = process_info()
            process[i_proc].process_GID = i_proc
            process[i_proc].region_GID = sd["n_x_regions"] * i_y_region + i_x_region
            process[i_proc].region_name = sd["region_names"][process[i_proc].region_GID]
            process[i_proc].inter_regional_process_pairs = [[] for i in range(len(rd_array[process[i_proc].region_GID]["inter_regional_connection"])) ]
            process[i_proc].inter_regional_connection_process = [[] for i in range(len(rd_array[process[i_proc].region_GID]["inter_regional_connection"])) ]
            process[i_proc].inter_regional_pre_post_address = [[] for i in range(len(rd_array[process[i_proc].region_GID]["inter_regional_connection"])) ]

            # vertices (x, y), (x+1, y+1)
            process[i_proc].spatial_extent = [ x * rd_array[process[i_proc].region_GID]["structure_info"]["xy_length_per_tile"] ,
                                               y *  rd_array[process[i_proc].region_GID]["structure_info"]["xy_length_per_tile"],
                                               rd_array[process[i_proc].region_GID]["structure_info"]["region_z_extent"][0],
                                               (x + 1) *  rd_array[process[i_proc].region_GID]["structure_info"]["xy_length_per_tile"],
                                               (y + 1) *  rd_array[process[i_proc].region_GID]["structure_info"]["xy_length_per_tile"],
                                               rd_array[process[i_proc].region_GID]["structure_info"]["region_z_extent"][1]]

            process[i_proc].center_position = [ (x + 0.5)*  rd_array[process[i_proc].region_GID]["structure_info"]["xy_length_per_tile"],
                                                (y + 0.5) *  rd_array[process[i_proc].region_GID]["structure_info"]["xy_length_per_tile"],
                                                (rd_array[process[i_proc].region_GID]["structure_info"]["region_z_extent"][0] + rd_array[process[i_proc].region_GID]["structure_info"]["region_z_extent"][1]) / 2 ]

            process[i_proc].xy_plane_4vertices = [ [ (x + 0) *  rd_array[process[i_proc].region_GID]["structure_info"]["xy_length_per_tile"],
                                                     (y + 0) *  rd_array[process[i_proc].region_GID]["structure_info"]["xy_length_per_tile"],
                                                     (rd_array[process[i_proc].region_GID]["structure_info"]["region_z_extent"][0] + rd_array[process[i_proc].region_GID]["structure_info"]["region_z_extent"][1]) / 2 ],
                                                   [ (x + 1) *  rd_array[process[i_proc].region_GID]["structure_info"]["xy_length_per_tile"],
                                                     (y + 0) *  rd_array[process[i_proc].region_GID]["structure_info"]["xy_length_per_tile"],
                                                     (rd_array[process[i_proc].region_GID]["structure_info"]["region_z_extent"][0] + rd_array[process[i_proc].region_GID]["structure_info"]["region_z_extent"][1]) / 2 ],
                                                   [ (x + 1) *  rd_array[process[i_proc].region_GID]["structure_info"]["xy_length_per_tile"],
                                                     (y + 1) *  rd_array[process[i_proc].region_GID]["structure_info"]["xy_length_per_tile"],
                                                     (rd_array[process[i_proc].region_GID]["structure_info"]["region_z_extent"][0] + rd_array[process[i_proc].region_GID]["structure_info"]["region_z_extent"][1]) / 2 ],
                                                   [ (x + 0) *  rd_array[process[i_proc].region_GID]["structure_info"]["xy_length_per_tile"],
                                                     (y + 1) *  rd_array[process[i_proc].region_GID]["structure_info"]["xy_length_per_tile"],
                                                     (rd_array[process[i_proc].region_GID]["structure_info"]["region_z_extent"][0] + rd_array[process[i_proc].region_GID]["structure_info"]["region_z_extent"][1]) / 2 ] ]

            # calculate TofuD coordinates based on the numbers of CMGs per CPU for each process
            # A64fx has 4 CMGs. 2 x 2 process are mapped on 2 by 2 xy plane.
            TOFUD_2D_coordinate_x = x / 2
            TOFUD_2D_coordinate_y = y / 2
            process[i_proc].TOFUD_2D_coordinate = [TOFUD_2D_coordinate_x, TOFUD_2D_coordinate_y]
                                                
            # cumulative depth of layer
            lt_sum=0
            initial_depth, final_depth = [], []
            for lt in rd_array[process[i_proc].region_GID]["structure_info"]["layer_thickness"]:
                initial_depth.append(lt_sum)
                lt_sum += lt
                final_depth.append(lt_sum)

            process[i_proc].subregion_positions = [ [process[i_proc].spatial_extent[0],
                                                     process[i_proc].spatial_extent[1],
                                                     process[i_proc].spatial_extent[2] + id,
                                                     process[i_proc].spatial_extent[3],
                                                     process[i_proc].spatial_extent[4],
                                                     process[i_proc].spatial_extent[2] + fd ] for id, fd in zip(initial_depth, final_depth) ]

            # when tile is larger than 32x32 tile, tile group is prepared
            if rd_array[process[i_proc].region_GID]["structure_info"]["n_tile_groups"] > 0:
                process[i_proc].i_x_tile_group = int(math.floor( x / float(rd_array[i_rn]["structure_info"]["n_points_on_a_side_of_tile_group"])))
                process[i_proc].i_y_tile_group = int(math.floor( y / float(rd_array[i_rn]["structure_info"]["n_points_on_a_side_of_tile_group"])))
                process[i_proc].i_tile_group = process[i_proc].i_y_tile_group * rd_array[i_rn]["structure_info"]["n_tile_groups_on_a_x_side"] + process[i_proc].i_x_tile_group
                rd_array[process[i_proc].region_GID]["structure_info"]["tile_group"][process[i_proc].i_tile_group].append( i_proc )

    #print "tile_group", rd_array[0]["structure_info"]["tile_group"], "len", len(rd_array[0]["structure_info"]["tile_group"])
    #for i_tg, tg in enumerate(rd_array[0]["structure_info"]["tile_group"]):
    #    print i_tg, tg, len(tg)
    
    # register position of regions
    for i_y_region in range(sd["n_y_regions"]):
        for i_x_region in range(sd["n_x_regions"]):
            i_region = sd["n_x_regions"] * i_y_region + i_x_region
            start_x_region =  rd_array[i_region]["structure_info"]["xy_length_per_tile"] * sd["region_x_points"] * i_x_region
            start_y_region =  rd_array[i_region]["structure_info"]["xy_length_per_tile"] * sd["region_y_points"] * i_y_region
            end_x_region =  rd_array[i_region]["structure_info"]["xy_length_per_tile"] * (i_x_region + 1)
            end_y_region =  rd_array[i_region]["structure_info"]["xy_length_per_tile"] * (i_y_region + 1)
            rd_array[i_region]["structure_info"]["region_center"]=[start_x_region +  rd_array[i_region]["structure_info"]["xy_length_per_tile"] * sd["region_x_points"] / 2,
                                                                    start_y_region +  rd_array[i_region]["structure_info"]["xy_length_per_tile"] * sd["region_y_points"] / 2,
                                                                    (rd_array[i_region]["structure_info"]["region_z_extent"][0] + rd_array[i_region]["structure_info"]["region_z_extent"][1]) / 2]
    
            
    return process

def set_process_info_human_cortical_sheet(sd, rd_array):
    process = [[] for x in range(sd["n_process"])]
    
    sheet_cortex_coordinates = []
    sheet_cortex_coordinates_x = []
    sheet_cortex_coordinates_y = []
    grid_points = []
    grid_points_x = []
    grid_points_y = []
    grid_span = int( rd_array[process[i].region_GID]["structure_info"]["xy_length_per_tile"]/1000)

    # read vertex data files
    f = open("S1vtx.tsv", "r")
    lines = f.readlines() 
    f.close()

    for line in lines:
        x, y, z = line.split("\t")
        if float(z) == 0.:
            sheet_cortex_coordinates.append([float(x), float(y)])
            sheet_cortex_coordinates_x.append(float(x))
            sheet_cortex_coordinates_y.append(float(y))

    # check spatial extent
    max_x = max(sheet_cortex_coordinates_x)
    min_x = min(sheet_cortex_coordinates_x)
    max_y = max(sheet_cortex_coordinates_y)
    min_y = min(sheet_cortex_coordinates_y)

    print(min_x, max_x, min_y, max_y)
    print("n_vertex", len(sheet_cortex_coordinates))

    # segament cortical sheet
    for i_scc, scc in enumerate(sheet_cortex_coordinates):
        x = int(scc[0]) - (int(scc[0])%grid_span)
        y = int(scc[1]) - (int(scc[1])%grid_span)
        if [ x, y ] not in grid_points:
            grid_points.append([x, y])
            grid_points_x.append(x)
            grid_points_y.append(y)
            
    # plot tiles
    plt.plot(sheet_cortex_coordinates_x, sheet_cortex_coordinates_y, "k.", ms=0.3)
    plt.plot(grid_points_x, grid_points_y, "r.")
    plt.show()

    # fill holes
    count_holes=0
    holes=[]
    holes_x=[]
    holes_y=[]
    for i, (x, y) in enumerate(grid_points):
        upper, lower, right, left, center = False, False, False, False, False
        if i % 10000 == 0:
            print(i, len(grid_points) )
            
        for other_x, other_y in grid_points:
            if x == other_x and y-2*grid_span == other_y:
                lower = True
            if x+1*grid_span == other_x and y-1*grid_span == other_y:
                right = True
            if x-1*grid_span == other_x and y-1*grid_span == other_y:
                left = True
            if x == other_x and y-1*grid_span == other_y:
                center = True
            
        if lower == True and right == True and left == True and center == False:
            holes.append([x, y-1*grid_span])
            holes_x.append(x)
            holes_y.append(y-1*grid_span)
            count_holes+=1
    grid_points += holes
    grid_points_x += holes_x
    grid_points_y += holes_y
    
    # make tiles from grid points
    tiles = []
    for i, (x, y) in enumerate(grid_points):
        ur, ll, lr = False, False, False
        
        for other_x, other_y in grid_points:
            if x+1*grid_span == other_x and y == other_y:
                ur = True
            if x == other_x and y-1*grid_span == other_y:
                ll = True
            if x+1*grid_span == other_x and y-1*grid_span == other_y:
                lr = True
            if ur == True and ll == True and lr == True:
                tiles.append([[x, y], [x+1*grid_span, y], [x, y-1*grid_span], [x+1*grid_span, y-1*grid_span]])
                break
    print( "n_tiles", len(tiles))

    # sort tilse using left upper angle
    tiles=sorted(tiles, key=lambda x: (x[0][0], x[0][1]) )
    
    # plot tiles
    plt.plot(grid_points_x, grid_points_y, "r.")
    for t in tiles:
        plt.monplot( (t[0][0]+t[3][0])/2., (t[0][1]+t[3][1])/2., "k.")
    plt.show()
    
    # assign tiles to processes
    for i_t, t in enumerate(tiles):
        # underconstruction
        if sd["n_process"] < i_t:
            break
        region_GID=0
        process[i_t] = process_info()
        process[i_t].process_GID = i_t
        process[i_t].region_GID = region_GID
        process[i_t].region_name = sd["region_names"][region_GID]
        process[i_t].spatial_extent = [ t[0][0],
                                        t[0][1],
                                        rd_array[region_GID]["structure_info"]["region_z_extent"][0],
                                        t[3][0],
                                        t[3][1],                                     
                                        rd_array[region_GID]["structure_info"]["region_z_extent"][1]]
    return process

def set_process_info_corticothalamic_circuit_for_cerebellum(sd, rd_array):
    process = [ [] for x in range(sd["n_process"])]
    
    for i in range(sd["n_process"]):
        process[i] = process_info()
        process[i].process_GID = i

    # set details by hand
    # M1
    i=0
    process[i].region_GID = 0
    process[i].region_name = "M1"
    process[i].region_GID = sd["region_name_to_region_GID"][process[i].region_name]
    x = 0
    y = 0
    process[i].spatial_extent =[ x *  rd_array[process[i].region_GID]["structure_info"]["xy_length_per_tile"] ,
                                 y *  rd_array[process[i].region_GID]["structure_info"]["xy_length_per_tile"],
                                 rd_array[process[i].region_GID]["structure_info"]["region_z_extent"][0],
                                 (x+1) *  rd_array[process[i].region_GID]["structure_info"]["xy_length_per_tile"],
                                 (y+1) *  rd_array[process[i].region_GID]["structure_info"]["xy_length_per_tile"],
                                 rd_array[process[i].region_GID]["structure_info"]["region_z_extent"][1]]
    
    process[i].center_position=[ (x+0.5)*  rd_array[process[i].region_GID]["structure_info"]["xy_length_per_tile"],
                                     (y+0.5) *  rd_array[process[i].region_GID]["structure_info"]["xy_length_per_tile"],
                                     (rd_array[process[i].region_GID]["structure_info"]["region_z_extent"][0] + rd_array[process[i].region_GID]["structure_info"]["region_z_extent"][1])/2 ]
    
    # cumulative depth of layer
    lt_sum=0
    initial_depth, final_depth = [], []
    for lt in rd_array[process[i].region_GID]["structure_info"]["layer_thickness"]:
        initial_depth.append(lt_sum)
        lt_sum += lt
        final_depth.append(lt_sum)

        process[i].subregion_positions = [ [process[i].spatial_extent[0],
                                               process[i].spatial_extent[1],
                                               process[i].spatial_extent[2] + id,
                                               process[i].spatial_extent[3],
                                               process[i].spatial_extent[4],
                                               process[i].spatial_extent[2] + fd ] for id, fd in zip(initial_depth, final_depth) ]
        
    # M1-Th
    i=1
    process[i].region_GID = 1
    process[i].region_name = "M1-Th"
    process[i].region_GID = sd["region_name_to_region_GID"][process[i].region_name]
    x = 0
    y = 0
    process[i].spatial_extent =[ x *  rd_array[process[i].region_GID]["structure_info"]["xy_length_per_tile"] ,
                                 y *  rd_array[process[i].region_GID]["structure_info"]["xy_length_per_tile"],
                                 rd_array[process[i].region_GID]["structure_info"]["region_z_extent"][0],
                                 (x+1) *  rd_array[process[i].region_GID]["structure_info"]["xy_length_per_tile"],
                                 (y+1) *  rd_array[process[i].region_GID]["structure_info"]["xy_length_per_tile"],
                                 rd_array[process[i].region_GID]["structure_info"]["region_z_extent"][1]]
    
    process[i].center_position=[ (x+0.5)*  rd_array[process[i].region_GID]["structure_info"]["xy_length_per_tile"],
                                     (y+0.5) *  rd_array[process[i].region_GID]["structure_info"]["xy_length_per_tile"],
                                     (rd_array[process[i].region_GID]["structure_info"]["region_z_extent"][0] + rd_array[process[i].region_GID]["structure_info"]["region_z_extent"][1])/2 ]
    
    # cumulative depth of layer
    lt_sum=0
    initial_depth, final_depth = [], []
    for lt in rd_array[process[i].region_GID]["structure_info"]["layer_thickness"]:
        initial_depth.append(lt_sum)
        lt_sum += lt
        final_depth.append(lt_sum)

        process[i].subregion_positions = [ [process[i].spatial_extent[0],
                                               process[i].spatial_extent[1],
                                               process[i].spatial_extent[2] + id,
                                               process[i].spatial_extent[3],
                                               process[i].spatial_extent[4],
                                               process[i].spatial_extent[2] + fd ] for id, fd in zip(initial_depth, final_depth) ]
    return process

def set_process_info_minimum_test(sd, rd_array):
    process = [ [] for x in range(sd["n_process"])]

    # set corresponding process IDs to each region data
    region_GID=0
    process_GID=0
    rd_array[region_GID]["process_ID"] = [process_GID]
    region_GID=1
    process_GID=1
    rd_array[region_GID]["process_ID"] = [process_GID]
    
    for i in range(sd["n_process"]):
        process[i] = process_info()
        process[i].process_GID = i

    # set details by hand
    # test1
    i=0
    process[i].region_GID = 0
    process[i].region_name = "test1"
    process[i].region_GID = sd["region_name_to_region_GID"][process[i].region_name]
    x = 0
    y = 0
    process[i].spatial_extent =[ x *  rd_array[process[i].region_GID]["structure_info"]["xy_length_per_tile"] ,
                                 y *  rd_array[process[i].region_GID]["structure_info"]["xy_length_per_tile"],
                                 rd_array[process[i].region_GID]["structure_info"]["region_z_extent"][0],
                                 (x+1) *  rd_array[process[i].region_GID]["structure_info"]["xy_length_per_tile"],
                                 (y+1) *  rd_array[process[i].region_GID]["structure_info"]["xy_length_per_tile"],
                                 rd_array[process[i].region_GID]["structure_info"]["region_z_extent"][1]]
    
    process[i].center_position=[ (x+0.5)*  rd_array[process[i].region_GID]["structure_info"]["xy_length_per_tile"],
                                     (y+0.5) *  rd_array[process[i].region_GID]["structure_info"]["xy_length_per_tile"],
                                     (rd_array[process[i].region_GID]["structure_info"]["region_z_extent"][0] + rd_array[process[i].region_GID]["structure_info"]["region_z_extent"][1])/2 ]
    
    # cumulative depth of layer
    lt_sum=0
    initial_depth, final_depth = [], []

    for lt in rd_array[process[i].region_GID]["structure_info"]["layer_thickness"]:
        initial_depth.append(lt_sum)
        lt_sum += lt
        final_depth.append(lt_sum)

        process[i].subregion_positions = [ [process[i].spatial_extent[0],
                                            process[i].spatial_extent[1],
                                            process[i].spatial_extent[2] + id,
                                            process[i].spatial_extent[3],
                                            process[i].spatial_extent[4],
                                            process[i].spatial_extent[2] + fd ] for id, fd in zip(initial_depth, final_depth) ]
        
    # test2
    i=1
    process[i].region_GID = 1
    process[i].region_name = "test2"
    
    process[i].region_GID = sd["region_name_to_region_GID"][process[i].region_name]
    x = 10
    y = 0
    process[i].spatial_extent =[ x *  rd_array[process[i].region_GID]["structure_info"]["xy_length_per_tile"] ,
                                 y *  rd_array[process[i].region_GID]["structure_info"]["xy_length_per_tile"],
                                 rd_array[process[i].region_GID]["structure_info"]["region_z_extent"][0],
                                 (x+1) *  rd_array[process[i].region_GID]["structure_info"]["xy_length_per_tile"],
                                 (y+1) *  rd_array[process[i].region_GID]["structure_info"]["xy_length_per_tile"],
                                 rd_array[process[i].region_GID]["structure_info"]["region_z_extent"][1]]
    
    process[i].center_position=[ (x+0.5)*  rd_array[process[i].region_GID]["structure_info"]["xy_length_per_tile"],
                                     (y+0.5) *  rd_array[process[i].region_GID]["structure_info"]["xy_length_per_tile"],
                                     (rd_array[process[i].region_GID]["structure_info"]["region_z_extent"][0] + rd_array[process[i].region_GID]["structure_info"]["region_z_extent"][1])/2 ]
    
    # cumulative depth of layer
    lt_sum=0
    initial_depth, final_depth = [], []
    for lt in rd_array[process[i].region_GID]["structure_info"]["layer_thickness"]:
        initial_depth.append(lt_sum)
        lt_sum += lt
        final_depth.append(lt_sum)

        process[i].subregion_positions = [ [process[i].spatial_extent[0],
                                               process[i].spatial_extent[1],
                                               process[i].spatial_extent[2] + id,
                                               process[i].spatial_extent[3],
                                               process[i].spatial_extent[4],
                                               process[i].spatial_extent[2] + fd ] for id, fd in zip(initial_depth, final_depth) ]
    return process
 
def set_process_info_cortico_thalamo_cerebello_sheets(sd, rd_array):
    process = [ [] for x in range(sd["n_process"])]

    # cortex
    count_process_GIDs = 0
    count_region_GIDs = 0

    for i_mrn, mrn in enumerate(sd["meta_region_names"]):
        
        x_region_width = sd[mrn]["x_points"] / sd[mrn]["n_x_regions"]
        y_region_width = sd[mrn]["y_points"] / sd[mrn]["n_y_regions"]
        
        for i_rn, rn in enumerate(sd[mrn]["region_names"]):
            # i_x_region per meta region. Therefore, i_rn is used here.
            i_x_region = i_rn % sd[mrn]["n_x_regions"]
            i_y_region = (i_rn / sd[mrn]["n_y_regions"]) % sd[mrn]["n_y_regions"]
            
            min_x_rank_tile = i_x_region * x_region_width
            max_x_rank_tile = (i_x_region + 1) * x_region_width
            min_y_rank_tile = i_y_region * y_region_width
            max_y_rank_tile = (i_y_region + 1) * y_region_width
            
            rd_array[i_rn + count_region_GIDs]["structure_info"]["n_tiles_on_a_x_side"] = x_region_width 
            rd_array[i_rn + count_region_GIDs]["structure_info"]["n_tiles_on_a_y_side"] = y_region_width 

            # tile group mode
            if max_x_rank_tile - min_x_rank_tile >= 32 and max_y_rank_tile - min_y_rank_tile >= 32:
                rd_array[i_rn + count_region_GIDs]["structure_info"]["n_points_on_a_side_of_tile_group"] = 2 * int(math.ceil( rd_array[i_rn + count_region_GIDs]["structure_info"]["tile_link_limit"]
                                                                                                                              / float(rd_array[i_rn + count_region_GIDs]["structure_info"]["xy_length_per_tile"]) )) + 1

                rd_array[i_rn + count_region_GIDs]["structure_info"]["n_tile_groups"] = int(math.ceil( rd_array[i_rn + count_region_GIDs]["structure_info"]["n_tiles_on_a_x_side"]
                                                                                                       / float(rd_array[i_rn + count_region_GIDs]["structure_info"]["n_points_on_a_side_of_tile_group"]))
                                                                                            * math.ceil( rd_array[i_rn + count_region_GIDs]["structure_info"]["n_tiles_on_a_y_side"]
                                                                                                         / float(rd_array[i_rn + count_region_GIDs]["structure_info"]["n_points_on_a_side_of_tile_group"] )) )

                rd_array[i_rn + count_region_GIDs]["structure_info"]["n_tile_groups_on_a_x_side"] = int(math.ceil( rd_array[i_rn + count_region_GIDs]["structure_info"]["n_tiles_on_a_x_side"]
                                                                                                                   / float(rd_array[i_rn + count_region_GIDs]["structure_info"]["n_points_on_a_side_of_tile_group"] )))
                rd_array[i_rn + count_region_GIDs]["structure_info"]["n_tile_groups_on_a_y_side"] = int(math.ceil( rd_array[i_rn + count_region_GIDs]["structure_info"]["n_tiles_on_a_y_side"]
                                                                                                                   / float(rd_array[i_rn + count_region_GIDs]["structure_info"]["n_points_on_a_side_of_tile_group"] )))
                rd_array[i_rn + count_region_GIDs]["structure_info"]["tile_group"] = [ [] for i in range(rd_array[i_rn + count_region_GIDs]["structure_info"]["n_tile_groups"])]
                rd_array[i_rn + count_region_GIDs]["structure_info"]["tile_group_pair"] = [ [] for i in range(rd_array[i_rn + count_region_GIDs]["structure_info"]["n_tile_groups"])]
                rd_array[i_rn + count_region_GIDs]["structure_info"]["tile_group_center_position"] = [ [] for i in range(rd_array[i_rn + count_region_GIDs]["structure_info"]["n_tile_groups"])]
                rd_array[i_rn + count_region_GIDs]["structure_info"]["xy_plane_4vertices"] = [ [ [] for j in range(4) ] for i in range(rd_array[i_rn + count_region_GIDs]["structure_info"]["n_tile_groups"]) ]

                # positions of tile groups
                for i_tg in range(rd_array[i_rn + count_region_GIDs]["structure_info"]["n_tile_groups"]):
                    # center position of tile group
                    rd_array[i_rn + count_region_GIDs]["structure_info"]["tile_group_center_position"][i_tg] = [ (i_tg % rd_array[i_rn + count_region_GIDs]["structure_info"]["n_tile_groups_on_a_x_side"])
                                                                                                                 * rd_array[i_rn + count_region_GIDs]["structure_info"]["n_points_on_a_side_of_tile_group"]
                                                                                                                 + rd_array[i_rn + count_region_GIDs]["structure_info"]["n_points_on_a_side_of_tile_group"]/2
                                                                                                                 + min_x_rank_tile * rd_array[i_rn + count_region_GIDs]["structure_info"]["xy_length_per_tile"],
                                                                                                                 i_tg / rd_array[i_rn + count_region_GIDs]["structure_info"]["n_tile_groups_on_a_x_side"]
                                                                                                                 * rd_array[i_rn + count_region_GIDs]["structure_info"]["n_points_on_a_side_of_tile_group"]
                                                                                                                 + rd_array[i_rn + count_region_GIDs]["structure_info"]["n_points_on_a_side_of_tile_group"]/2
                                                                                                                 + min_y_rank_tile * rd_array[i_rn + count_region_GIDs]["structure_info"]["xy_length_per_tile"] ]
                
                    # four edge vetices of tile group
                    rd_array[i_rn + count_region_GIDs]["structure_info"]["xy_plane_4vertices"][i_tg][0] = [ (( 0 + (i_tg % rd_array[i_rn + count_region_GIDs]["structure_info"]["n_tile_groups_on_a_x_side"]))
                                                                                                             * rd_array[i_rn + count_region_GIDs]["structure_info"]["n_points_on_a_side_of_tile_group"]
                                                                                                             + min_x_rank_tile) * rd_array[i_rn + count_region_GIDs]["structure_info"]["xy_length_per_tile"],
                                                                                                            (( 0 + (i_tg / rd_array[i_rn + count_region_GIDs]["structure_info"]["n_tile_groups_on_a_x_side"]) )
                                                                                                             * rd_array[i_rn + count_region_GIDs]["structure_info"]["n_points_on_a_side_of_tile_group"]
                                                                                                             + min_y_rank_tile) * rd_array[i_rn + count_region_GIDs]["structure_info"]["xy_length_per_tile"] ]
                    rd_array[i_rn + count_region_GIDs]["structure_info"]["xy_plane_4vertices"][i_tg][1] = [ (( 1 + (i_tg % rd_array[i_rn + count_region_GIDs]["structure_info"]["n_tile_groups_on_a_x_side"]))
                                                                                                             * rd_array[i_rn + count_region_GIDs]["structure_info"]["n_points_on_a_side_of_tile_group"]
                                                                                                             + min_x_rank_tile) * rd_array[i_rn + count_region_GIDs]["structure_info"]["xy_length_per_tile"],
                                                                                                            (( 0 + (i_tg / rd_array[i_rn + count_region_GIDs]["structure_info"]["n_tile_groups_on_a_x_side"]) )
                                                                                                             * rd_array[i_rn + count_region_GIDs]["structure_info"]["n_points_on_a_side_of_tile_group"]
                                                                                                             + min_y_rank_tile) * rd_array[i_rn + count_region_GIDs]["structure_info"]["xy_length_per_tile"] ]
                    rd_array[i_rn + count_region_GIDs]["structure_info"]["xy_plane_4vertices"][i_tg][2] = [ (( 1 + (i_tg % rd_array[i_rn + count_region_GIDs]["structure_info"]["n_tile_groups_on_a_x_side"]))
                                                                                                             * rd_array[i_rn + count_region_GIDs]["structure_info"]["n_points_on_a_side_of_tile_group"]
                                                                                                             + min_x_rank_tile) * rd_array[i_rn + count_region_GIDs]["structure_info"]["xy_length_per_tile"],
                                                                                                            (( 1 + (i_tg / rd_array[i_rn + count_region_GIDs]["structure_info"]["n_tile_groups_on_a_x_side"]) )
                                                                                                             * rd_array[i_rn + count_region_GIDs]["structure_info"]["n_points_on_a_side_of_tile_group"]
                                                                                                             + min_y_rank_tile) * rd_array[i_rn + count_region_GIDs]["structure_info"]["xy_length_per_tile"] ]
                    rd_array[i_rn + count_region_GIDs]["structure_info"]["xy_plane_4vertices"][i_tg][3] = [ (( 0 + (i_tg % rd_array[i_rn + count_region_GIDs]["structure_info"]["n_tile_groups_on_a_x_side"]))
                                                                                                             * rd_array[i_rn + count_region_GIDs]["structure_info"]["n_points_on_a_side_of_tile_group"]
                                                                                                             + min_x_rank_tile) * rd_array[i_rn + count_region_GIDs]["structure_info"]["xy_length_per_tile"],
                                                                                                            (( 1 + (i_tg / rd_array[i_rn + count_region_GIDs]["structure_info"]["n_tile_groups_on_a_x_side"]) )
                                                                                                             * rd_array[i_rn + count_region_GIDs]["structure_info"]["n_points_on_a_side_of_tile_group"]
                                                                                                             + min_y_rank_tile) * rd_array[i_rn + count_region_GIDs]["structure_info"]["xy_length_per_tile"] ]
                
            # no tile group mode
            else:
                rd_array[i_rn + count_region_GIDs]["structure_info"]["n_tile_groups"] = 0
                
            # match x y coorditates of tiles and spatial extent of region, and register corresponding process ID   
            rd_array[i_rn + count_region_GIDs]["process_ID"] = []
            
            for x in range(sd[mrn]["x_points"]):
                for y in range(sd[mrn]["y_points"]):
                    if (min_x_rank_tile <= x < max_x_rank_tile) and (min_y_rank_tile <= y < max_y_rank_tile):
                        process_ID = y * sd[mrn]["x_points"] + x + count_process_GIDs
                        rd_array[i_rn + count_region_GIDs]["process_ID"].append(process_ID)

            rd_array[i_rn + count_region_GIDs]["process_ID"].sort()    
            
        # assign region to process and set positions of process in cortical sheet
        sd[mrn]["region_y_points"] = sd[mrn]["y_points"] / sd[mrn]["n_y_regions"]
        sd[mrn]["region_x_points"] = sd[mrn]["x_points"] / sd[mrn]["n_x_regions"]

        # set process spatial info
        for y in range(sd[mrn]["y_points"]):
            i_y_region = int(y / sd[mrn]["region_y_points"])
                
            for x in range(sd[mrn]["x_points"]):
                                        
                i_x_region = int(x / sd[mrn]["region_x_points"])

                i_proc = y * sd[mrn]["x_points"] + x + count_process_GIDs
                
                process[i_proc] = process_info()
                process[i_proc].process_GID = i_proc
                process[i_proc].region_GID = sd[mrn]["n_x_regions"] * i_y_region + i_x_region + count_region_GIDs
                process[i_proc].region_name = sd["region_names"][process[i_proc].region_GID]
                process[i_proc].inter_regional_process_pairs = [[] for i in range(len(rd_array[process[i_proc].region_GID]["inter_regional_connection"])) ]
                process[i_proc].inter_regional_connection_process = [[] for i in range(len(rd_array[process[i_proc].region_GID]["inter_regional_connection"])) ]
                process[i_proc].inter_regional_pre_post_address = [[] for i in range(len(rd_array[process[i_proc].region_GID]["inter_regional_connection"])) ]

                process[i_proc].spatial_extent =[ x * rd_array[process[i_proc].region_GID]["structure_info"]["xy_length_per_tile"] ,
                                                  y * rd_array[process[i_proc].region_GID]["structure_info"]["xy_length_per_tile"],
                                                  rd_array[process[i_proc].region_GID]["structure_info"]["region_z_extent"][0],
                                                  (x + 1) * rd_array[process[i_proc].region_GID]["structure_info"]["xy_length_per_tile"],
                                                  (y + 1) * rd_array[process[i_proc].region_GID]["structure_info"]["xy_length_per_tile"],
                                                  rd_array[process[i_proc].region_GID]["structure_info"]["region_z_extent"][1]]
                    
                process[i_proc].center_position=[ (x + 0.5) * rd_array[process[i_proc].region_GID]["structure_info"]["xy_length_per_tile"],
                                                  (y + 0.5) * rd_array[process[i_proc].region_GID]["structure_info"]["xy_length_per_tile"],
                                                  (rd_array[process[i_proc].region_GID]["structure_info"]["region_z_extent"][0] + rd_array[process[i_proc].region_GID]["structure_info"]["region_z_extent"][1]) / 2 ]

                process[i_proc].xy_plane_4vertices = [ [ (x + 0) *  rd_array[process[i_proc].region_GID]["structure_info"]["xy_length_per_tile"],
                                                         (y + 0) *  rd_array[process[i_proc].region_GID]["structure_info"]["xy_length_per_tile"],
                                                         (rd_array[process[i_proc].region_GID]["structure_info"]["region_z_extent"][0] + rd_array[process[i_proc].region_GID]["structure_info"]["region_z_extent"][1]) / 2 ],
                                                       [ (x + 1) *  rd_array[process[i_proc].region_GID]["structure_info"]["xy_length_per_tile"],
                                                         (y + 0) *  rd_array[process[i_proc].region_GID]["structure_info"]["xy_length_per_tile"],
                                                         (rd_array[process[i_proc].region_GID]["structure_info"]["region_z_extent"][0] + rd_array[process[i_proc].region_GID]["structure_info"]["region_z_extent"][1]) / 2 ],
                                                       [ (x + 1) *  rd_array[process[i_proc].region_GID]["structure_info"]["xy_length_per_tile"],
                                                         (y + 1) *  rd_array[process[i_proc].region_GID]["structure_info"]["xy_length_per_tile"],
                                                         (rd_array[process[i_proc].region_GID]["structure_info"]["region_z_extent"][0] + rd_array[process[i_proc].region_GID]["structure_info"]["region_z_extent"][1]) / 2 ],
                                                       [ (x + 0) *  rd_array[process[i_proc].region_GID]["structure_info"]["xy_length_per_tile"],
                                                         (y + 1) *  rd_array[process[i_proc].region_GID]["structure_info"]["xy_length_per_tile"],
                                                         (rd_array[process[i_proc].region_GID]["structure_info"]["region_z_extent"][0] + rd_array[process[i_proc].region_GID]["structure_info"]["region_z_extent"][1]) / 2 ] ]
                     
                # calculate TofuD coordinates based on the numbers of CMGs per CPU for each process
                # A64fx has 4 CMGs. 2 x 2 process are mapped on 2 by 2 xy plane.
                TOFUD_2D_coordinate_x = x / 2
                TOFUD_2D_coordinate_y = y / 2
                process[i_proc].TOFUD_2D_coordinate = [TOFUD_2D_coordinate_x, TOFUD_2D_coordinate_y]
            
                # cumulative depth of layer
                lt_sum=0
                initial_depth, final_depth = [], []
                for lt in rd_array[process[i_proc].region_GID]["structure_info"]["layer_thickness"]:
                    initial_depth.append(lt_sum)
                    lt_sum += lt
                    final_depth.append(lt_sum)

                process[i_proc].subregion_positions = [ [process[i_proc].spatial_extent[0],
                                                         process[i_proc].spatial_extent[1],
                                                         process[i_proc].spatial_extent[2] + id,
                                                         process[i_proc].spatial_extent[3],
                                                         process[i_proc].spatial_extent[4],
                                                         process[i_proc].spatial_extent[2] + fd ] for id, fd in zip(initial_depth, final_depth) ]

                # when tile is larger than 32x32 tile, tile group is prepared
                if rd_array[process[i_proc].region_GID]["structure_info"]["n_tile_groups"] > 0:
                    process[i_proc].i_x_tile_group = int(math.floor( x / float(rd_array[i_rn]["structure_info"]["n_points_on_a_side_of_tile_group"])))
                    process[i_proc].i_y_tile_group = int(math.floor( y / float(rd_array[i_rn]["structure_info"]["n_points_on_a_side_of_tile_group"])))
                    process[i_proc].i_tile_group = process[i_proc].i_y_tile_group * rd_array[i_rn]["structure_info"]["n_tile_groups_on_a_x_side"] + process[i_proc].i_x_tile_group
                    rd_array[process[i_proc].region_GID]["structure_info"]["tile_group"][process[i_proc].i_tile_group].append( i_proc )

                # cluster
                # hexagonal_cluster(process[i_proc], rd_array, x, y)
                            
        for i_y_region in range(sd[mrn]["n_y_regions"]):
            for i_x_region in range(sd[mrn]["n_x_regions"]):
                i_region = sd[mrn]["n_x_regions"] * i_y_region + i_x_region
                start_x_region =  rd_array[i_region]["structure_info"]["xy_length_per_tile"] * sd[mrn]["region_x_points"] * i_x_region
                start_y_region =  rd_array[i_region]["structure_info"]["xy_length_per_tile"] * sd[mrn]["region_y_points"] * i_y_region
                end_x_region =  rd_array[i_region]["structure_info"]["xy_length_per_tile"] * (i_x_region + 1)
                end_y_region =  rd_array[i_region]["structure_info"]["xy_length_per_tile"] * (i_y_region + 1)
                rd_array[i_region]["structure_info"]["region_center"]=[start_x_region +  rd_array[i_region]["structure_info"]["xy_length_per_tile"] * sd[mrn]["region_x_points"] / 2,
                                                                       start_y_region +  rd_array[i_region]["structure_info"]["xy_length_per_tile"] * sd[mrn]["region_y_points"] / 2,
                                                                       (rd_array[i_region]["structure_info"]["region_z_extent"][0] + rd_array[i_region]["structure_info"]["region_z_extent"][1]) / 2]


        """
        # cluster
        # horizontal grid
        if rd["structure_info"]["cluster_shape"] == "hexagonal_grid":
            upper_plane_vertex = []
            lower_plane_vertex = []
                        
        # regular grid
        elif rd["structure_info"]["cluster_shape"] == "regular_grid":
            
            for i_x in range(rd["structure_info"]["n_clusters_on_a_side"] * sd[mrn]["x_points"] ):
                for i_y in range(rd["structure_info"]["n_clusters_on_a_side"] * sd[mrn]["y_points"]):
                    center_x = rd["structure_info"]["cluster_size"] * 0.5 + i_x * rd["structure_info"]["cluster_size"]
                    center_y = rd["structure_info"]["cluster_size"] * 0.5 + i_y * rd["structure_info"]["cluster_size"]
                    center_z = (rd["structure_info"]["region_z_extent"][0] + rd["structure_info"]["region_z_extent"][1]) / 2. 
                        
        # exceptional process
        else:
            print "There is no registered name in cluster_shape:", rd["structure_info"]["cluster_shape"]
            exit()
        """
                    
        count_process_GIDs += sd[mrn]["n_process"]
        count_region_GIDs += sd[mrn]["n_regions"]
        
    return process

"""
def hexagonal_cluster(process, rd_array, x, y):
    # spatial extent of process
    x_min = x * rd_array[process.region_GID]["structure_info"]["xy_length_per_tile"]
    x_max = (x+1) * rd_array[process.region_GID]["structure_info"]["xy_length_per_tile"]
    y_min = y * rd_array[process.region_GID]["structure_info"]["xy_length_per_tile"]
    y_max = (y+1) * rd_array[process.region_GID]["structure_info"]["xy_length_per_tile"]

    #y_min = int(math.ceil(y_min /  (1.5 * rd["structure_info"]["cluster_size"])) * (1.5 * rd["structure_info"]["cluster_size"]))
    #y_max = int(math.floor(y_max /  (1.5 * rd["structure_info"]["cluster_size"])) * (1.5 * rd["structure_info"]["cluster_size"]))
   
    # region directory
    rd = rd_array[process.region_GID]

    # center of hexaponal grid  
    for i_y_center, y_center in enumerate(range(y_min, y_max, int(1.5 * rd["structure_info"]["cluster_size"]))):
        # determine x extent
        x_min_center = int( math.ceil(x_min /  (np.sqrt(3) * rd["structure_info"]["cluster_size"])) * np.sqrt(3) * rd["structure_info"]["cluster_size"] + (i_y_center%2) * np.sqrt(3) * rd["structure_info"]["cluster_size"] * 0.5)
        x_max_center =  int( math.ceil(x_max /  (np.sqrt(3) * rd["structure_info"]["cluster_size"])) * np.sqrt(3) * rd["structure_info"]["cluster_size"] + (i_y_center%2) * np.sqrt(3) * rd["structure_info"]["cluster_size"] * 0.5)

        for x_center in range(x_min_center, x_max_center, int(np.sqrt(3) * rd["structure_info"]["cluster_size"])):
            center_x = x_center
            center_y = y_center
            center_z = (rd["structure_info"]["region_z_extent"][0] + rd["structure_info"]["region_z_extent"][1]) / 2.

            shift_x = np.sqrt(3) * rd["structure_info"]["cluster_size"] * 0.5
            shift_y = rd["structure_info"]["cluster_size"] * 0.5
            
            upper_plane_vertex = []
            lower_plane_vertex = []
            
            upper_plane_vertex.append([
                center_x - shift_x,
                center_y - shift_y,
                rd["structure_info"]["region_z_extent"][1]
            ])
            upper_plane_vertex.append([
                center_x - shift_x,
                center_y + shift_y,
                rd["structure_info"]["region_z_extent"][1]
            ])
            upper_plane_vertex.append([
                center_x,
                center_y + 2 * shift_y,
                rd["structure_info"]["region_z_extent"][1]
            ])
            upper_plane_vertex.append([
                center_x + shift_x,
                center_y + shift_y,
            rd["structure_info"]["region_z_extent"][1]
            ])
            upper_plane_vertex.append([
                center_x + shift_x,
                center_y - shift_y,
                rd["structure_info"]["region_z_extent"][1]
            ])
            upper_plane_vertex.append([
                center_x,
                center_y - 2 * shift_y,
            rd["structure_info"]["region_z_extent"][1]
            ])
            
            lower_plane_vertex.append([
            center_x - shift_x,
                center_y - shift_y,
                rd["structure_info"]["region_z_extent"][0]
            ])
            lower_plane_vertex.append([
                center_x - shift_x,
                center_y + shift_y,
                rd["structure_info"]["region_z_extent"][0]
            ])
            lower_plane_vertex.append([
                center_x,
                center_y + 2 * shift_y,
                rd["structure_info"]["region_z_extent"][0]
            ])
            lower_plane_vertex.append([
                center_x + shift_x,
                center_y + shift_y,
                rd["structure_info"]["region_z_extent"][0]
            ])
            lower_plane_vertex.append([
                center_x + shift_x,
                center_y - shift_y,
                rd["structure_info"]["region_z_extent"][0]
            ])
            lower_plane_vertex.append([
                center_x,
                center_y - 2 * shift_y,
                rd["structure_info"]["region_z_extent"][0]
            ])
            
            process.cluster_center.append([center_x, center_y, center_z])
            process.cluster_upper_plane.append([upper_plane_vertex])
            process.cluster_lower_plane.append([lower_plane_vertex])
"""

                            
def main():

    #
    """
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('integers', metavar='N', type=int, nargs='+',
                        help='an integer for the accumulator')
    parser.add_argument('--sum', dest='accumulate', action='store_const',
                        const=sum, default=max,
                        help='sum the integers (default: find the max)')
    
    args = parser.parse_args()
    print(args.accumulate(args.integers))
    """
    
    # read system parameters
    start = time.time()
    sd = mone_in.read_system_parameters()
    elapsed_time = time.time() - start
    print("read_system_parameters:{0}".format(elapsed_time) + "[sec]")
    
    # read region parameters
    start = time.time()
    rd_array = mone_in.read_region_parameters(sd)
    elapsed_time = time.time() - start
    print("read_region_parameters:{0}".format(elapsed_time) + "[sec]")

    # set estimated numbers of connections 
    start = time.time()
    mone_in.estimate_n_connections(rd_array, sd)
    elapsed_time = time.time() - start
    print("estimate_n_connections:{0}".format(elapsed_time) + "[sec]")

    # set process information
    start = time.time()
    process = set_process_info(sd, rd_array)
    elapsed_time = time.time() - start
    print("set_process_info: {0}".format(elapsed_time) + "[sec]")
    
    # set connection information
    start = time.time()
    mone_conn.set_region_connection(sd, rd_array, process)
    elapsed_time = time.time() - start
    print("set_region_connection: {0}".format(elapsed_time) + "[sec]")

    """
    # shift degree for commissural connection caculation
    count_process_GIDs = 0
    count_region_GIDs = 0
    for i_mrn, mrn in enumerate(sd["meta_region_names"]):
        x_region_width = sd[mrn]["x_points"] / sd[mrn]["n_x_regions"]
        y_region_width = sd[mrn]["y_points"] / sd[mrn]["n_y_regions"]
        #print(x_region_width, y_region_width)
        
        if mrn == "cortex":
            print(mrn)
            for i_rn, rn in enumerate(sd[mrn]["region_names"]):
                #print(-x_region_width * rd_array[i_rn + count_region_GIDs]["structure_info"]["xy_length_per_tile"])
                #print(x_region_width * rd_array[i_rn + count_region_GIDs]["structure_info"]["xy_length_per_tile"])
                #print(rn)
                for bn in rd_array[i_rn + count_region_GIDs]["bdl"]:
                    #print(bn)
                    if bn == "bundle_M1_to_M2" and rn == "M2":
                        print("before", bn, rn, rd_array[i_rn + count_region_GIDs]["bdl"][bn].shift[0])
                        rd_array[i_rn + count_region_GIDs]["bdl"][bn].shift[0] = - x_region_width * rd_array[i_rn + count_region_GIDs]["structure_info"]["xy_length_per_tile"]
                        print("after", bn, rn, rd_array[i_rn + count_region_GIDs]["bdl"][bn].shift[0])

                    if bn == "bundle_M2_to_M1" and rn == "M2":
                        print("before", bn, rn, rd_array[i_rn + count_region_GIDs]["bdl"][bn].shift[0])
                        rd_array[i_rn + count_region_GIDs]["bdl"][bn].shift[0] =  x_region_width * rd_array[i_rn + count_region_GIDs]["structure_info"]["xy_length_per_tile"]
                        print("after", bn, rn, rd_array[i_rn + count_region_GIDs]["bdl"][bn].shift[0])

                    elif bn == "bundle_M2_to_M1" and rn == "M1":
                        print("before", bn, rn, rd_array[i_rn + count_region_GIDs]["bdl"][bn].shift[0])
                        rd_array[i_rn + count_region_GIDs]["bdl"][bn].shift[0] = x_region_width * rd_array[i_rn + count_region_GIDs]["structure_info"]["xy_length_per_tile"]
                        print("after", bn, rn, rd_array[i_rn + count_region_GIDs]["bdl"][bn].shift[0])
                        
                    elif bn == "bundle_M1_to_M2" and rn == "M1":
                        print("before", bn, rn, rd_array[i_rn + count_region_GIDs]["bdl"][bn].shift[0])
                        rd_array[i_rn + count_region_GIDs]["bdl"][bn].shift[0] = - x_region_width * rd_array[i_rn + count_region_GIDs]["structure_info"]["xy_length_per_tile"]
                        print("after", bn, rn, rd_array[i_rn + count_region_GIDs]["bdl"][bn].shift[0])
                        

        count_region_GIDs += sd[mrn]["n_regions"]
    """
    
    # DTI
    if sd["use_DTI_data"] == True:
        start = time.time()
        mone_DTI.read_DTI_data(sd, rd_array, process)    
        elapsed_time = time.time() - start
        print("DTI_connection: {0}".format(elapsed_time) + "[sec]")

    # wirte process configuration files
    start = time.time()
    
    """
    # wirte process configuration files
    start = time.time()

    n_mp = 8
    procs_per_mp = [[]for i in range(n_mp)]
    n_proc = sd["n_process"]
    n_proc_mp = int(n_proc / n_mp)
    
    for i in range(n_mp):
        start_proc, end_proc = i * n_proc_mp, (i+1) * n_proc_mp
        procs_per_mp[i] = list(range(start_proc, end_proc))

    modulo = n_proc % n_mp
    for i in range(modulo):
        procs_per_mp[i].append(i + n_proc_mp * n_mp)

    jobs=[]
    for i in range(n_mp):
        #start_proc, end_proc = i*n_proc_mp, (i+1)*n_proc_mp
        #p = Process(target=mone_out.write_process_configuration_files, args=(sd, rd_array, process, range(start_proc, end_proc)))
        p = Process(target=mone_out.write_process_configuration_files, args=(sd, rd_array, process, procs_per_mp[i]))
        jobs.append(p)
        p.start()
    for job in jobs:
        job.join()
    elapsed_time = time.time() - start
    """
    
    # write rank node map file
    mone_out.write_rank_node_map(sd, process)
    mone_out.write_global_shared_info_file(sd, rd_array)
    mone_out.write_region_definition_file(sd, rd_array)
    mone_out.write_process_info_file(sd, process)
    #mone_out.write_job_script_fugaku_detail_profile(sd, rd_array)
    mone_out.write_job_script_fugaku(sd, rd_array)
    elapsed_time = time.time() - start
    print("write_process_configuration_files: {0}".format(elapsed_time) + "[sec]")

    # draw figures of network
    if switch_vtk_availablity == True and sd["network_configuration_visualization"] == True:
        mone_vis.visualize_network(process, sd)
    elif switch_vtk_availablity == False and sd["network_configuration_visualization"] == True:
        print('vtk is not installed. Please install to use visualization')

if __name__ == "__main__":
    main()
