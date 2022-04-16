#!/bin/env python
import math
import copy
import time
import multiprocessing as mp
from multiprocessing import Process
from multiprocessing import Pool
import sys
#import matplotlib.pyplot as plt
import itertools
#import numpy as np

class connection_type:
    """ connection_type info class"""
    def __init__(self):
        self.connection_type_SerLID = []
        self.connection_type_full_name=""
        self.connection_type_LID_in_nt=[]
        self.pre_rg_GID = []
        self.post_rg_GID = []
        self.pre_nt_SerLID = []
        self.post_nt_SerLID = []
        self.pre_nt_LID = []
        self.pre_sr_LID = []
        self.post_nt_LID = []
        self.post_sr_LID = []
        self.LID_in_nt = []
        self.pre_rg_name = []
        self.post_rg_name = []
        self.pre_sr_name = []
        self.post_sr_name = []
        self.pre_nt_name = []
        self.post_nt_name = []
        
        self.inter_or_intra = []
        self.bundle_name = []
        self.bundle_LID = []
        self.connection_name = []
        self.connection_type_LID_in_bundle = []
        
        self.topology = []
        self.spatial_pattern = []
        self.synaptic_channel = []
        self.weight_distribution = []
        self.delay = []
        self.LTP = []

class bundle_class:
    """ bundle class"""
    def __init__(self):
        self.n_indegree_inter_regional_connection_types_per_bundle=0
        self.bundle_name=""
        self.bundle_LID = []
        self.inter_regional_pre_process_per_bundle = []
        self.inter_regional_post_process_per_bundle = []
        self.inter_regional_pre_post_address_per_bundle = []
        self.inter_regional_post_neuron_types_per_pre_neuron_type_per_bundle = []
        self.new_inter_regional_post_neuron_types_per_pre_neuron_type_at_bundle = []
        self.n_pre_nt_at_bundle = 0
        
        
class connection_in_bundle:
    """connection in bundle class """
    def __int__(self):
        self.connection_type_SerLID = []
        self.pre_rg_GID = []
        self.post_rg_GID = []
        self.pre_nt_SerLID = []
        self.post_nt_SerLID = []
        self.pre_nt_LID = []
        self.pre_sr_LID = []
        self.post_nt_LID = []
        self.post_sr_LID = []
        self.LID_in_nt = []
        self.pre_rg_name = []
        self.post_rg_name = []
        self.pre_sr_name = []
        self.post_sr_name = []
        self.pre_nt_name = []
        self.post_nt_name = []
        self.connection_type_LID_in_bundle = []
        self.pre_post_colocalization = []
        self.n_post_nt = 0

def set_region_connection(sd, rd_array, process):

    start = time.time()
    #####################################################
    # set address based on process IDs and neuronal IDs #
    #####################################################
    # 1. first, prepare intra connection pair data
    for i_rd, rd in enumerate(rd_array):
        temp = []
        for i_post_nt, cpm in enumerate(rd["intra_regional_connection"]["connection_parameter_matrix"]):
            for i_pre_nt, cpm_row in enumerate(cpm):
                if cpm_row["spatial_pattern"] != "None":
                        
                    temp.append([
                        rd["structure_info"]["neuron_type_SerLID_to_subregion_LID"][i_pre_nt],
                        rd["structure_info"]["neuron_type_SerLID_to_neuron_type_LID"][i_pre_nt],
                        rd["structure_info"]["neuron_type_SerLID_to_subregion_LID"][i_post_nt],
                        rd["structure_info"]["neuron_type_SerLID_to_neuron_type_LID"][i_post_nt]
                    ])
                    
        rd_array[i_rd]["intra_regional_connection"]["pre_post_neuron_type_LID_pairs"] = temp
    
    elapsed_time = time.time() - start
    print ("set_region_connection1.0: {0}".format(elapsed_time) + "[sec]")
    start = time.time()
    """
    # 1.1 tile group pair
    for i_rd, rd in enumerate(rd_array):
        if rd["structure_info"]["n_tile_groups"] > 0:
            for i_ctg, center_tg in enumerate(rd["structure_info"]["xy_plane_4vertices"]):
                for i_stg, surround_tg in enumerate(rd["structure_info"]["xy_plane_4vertices"]):
                    distance = minimum_distance_between_4vertices(surround_tg, center_tg)

                    if distance < rd["structure_info"]["tile_link_limit"]:
                        rd["structure_info"]["tile_group_pair"][i_ctg].append(i_stg)
            
            #print "tg", rd["structure_info"]["tile_group"]
            #print "tg pos", rd["structure_info"]["tile_group_center_position"]
            #print "tg pair", rd["structure_info"]["tile_group_pair"]
            #print len(rd["structure_info"]["tile_group_pair"])
            
    # 2. second, search intra regional process pairs
    for i_rd, rd in enumerate(rd_array):
        # no tile group mode
        if rd["structure_info"]["n_tile_groups"] == 0:
            for i_post_proc in rd["process_ID"]:
                for i_pre_proc in rd["process_ID"]:
                    
                    # calculate distance between process
                    distance = minimum_distance_between_4vertices(process[i_pre_proc].xy_plane_4vertices, process[i_post_proc].xy_plane_4vertices)
                    
                    if (distance < rd_array[process[i_post_proc].region_GID]["structure_info"]["tile_link_limit"]):
                        process[i_post_proc].intra_regional_process_pairs.append([process[i_pre_proc].process_GID, process[i_post_proc].process_GID])
                        process[i_post_proc].intra_regional_connection_process.append(process[i_pre_proc].process_GID)

                process[i_post_proc].intra_regional_process_pairs.sort(key=lambda x:x[0])

        # tile group mode
        elif rd["structure_info"]["n_tile_groups"] > 0:
            for i_post_proc in rd["process_ID"]:
                for i_pre_tg in rd["structure_info"]["tile_group_pair"][process[i_post_proc].i_tile_group]:
                    for i_pre_proc in rd["structure_info"]["tile_group"][i_pre_tg]:
                        distance = minimum_distance_between_4vertices(process[i_pre_proc].xy_plane_4vertices, process[i_post_proc].xy_plane_4vertices)

                        if (distance < rd["structure_info"]["tile_link_limit"]):
                            process[i_post_proc].intra_regional_process_pairs.append([process[i_pre_proc].process_GID, process[i_post_proc].process_GID])
                            process[i_post_proc].intra_regional_connection_process.append(process[i_pre_proc].process_GID)

                process[i_post_proc].intra_regional_process_pairs.sort(key=lambda x:x[0])
    """

    # preparation for conn pre_nt post_nt pairs
    #print(len(rd["intra_regional_connection"]["connection_parameter_matrix"]))
    #print(rd["intra_regional_connection"]["connection_parameter_matrix"])
    n_pre = len(rd["intra_regional_connection"]["connection_parameter_matrix"])
    conn_post_nt_SIDs_per_pre_nt_SID = [ [] for i in range(n_pre) ]
    conn_post_nts_per_pre_nt = [ [ [] for i in range(n_pre) ] for j in range(n_pre) ]
    n_conn_pairs = 0
    n_no_conn_pairs = 0
    #print(n_pre)
    #print(conn_post_nts_per_pre_nt )
    
    
    elapsed_time = time.time() - start
    print ("set_region_connection1.1: {0}".format(elapsed_time) + "[sec]")
    start = time.time()
    """
    # 3. finally, prepare address data from pre to post with inforamtion of process GID and neuronal LIDs
    for i_post_proc in range(sd["n_process"]):
        region_GID = process[i_post_proc].region_GID

        for irpp in process[i_post_proc].intra_regional_process_pairs:
            for nlp in rd_array[region_GID]["intra_regional_connection"]["pre_post_neuron_type_LID_pairs"]:
                    process[i_post_proc].intra_regional_pre_post_address.append([irpp[0], region_GID, nlp[0], nlp[1], irpp[1], region_GID, nlp[2], nlp[3] ])
    """
    elapsed_time = time.time() - start
    print ("set_region_connection1.2: {0}".format(elapsed_time) + "[sec]")


    ### inter
    start = time.time()

    # 1. first, prepare inter connection pair data
    for i_rd, rd in enumerate(rd_array):
        temp = []
        for i_bundle, bundle in enumerate(rd["inter_regional_connection"]):
            pre_region_GID = sd["region_name_to_region_GID"][bundle["pre"]["region"]]
            post_region_GID = sd["region_name_to_region_GID"][bundle["post"]["region"]]
           
            for i_ct, ct in enumerate(bundle["connection_type"]):
                pre_subregion_LID = rd_array[pre_region_GID]["structure_info"]["subregion_name_to_subregion_LID"][ct["pre"]["subregion"]],
                pre_neuron_type_LID = rd_array[pre_region_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][ct["pre"]["subregion"]][ct["pre"]["neuron_type"]],
                post_subregion_LID = rd_array[post_region_GID]["structure_info"]["subregion_name_to_subregion_LID"][ct["post"]["subregion"]],
                post_neuron_type_LID = rd_array[post_region_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][ct["post"]["subregion"]][ct["post"]["neuron_type"]],
                if ct["spatial_pattern"] != "None":
                        
                    temp.append([
                        pre_subregion_LID,
                        pre_neuron_type_LID,
                        post_subregion_LID,
                        post_neuron_type_LID
                    ])

            rd_array[i_rd]["inter_regional_connection"][i_bundle]["pre_post_neuron_type_LID_pairs"] = temp
        #print rd_array[i_rd]["inter_regional_connection"]
    """
    # 2. second, search inter regional process pairs
    for i_rd, rd in enumerate(rd_array):
        for i_bundle, bundle in enumerate(rd["inter_regional_connection"]):

            pre_region_GID = sd["region_name_to_region_GID"][bundle["pre"]["region"]]
            post_region_GID = sd["region_name_to_region_GID"][bundle["post"]["region"]]
            #print ("process_ID", rd_array[post_region_GID]["process_ID"] )
            
            for i_post_proc in rd_array[post_region_GID]["process_ID"]:
                for i_pre_proc in rd_array[pre_region_GID]["process_ID"]:
                    # calculate distance between process
                    distance = minimum_distance_between_4vertices(process[i_pre_proc].xy_plane_4vertices, process[i_post_proc].xy_plane_4vertices)
                    #print ("dist", distance)
                    if (distance < rd_array[process[i_post_proc].region_GID]["structure_info"]["tile_link_limit"]):
                        process[i_post_proc].inter_regional_process_pairs[i_bundle].append([process[i_pre_proc].process_GID, process[i_post_proc].process_GID])
                        process[i_post_proc].inter_regional_connection_process[i_bundle].append(process[i_pre_proc].process_GID)

                # might need sepreated to bundle
                #print (i_rd, i_post_proc, process[i_post_proc].inter_regional_process_pairs)
                #process[i_post_proc].inter_regional_process_pairs.sort(key=lambda x:x[0])
        
    # 3. finally, prepare address data of from pre to post with inforamtion of process GID and neuronal LIDs
    for i_post_proc in range(sd["n_process"]):
        post_region_GID = process[i_post_proc].region_GID
        for i_bundle, irpps in enumerate(process[i_post_proc].inter_regional_process_pairs):
            for irpp in irpps:
                for ppntlp in rd_array[post_region_GID]["inter_regional_connection"][i_bundle]["pre_post_neuron_type_LID_pairs"]:
                    process[i_post_proc].inter_regional_pre_post_address[i_bundle].append([irpp[0], process[irpp[0]].region_GID, ppntlp[0], ppntlp[1], irpp[1], post_region_GID, ppntlp[2], ppntlp[3] ])
    """
    
    elapsed_time = time.time() - start
    print ("set inter regional process connect: {0}".format(elapsed_time) + "[sec]")

    start = time.time()
    #########################################################
    # create connection type and bundle type data structure #
    #########################################################
    # loop post regions
    for i_rd, rd in enumerate(rd_array):
        ct = []
        count_ct=0
        connection_type_SerLIDs_per_post = [ 0 for j in range(rd["structure_info"]["n_neuron_types_per_process"] ) ]
        n_connection_types_per_post = [ 0 for j in range(rd["structure_info"]["n_neuron_types_per_process"] ) ]
        n_indegree_intra_regional_connection_types_per_process = 0
        n_indegree_inter_regional_connection_types_per_process = 0
        intra_regional_connection_nt_pairs = []
        inter_regional_connection_nt_pairs = []
        bdl = {}

        # count numbers of connection types per neurons in advance
        # intra connections
        for post_nt_SerLID, row in enumerate(rd["intra_regional_connection"]["connection_parameter_matrix"]):
            for pre_nt_SerLID, cpe in enumerate(row):
                for i_sc, sc in enumerate(cpe["synaptic_channel"]):
                    n_connection_types_per_post[post_nt_SerLID] += 1
        """
        # intra ID connections
        if "ID_connect" in rd["intra_regional_connection"].keys():
            for i_ID_conn, ID_conn in enumerate(rd["intra_regional_connection"]["ID_connect"]):
                for i_sc, sc in enumerate(ID_conn["synaptic_channel"]):
                    n_connection_types_per_post[post_nt_SerLID] += 1
        """
        
        # inter connections
        for i_b, bundle in enumerate(rd["inter_regional_connection"]):
            for cpe in bundle["connection_type"]:
                for i_sc, sc in enumerate(cpe["synaptic_channel"]):
                    post_sr_LID = rd_array[i_rd]["structure_info"]["subregion_name_to_subregion_LID"][cpe["post"]["subregion"]] 
                    post_nt_LID = rd_array[i_rd]["structure_info"]["neuron_type_name_to_neuron_type_LID"][cpe["post"]["subregion"]][cpe["post"]["neuron_type"]]
                    post_nt_SerLID = rd_array[i_rd]["structure_info"]["LID_to_neuron_type_SerLID"][post_sr_LID][post_nt_LID]
                    n_connection_types_per_post[post_nt_SerLID] += 1
                
        # padding dummy connection type
        n_connection_types_per_post_with_padding = copy.deepcopy(n_connection_types_per_post)
        
        max_n_connection_types_per_post = max(n_connection_types_per_post)
        for i, n_ct in enumerate(n_connection_types_per_post):
            if n_ct != max_n_connection_types_per_post:
                n_connection_types_per_post_with_padding[i] = max_n_connection_types_per_post
        
        #print n_connection_types_per_post
        #print n_connection_types_per_post_with_padding
        
        rd_array[i_rd]["n_connection_types_per_post_with_padding"] = n_connection_types_per_post_with_padding
        rd_array[i_rd]["intra_regional_post_neuron_types_per_pre_neuron_type"] = [ [] for i in range(rd["structure_info"]["n_neuron_types_per_process"])]
        
        # connection type of intra regianl connection (loop in columns and rows in conection matrix)
        for post_nt_SerLID, row in enumerate(rd["intra_regional_connection"]["connection_parameter_matrix"]):
            for pre_nt_SerLID, cpe in enumerate(row):
                
                if len(cpe["synaptic_channel"]) == 1:
                    rd["intra_regional_post_neuron_types_per_pre_neuron_type"][pre_nt_SerLID].append([
                        0,
                        count_ct,
                        -1,
                        sd["region_name_to_region_GID"][cpe["pre"]["region"]],
                        rd["structure_info"]["subregion_name_to_subregion_LID"][cpe["pre"]["subregion"]],
                        rd["structure_info"]["neuron_type_name_to_neuron_type_LID"][cpe["pre"]["subregion"]][cpe["pre"]["neuron_type"]],
                        pre_nt_SerLID,
                        sd["region_name_to_region_GID"][cpe["post"]["region"]],
                        rd["structure_info"]["subregion_name_to_subregion_LID"][cpe["post"]["subregion"]],
                        rd["structure_info"]["neuron_type_name_to_neuron_type_LID"][cpe["post"]["subregion"]][cpe["post"]["neuron_type"]],
                        post_nt_SerLID
                    ])
                    
                elif len(cpe["synaptic_channel"]) == 2:
                    rd["intra_regional_post_neuron_types_per_pre_neuron_type"][pre_nt_SerLID].append([
                        1,
                        count_ct,
                        count_ct+1,
                        sd["region_name_to_region_GID"][cpe["pre"]["region"]],
                        rd["structure_info"]["subregion_name_to_subregion_LID"][cpe["pre"]["subregion"]],
                        rd["structure_info"]["neuron_type_name_to_neuron_type_LID"][cpe["pre"]["subregion"]][cpe["pre"]["neuron_type"]],
                        pre_nt_SerLID,
                        sd["region_name_to_region_GID"][cpe["post"]["region"]],
                        rd["structure_info"]["subregion_name_to_subregion_LID"][cpe["post"]["subregion"]],
                        rd["structure_info"]["neuron_type_name_to_neuron_type_LID"][cpe["post"]["subregion"]][cpe["post"]["neuron_type"]],
                        post_nt_SerLID
                    ])
                else:
                    print ("number of channel was wrong", len(cpe["synaptic_channel"]))
                
 
                # pre, post colocalization, and connection_type_SerLID
                if len(cpe["synaptic_channel"]) == 1:
                    rd["intra_regional_connection"]["connection_parameter_matrix"][post_nt_SerLID][pre_nt_SerLID]["pre_post_colocalization"] = [
                        sd["region_name_to_region_GID"][cpe["pre"]["region"]],
                        rd["structure_info"]["subregion_name_to_subregion_LID"][cpe["pre"]["subregion"]],
                        rd["structure_info"]["neuron_type_name_to_neuron_type_LID"][cpe["pre"]["subregion"]][cpe["pre"]["neuron_type"]],
                        sd["region_name_to_region_GID"][cpe["post"]["region"]],
                        rd["structure_info"]["subregion_name_to_subregion_LID"][cpe["post"]["subregion"]],
                        rd["structure_info"]["neuron_type_name_to_neuron_type_LID"][cpe["post"]["subregion"]][cpe["post"]["neuron_type"]],
                        0,
                        count_ct,
                        -1,
                    ]
                elif len(cpe["synaptic_channel"]) == 2:
                    rd["intra_regional_connection"]["connection_parameter_matrix"][post_nt_SerLID][pre_nt_SerLID]["pre_post_colocalization"] = [
                        sd["region_name_to_region_GID"][cpe["pre"]["region"]],
                        rd["structure_info"]["subregion_name_to_subregion_LID"][cpe["pre"]["subregion"]],
                        rd["structure_info"]["neuron_type_name_to_neuron_type_LID"][cpe["pre"]["subregion"]][cpe["pre"]["neuron_type"]],
                        sd["region_name_to_region_GID"][cpe["post"]["region"]],
                        rd["structure_info"]["subregion_name_to_subregion_LID"][cpe["post"]["subregion"]],
                        rd["structure_info"]["neuron_type_name_to_neuron_type_LID"][cpe["post"]["subregion"]][cpe["post"]["neuron_type"]],
                        1,
                        count_ct,
                        count_ct+1
                    ]
                else:
                    print ("number of channel was wrong", len(cpe["synaptic_channel"]))
                
                
                for i_sc, sc in enumerate(cpe["synaptic_channel"]):

                    ct.append(connection_type())
                    ct[-1].connection_type_SerLID = count_ct
                    ct[-1].connection_type_full_name = str( cpe["pre"]["region"] + "_" + cpe["pre"]["subregion"] + "_" + cpe["pre"]["neuron_type"]
                                                            + "_to_" + cpe["post"]["region"] + "_" + cpe["post"]["subregion"] + "_" + cpe["post"]["neuron_type"] )
                    # pre, post names 
                    ct[-1].pre_rg_name = cpe["pre"]["region"] 
                    ct[-1].pre_sr_name = cpe["pre"]["subregion"] 
                    ct[-1].pre_nt_name = cpe["pre"]["neuron_type"] 
                    ct[-1].post_rg_name = cpe["post"]["region"] 
                    ct[-1].post_sr_name = cpe["post"]["subregion"] 
                    ct[-1].post_nt_name = cpe["post"]["neuron_type"]

                    # pre, post IDs
                    ct[-1].pre_rg_GID = i_rd
                    ct[-1].pre_sr_LID = rd_array[ct[-1].pre_rg_GID]["structure_info"]["subregion_name_to_subregion_LID"][cpe["pre"]["subregion"]]
                    ct[-1].pre_nt_LID = rd_array[ct[-1].pre_rg_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][cpe["pre"]["subregion"]][cpe["pre"]["neuron_type"]]
                    ct[-1].pre_nt_SerLID = rd_array[ct[-1].pre_rg_GID]["structure_info"]["LID_to_neuron_type_SerLID"][ct[-1].pre_sr_LID][ct[-1].pre_nt_LID]
                    ct[-1].post_rg_GID = i_rd
                    ct[-1].post_sr_LID = rd_array[ct[-1].post_rg_GID]["structure_info"]["subregion_name_to_subregion_LID"][cpe["post"]["subregion"]] 
                    ct[-1].post_nt_LID = rd_array[ct[-1].post_rg_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][cpe["post"]["subregion"]][cpe["post"]["neuron_type"]]

                    ct[-1].post_nt_SerLID = rd_array[ct[-1].post_rg_GID]["structure_info"]["LID_to_neuron_type_SerLID"][ct[-1].post_sr_LID][ct[-1].post_nt_LID]

                    # LID per neuron_type
                    ct[-1].connection_type_LID_in_nt = connection_type_SerLIDs_per_post[ct[-1].post_nt_SerLID]

                    # intra or inter
                    ct[-1].inter_or_intra = "intra"

                    # connection properties                    
                    ct[-1].spatial_pattern = cpe["spatial_pattern"]
                    ct[-1].synaptic_channel = sc
                    ct[-1].colocalization = 1 if len(cpe["synaptic_channel"]) == 2 else 0
                    ct[-1].weight_distribution = sc[str(*sc.keys())]["weight_distribution"]
                    ct[-1].delay = cpe["delay"]
                    ct[-1].LTP = cpe["LTP"]

                    # make intra pre post pair
                    if [ct[-1].pre_nt_SerLID, ct[-1].pre_rg_GID, ct[-1].pre_sr_LID, ct[-1].pre_nt_LID,
                        ct[-1].post_nt_SerLID, ct[-1].post_rg_GID, ct[-1].post_sr_LID, ct[-1].post_nt_LID] not in intra_regional_connection_nt_pairs:
                        intra_regional_connection_nt_pairs.append([ct[-1].pre_nt_SerLID, ct[-1].pre_rg_GID, ct[-1].pre_sr_LID, ct[-1].pre_nt_LID,
                                                                   ct[-1].post_nt_SerLID, ct[-1].post_rg_GID, ct[-1].post_sr_LID, ct[-1].post_nt_LID])
                        
                    # count numbers
                    connection_type_SerLIDs_per_post[ct[-1].post_nt_SerLID] += 1
                    count_ct+=1
                    n_indegree_intra_regional_connection_types_per_process += 1
                    
                    if str(*cpe["spatial_pattern"].keys()) == "two_dimensional_gaussian":
                        if isinstance(cpe["spatial_pattern"]["two_dimensional_gaussian"]["mu"], str):
                            temp_coefficient_name = cpe["spatial_pattern"]["two_dimensional_gaussian"]["mu"]
                            temp_coefficient = rd["intra_regional_connection"]["coefficients"][temp_coefficient_name]
                            if temp_coefficient > 0.:
                                conn_post_nts_per_pre_nt[ct[-1].pre_nt_SerLID][ct[-1].post_nt_SerLID].append(1)
                                conn_post_nt_SIDs_per_pre_nt_SID[ct[-1].pre_nt_SerLID].append(ct[-1].post_nt_SerLID)
                                n_conn_pairs += 1
                                
                            else:
                                n_no_conn_pairs += 1

                        else:
                            if cpe["spatial_pattern"]["two_dimensional_gaussian"]["mu"] > 0.:

                                n_conn_pairs += 1
                                #n_conn_posts_per_pre[ct[-1].pre_nt_SerLID] += 1
                                conn_post_nts_per_pre_nt[ct[-1].pre_nt_SerLID][ct[-1].post_nt_SerLID].append(1)
                                """
                                ct[-1].post_nt_SerLID, ct[-1].post_nt_LID
                                """

                            else:
                                n_no_conn_pairs += 1

                    elif str(*cpe["spatial_pattern"].keys()) == "orthogornal_cross":
                        if cpe["spatial_pattern"]["orthogornal_cross"]["probability"] > 0.:
                            print("B")
                    elif str(*cpe["spatial_pattern"].keys()) == "circular":
                        if cpe["spatial_pattern"]["circular"]["probability"] > 0.:
                            print("B")
                    elif str(*cpe["spatial_pattern"].keys()) == "orthogornal_cross":
                        if cpe["spatial_pattern"]["square"]["probability"] > 0.:
                            print("B")
                    elif str(*cpe["spatial_pattern"].keys()) == "ID_connect":
                        if cpe["spatial_pattern"]["ID_connect"]["probability"] > 0.:
                            print("B")
                            
                
        #print("n_conn_/no_conn pairs", n_conn_pairs, n_no_conn_pairs )
        #print(len(conn_post_nts_per_pre_nt))
        #print(conn_post_nts_per_pre_nt)
        #print(conn_post_nt_SIDs_per_pre_nt_SID)

        """
        # intra ID connections
        if "ID_connect" in rd["intra_regional_connection"].keys():
            for i_ID_conn, ID_conn in enumerate(rd["intra_regional_connection"]["ID_connect"]):
            
                # pre, post colocalization, and connection_type_SerLID
                if len(ID_conn["synaptic_channel"]) == 1:
                    rd["intra_regional_connection"]["ID_connect"][i_ID_conn]["pre_post_colocalization"] = [
                        sd["region_name_to_region_GID"][ID_conn["pre"]["region"]],
                        rd["structure_info"]["subregion_name_to_subregion_LID"][ID_conn["pre"]["subregion"]],
                        rd["structure_info"]["neuron_type_name_to_neuron_type_LID"][ID_conn["pre"]["subregion"]][ID_conn["pre"]["neuron_type"]],
                        sd["region_name_to_region_GID"][ID_conn["post"]["region"]],
                        rd["structure_info"]["subregion_name_to_subregion_LID"][ID_conn["post"]["subregion"]],
                        rd["structure_info"]["neuron_type_name_to_neuron_type_LID"][ID_conn["post"]["subregion"]][ID_conn["post"]["neuron_type"]],
                        0,
                        count_ct,
                        -1,
                    ]
                elif len(ID_conn["synaptic_channel"]) == 2:
                    rd["intra_regional_connection"]["ID_connect"][i_ID_conn]["pre_post_colocalization"] = [
                        sd["region_name_to_region_GID"][ID_conn["pre"]["region"]],
                        rd["structure_info"]["subregion_name_to_subregion_LID"][ID_conn["pre"]["subregion"]],
                        rd["structure_info"]["neuron_type_name_to_neuron_type_LID"][ID_conn["pre"]["subregion"]][ID_conn["pre"]["neuron_type"]],
                        sd["region_name_to_region_GID"][ID_conn["post"]["region"]],
                        rd["structure_info"]["subregion_name_to_subregion_LID"][ID_conn["post"]["subregion"]],
                        rd["structure_info"]["neuron_type_name_to_neuron_type_LID"][ID_conn["post"]["subregion"]][ID_conn["post"]["neuron_type"]],
                        1,
                        count_ct,
                        count_ct+1
                    ]
                else:
                    print ("number of channel was wrong", len(ID_conn["synaptic_channel"]))


                for i_sc, sc in enumerate(ID_conn["synaptic_channel"]):

                    ct.append(connection_type())
                    ct[-1].connection_type_SerLID = count_ct
                    ct[-1].connection_type_full_name = str( ID_conn["pre"]["region"] + "_" + ID_conn["pre"]["subregion"] + "_" + ID_conn["pre"]["neuron_type"]
                                                            + "_to_" + ID_conn["post"]["region"] + "_" + ID_conn["post"]["subregion"] + "_" + ID_conn["post"]["neuron_type"] )
                    # pre, post names 
                    ct[-1].pre_rg_name = ID_conn["pre"]["region"] 
                    ct[-1].pre_sr_name = ID_conn["pre"]["subregion"] 
                    ct[-1].pre_nt_name = ID_conn["pre"]["neuron_type"] 
                    ct[-1].post_rg_name = ID_conn["post"]["region"] 
                    ct[-1].post_sr_name = ID_conn["post"]["subregion"] 
                    ct[-1].post_nt_name = ID_conn["post"]["neuron_type"]

                    # pre, post IDs
                    ct[-1].pre_rg_GID = i_rd
                    ct[-1].pre_sr_LID = rd_array[ct[-1].pre_rg_GID]["structure_info"]["subregion_name_to_subregion_LID"][ID_conn["pre"]["subregion"]]
                    ct[-1].pre_nt_LID = rd_array[ct[-1].pre_rg_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][ID_conn["pre"]["subregion"]][ID_conn["pre"]["neuron_type"]]
                    ct[-1].pre_nt_SerLID = rd_array[ct[-1].pre_rg_GID]["structure_info"]["LID_to_neuron_type_SerLID"][ct[-1].pre_sr_LID][ct[-1].pre_nt_LID]
                    ct[-1].post_rg_GID = i_rd
                    ct[-1].post_sr_LID = rd_array[ct[-1].post_rg_GID]["structure_info"]["subregion_name_to_subregion_LID"][ID_conn["post"]["subregion"]] 
                    ct[-1].post_nt_LID = rd_array[ct[-1].post_rg_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][ID_conn["post"]["subregion"]][ID_conn["post"]["neuron_type"]]
                    ct[-1].post_nt_SerLID = rd_array[ct[-1].post_rg_GID]["structure_info"]["LID_to_neuron_type_SerLID"][ct[-1].post_sr_LID][ct[-1].post_nt_LID]

                    # LID per neuron_type
                    ct[-1].connection_type_LID_in_nt = connection_type_SerLIDs_per_post[ct[-1].post_nt_SerLID]

                    # intra or inter
                    ct[-1].inter_or_intra = "intra"

                    # connection properties                    
                    ct[-1].spatial_pattern = ID_conn["spatial_pattern"]
                    ct[-1].synaptic_channel = sc
                    ct[-1].colocalization = 1 if len(ID_conn["synaptic_channel"])==2 else 0
                    ct[-1].weight_distribution = sc[str(*sc.keys())]["weight_distribution"]
                    ct[-1].delay = ID_conn["delay"]
                    ct[-1].LTP = ID_conn["LTP"]

                    # make intra pre post pair
                    if [ct[-1].pre_nt_SerLID, ct[-1].pre_rg_GID, ct[-1].pre_sr_LID, ct[-1].pre_nt_LID,
                        ct[-1].post_nt_SerLID, ct[-1].post_rg_GID, ct[-1].post_sr_LID, ct[-1].post_nt_LID] not in intra_regional_connection_nt_pairs:
                        intra_regional_connection_nt_pairs.append([ct[-1].pre_nt_SerLID, ct[-1].pre_rg_GID, ct[-1].pre_sr_LID, ct[-1].pre_nt_LID,
                                                                   ct[-1].post_nt_SerLID, ct[-1].post_rg_GID, ct[-1].post_sr_LID, ct[-1].post_nt_LID])
                        
                    # count numbers
                    connection_type_SerLIDs_per_post[ct[-1].post_nt_SerLID]+=1
                    count_ct+=1
                    n_indegree_intra_regional_connection_types_per_process+=1
        """
        

        # connection type of inter regional connection (loop bundles and connection types)
        for i_b, bundle in enumerate(rd["inter_regional_connection"]):
            count_connection_type_LID_in_bundle = 0
            
            #bundle_name = "bundle_" + bundle[i_b]["pre"]["region"] + "_to_" + bundle[i_b]["post"]["region"]
            bundle_name = "bundle_" + bundle["pre"]["region"] + "_to_" + bundle["post"]["region"]
            bdl[bundle_name] = bundle_class()
            bdl[bundle_name].bundle_name = bundle_name
            bdl[bundle_name].conn = []
            bdl[bundle_name].post_nt_per_pre_nt_at_bundle = []
            bdl[bundle_name].n_pre_nt_at_bundle = 0
            bdl[bundle_name].pre_nt_SIDs_at_bundle = []
            bdl[bundle_name].post_nt_per_pre_nt_at_bundle = []
            bdl[bundle_name].pre_region_GID = sd["region_name_to_region_GID"][bundle["pre"]["region"]]
            bdl[bundle_name].post_region_GID = sd["region_name_to_region_GID"][bundle["post"]["region"]]
            bdl[bundle_name].shift = bundle["topology"]["shift"]
            bdl[bundle_name].scale = bundle["topology"]["scale"]
            bdl[bundle_name].tile_link_limit = bundle["tile_link_limit"]
            
            #for i_ct, ct in enumerate(bundle["connection_type"]):            
            for cpe in bundle["connection_type"]:
                # regsiter IDs to connection (redundunt expresission)
                bdl[bundle_name].conn.append(connection_in_bundle())
                bdl[bundle_name].conn[-1].pre_rg_GID = sd["region_name_to_region_GID"][cpe["pre"]["region"]]
                bdl[bundle_name].conn[-1].pre_sr_LID = rd_array[bdl[bundle_name].conn[-1].pre_rg_GID]["structure_info"]["subregion_name_to_subregion_LID"][cpe["pre"]["subregion"]]
                bdl[bundle_name].conn[-1].pre_nt_LID = rd_array[bdl[bundle_name].conn[-1].pre_rg_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][cpe["pre"]["subregion"]][cpe["pre"]["neuron_type"]]
                bdl[bundle_name].conn[-1].post_rg_GID = sd["region_name_to_region_GID"][cpe["post"]["region"]]
                bdl[bundle_name].conn[-1].post_sr_LID = rd_array[bdl[bundle_name].conn[-1].post_rg_GID]["structure_info"]["subregion_name_to_subregion_LID"][cpe["post"]["subregion"]]
                bdl[bundle_name].conn[-1].post_nt_LID = rd_array[bdl[bundle_name].conn[-1].post_rg_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][cpe["post"]["subregion"]][cpe["post"]["neuron_type"]]
                
                # pre, post colocalization, and connection_type_SerLID
                if len(cpe["synaptic_channel"]) == 1:
                    bdl[bundle_name].conn[-1].pre_post_colocalization = [
                        sd["region_name_to_region_GID"][cpe["pre"]["region"]],
                        rd_array[bdl[bundle_name].conn[-1].pre_rg_GID]["structure_info"]["subregion_name_to_subregion_LID"][cpe["pre"]["subregion"]],
                        rd_array[bdl[bundle_name].conn[-1].pre_rg_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][cpe["pre"]["subregion"]][cpe["pre"]["neuron_type"]],
                        sd["region_name_to_region_GID"][cpe["post"]["region"]],
                        rd_array[bdl[bundle_name].conn[-1].post_rg_GID]["structure_info"]["subregion_name_to_subregion_LID"][cpe["post"]["subregion"]],
                        rd_array[bdl[bundle_name].conn[-1].post_rg_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][cpe["post"]["subregion"]][cpe["post"]["neuron_type"]],
                        0,
                        count_ct,
                        -1,
                    ]
                elif len(cpe["synaptic_channel"]) == 2:
                    bdl[bundle_name].conn[-1].pre_post_colocalization = [
                        sd["region_name_to_region_GID"][cpe["pre"]["region"]],
                        rd_array[ct[-1].pre_rg_GID]["structure_info"]["subregion_name_to_subregion_LID"][cpe["pre"]["subregion"]],
                        rd_array[ct[-1].pre_rg_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][cpe["pre"]["subregion"]][cpe["pre"]["neuron_type"]],
                        sd["region_name_to_region_GID"][cpe["post"]["region"]],
                        rd_array[ct[-1].post_rg_GID]["structure_info"]["subregion_name_to_subregion_LID"][cpe["post"]["subregion"]],
                        rd_array[ct[-1].post_rg_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][cpe["post"]["subregion"]][cpe["post"]["neuron_type"]],
                        1,
                        count_ct,
                        count_ct+1
                    ]
                else:
                    print ("number of channel was wrong", len(cpe["synaptic_channel"]))
                    
                for i_sc, sc in enumerate(cpe["synaptic_channel"]):
                    ct.append(connection_type())
                    ct[-1].connection_type_SerLID = count_ct
                    ct[-1].connection_type_full_name = str( cpe["pre"]["region"] + "_" + cpe["pre"]["subregion"] + "_" + cpe["pre"]["neuron_type"]
                                                            + "_to_" + cpe["post"]["region"] + "_" + cpe["post"]["subregion"] + "_" + cpe["post"]["neuron_type"] )
                    # pre, post names 
                    ct[-1].pre_rg_name = cpe["pre"]["region"] 
                    ct[-1].pre_sr_name = cpe["pre"]["subregion"] 
                    ct[-1].pre_nt_name = cpe["pre"]["neuron_type"] 
                    ct[-1].post_rg_name = cpe["post"]["region"] 
                    ct[-1].post_sr_name = cpe["post"]["subregion"] 
                    ct[-1].post_nt_name = cpe["post"]["neuron_type"]
                    
                    # pre, post IDs
                    ct[-1].pre_rg_GID = sd["region_name_to_region_GID"][cpe["pre"]["region"]]
                    ct[-1].pre_sr_LID = rd_array[ct[-1].pre_rg_GID]["structure_info"]["subregion_name_to_subregion_LID"][cpe["pre"]["subregion"]]
                    ct[-1].pre_nt_LID = rd_array[ct[-1].pre_rg_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][cpe["pre"]["subregion"]][cpe["pre"]["neuron_type"]]
                    ct[-1].pre_nt_SerLID = rd_array[ct[-1].pre_rg_GID]["structure_info"]["LID_to_neuron_type_SerLID"][ct[-1].pre_sr_LID][ct[-1].pre_nt_LID]
                    ct[-1].post_rg_GID = sd["region_name_to_region_GID"][cpe["post"]["region"]]
                    ct[-1].post_sr_LID = rd_array[ct[-1].post_rg_GID]["structure_info"]["subregion_name_to_subregion_LID"][cpe["post"]["subregion"]] 
                    ct[-1].post_nt_LID = rd_array[ct[-1].post_rg_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][cpe["post"]["subregion"]][cpe["post"]["neuron_type"]]
                    ct[-1].post_nt_SerLID = rd_array[ct[-1].post_rg_GID]["structure_info"]["LID_to_neuron_type_SerLID"][ct[-1].post_sr_LID][ct[-1].post_nt_LID]

                    # LID per neuron_type 
                    ct[-1].connection_type_LID_in_nt = connection_type_SerLIDs_per_post[ct[-1].post_nt_SerLID]

                    # intra or inter
                    ct[-1].inter_or_intra = "inter"

                    # bundle information
                    ct[-1].bundle_name = "bundle_" + cpe["pre"]["region"] + "_to_" + cpe["post"]["region"] 
                    ct[-1].connection_type_LID_in_bundle = count_connection_type_LID_in_bundle

                    # connection properties
                    ct[-1].spatial_pattern = cpe["spatial_pattern"]
                    ct[-1].synaptic_channel = sc
                    ct[-1].colocalization = 1 if len(cpe["synaptic_channel"]) == 2 else 0
                    ct[-1].weight_distribution = sc[str(*sc.keys())]["weight_distribution"]
                    ct[-1].delay = cpe["delay"]
                    ct[-1].LTP = cpe["LTP"]

                    # make inter pre post pair, not connection type pair including colocalization channels
                    if  [ct[-1].pre_nt_SerLID, ct[-1].pre_rg_GID, ct[-1].pre_sr_LID, ct[-1].pre_nt_LID,
                         ct[-1].post_nt_SerLID, ct[-1].post_rg_GID, ct[-1].post_sr_LID, ct[-1].post_nt_LID] not in inter_regional_connection_nt_pairs:
                        inter_regional_connection_nt_pairs.append([ct[-1].pre_nt_SerLID, ct[-1].pre_rg_GID, ct[-1].pre_sr_LID, ct[-1].pre_nt_LID,
                                                                   ct[-1].post_nt_SerLID, ct[-1].post_rg_GID, ct[-1].post_sr_LID, ct[-1].post_nt_LID])

                    # count numbers of connection type per post neuron type
                    connection_type_SerLIDs_per_post[ct[-1].post_nt_SerLID] += 1

                    # count interregional connection types per bundle
                    count_connection_type_LID_in_bundle+=1

                    # count interregional connection types per process
                    n_indegree_inter_regional_connection_types_per_process+=1

                    # count interregional connection types per bundle (for preservation in bdl)
                    bdl[bundle_name].n_indegree_inter_regional_connection_types_per_bundle += 1

                    # count all connection types in intra, inter regional connection
                    count_ct+=1


        # dummy connection type for padding
        
        # prepare connection_type_SerLIDs_per_post
        connection_type_SerLIDs_per_post = [ [] for i in range(rd["structure_info"]["n_neuron_types_per_process"]) ]
        
        for i_ins_ct, ins_ct in enumerate(ct):
            connection_type_SerLIDs_per_post[ins_ct.post_nt_SerLID].append([ins_ct.connection_type_LID_in_nt, ins_ct.connection_type_SerLID])

        ct.append(connection_type())
        ct[-1].connection_type_SerLID = count_ct
        ct[-1].connection_type_full_name = "dummy"
        # pre, post names 
        ct[-1].pre_rg_name = "dummy" 
        ct[-1].pre_sr_name = "dummy" 
        ct[-1].pre_nt_name = "dummy" 
        ct[-1].post_rg_name = "dummy"
        ct[-1].post_sr_name = "dummy"
        ct[-1].post_nt_name = "dummy"
        
        # dummy use #0 of pre, post IDs
        ct[-1].pre_rg_GID = 0
        ct[-1].pre_sr_LID = 0
        ct[-1].pre_nt_LID = 0
        ct[-1].pre_nt_SerLID = 0
        ct[-1].post_rg_GID = 0
        ct[-1].post_sr_LID = 0
        ct[-1].post_nt_LID = 0
        ct[-1].post_nt_SerLID = 0
        
        # LID per neuron_type
        ct[-1].connection_type_LID_in_nt = 0
        
        # intra or inter
        ct[-1].inter_or_intra = "dummy" 
        
        # connection properties                    
        ct[-1].spatial_pattern = { "two_dimensional_gaussian": { "mu":0, "sd":200 } }
        ct[-1].synaptic_channel = { "AMPA": { "dynamics": { "alpha": { "tau":2 } }, "reversal_potential":0, "weight_distribution": { "uniform_dist": {"uniform_value": 0} } } }
        ct[-1].colocalization = 0
        ct[-1].weight_distribution = { "uniform_dist": {"uniform_value": 0} }
        ct[-1].delay = 1.5
        ct[-1].LTP = {}

        rd_array[i_rd]["n_indegree_dummy_connection_types_per_process"] = 1
        # count all connection types in intra, inter regional connection
        count_ct+=1
                    
        # register bundles and connection types to rd_array
        rd["bdl"] = copy.deepcopy(bdl)
        rd["ct"] = copy.deepcopy(ct)
 
        # padding dummy if n_ct_per_post is short
        n_ct_for_padding = len(rd["ct"])
        for i_cstspp, ctspp in enumerate(connection_type_SerLIDs_per_post):
            if len(ctspp) < rd_array[i_rd]["n_connection_types_per_post_with_padding"][i_cstspp]:
                connection_type_SerLIDs_per_post[i_cstspp].append([ctspp[-1][0]+1, n_ct_for_padding-1])
            
       
        # finally, register various count of connection type
        rd["connection_type_SerLIDs_per_post"] = copy.deepcopy(connection_type_SerLIDs_per_post)
        rd["n_indegree_intra_regional_connection_types_per_process"] = n_indegree_intra_regional_connection_types_per_process
        rd["n_indegree_inter_regional_connection_types_per_process"] = n_indegree_inter_regional_connection_types_per_process
        rd["intra_regional_connection_nt_pairs"] = intra_regional_connection_nt_pairs
        rd["inter_regional_connection_nt_pairs"] = inter_regional_connection_nt_pairs
        rd["count_ct"] = count_ct
        
    # set outdegree bundle information for mpi communication send ????
    for i_rd, rd in enumerate(rd_array):
        for i_b, bundle in enumerate(rd["inter_regional_connection"]):
            pre_region_GID = sd["region_name_to_region_GID"][bundle["pre"]["region"]]
            post_region_GID = sd["region_name_to_region_GID"][bundle["post"]["region"]]
            bundle_name = "bundle_" + bundle["pre"]["region"] + "_to_" + bundle["post"]["region"]
            rd_array[pre_region_GID]["bdl"][bundle_name] = bundle_class()
            rd_array[pre_region_GID]["bdl"][bundle_name].bundle_name = bundle_name
            rd_array[post_region_GID]["bdl"][bundle_name].post_nt_per_pre_nt_at_bundle = [ [] for i_pns in range(rd_array[pre_region_GID]["structure_info"]["n_neuron_types_per_process"])]
            rd_array[post_region_GID]["bdl"][bundle_name].topology = rd_array[post_region_GID]["inter_regional_connection"][i_b]["topology"]
            rd_array[pre_region_GID]["bdl"][bundle_name].pre_region_GID = pre_region_GID 
            rd_array[pre_region_GID]["bdl"][bundle_name].post_region_GID = post_region_GID 
            rd_array[pre_region_GID]["bdl"][bundle_name].shift = bundle["topology"]["shift"]
            rd_array[pre_region_GID]["bdl"][bundle_name].scale = bundle["topology"]["scale"]
            rd_array[pre_region_GID]["bdl"][bundle_name].tile_link_limit = bundle["tile_link_limit"]
             #print(rd_array[post_region_GID]["inter_regional_connection"][i_b]["topology"])
            
            # make post_nt_per_per_nt_at_bundle
            for i_ct, conn in enumerate(bundle["connection_type"]):
                pre_sr_LID = rd_array[pre_region_GID]["structure_info"]["subregion_name_to_subregion_LID"][conn["pre"]["subregion"]] 
                pre_nt_LID = rd_array[pre_region_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][conn["pre"]["subregion"]][conn["pre"]["neuron_type"]]
                pre_nt_SerLID = rd_array[pre_region_GID]["structure_info"]["LID_to_neuron_type_SerLID"][pre_sr_LID][pre_nt_LID]
                post_sr_LID = rd_array[post_region_GID]["structure_info"]["subregion_name_to_subregion_LID"][conn["post"]["subregion"]] 
                post_nt_LID = rd_array[post_region_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][conn["post"]["subregion"]][conn["post"]["neuron_type"]]
                post_nt_SerLID = rd_array[post_region_GID]["structure_info"]["LID_to_neuron_type_SerLID"][post_sr_LID][post_nt_LID]
                rd_array[post_region_GID]["bdl"][bundle_name].post_nt_per_pre_nt_at_bundle[pre_nt_SerLID].append([post_sr_LID, post_nt_LID, post_nt_SerLID])

            # make pre_nt_SID_to_pre_nt_LID_at_bundle and count n_pre_nt_at_bundle 
            count_pre_nt_LID_at_bundle = 0
            for i_post_nt_group, post_nt_group in enumerate(rd_array[post_region_GID]["bdl"][bundle_name].post_nt_per_pre_nt_at_bundle):
                if post_nt_group != []:
                    #pre_nt_SID_to_pre_nt_LID_at_bundle[i_post_nt_group] = count_pre_nt_LID_at_bundle
                    pre_nt_LID_at_bundle = count_pre_nt_LID_at_bundle
                    rd_array[post_region_GID]["bdl"][bundle_name].pre_nt_SIDs_at_bundle.append(i_post_nt_group)
                    count_pre_nt_LID_at_bundle += 1
            rd_array[post_region_GID]["bdl"][bundle_name].n_pre_nt_at_bundle = count_pre_nt_LID_at_bundle

            
    for i_rd, rd in enumerate(rd_array):
        for i_bn, bundle_name in enumerate(rd_array[i_rd]["bdl"].keys()):
            rd_array[i_rd]["bdl"][bundle_name].bundle_LID = i_bn
            #print("a", rd_array[i_rd]["structure_info"]["region_name"], bundle_name, rd_array[i_rd]["bdl"][bundle_name].bundle_LID)
            
            for i_ct, ct in enumerate(rd_array[i_rd]["ct"]):
                if ct.inter_or_intra == "inter":
                    if "bundle_" + ct.pre_rg_name + "_to_" + ct.post_rg_name == bundle_name:
                        rd_array[i_rd]["ct"][i_ct].bundle_LID = i_bn
                        #print("b", rd_array[i_rd]["structure_info"]["region_name"], bundle_name, rd_array[i_rd]["ct"][i_ct].pre_rg_name, rd_array[i_rd]["ct"][i_ct].post_rg_name, rd_array[i_rd]["ct"][i_ct].bundle_LID, rd_array[i_rd]["ct"][i_ct].post_rg_GID, "i_ct:", i_ct)
            
     
    elapsed_time = time.time() - start
    print ("set_region_connection2: {0}".format(elapsed_time) + "[sec]")
            
    #####################################
    # set region informaiton to process #
    #####################################

    start = time.time()
    """
    for i_rd, rd in enumerate(rd_array):
        for i_post_proc in rd["process_ID"]:
            process[i_post_proc].n_indegree_intra_regional_connection_types_per_process = rd["n_indegree_intra_regional_connection_types_per_process"]
            process[i_post_proc].bdl = copy.deepcopy(rd["bdl"])
    """
    elapsed_time = time.time() - start
    print ("set_region_connection3.1: {0}".format(elapsed_time) + "[sec]")

    # make intra_regional_process_pairs
    start = time.time()
    
    for i_rd, rd in enumerate(rd_array):
        """
        # no tile group mode
        if rd["structure_info"]["n_tile_groups"] == 0:

            for i_post_proc in rd["process_ID"]:
                process[i_post_proc].intra_regional_pre_post_address3 = [ [ [] for i in rd["intra_regional_connection"]["connection_parameter_matrix"] ] for j in range(len(process[i_post_proc].intra_regional_process_pairs)) ]
                for i_pcp, pcp in enumerate(process[i_post_proc].intra_regional_process_pairs):
                    for post_nt_SerLID, row in enumerate(rd["intra_regional_connection"]["connection_parameter_matrix"]):
                        for pre_nt_SerLID, cpe in enumerate(row):
                            process[i_post_proc].intra_regional_pre_post_address3[i_pcp][pre_nt_SerLID].append([
                                cpe["pre_post_colocalization"][6],
                                cpe["pre_post_colocalization"][7],
                                cpe["pre_post_colocalization"][8],
                                pcp[0],
                                cpe["pre_post_colocalization"][0],
                                cpe["pre_post_colocalization"][1],
                                cpe["pre_post_colocalization"][2],
                                process[i_post_proc].process_GID,
                                cpe["pre_post_colocalization"][3],
                                cpe["pre_post_colocalization"][4],
                                cpe["pre_post_colocalization"][5],
                            ])
                            
        # tile group mode 
        elif rd["structure_info"]["n_tile_groups"] > 0:
            for i_post_proc in rd["process_ID"]:
                process[i_post_proc].intra_regional_pre_post_address3 = [ [ [] for i in rd["intra_regional_connection"]["connection_parameter_matrix"] ] for j in range(len(process[i_post_proc].intra_regional_process_pairs)) ]

                for i_pcp, pcp in enumerate(process[i_post_proc].intra_regional_process_pairs):
                    for post_nt_SerLID, row in enumerate(rd["intra_regional_connection"]["connection_parameter_matrix"]):
                        for pre_nt_SerLID, cpe in enumerate(row):
                            process[i_post_proc].intra_regional_pre_post_address3[i_pcp][pre_nt_SerLID].append([
                                cpe["pre_post_colocalization"][6],
                                cpe["pre_post_colocalization"][7],
                                cpe["pre_post_colocalization"][8],
                                pcp[0],
                                cpe["pre_post_colocalization"][0],
                                cpe["pre_post_colocalization"][1],
                                cpe["pre_post_colocalization"][2],
                                process[i_post_proc].process_GID,
                                cpe["pre_post_colocalization"][3],
                                cpe["pre_post_colocalization"][4],
                                cpe["pre_post_colocalization"][5],
                            ])
                            
        """
        
        # set inter regional connections to process
        for i_b, bundle in enumerate(rd["inter_regional_connection"]):
            post_region_GID = sd["region_name_to_region_GID"][bundle["post"]["region"]]
            pre_region_GID = sd["region_name_to_region_GID"][bundle["pre"]["region"]]
            bundle_name = "bundle_" + bundle["pre"]["region"] + "_to_" + bundle["post"]["region"]
            rd_array[post_region_GID]["bdl"][bundle_name].post_nt_per_pre_nt_at_bundle = [ [] for i_pns in range(rd_array[pre_region_GID]["structure_info"]["n_neuron_types_per_process"])]
            pre_nt_SID_to_pre_nt_LID_at_bundle = [ 0 for i_npns in range(rd_array[pre_region_GID]["structure_info"]["n_neuron_types_per_process"])]

            for i_conn, conn in enumerate(bundle["connection_type"]):
                pre_sr_LID = rd_array[pre_region_GID]["structure_info"]["subregion_name_to_subregion_LID"][conn["pre"]["subregion"]] 
                pre_nt_LID = rd_array[pre_region_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][conn["pre"]["subregion"]][conn["pre"]["neuron_type"]]
                pre_nt_SerLID = rd_array[pre_region_GID]["structure_info"]["LID_to_neuron_type_SerLID"][pre_sr_LID][pre_nt_LID]
                post_sr_LID = rd_array[post_region_GID]["structure_info"]["subregion_name_to_subregion_LID"][conn["post"]["subregion"]] 
                post_nt_LID = rd_array[post_region_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][conn["post"]["subregion"]][conn["post"]["neuron_type"]]
                post_nt_SerLID = rd_array[post_region_GID]["structure_info"]["LID_to_neuron_type_SerLID"][post_sr_LID][post_nt_LID]

                rd_array[post_region_GID]["bdl"][bundle_name].post_nt_per_pre_nt_at_bundle[pre_nt_SerLID].append([post_sr_LID, post_nt_LID, post_nt_SerLID])
            #print(rd_array[post_region_GID]["bdl"][bundle_name].post_nt_per_pre_nt_at_bundle)
            
            count_pre_nt_LID_at_bundle = 0
            for i_post_nt_group, post_nt_group in enumerate(rd_array[post_region_GID]["bdl"][bundle_name].post_nt_per_pre_nt_at_bundle):
                if post_nt_group != []:
                    pre_nt_SID_to_pre_nt_LID_at_bundle[i_post_nt_group] = count_pre_nt_LID_at_bundle

                    rd_array[post_region_GID]["bdl"][bundle_name].pre_nt_SIDs_at_bundle.append(i_post_nt_group)
                    count_pre_nt_LID_at_bundle += 1
            rd_array[post_region_GID]["bdl"][bundle_name].n_pre_nt_at_bundle = count_pre_nt_LID_at_bundle

            #print(pre_nt_SID_to_pre_nt_LID_at_bundle)
            #print(rd_array[post_region_GID]["bdl"][bundle_name].pre_nt_SIDs_at_bundle)
            #print(rd_array[post_region_GID]["bdl"][bundle_name].n_pre_nt_at_bundle)

            rd_array[post_region_GID]["bdl"][bundle_name].new_inter_regional_post_neuron_types_per_pre_neuron_type_at_bundle = [ [] for i_pns in range(rd_array[post_region_GID]["bdl"][bundle_name].n_pre_nt_at_bundle)]
            for i_conn, conn in enumerate(bundle["connection_type"]):
                pre_sr_LID = rd_array[pre_region_GID]["structure_info"]["subregion_name_to_subregion_LID"][conn["pre"]["subregion"]] 
                pre_nt_LID = rd_array[pre_region_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][conn["pre"]["subregion"]][conn["pre"]["neuron_type"]]
                pre_nt_SerLID = rd_array[pre_region_GID]["structure_info"]["LID_to_neuron_type_SerLID"][pre_sr_LID][pre_nt_LID]
                post_sr_LID = rd_array[post_region_GID]["structure_info"]["subregion_name_to_subregion_LID"][conn["post"]["subregion"]] 
                post_nt_LID = rd_array[post_region_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][conn["post"]["subregion"]][conn["post"]["neuron_type"]]
                post_nt_SerLID = rd_array[post_region_GID]["structure_info"]["LID_to_neuron_type_SerLID"][post_sr_LID][post_nt_LID]
            
                pre_nt_LID_at_bundle = pre_nt_SID_to_pre_nt_LID_at_bundle[pre_nt_SerLID]
                #print(pre_nt_LID_at_bundle, rd_array[post_region_GID]["bdl"][bundle_name].new_inter_regional_post_neuron_types_per_pre_neuron_type_at_bundle)
                rd_array[post_region_GID]["bdl"][bundle_name].new_inter_regional_post_neuron_types_per_pre_neuron_type_at_bundle[pre_nt_LID_at_bundle].append([
                    rd["bdl"][bundle_name].conn[i_conn].pre_post_colocalization[6],
                    rd["bdl"][bundle_name].conn[i_conn].pre_post_colocalization[7],
                    rd["bdl"][bundle_name].conn[i_conn].pre_post_colocalization[8],
                    pre_region_GID,
                    pre_sr_LID,
                    pre_nt_LID,
                    pre_nt_SerLID,
                    post_region_GID, 
                    post_sr_LID,
                    post_nt_LID,
                    post_nt_SerLID,
                ])
                
            #print(rd_array[post_region_GID]["bdl"][bundle_name].new_inter_regional_post_neuron_types_per_pre_neuron_type_at_bundle)
            """    
            # process to process
            for i_post_proc in rd_array[post_region_GID]["process_ID"]:
                for i_pre_proc in rd_array[pre_region_GID]["process_ID"]:
                    # calculate distance between process
                    
                    distance = minimum_distance_between_4vertices_with_shift(process[i_pre_proc].xy_plane_4vertices, process[i_post_proc].xy_plane_4vertices, bundle["topology"]["shift"])
                    if (distance < bundle["tile_link_limit"]):

                        # add pre process GID to post synaptic process
                        if process[i_pre_proc].process_GID not in process[i_post_proc].bdl[bundle_name].inter_regional_pre_process_per_bundle:
                            process[i_post_proc].bdl[bundle_name].inter_regional_pre_process_per_bundle.append(process[i_pre_proc].process_GID)

                        # add post process GID to pre synaptic process
                        if process[i_post_proc].process_GID not in process[i_pre_proc].bdl[bundle_name].inter_regional_post_process_per_bundle:
                            process[i_pre_proc].bdl[bundle_name].inter_regional_post_process_per_bundle.append(process[i_post_proc].process_GID)
            
            # connection type 
            for i_ct, conn in enumerate(bundle["connection_type"]):
                for i_post_proc in rd_array[post_region_GID]["process_ID"]:
                    for i_pre_proc in rd_array[pre_region_GID]["process_ID"]:

                        distance = minimum_distance_between_4vertices_with_shift(process[i_pre_proc].xy_plane_4vertices, process[i_post_proc].xy_plane_4vertices, bundle["topology"]["shift"])

                        if (distance < bundle["tile_link_limit"]):
                            process[i_post_proc].bdl[bundle_name].inter_regional_pre_post_address_per_bundle.append([
                                rd["bdl"][bundle_name].conn[i_ct].pre_post_colocalization[6],
                                rd["bdl"][bundle_name].conn[i_ct].pre_post_colocalization[7],
                                rd["bdl"][bundle_name].conn[i_ct].pre_post_colocalization[8],
                                i_pre_proc,
                                pre_region_GID,
                                rd_array[pre_region_GID]["structure_info"]["subregion_name_to_subregion_LID"][conn["pre"]["subregion"]],
                                rd_array[pre_region_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][conn["pre"]["subregion"]][conn["pre"]["neuron_type"]],
                                i_post_proc,
                                post_region_GID, 
                                rd_array[post_region_GID]["structure_info"]["subregion_name_to_subregion_LID"][conn["post"]["subregion"]],
                                rd_array[post_region_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][conn["post"]["subregion"]][conn["post"]["neuron_type"]]
                            ])

            """
            
            # calculate post neuron types per pre neuron types at each process-to-process combination
            pre_nt_SID_to_pre_nt_LID_at_bundle = [ 0 for i_npns in range(rd_array[pre_region_GID]["structure_info"]["n_neuron_types_per_process"])]
            rd_array[post_region_GID]["bdl"][bundle_name].post_nt_per_pre_nt_at_bundle = [ [] for i_pns in range(rd_array[pre_region_GID]["structure_info"]["n_neuron_types_per_process"])]
            
            # make post_nt_per_per_nt_at_bundle
            for i_ct, conn in enumerate(bundle["connection_type"]):
                pre_sr_LID = rd_array[pre_region_GID]["structure_info"]["subregion_name_to_subregion_LID"][conn["pre"]["subregion"]] 
                pre_nt_LID = rd_array[pre_region_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][conn["pre"]["subregion"]][conn["pre"]["neuron_type"]]
                pre_nt_SerLID = rd_array[pre_region_GID]["structure_info"]["LID_to_neuron_type_SerLID"][pre_sr_LID][pre_nt_LID]
                post_sr_LID = rd_array[post_region_GID]["structure_info"]["subregion_name_to_subregion_LID"][conn["post"]["subregion"]] 
                post_nt_LID = rd_array[post_region_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][conn["post"]["subregion"]][conn["post"]["neuron_type"]]
                post_nt_SerLID = rd_array[post_region_GID]["structure_info"]["LID_to_neuron_type_SerLID"][post_sr_LID][post_nt_LID]
                rd_array[post_region_GID]["bdl"][bundle_name].post_nt_per_pre_nt_at_bundle[pre_nt_SerLID].append([post_sr_LID, post_nt_LID, post_nt_SerLID])

            # make pre_nt_SID_to_pre_nt_LID_at_bundle and count n_pre_nt_at_bundle 
            count_pre_nt_LID_at_bundle = 0
            for i_post_nt_group, post_nt_group in enumerate(rd_array[post_region_GID]["bdl"][bundle_name].post_nt_per_pre_nt_at_bundle):
                if post_nt_group != []:
                    pre_nt_SID_to_pre_nt_LID_at_bundle[i_post_nt_group] = count_pre_nt_LID_at_bundle
                    pre_nt_LID_at_bundle = count_pre_nt_LID_at_bundle
                    rd_array[post_region_GID]["bdl"][bundle_name].pre_nt_SIDs_at_bundle.append(i_post_nt_group)
                    count_pre_nt_LID_at_bundle += 1
            rd_array[post_region_GID]["bdl"][bundle_name].n_pre_nt_at_bundle = count_pre_nt_LID_at_bundle

            """
            # prepare list, inter_regional_post_neuron_types_per_pre_neuron_type_per_bundle
            for post_proc_GID in rd_array[post_region_GID]["process_ID"]:
                n_pre_process = len(process[post_proc_GID].bdl[bundle_name].inter_regional_pre_process_per_bundle)
                process[post_proc_GID].bdl[bundle_name].inter_regional_post_neuron_types_per_pre_neuron_type_per_bundle = [[[] for j in range(rd_array[post_region_GID]["bdl"][bundle_name].n_pre_nt_at_bundle)] for i in range(n_pre_process)]

                
            # conn : nt x nt
            for i_ct, conn in enumerate(bundle["connection_type"]):
                pre_sr_LID = rd_array[pre_region_GID]["structure_info"]["subregion_name_to_subregion_LID"][conn["pre"]["subregion"]] 
                pre_nt_LID = rd_array[pre_region_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][conn["pre"]["subregion"]][conn["pre"]["neuron_type"]]
                pre_nt_SerLID = rd_array[pre_region_GID]["structure_info"]["LID_to_neuron_type_SerLID"][pre_sr_LID][pre_nt_LID]
                post_sr_LID = rd_array[post_region_GID]["structure_info"]["subregion_name_to_subregion_LID"][conn["post"]["subregion"]] 
                post_nt_LID = rd_array[post_region_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][conn["post"]["subregion"]][conn["post"]["neuron_type"]]
                post_nt_SerLID = rd_array[post_region_GID]["structure_info"]["LID_to_neuron_type_SerLID"][post_sr_LID][post_nt_LID]

                # process x process
                for post_proc_GID in rd_array[post_region_GID]["process_ID"]:
                    for i_pre_proc, pre_proc_GID in enumerate(process[post_proc_GID].bdl[bundle_name].inter_regional_pre_process_per_bundle):
                        pre_nt_LID_at_bundle = pre_nt_SID_to_pre_nt_LID_at_bundle[pre_nt_SerLID]
                        process[post_proc_GID].bdl[bundle_name].inter_regional_post_neuron_types_per_pre_neuron_type_per_bundle[i_pre_proc][pre_nt_LID_at_bundle].append([
                            rd["bdl"][bundle_name].conn[i_ct].pre_post_colocalization[6],
                            rd["bdl"][bundle_name].conn[i_ct].pre_post_colocalization[7],
                            rd["bdl"][bundle_name].conn[i_ct].pre_post_colocalization[8],
                            pre_proc_GID,
                            pre_region_GID,
                            pre_sr_LID,
                            pre_nt_LID,
                            post_proc_GID,
                            post_region_GID, 
                            post_sr_LID,
                            post_nt_LID,
                        ])
            """
            """            
            for post_proc_GID in rd_array[post_region_GID]["process_ID"]:
                print(bundle_name, post_region_GID, post_proc_GID, process[post_proc_GID].bdl[bundle_name].inter_regional_post_neuron_types_per_pre_neuron_type_per_bundle)
            """

    
    elapsed_time = time.time() - start
    print ("set_region_connection3.2: {0}".format(elapsed_time) + "[sec]")
    
    start = time.time()
    # set numbers of bundles to process
    for i_post_proc in range(sd["n_process"]):
        process[i_post_proc].n_bundles = len(rd_array[process[i_post_proc].region_GID]["bdl"])

    elapsed_time = time.time() - start
    print ("set_region_connection4: {0}".format(elapsed_time) + "[sec]")

    
    
def minimum_distance_between_4vertices(pre_4vertices, post_4vertices):
    # calculate distance of all combination of 4 vertices of pre and post tiles and judge
    distance = sys.float_info.max
    for pre_vs in pre_4vertices:
        for post_vs in post_4vertices:
            temp = math.sqrt( (pre_vs[0] - post_vs[0] )**2 + (pre_vs[1] - post_vs[1] )**2)
            if temp < distance:
                distance = temp
        
    return distance

def minimum_distance_between_4vertices_with_shift(pre_4vertices, post_4vertices, shift):
    # calculate distance of all combination of 4 vertices of pre and post tiles and judge
    distance = sys.float_info.max
    for pre_vs in pre_4vertices:
        for post_vs in post_4vertices:
            temp = math.sqrt( ( post_vs[0] - pre_vs[0] - shift[0] )**2 + ( post_vs[1] - pre_vs[1] - shift[1] )**2)
            if temp < distance:
                distance = temp
        
    return distance
