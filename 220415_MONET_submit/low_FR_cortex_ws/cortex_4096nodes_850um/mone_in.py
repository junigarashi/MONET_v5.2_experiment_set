#!/bin/env python
import json
#import numpy as np
import math
import time
import subprocess
import os

# read system parameter
def read_system_parameters():
    f = open("system.json", "r")

    #sd: system data
    sd_candidates = json.load(f)
    f.close()
    sd=[]

    print("mode:", sd_candidates["selected_mode"])
    
    # "cortical_sheet_type": "artificial_square"
    if sd_candidates["selected_mode"] == "artificial_square":
        sd = sd_candidates["artificial_square"]

        # count numbers of items
        sd["n_process"] = sd["x_points"]*sd["y_points"]
        sd["n_regions"] = sd["n_x_regions"]*sd["n_y_regions"]
        if sd["n_regions"] != len(sd["region_names"]):
            print( "Error: number of region names:" + str(len(sd["region_names"])) + " was smaller than n_x_regions x n_yregions:" + str(sd["n_x_regions"] * sd["n_y_regions"]))
            exit()
    elif sd_candidates["selected_mode"] == "minimum_artificial_square":
        sd = sd_candidates["minimum_artificial_square"]

        # count numbers of items
        sd["n_process"] = sd["x_points"]*sd["y_points"]
        sd["n_regions"] = sd["n_x_regions"]*sd["n_y_regions"]
        if sd["n_regions"] != len(sd["region_names"]):
            print( "Error: number of region names:" + str(len(sd["region_names"])) + " was smaller than n_x_regions x n_yregions:" + str(sd["n_x_regions"] * sd["n_y_regions"]))
            exit()
            
    elif sd_candidates["selected_mode"] == "example_minimum_thalamic_neuron":
        sd = sd_candidates["example_minimum_thalamic_neuron"]
        
        # count numbers of items
        sd["n_process"] = sd["x_points"]*sd["y_points"]
        sd["n_regions"] = sd["n_x_regions"]*sd["n_y_regions"]
        if sd["n_regions"] != len(sd["region_names"]):
            print( "Error: number of region names:" + str(len(sd["region_names"])) + " was smaller than n_x_regions x n_yregions:" + str(sd["n_x_regions"] * sd["n_y_regions"]))
            exit()

   # human cortical surface
    elif sd_candidates["selected_mode"] == "human_cortical_sheet":
        sd = sd_candidates["human_cortical_sheet"]
        if sd["n_regions"] != len(sd["region_names"]):
            print( "Error: number of region names:" + str(len(sd["region_name"])) + " was smaller than n_region_human", str(sd["n_region"] ))
            exit()
    
    # corticothalamic_circuit_for_cerebellum
    elif sd_candidates["selected_mode"] == "corticothalamic_circuit_for_cerebellum":
        sd = sd_candidates["corticothalamic_circuit_for_cerebellum"]
        
        if sd["n_regions"] != len(sd["region_names"]):
            print( "Error: number of region names:" + str(len(sd["region_name"])) + " was smaller than n_region_human", str(sd["n_region"] ))
            exit()

    # cerebellum 
    elif sd_candidates["selected_mode"] == "cerebellum":
        sd = sd_candidates["cerebellum"]
        
        # count numbers of items
        sd["n_process"] = sd["x_points"]*sd["y_points"]
        sd["n_regions"] = sd["n_x_regions"]*sd["n_y_regions"]
        if sd["n_regions"] != len(sd["region_names"]):
            print( "Error: number of region names:" + str(len(sd["region_names"])) + " was smaller than n_x_regions x n_yregions:" + str(sd["n_x_regions"] * sd["n_y_regions"]))
            exit()
            
    # PFCFP test
    elif sd_candidates["selected_mode"] == "PFCFP_test":
        sd = sd_candidates["PFCFP_test"]
        
        # count numbers of items
        sd["n_process"] = sd["x_points"]*sd["y_points"]
        sd["n_regions"] = sd["n_x_regions"]*sd["n_y_regions"]
        if sd["n_regions"] != len(sd["region_names"]):
            print( "Error: number of region names:" + str(len(sd["region_names"])) + " was smaller than n_x_regions x n_yregions:" + str(sd["n_x_regions"] * sd["n_y_regions"]))
            exit()
    
    # minimum test
    elif sd_candidates["selected_mode"] == "minimum_test":
        sd = sd_candidates["minimum_test"]

    elif sd_candidates["selected_mode"] == "example_colocalization_channels":
        sd = sd_candidates["example_colocalization_channels"]
        
    elif sd_candidates["selected_mode"] == "example_STDP":
        sd = sd_candidates["example_STDP"]
        
    elif sd_candidates["selected_mode"] == "example_HTC":
        sd = sd_candidates["example_HTC"]
        
    # cortico_thalamo_cerebello_sheets
    elif sd_candidates["selected_mode"] == "cortico_thalamo_cerebello_sheets":
        sd = sd_candidates["cortico_thalamo_cerebello_sheets"]
        sd["n_process"] = 0
        sd["n_regions"] = 0
        sd["region_names"] = []
        sd["x_points_array"] = []
        sd["y_points_array"] = []
        
        # count numbers of items
        for i_mrn, mrn in enumerate(sd["meta_region_names"]):
            sd[mrn]["n_process"] = sd[mrn]["x_points"] * sd[mrn]["y_points"]
            sd[mrn]["n_regions"] = sd[mrn]["n_x_regions"] * sd[mrn]["n_y_regions"]
            sd["n_process"] += sd[mrn]["x_points"] * sd[mrn]["y_points"]
            sd["n_regions"] += sd[mrn]["n_x_regions"]*sd[mrn]["n_y_regions"]
            sd["region_names"] += sd[mrn]["region_names"]
            
            for i_rn, rn in enumerate(sd[mrn]["region_names"]):
                sd["x_points_array"].append( (sd[mrn]["x_points"] / sd[mrn]["n_x_regions"]) )
                sd["y_points_array"].append( (sd[mrn]["y_points"] / sd[mrn]["n_y_regions"]) )

    # "cortical_sheet_type": "artificial_square"
    elif sd_candidates["selected_mode"] == "DTI_test":
        sd = sd_candidates["DTI_test"]

        # count numbers of items
        sd["n_process"] = sd["x_points"]*sd["y_points"]
        sd["n_regions"] = sd["n_x_regions"]*sd["n_y_regions"]
        if sd["n_regions"] != len(sd["region_names"]):
            print( "Error: number of region names:" + str(len(sd["region_names"])) + " was smaller than n_x_regions x n_yregions:" + str(sd["n_x_regions"] * sd["n_y_regions"]))
            exit()

    elif sd_candidates["selected_mode"] == "cluster_test":
        sd = sd_candidates["cluster_test"]

        # count numbers of items
        sd["n_process"] = sd["x_points"]*sd["y_points"]
        sd["n_regions"] = sd["n_x_regions"]*sd["n_y_regions"]
        if sd["n_regions"] != len(sd["region_names"]):
            print( "Error: number of region names:" + str(len(sd["region_names"])) + " was smaller than n_x_regions x n_yregions:" + str(sd["n_x_regions"] * sd["n_y_regions"]))
            exit()
            
    elif sd_candidates["selected_mode"] == "ID_connect_example":
        sd = sd_candidates["ID_connect_example"]

        # count numbers of items
        sd["n_process"] = sd["x_points"]*sd["y_points"]
        sd["n_regions"] = sd["n_x_regions"]*sd["n_y_regions"]
        if sd["n_regions"] != len(sd["region_names"]):
            print( "Error: number of region names:" + str(len(sd["region_names"])) + " was smaller than n_x_regions x n_yregions:" + str(sd["n_x_regions"] * sd["n_y_regions"]))
            exit()
            
    else:
        print( "Can't find suitable mode for input mode name:", sd["selected_mode"] )
        exit()

    # set region GID
    sd["region_name_to_region_GID"]={}
    for i_rn, rn in enumerate(sd["region_names"]):
        sd["region_name_to_region_GID"][rn]=i_rn

    # PRNG seed preparation
    if isinstance(sd["PRNG_seed"], str):
        # case to use time
        if sd["PRNG_seed"] == "time":
            sd["PRNG_seed"] = str(int(time.time()))
        
        else:
            print( "Error: Invalid value in PRNG_seed, ", sd["PRNG_seed"])
            exit()

    return sd
    
# read region files
def read_region_parameters(sd):
    rd_array = [ [] for rn in sd["region_names"] ]
    
    for i_rn, rn in enumerate(sd["region_names"]):
        filename = "./mode_setting/" + sd["mode"] + "/region_" + rn + ".json"
        print( "read", filename)
        f = open(filename, "r")

        # rd: region data
        rd = json.load(f)
        f.close()

        # prepare structure information
        # subregion info
        rd["structure_info"]["n_subregions"] = len(rd["neuron_info"][rn])
        
        rd["structure_info"]["subregion_name_to_subregion_LID"]={}
        for i_sr, sr in enumerate(rd["neuron_info"][rn]):
            sr_name = str(*sr.keys())
            rd["structure_info"]["subregion_name_to_subregion_LID"][sr_name] =  i_sr
            
        rd["structure_info"]["n_neuron_types_per_subregion"] = [ len(sr[str(*sr.keys())]) for sr in rd["neuron_info"][rn] ]
        
        # conversion table of neuron_type_LID_to_subregion_LID
        rd["structure_info"]["neuron_type_SerLID_to_subregion_LID"] = []
        for i_nntps, nntps in enumerate(rd["structure_info"]["n_neuron_types_per_subregion"]):
            for j in range(nntps):
                rd["structure_info"]["neuron_type_SerLID_to_subregion_LID"].append(i_nntps)

        # converison table of LID to SerLID
        count=0
        rd["structure_info"]["LID_to_neuron_type_SerLID"]=[]
        for i_sr, sr in enumerate(rd["structure_info"]["n_neuron_types_per_subregion"]):
            list_sr = list(range(sr))
            rd["structure_info"]["LID_to_neuron_type_SerLID"].append(list_sr)
            for x in range(sr):
                rd["structure_info"]["LID_to_neuron_type_SerLID"][i_sr][x]=count
                count+=1
        
        # neuron type info
        rd["structure_info"]["n_neuron_types_per_process"] = sum(rd["structure_info"]["n_neuron_types_per_subregion"])
        
        n_neuron_types_per_region = 0
        rd["structure_info"]["neuron_type_names"] = []
        for sr in rd["neuron_info"][rn]:
            sn=str(*sr.keys())
            temp=[]
            for nts in sr[str(*sr.keys())]:
                temp.append( str(*nts.keys()))
                n_neuron_types_per_region += 1
            rd["structure_info"]["neuron_type_names"].append(temp)

        rd["structure_info"]["n_neuron_types_per_region"] = n_neuron_types_per_region

        rd["structure_info"]["neuron_type_name_to_neuron_type_LID"] = {}
        for i_sn, ntns  in enumerate(rd["structure_info"]["neuron_type_names"]):
            temp={}
            for i_ntn, ntn in enumerate(ntns):
                temp[ntn]=i_ntn
            sn = str(*rd["neuron_info"][rn][i_sn].keys())
            rd["structure_info"]["neuron_type_name_to_neuron_type_LID"][sn] = temp
        
        i_SerLID = 0
        rd["structure_info"]["neuron_type_name_to_neuron_type_SerLID"] = {}
        for i_sn, ntns  in enumerate(rd["structure_info"]["neuron_type_names"]):
            temp={}
            for i_ntn, ntn in enumerate(ntns):
                temp[ntn] = i_SerLID
                i_SerLID+=1
            sn = str(*rd["neuron_info"][rn][i_sn].keys())
            rd["structure_info"]["neuron_type_name_to_neuron_type_SerLID"][sn] = temp
            
        # conversion table of neuron_type_SerLID to neuron_type_LID
        rd["structure_info"]["neuron_type_SerLID_to_neuron_type_LID"] = []
        for i, dumy in enumerate(rd["structure_info"]["neuron_type_names"]):
            for j, dumy2 in enumerate(dumy):
                rd["structure_info"]["neuron_type_SerLID_to_neuron_type_LID"].append(j)
        
        # set n_neuron per neuron type
        rd["structure_info"]["n_neurons_per_neuron_type"] = []
        for i_sr, sr in enumerate( rd["neuron_info"][rn]):
            sn=str(*sr.keys())
            for i_nt, ntn in enumerate(rd["structure_info"]["neuron_type_names"][i_sr]):
                
                rd["structure_info"]["n_neurons_per_neuron_type"].append(int(rd["neuron_info"][rn][i_sr][sn][i_nt][ntn]["n_cells_per_mm2"]* rd["structure_info"]["xy_length_per_tile"] / 1000. * rd["structure_info"]["xy_length_per_tile"] / 1000. * rd["structure_info"]["n_neurons_factor"]))
                
        # check the number of neuron
        n_neurons_per_process = sum(rd["structure_info"]["n_neurons_per_neuron_type"])

        # set neuron name array
        neuron_type_SerLID=0
        rd["structure_info"]["neuron_type_name_array"] = []
        rd["structure_info"]["full_neuron_name_array"] = []
        for i_sr, sr in enumerate( rd["neuron_info"][rn]):
            for i_nt, nt in enumerate(sr[str(*sr.keys())]):
                rd["structure_info"]["neuron_type_name_array"].append(str(*nt.keys()))
                rd["structure_info"]["full_neuron_name_array"].append(str(rn) + "_" + str(*sr.keys()) + "_" + str(*nt.keys()))
                neuron_type_SerLID+=1

        # set neuron_model array
        neuron_type_SerLID=0
        rd["structure_info"]["neuron_model_array"] = []
        for i_sr, sr in enumerate( rd["neuron_info"][rn]):
            temp={}
            for i_nt, nt in enumerate(sr[str(*sr.keys())]):
                temp[str(neuron_type_SerLID)]=str(nt[str(*nt.keys())]["neuron_model"])
                neuron_type_SerLID+=1
            rd["structure_info"]["neuron_model_array"].append([ n_cell for (rnt, n_cell) in sorted(temp.items(), key=lambda x:x[0])])

        # set E_or_I array
        neuron_type_SerLID=0
        rd["structure_info"]["E_or_I_array"] = []
        for i_sr, sr in enumerate( rd["neuron_info"][rn]):
            temp={}
            for i_nt, nt in enumerate(sr[str(*sr.keys())]):
                temp[str(neuron_type_SerLID)]=str(nt[str(*nt.keys())]["EorI"])
                neuron_type_SerLID+=1
            rd["structure_info"]["E_or_I_array"].append([ n_cell for (rnt, n_cell) in sorted(temp.items(), key=lambda x:x[0])])
            
        # set membrane_time_constant_array
        neuron_type_SerLID=0
        rd["structure_info"]["membrane_time_constant_array"] = []
        for i_sr, sr in enumerate( rd["neuron_info"][rn]):
            temp={}
            for i_nt, nt in enumerate(sr[str(*sr.keys())]):
                temp[str(neuron_type_SerLID)]=str(nt[str(*nt.keys())]["membrane_time_constant"])
                neuron_type_SerLID+=1            
            rd["structure_info"]["membrane_time_constant_array"].append([ n_cell for (rnt, n_cell) in sorted(temp.items(), key=lambda x:x[0])])
        
        # set spike_threshold_array
        neuron_type_SerLID=0
        rd["structure_info"]["spike_threshold_array"] = []
        for i_sr, sr in enumerate( rd["neuron_info"][rn]):
            temp={}
            for i_nt, nt in enumerate(sr[str(*sr.keys())]):
                temp[str(neuron_type_SerLID)]=str(nt[str(*nt.keys())]["spike_threshold"])
                neuron_type_SerLID+=1
            rd["structure_info"]["spike_threshold_array"].append([ n_cell for (rnt, n_cell) in sorted(temp.items(), key=lambda x:x[0])])
        
        # set reset_value_array
        neuron_type_SerLID=0
        rd["structure_info"]["reset_value_array"] = []
        for i_sr, sr in enumerate( rd["neuron_info"][rn]):
            temp={}
            for i_nt, nt in enumerate(sr[str(*sr.keys())]):
                temp[str(neuron_type_SerLID)]=str(nt[str(*nt.keys())]["reset_value"])
                neuron_type_SerLID+=1
            rd["structure_info"]["reset_value_array"].append([ n_cell for (rnt, n_cell) in sorted(temp.items(), key=lambda x:x[0])])
        
        # set E_rest_array
        neuron_type_SerLID=0
        rd["structure_info"]["E_rest_array"] = []
        for i_sr, sr in enumerate( rd["neuron_info"][rn]):
            temp={}
            for i_nt, nt in enumerate(sr[str(*sr.keys())]):
                temp[str(neuron_type_SerLID)]=str(nt[str(*nt.keys())]["E_rest"])
                neuron_type_SerLID+=1            
            rd["structure_info"]["E_rest_array"].append([ n_cell for (rnt, n_cell) in sorted(temp.items(), key=lambda x:x[0])])

        # set I_ex_array
        neuron_type_SerLID=0
        rd["structure_info"]["I_ex_array"] = []
        for i_sr, sr in enumerate( rd["neuron_info"][rn]):
            temp={}
            for i_nt, nt in enumerate(sr[str(*sr.keys())]):
                temp[str(neuron_type_SerLID)]=str(nt[str(*nt.keys())]["I_ex"])
                neuron_type_SerLID+=1            
            rd["structure_info"]["I_ex_array"].append([ n_cell for (rnt, n_cell) in sorted(temp.items(), key=lambda x:x[0])])
    
        # set absolute_refractory_period_array
        neuron_type_SerLID=0
        rd["structure_info"]["absolute_refractory_period_array"] = []
        for i_sr, sr in enumerate( rd["neuron_info"][rn]):
            temp={}
            for i_nt, nt in enumerate(sr[str(*sr.keys())]):
                temp[str(neuron_type_SerLID)]=str(nt[str(*nt.keys())]["absolute_refractory_period"])
                neuron_type_SerLID+=1
            rd["structure_info"]["absolute_refractory_period_array"].append([ n_cell for (rnt, n_cell) in sorted(temp.items(), key=lambda x:x[0])])

        # check the number of neurons whether the numbers are not zero

        # estimate numbers of connection
        #estimate_n_connections(rd_array, sd, i_rn)
        
        # set region directory
        rd_array[i_rn] = rd
                    
    return rd_array


def estimate_n_connections(rd_array, sd):

    # region loop 
    for i_rd, rd in enumerate(rd_array):

        intra_sum_per_pre = [ 0 for nt in range(rd["structure_info"]["n_neuron_types_per_region"]) ]
        intra_sum_per_post = [ 0 for nt in range(rd["structure_info"]["n_neuron_types_per_region"]) ]

        inter_sum_per_post = [ 0 for nt in range(rd["structure_info"]["n_neuron_types_per_region"]) ]
        inter_sum_per_pre = [ [] for bundle in rd["inter_regional_connection"] ]
        for i_bundle, bundle in enumerate(rd["inter_regional_connection"]):
            inter_sum_per_pre[i_bundle] = [ 0 for ct in bundle["connection_type"] ]
        
        # intra regional connection
        for i_post, conns in enumerate(rd["intra_regional_connection"]["connection_parameter_matrix"]):
            for i_pre, conn in enumerate(conns):
                # two dimensional gaussian
                if str(*conn["spatial_pattern"].keys()) == "two_dimensional_gaussian":
                    temp_conn_prob = 0.
                    temp_conn_sd = 0.
                
                    if isinstance(conn["spatial_pattern"]["two_dimensional_gaussian"]["mu"], str):
                        temp_conn_prob = rd["intra_regional_connection"]["coefficients"][ conn["spatial_pattern"]["two_dimensional_gaussian"]["mu"] ]
                        temp_conn_sd = rd["intra_regional_connection"]["coefficients"][ conn["spatial_pattern"]["two_dimensional_gaussian"]["sd"] ]
                    else:
                        temp_conn_prob = conn["spatial_pattern"]["two_dimensional_gaussian"]["mu"] 
                        temp_conn_sd = conn["spatial_pattern"]["two_dimensional_gaussian"]["sd"]
                
                    if temp_conn_prob > 0.:
                        pre_neuron_num = rd["structure_info"]["n_neurons_per_neuron_type"][i_pre]
                        post_neuron_num = rd["structure_info"]["n_neurons_per_neuron_type"][i_post]
                        full_num_conn = pre_neuron_num * post_neuron_num
                        gf_conn_area = math.sqrt(math.pi * (2*temp_conn_sd**2)) * temp_conn_prob
                        full_conn_area = rd["structure_info"]["xy_length_per_tile"] * 1.

                        intra_sum_per_pre[i_pre] += full_num_conn * gf_conn_area / full_conn_area
                        intra_sum_per_post[i_post] += full_num_conn * gf_conn_area / full_conn_area

                elif str(*conn["spatial_pattern"].keys()) == "orthogornal_cross":
                    if conn["spatial_pattern"]["orthogornal_cross"]["probability"] > 0.:
                        pre_neuron_num = rd["structure_info"]["n_neurons_per_neuron_type"][i_pre]
                        post_neuron_num = rd["structure_info"]["n_neurons_per_neuron_type"][i_post]
                        area_fraction = 2 * conn["spatial_pattern"]["orthogornal_cross"]["pre_width"] / float(rd["structure_info"]["xy_length_per_tile"]) * 2 * conn["spatial_pattern"]["orthogornal_cross"]["post_width"] / float(rd["structure_info"]["xy_length_per_tile"]) * 1.2
                        intra_sum_per_pre[i_pre] += area_fraction * pre_neuron_num * post_neuron_num * conn["spatial_pattern"]["orthogornal_cross"]["probability"] 
                        intra_sum_per_post[i_post] += area_fraction * pre_neuron_num * post_neuron_num * conn["spatial_pattern"]["orthogornal_cross"]["probability"]
                    
                elif str(*conn["spatial_pattern"].keys()) == "circular":
                    print( str(*conn["spatial_pattern"].keys()), "not yet")
                elif str(*conn["spatial_pattern"].keys()) == "square":
                    print( str(*conn["spatial_pattern"].keys()), "not yet")
                elif str(*conn["spatial_pattern"].keys()) == "ID_connect":
                    intra_sum_per_pre[i_pre] += len(conn["spatial_pattern"]["ID_connect"]["neuron_ID_pairs"])
                    intra_sum_per_pre[i_post] += len(conn["spatial_pattern"]["ID_connect"]["neuron_ID_pairs"])
                else:
                    print( str(*conn["spatial_pattern"].keys()), "is not supported" )
                    
        # intra ID connection
        if "ID_connection" in rd["intra_regional_connection"].keys():
            for ID_conn in rd["intra_regional_connection"]["ID_connection"]:

                # read numbers of neuron ID pairs from file
                neuron_ID_file_name = ID_conn["spatial_pattern"]["ID_connect"]["neuron_ID_file_name"]
                file_names = os.listdir(path='.')
                if neuron_ID_file_name in file_names:
                    line_count = int(subprocess.check_output(['wc', '-l', neuron_ID_file_name]).decode().split(' ')[0])
                    # underconstruction
                    os.exit()
                    
                else:
                    print("There is not", neuron_ID_file_name, "at current directory")

                # set from region.json file
                if ID_conn["pair_IDs"] != []:
                    n_ID_pairs = len(ID_conn["pair_IDs"])
                    pre_neuron_type_SerLID = rd["structure_info"]["neuron_type_name_to_neuron_type_SerLID"][ID_conn["pre"]["subregion"]][ID_conn["pre"]["neuron_type"]]
                    post_neuron_type_SerLID = rd["structure_info"]["neuron_type_name_to_neuron_type_SerLID"][ID_conn["post"]["subregion"]][ID_conn["post"]["neuron_type"]]

                    intra_sum_per_pre[pre_neuron_type_SerLID] += 1
                    intra_sum_per_post[post_neuron_type_SerLID] += 1
        
        # inter regional connection
        for i_bundle, bundle in enumerate(rd["inter_regional_connection"]):
            
            for i_ct, ct in enumerate(bundle["connection_type"]):
                # two dimensional gaussian
                if str(*ct["spatial_pattern"].keys()) == "two_dimensional_gaussian":
                    temp_conn_prob = 0.
                    temp_conn_sd = 0.
                    
                    temp_conn_prob = ct["spatial_pattern"]["two_dimensional_gaussian"]["mu"] 
                    temp_conn_sd = ct["spatial_pattern"]["two_dimensional_gaussian"]["sd"]
                
                    if temp_conn_prob > 0.:
                        pre_region_GID = sd["region_name_to_region_GID"][ct["pre"]["region"]]
                        post_region_GID = sd["region_name_to_region_GID"][ct["post"]["region"]]

                        pre_neuron_type_SerLID = rd_array[pre_region_GID]["structure_info"]["neuron_type_name_to_neuron_type_SerLID"][ct["pre"]["subregion"]][ct["pre"]["neuron_type"]]
                        post_neuron_type_SerLID = rd_array[post_region_GID]["structure_info"]["neuron_type_name_to_neuron_type_SerLID"][ct["post"]["subregion"]][ct["post"]["neuron_type"]]

                        pre_neuron_num = rd_array[pre_region_GID]["structure_info"]["n_neurons_per_neuron_type"][pre_neuron_type_SerLID]
                        post_neuron_num = rd_array[post_region_GID]["structure_info"]["n_neurons_per_neuron_type"][post_neuron_type_SerLID]
                        full_num_conn = pre_neuron_num * post_neuron_num
                        gf_conn_area = math.sqrt(math.pi * (2*temp_conn_sd**2)) * temp_conn_prob
                        full_conn_area = rd_array[post_region_GID]["structure_info"]["xy_length_per_tile"] * 1.

                        # need to be finalize
                        inter_sum_per_post[post_neuron_type_SerLID] += full_num_conn * gf_conn_area / full_conn_area
                        inter_sum_per_pre[i_bundle][i_ct] += full_num_conn * gf_conn_area / full_conn_area
                    
                elif str(*ct["spatial_pattern"].keys()) == "circular":
                    print( str(*ct["spatial_pattern"].keys()), "not yet")
                    os.exit()
                elif str(*conn["spatial_pattern"].keys()) == "square":
                    print( str(*ct["spatial_pattern"].keys()), "not yet")
                    os.exit()
                    
            """
            # inter ID connection
            print(rd["inter_regional_connection"])
            if "ID_connection" in bundle["inter_regional_connection"].keys():
                for ID_conn in rd["inter_regional_connection"]["ID_connection"]:

                    # read numbers of neuron ID pairs from file
                    neuron_ID_file_name = ID_conn["spatial_pattern"]["ID_connect"]["neuron_ID_file_name"]
                    file_names = os.listdir(path='.')
                    if neuron_ID_file_name in file_names:
                        line_count = int(subprocess.check_output(['wc', '-l', neuron_ID_file_name]).decode().split(' ')[0])
                        print( "neuron_ID_file_name is not yet implemented")
                        # underconstruction
                        os.exit()
                    else:
                        print("There is not", neuron_ID_file_name, "at current directory")

                        # set from region.json file
                    if ID_conn["pair_IDs"] != []:
                        n_ID_pairs = len(ID_conn["pair_IDs"])
                        pre_region_GID = sd["region_name_to_region_GID"][ct["pre"]["region"]]
                        post_region_GID = sd["region_name_to_region_GID"][ct["post"]["region"]]
                       
                        pre_neuron_type_SerLID = rd_array[pre_region_GID]["structure_info"]["neuron_type_name_to_neuron_type_SerLID"][ID_conn["pre"]["subregion"]][ID_conn["pre"]["neuron_type"]]
                        post_neuron_type_SerLID = rd_array[post_region_GID]["structure_info"]["neuron_type_name_to_neuron_type_SerLID"][ID_conn["post"]["subregion"]][ID_conn["post"]["neuron_type"]]
                        
                        inter_sum_per_pre[i_bundle][pre_neuron_type_SerLID] += len(ID_conn["pair_IDs"])
                        inter_sum_per_post[post_neuron_type_SerLID] += len(ID_conn["pair_IDs"])
            """
        
        # num of intra outdegree connections of presynaptic neuron type
        rd_array[i_rd]["structure_info"]["n_estimated_intra_outdegree_connections_per_presynaptic_neuron_type"] = [ 0 for nt in range(rd_array[i_rd]["structure_info"]["n_neuron_types_per_region"]) ]
        rd_array[i_rd]["structure_info"]["n_estimated_intra_outdegree_connections_per_one_presynaptic_neuron"] =  [ 0 for nt in range(rd_array[i_rd]["structure_info"]["n_neuron_types_per_region"]) ]
        for i_pre, intra_sum_pre in enumerate(intra_sum_per_pre):
            rd_array[i_rd]["structure_info"]["n_estimated_intra_outdegree_connections_per_presynaptic_neuron_type"][i_pre] = math.ceil(intra_sum_per_pre[i_pre]) + 10
            rd_array[i_rd]["structure_info"]["n_estimated_intra_outdegree_connections_per_one_presynaptic_neuron"][i_pre] = math.ceil(intra_sum_per_pre[i_pre] / float(rd_array[i_rd]["structure_info"]["n_neurons_per_neuron_type"][i_pre]))+10

        # total intra outdegree connections
        rd_array[i_rd]["structure_info"]["n_estimated_intra_outdegree_connections"] = sum(rd_array[i_rd]["structure_info"]["n_estimated_intra_outdegree_connections_per_presynaptic_neuron_type"])
    
        # num of intra indegree connections of postsynaptic neuron type
        rd_array[i_rd]["structure_info"]["n_estimated_intra_indegree_connections_per_postsynaptic_neuron_type"] = [ 0 for nt in range(rd_array[i_rd]["structure_info"]["n_neuron_types_per_region"]) ]
        rd_array[i_rd]["structure_info"]["n_estimated_intra_indegree_connections_per_one_postsynaptic_neuron"] =  [ 0 for nt in range(rd_array[i_rd]["structure_info"]["n_neuron_types_per_region"]) ]    
        for i_post, intra_sum_post in enumerate(intra_sum_per_post):
            rd_array[i_rd]["structure_info"]["n_estimated_intra_indegree_connections_per_postsynaptic_neuron_type"][i_post] = math.ceil(intra_sum_per_post[i_post]) + 10
            rd_array[i_rd]["structure_info"]["n_estimated_intra_indegree_connections_per_one_postsynaptic_neuron"][i_post] = math.ceil(intra_sum_per_post[i_post] / float(rd_array[i_rd]["structure_info"]["n_neurons_per_neuron_type"][i_post]))+10

        # total intra indegree connections
        rd_array[i_rd]["structure_info"]["n_estimated_intra_indegree_connections"] = sum(rd_array[i_rd]["structure_info"]["n_estimated_intra_indegree_connections_per_postsynaptic_neuron_type"])

        # num of inter outdegree connections of postsynaptic neuron type
        #pre_region_GID = sd["region_name_to_region_GID"][ct["pre"]["region"]]
        rd_array[i_rd]["structure_info"]["n_estimated_inter_outdegree_connections_per_presynaptic_neuron_type"] = [ [] for bundle in rd_array[i_rd]["inter_regional_connection"] ]
        rd_array[i_rd]["structure_info"]["n_estimated_inter_outdegree_connections_per_one_presynaptic_neuron"] =  [ [] for bundle in rd_array[i_rd]["inter_regional_connection"] ]

        # bundle loop
        for i_b, bundle in enumerate(rd["inter_regional_connection"]):
            pre_region_GID = sd["region_name_to_region_GID"][bundle["pre"]["region"]]
            post_region_GID = sd["region_name_to_region_GID"][bundle["post"]["region"]]
            bundle_name = "bundle_" + bundle["pre"]["region"] + "_to_" + bundle["post"]["region"]
            pre_post_nt_pair_per_pre_nt_at_bundle = [ [] for i_pns in range(rd_array[pre_region_GID]["structure_info"]["n_neuron_types_per_process"])]
            pre_nt_SID_to_pre_nt_LID_at_bundle = [ 0 for i_npns in range(rd_array[pre_region_GID]["structure_info"]["n_neuron_types_per_process"])]

            # make pre post conn array per pre nt 
            for i_ct, conn in enumerate(bundle["connection_type"]):
                pre_sr_LID = rd_array[pre_region_GID]["structure_info"]["subregion_name_to_subregion_LID"][conn["pre"]["subregion"]] 
                pre_nt_LID = rd_array[pre_region_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][conn["pre"]["subregion"]][conn["pre"]["neuron_type"]]
                pre_nt_SerLID = rd_array[pre_region_GID]["structure_info"]["LID_to_neuron_type_SerLID"][pre_sr_LID][pre_nt_LID]
                post_sr_LID = rd_array[post_region_GID]["structure_info"]["subregion_name_to_subregion_LID"][conn["post"]["subregion"]] 
                post_nt_LID = rd_array[post_region_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][conn["post"]["subregion"]][conn["post"]["neuron_type"]]
                post_nt_SerLID = rd_array[post_region_GID]["structure_info"]["LID_to_neuron_type_SerLID"][post_sr_LID][post_nt_LID]
                pre_post_nt_pair_per_pre_nt_at_bundle[pre_nt_SerLID].append([pre_sr_LID, pre_nt_LID, pre_nt_SerLID, post_sr_LID, post_nt_LID, post_nt_SerLID, conn])
            #print(pre_post_nt_pair_per_pre_nt_at_bundle)
            
            # count pre nt with connections and get the SIDs
            count_pre_nt_LID_at_bundle = 0
            pre_nt_SIDs_at_bundle = []
            for i_ppntpg, pre_post_nt_pair_group in enumerate(pre_post_nt_pair_per_pre_nt_at_bundle):
                if pre_post_nt_pair_group != []:
                    pre_nt_SID_to_pre_nt_LID_at_bundle[i_ppntpg] = count_pre_nt_LID_at_bundle
                    pre_nt_LID_at_bundle = count_pre_nt_LID_at_bundle
                    pre_nt_SIDs_at_bundle.append(i_ppntpg)
                    count_pre_nt_LID_at_bundle += 1
            n_pre_nt_at_bundle = count_pre_nt_LID_at_bundle
            #print("n_pre_nt_at_bundle", n_pre_nt_at_bundle)
            
            # make data array of inter outdgeree connections per pre nt
            rd_array[i_rd]["structure_info"]["n_estimated_inter_outdegree_connections_per_presynaptic_neuron_type"][i_b] = [ 0 for i_pre in range(n_pre_nt_at_bundle) ]

            #print(pre_post_nt_pair_per_pre_nt_at_bundle)
            # sum estimated connections per pre nt
            for i_ppntpg, pre_post_nt_pair_group in enumerate(pre_post_nt_pair_per_pre_nt_at_bundle):
                #print(i_ppntpg, pre_post_nt_pair_group)
                for i_ppntp, pre_post_nt_pair in enumerate(pre_post_nt_pair_group):
                    if pre_post_nt_pair != []:
                        temp_conn_prob = pre_post_nt_pair[6]["spatial_pattern"]["two_dimensional_gaussian"]["mu"] 
                        temp_conn_sd = pre_post_nt_pair[6]["spatial_pattern"]["two_dimensional_gaussian"]["sd"]
                        if temp_conn_prob > 0.:
                            pre_neuron_num = rd_array[pre_region_GID]["structure_info"]["n_neurons_per_neuron_type"][pre_post_nt_pair[2] ]
                            post_neuron_num = rd_array[post_region_GID]["structure_info"]["n_neurons_per_neuron_type"][pre_post_nt_pair[5] ]
                            full_num_conn = pre_neuron_num * post_neuron_num
                            gf_conn_area = math.sqrt(math.pi * (2*temp_conn_sd**2)) * temp_conn_prob
                            full_conn_area = rd_array[post_region_GID]["structure_info"]["xy_length_per_tile"] * 1.
                            rd_array[i_rd]["structure_info"]["n_estimated_inter_outdegree_connections_per_presynaptic_neuron_type"][i_b][pre_nt_SID_to_pre_nt_LID_at_bundle[i_ppntpg]] += full_num_conn * gf_conn_area / full_conn_area

        
        # num of inter indegree connections of postsynaptic neuron type
        rd_array[i_rd]["structure_info"]["n_estimated_inter_indegree_connections_per_postsynaptic_neuron_type"] = [ 0 for nt in range(rd_array[i_rd]["structure_info"]["n_neuron_types_per_region"]) ]
        rd_array[i_rd]["structure_info"]["n_estimated_inter_indegree_connections_per_one_postsynaptic_neuron"] =  [ 0 for nt in range(rd_array[i_rd]["structure_info"]["n_neuron_types_per_region"]) ]    
        for i_post, inter_sum_post in enumerate(inter_sum_per_post):
            rd_array[i_rd]["structure_info"]["n_estimated_inter_indegree_connections_per_postsynaptic_neuron_type"][i_post] = math.ceil(inter_sum_per_post[i_post]) + 10
            rd_array[i_rd]["structure_info"]["n_estimated_inter_indegree_connections_per_one_postsynaptic_neuron"][i_post] = math.ceil(inter_sum_per_post[i_post] / float(rd_array[i_rd]["structure_info"]["n_neurons_per_neuron_type"][i_post])) + 10

        # total inter outdegree connections
        rd_array[i_rd]["structure_info"]["n_estimated_inter_indegree_connections"] = sum(rd_array[i_rd]["structure_info"]["n_estimated_inter_indegree_connections_per_postsynaptic_neuron_type"])

        
        # one connection consume 18 B in simulation
        memory_consumption = rd_array[i_rd]["structure_info"]["n_estimated_intra_indegree_connections"] * 18

        
        if sd["memory_per_node[GB]"] * 10**9 < memory_consumption:
            print( "#############################################################################" )
            print( "# Warning: estimation of memory consumption exceeded limit of system memory" )
            print( "# System memory (compute node):", sd["memory_per_node[GB]"], " GB" )
            print( "# Estimated memory consumption:", memory_consumption / (10**9), " GB" )
            print( "# Num neurons:", n_neurons_per_process, "Estimated Num connections:", n_neurons_per_process * n_connections_per_neuron )
            print( "# (assumption:", n_connections_per_neuron, "connections per neuron)" )
            print( "# You should reduce the area per tile" )
            print( "#############################################################################" )
            print( "\n" )

        else:
            print( memory_consumption / (float)(10**9), "GB per compute node will be consumed for connections in simulation" )
        
    
