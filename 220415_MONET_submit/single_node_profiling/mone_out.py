#!/bin/env python
import sys
import math
import time, os

def write_rank_node_map(sd, process):
    rank_node_map = ''
    for i_proc in range(sd["n_process"]):
        rank_node_map += "(" + str(process[i_proc].TOFUD_2D_coordinate[0]) + ", " + str(process[i_proc].TOFUD_2D_coordinate[1]) + ")\n" 

    with open('./rank_node_map.dat', mode='w') as f:
        f.write(rank_node_map)
    f.close()

def write_process_configuration_files(sd, rd_array, process, process_IDs):

    # convert table
    region_subregion_neuron_type_names_to_neuron_type_GID = {}
    count_nt_GID = 0
    for i_rd, rd in enumerate(rd_array):
        rn = rd["structure_info"]["region_name"]
        region_subregion_neuron_type_names_to_neuron_type_GID[rn] = {}

        for i_sr, sr in enumerate(rd["neuron_info"][rn]):
            srn = str(*sr.keys())
            region_subregion_neuron_type_names_to_neuron_type_GID[rn][srn]={}
            for i_nt, nt in enumerate(rd[ "neuron_info"][rn][i_sr][srn]):
                ntn = str(*nt.keys())
                region_subregion_neuron_type_names_to_neuron_type_GID[rn][srn][ntn] = count_nt_GID
                count_nt_GID+=1

    # preparetaion of gsi strings
    gsi_strings = make_global_shared_information_strings(sd, rd_array)

    # preparation of switch strings
    switch_strings = make_switch_strings(sd)

    # preparation of hardware info strings
    hardware_info_strings = make_hardware_info_strings(sd)
    
    # preparation of region ct
    region_ct_strings = [ [] for rd in rd_array ]
    for i_rd, rd in enumerate(rd_array):
        region_ct_strings[i_rd].extend(make_connection_type_strings(rd))
        
    # start preparation of input data
    #for i_proc in range(sd["n_process"]):
    for i_proc in process_IDs:
        region_GID = process[i_proc].region_GID

        #f = open("./process_"+str(process[i_proc].process_GID)+".dat","w")
        output="""simulation_time, %(simulation_time)d
dt, %(dt)f
memory_per_node, %(memory_per_node[GB])d
""" % sd

        # shared global information that do not depends on process information
        L = []
        L.extend(gsi_strings)


        L.append("PRNG_seed, ")
        L.append(str(sd["PRNG_seed"]))
        L.append("\n")

        L.append("process_GID, ")
        L.append(str(process[i_proc].process_GID))
        L.append("\n")
        
        L.extend(hardware_info_strings)

        L.append("region_GID, ")
        L.append(str(process[i_proc].region_GID))
        L.append("\n")
        L.append("region_name, ")
        L.append(str(process[i_proc].region_name))
        L.append("\n")
        L.append("xy_length_per_tile, " +str(rd_array[region_GID]["structure_info"]["xy_length_per_tile"]))
        L.append("\n")
        
        L.append("region_position")
        for i in process[i_proc].spatial_extent:
            L.append(", ")
            L.append(str(i))
        L.append("\n")

        L.append("n_same_region_process_GIDs, ")
        L.append(str(len(rd_array[process[i_proc].region_GID]["process_ID"])))
        L.append("\n")

        
        L.append("same_region_process_GIDs")
        for i in rd_array[process[i_proc].region_GID]["process_ID"]:
            L.append(", ")
            L.append(str(i))
        L.append("\n")
        

        # structure info
        L.append("n_subregions, ")
        L.append(str(rd_array[region_GID]["structure_info"]["n_subregions"]))
        L.append("\n")

        L.append("subregion_names")
        for i_sr, sr in enumerate(rd_array[region_GID]["neuron_info"][process[i_proc].region_name]):
            L.append(", ")
            L.append(str(*sr.keys()))
        L.append("\n")

        L.append("n_neuron_types_per_subregion")
        for i in rd_array[region_GID]["structure_info"]["n_neuron_types_per_subregion"]:
            L.append(", ")
            L.append(str(i))
        L.append("\n")

        L.append("subregion_positions")
        for i in process[i_proc].subregion_positions:
            for j in i:
                L.append(", ")
                L.append(str(j))
        L.append("\n")

        L.append("neuron_type_names")
        for i in rd_array[region_GID]["structure_info"]["neuron_type_names"]:
            for j in i:
                L.append(", ")
                L.append(str(j))
        L.append("\n")

        L.append("n_neurons_per_neuron_type")
        for i in rd_array[region_GID]["structure_info"]["n_neurons_per_neuron_type"]:
            L.append(", ")
            L.append(str(i))
        L.append("\n")

        # neuron info
        L.append("neuron_model")
        for i in rd_array[region_GID]["structure_info"]["neuron_model_array"]:
            for j in i:
                L.append(", ")
                L.append(str(j))
        L.append("\n")
        
        L.append("E_or_I")
        for i in rd_array[region_GID]["structure_info"]["E_or_I_array"]:
            for j in i:
                L.append(", ")
                L.append(str(j))
        L.append("\n")
         
        L.append("membrane_time_constant")
        for i in rd_array[region_GID]["structure_info"]["membrane_time_constant_array"]:
            for j in i:
                L.append(", ")
                L.append(str(j))
        L.append("\n")
 
        L.append("spike_threshold")
        for i in rd_array[region_GID]["structure_info"]["spike_threshold_array"]:
            for j in i:
                L.append(", ")
                L.append(str(j))
        L.append("\n")
 
        L.append("reset_value")
        for i in rd_array[region_GID]["structure_info"]["reset_value_array"]:
            for j in i:
                L.append(", ")
                L.append(str(j))
        L.append("\n")
            
        L.append("E_rest")
        for i in rd_array[region_GID]["structure_info"]["E_rest_array"]:
            for j in i:
                L.append(", ")
                L.append(str(j))
        L.append("\n")

        """
        L.append("I_ex")
        for i in rd_array[region_GID]["structure_info"]["I_ex_array"]:
            for j in i:
                L.append(", ")
                L.append(str(j))
        L.append("\n")
        """
        
        L.append("absolute_refractory_period")
        for i in rd_array[region_GID]["structure_info"]["absolute_refractory_period_array"]:
            for j in i:
                L.append(", ")
                L.append(str(j))
        L.append("\n")

        """
        # I_ex
        region_name = rd_array[region_GID]["neuron_info"].keys()[0]
        for i in rd_array[region_GID]["neuron_info"][region_name]:
            subregion_name = i.keys()[0]
            for j in i[subregion_name]:
                neuron_type_name = j.keys()[0]
                L.append("I_ex")
                L.append(", ")
                L.append(str(sd["region_name_to_region_GID"][region_name]))
                L.append(", ")
                L.append(str(rd_array[region_GID]["structure_info"]["subregion_name_to_subregion_LID"][subregion_name]))
                L.append(", ")
                L.append(str(rd_array[region_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][subregion_name][neuron_type_name]))
                L.append(", ")
                L.append(j[neuron_type_name]["I_ex"]["initial_setting_type"])

                if j[neuron_type_name]["I_ex"]["initial_setting_type"] == "homogeneous":
                    L.append(", value, ")
                    L.append(str(j[neuron_type_name]["I_ex"]["value"]))

                elif j[neuron_type_name]["I_ex"]["initial_setting_type"] == "gaussian":
                    L.append(", mean, ")
                    L.append(str(j[neuron_type_name]["I_ex"]["mean"]))
                    L.append(", sd, ")
                    L.append(str(j[neuron_type_name]["I_ex"]["sd"]))
                else:
                    "Error: There is no classification for I_ex:", j[neuron_type_name]["initial_setting_type"]["I_ex"]
                L.append("\n")
        """
        
        for neuron_variable in ["I_ex", "membrane_potential"]:
            L += write_neuron_parameter(sd, rd_array, region_GID, neuron_variable)
        
        # extended neuron position setting (this extended way arrowing hierachical description will be further extended to describe neuron information)
        region_name = str(*rd_array[region_GID]["neuron_info"].keys())
        for i in rd_array[region_GID]["neuron_info"][region_name]:
            subregion_name = str(*i.keys())
            for j in i[subregion_name]:
                neuron_type_name = str(*j.keys())
                L.append("position_type")
                L.append(", ")
                L.append(str(sd["region_name_to_region_GID"][region_name]))
                L.append(", ")
                L.append(str(rd_array[region_GID]["structure_info"]["subregion_name_to_subregion_LID"][subregion_name]))
                L.append(", ")
                L.append(str(rd_array[region_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][subregion_name][neuron_type_name]))
                L.append(", ")
                L.append(j[neuron_type_name]["position_type"]["position_type_name"])

                # settings by different position types
                # random
                if j[neuron_type_name]["position_type"]["position_type_name"] == "random":
                    pass
                    
                # gird square
                elif j[neuron_type_name]["position_type"]["position_type_name"] == "grid_square":
                    # check the parameter
                    n_neurons_from_grid_square = (j[neuron_type_name]["position_type"]["n_neurons_on_a_side"] **2) *j[neuron_type_name]["position_type"]["n_sheets"]
                    n_neurons_from_density =  j[neuron_type_name]["n_cells_per_mm2"] * ((rd_array[region_GID]["structure_info"]["xy_length_per_tile"] / 1000.) ** 2) * rd_array[region_GID]["structure_info"]["n_neurons_factor"]
                    
                    if (n_neurons_from_grid_square != n_neurons_from_density):
                        print("Error: mone_out:", region_name, subregion_name, neuron_type_name)
                        print( "n_neurons_from_grid_square", n_neurons_from_grid_square, "n_neurons_from_density", n_neurons_from_density)
                        sys.exit()
                    
                    L.append(", n_sheets, ")
                    L.append(str(j[neuron_type_name]["position_type"]["n_sheets"]))
                    L.append(", n_neurons_on_a_side, ")
                    L.append(str(int(j[neuron_type_name]["position_type"]["n_neurons_on_a_side"] * math.sqrt(rd_array[region_GID]["structure_info"]["n_neurons_factor"]))))

                # cluster (underdevelopment)
                elif j[neuron_type_name]["position_type"]["position_type_name"] == "cluster":
                    L.append(", n_clusters, ")
                    L.append(str(len(process[i_proc].cluster_center)))

                    for cc in process[i_proc].cluster_center:
                        L.append(str(cc))
                    pass

                # manual (underdevelopment)
                elif j[neuron_type_name]["position_type"]["position_type_name"] == "manual":
                    L.append("")

                else:
                    "Error: There is no classification for position type:", j[neuron_type_name]["position_type"]["position_type_name"]
                
                L.append("\n")
                
        # prepare counting of pre inputs per post
        intra_regional_pre_inputs_per_post = []
        n_indegree_intra_regional_connection_types_per_process = 0

        # count number of post neurons and initialize the arrays
        for i_post, row in enumerate(range(rd_array[region_GID]["structure_info"]["n_neuron_types_per_process"])):
            intra_regional_pre_inputs_per_post.append([])

        # count number of pre inputs per post
        i_connected_post_pre_neuron_type_pairs = 0
        connected_post_pre_neuron_type_pair_SerLIDs = [ [ [] for j in range(rd_array[region_GID]["structure_info"]["n_neuron_types_per_process"] ) ] for i in range(rd_array[region_GID]["structure_info"]["n_neuron_types_per_process"] ) ]
        
        for i_post, row in enumerate(rd_array[region_GID]["intra_regional_connection"]["connection_parameter_matrix"]):
            rank_synapse_in_post = 0
            for i_pre, cpm in enumerate(row):
                if len(cpm["synaptic_channel"]) > 0 and cpm["spatial_pattern"] != "None":
                
                    # count number of synaptic channels per post
                    n_indegree_intra_regional_connection_types_per_process += len(cpm["synaptic_channel"])

                    # set running number of post pre of connections
                    
                    connected_post_pre_neuron_type_pair_SerLIDs[i_post][i_pre] = [ i_connected_post_pre_neuron_type_pairs, len(cpm["synaptic_channel"]), str(*cpm["synaptic_channel"][0].keys()), str(*cpm["synaptic_channel"][0][str(*cpm["synaptic_channel"][0].keys())]['dynamics'].keys()),  rank_synapse_in_post]

                    i_connected_post_pre_neuron_type_pairs+=1
                    rank_synapse_in_post+=1
                    if len(cpm["synaptic_channel"]) == 1:
                        # set running number of pre according to number of channels
                        intra_regional_pre_inputs_per_post[i_post].append(i_pre)
                        
                    elif len(cpm["synaptic_channel"]) == 2:
                        # set running number of pre according to number of channels
                        intra_regional_pre_inputs_per_post[i_post].append(i_pre)
                        intra_regional_pre_inputs_per_post[i_post].append(i_pre)

        n_connected_post_pre_neuron_type_pairs=i_connected_post_pre_neuron_type_pairs

        # wirte information of pre inputs per post
        L.append("n_connection_type_SerLIDs_per_post")
        #for ctspp in rd_array[process[i_proc].region_GID]["connection_type_SerLIDs_per_post"]:
        for nctppwp in rd_array[process[i_proc].region_GID]["n_connection_types_per_post_with_padding"]:
            L.append(", ")
            L.append(str(nctppwp))
        L.append("\n")

        #print rd_array[process[i_proc].region_GID]["connection_type_SerLIDs_per_post"]
        for i_post, ct_SerLIDs in enumerate(rd_array[process[i_proc].region_GID]["connection_type_SerLIDs_per_post"]):
            L.append("connection_type_SerLIDs_per_post, ")
            L.append(str(i_post))
            for ct_SerLID in ct_SerLIDs:
                L.append(", ")
                L.append(str(ct_SerLID[1]))
            L.append("\n")
                
        for i_pi, pre_inputs in enumerate(intra_regional_pre_inputs_per_post):
            L.append("intra_regional_pre_inputs_per_post, ")
            L.append(str(i_pi))
            for pre_input in pre_inputs:
                L.append(", ")
                L.append(str(pre_input))
            L.append("\n")

        L.append("n_connected_post_pre_neuron_type_pairs, ")
        L.append(str(n_connected_post_pre_neuron_type_pairs) )
        L.append("\n")

        # connection_name will be removed and use pre post name. The followings will be updated to the new one.
        """ connection_name 
        # check the setting file
        for i_post, row in enumerate(rd_array[region_GID]["intra_regional_connection"]["connection_parameter_matrix"]):
            for i_pre, cpm in enumerate(row):
                if cpm["spatial_pattern"] != "None":
                    connection_words = cpm["connection_name"].split(" ")
                
                    # presynaptic cell is inhibitory cell. The connection should be inhibitory
                    if "FS" in connection_words[1] or "LTS" in connection_words[1] or "SBC" in connection_words[1] or "ENGC" in connection_words[1]:
                        if "AMPA" in cpm["synaptic_channel"] or "NMDA" in cpm["synaptic_channel"]:
                            print "Error", i_post, i_pre, connection_words[1], cpm["synaptic_channel"]
                        if "log_normal_dist" in cpm["weight_distribution"]:
                            print "Error", i_post, i_pre, connection_words[1], cpm["weight_distribution"]
                        if "gaussian_dist" in cpm["weight_distribution"]:
                            if cpm["weight_distribution"]["gaussian_dist"]["mean"] > 0:
                                print "Error", i_post, i_pre, connection_words[1], cpm["weight_distribution"]

                    # presynaptic cell is excitatory cell. The connection should be excitatory
                    if "CC" in connection_words[1] or "CS" in connection_words[1] or "PT" in connection_words[1] or "CT" in connection_words[1]:
                        if "GABAA" in cpm["synaptic_channel"] or "GABAB" in cpm["synaptic_channel"]:
                            print "Error", i_post, i_pre, connection_words[1], cpm["synaptic_channel"]
                        if "gaussian_dist" in cpm["weight_distribution"]:
                            if cpm["weight_distribution"]["gaussian_dist"]["mean"] < 0:
                                print "Error", i_post, i_pre, connection_words[1], cpm["weight_distribution"]
        """
              
        # only post
        L.append("intra_regional_process")
        for i in process[i_proc].intra_regional_process_pairs:
            L.append(", ")
            L.append(str(i[0]))
        L.append("\n")

        # intra_regional_pre_neuron_types_in_pre_process_per_post_neuron_type
        for i_pre_process, pre_process in enumerate(process[i_proc].intra_regional_pre_post_address3):
            for pre_nt_SerLID, infos_per_pre_nt in enumerate(pre_process):
                L.append("intra_regional_post_neuron_types_per_pre_process_per_pre_neuron_type, i_pre_process, " + str(i_pre_process) )
                for post_nt_SerLID, info in enumerate(infos_per_pre_nt):
                    L.append(", pre_nt_SerLID, " + str(pre_nt_SerLID) )
                    L.append(", " + str(info[0]) + ", " + str(info[1]) + ", " + str(info[2]) + ", " + str(info[3]) + ", " + str(info[4]) + ", " + str(info[5]) + ", " + str(info[6]) + ", " + str(info[7]) + ", " + str(info[8]) + ", " +str(info[9]) + ", " + str(info[10]) )
                L.append("\n")

        # inter_regional_pre_neuron_types_in_pre_process_per_post_neuron_type
        for bundle_LID_at_process, bundle in enumerate(rd_array[process[i_proc].region_GID]["inter_regional_connection"]):
            bundle_name = "bundle_" + bundle["pre"]["region"] + "_to_" + bundle["post"]["region"]
            
            for pre_proc_LID_at_bundle, multiple_pre_addresses in enumerate(process[i_proc].bdl[bundle_name].inter_regional_post_neuron_types_per_pre_neuron_type_per_bundle):
                for pre_nt_LID_at_bundle, single_pre_addresses in enumerate(multiple_pre_addresses):
                    L.append("inter_regional_post_neuron_types_per_pre_process_per_pre_neuron_type, bundle_LID_at_process, " + str(bundle_LID_at_process) + ", pre_proc_LID_at_bundle, " + str(pre_proc_LID_at_bundle) + ", pre_nt_LID_at_bundle, " + str(pre_nt_LID_at_bundle) )
                    #print(rd_array[process[i_proc].region_GID]["bdl"][bundle_name].post_nt_per_pre_nt_at_bundle, pre_nt_LID_at_bundle)
                    #print(process[i_proc].bdl[bundle_name].inter_regional_post_neuron_types_per_pre_neuron_type_per_bundle)
                    #print(rd_array[process[i_proc].region_GID]["bdl"][bundle_name].pre_nt_SIDs_at_bundle[pre_nt_LID_at_bundle])
                    pre_nt_SID_at_bundle = rd_array[process[i_proc].region_GID]["bdl"][bundle_name].pre_nt_SIDs_at_bundle[pre_nt_LID_at_bundle]
                    #print(str(len(rd_array[process[i_proc].region_GID]["bdl"][bundle_name].post_nt_per_pre_nt_at_bundle[ pre_nt_SID_at_bundle ] )))
                    L.append(", n_post_nt_per_pre_nt, " + str(len(rd_array[process[i_proc].region_GID]["bdl"][bundle_name].post_nt_per_pre_nt_at_bundle[ pre_nt_SID_at_bundle ] )) )
                    for pre_post_pair_address in single_pre_addresses:
                        L.append(", " + str(pre_post_pair_address[0]) + ", " + str(pre_post_pair_address[1]) + ", " + str(pre_post_pair_address[2]) + ", " + str(pre_post_pair_address[3]) + ", " + str(pre_post_pair_address[4]) + ", " + str(pre_post_pair_address[5]) + ", " + str(pre_post_pair_address[6]) + ", " + str(pre_post_pair_address[7]) + ", " + str(pre_post_pair_address[8]) + ", " +str(pre_post_pair_address[9]) + ", " + str(pre_post_pair_address[10]) )
                    L.append("\n")

        # numbers of estimated_intra outdegree connections per one presynaptic neuron type
        L.append("n_estimated_intra_outdegree_connections_per_presynaptic_neuron_type")
        for neocppnt in rd_array[process[i_proc].region_GID]["structure_info"]["n_estimated_intra_outdegree_connections_per_presynaptic_neuron_type"]:
            L.append(", " + str(int(neocppnt)))
        L.append("\n")

        # numbers of estimated_intra indegree connections per one postsynaptic neuron type
        L.append("n_estimated_intra_indegree_connections_per_postsynaptic_neuron_type")
        for neicppnt in rd_array[process[i_proc].region_GID]["structure_info"]["n_estimated_intra_indegree_connections_per_postsynaptic_neuron_type"]:
            L.append(", " + str(int(neicppnt)))
        L.append("\n")
        
        # numbers of estimated_intra outdegree connections per one presynaptic neuron type
        L.append("n_estimated_intra_outdegree_connections_per_one_presynaptic_neuron")
        for neocppn in rd_array[process[i_proc].region_GID]["structure_info"]["n_estimated_intra_outdegree_connections_per_one_presynaptic_neuron"]:
            L.append(", " + str(int(neocppn)))
        L.append("\n")

        # numbers of estimated_intra indegree connections per one postsynaptic neuron type
        L.append("n_estimated_intra_indegree_connections_per_one_postsynaptic_neuron")
        for neicppn in rd_array[process[i_proc].region_GID]["structure_info"]["n_estimated_intra_indegree_connections_per_one_postsynaptic_neuron"]:
            L.append(", " + str(int(neicppn)))
        L.append("\n")

        # numbers of estimated_inter indegree connections per one postsynaptic neuron type
        L.append("n_estimated_inter_indegree_connections_per_postsynaptic_neuron_type")
        for neicppnt in rd_array[process[i_proc].region_GID]["structure_info"]["n_estimated_inter_indegree_connections_per_postsynaptic_neuron_type"]:
            L.append(", " + str(int(neicppnt)))
        L.append("\n")

        # numbers of estimated_inter indegree connections per one postsynaptic neuron type
        L.append("n_estimated_inter_indegree_connections_per_one_postsynaptic_neuron")
        for neicppn in rd_array[process[i_proc].region_GID]["structure_info"]["n_estimated_inter_indegree_connections_per_one_postsynaptic_neuron"]:
            L.append(", " + str(int(neicppn)))
        L.append("\n")

        # numbers of estimated_inter outdegree connections per one presynaptic neuron type
        L.append("n_estimated_inter_outdegree_connections_per_presynaptic_neuron_type")
        for i_neiocppnt, neiocppnt in enumerate(rd_array[process[i_proc].region_GID]["structure_info"]["n_estimated_inter_outdegree_connections_per_presynaptic_neuron_type"]):
             L.append(", bundle_LID, " + str(int(i_neiocppnt)))
             for i_neiocppntb, neiocppntb  in enumerate(neiocppnt):
                 L.append(", pre_nt_LID_at_bundle, " + str(i_neiocppntb) + ", " + str(int(neiocppntb)))
        L.append("\n")
        

        # inter regional connection
        L.append("n_bundles, ")
        L.append(str(process[i_proc].n_bundles))
        L.append("\n")

        # output of connection type 
        L.extend(region_ct_strings[process[i_proc].region_GID])
        
        L.append("n_intra_regional_connection_process, ")
        L.append(str(len(process[i_proc].intra_regional_connection_process)))
        L.append("\n")
        
        L.append("intra_regional_connection_process")
        for i in process[i_proc].intra_regional_connection_process:
            L.append(", ")
            L.append(str(i))
        L.append("\n")
        
        L.append("n_connection_type_per_process, ")
        L.append(str( rd_array[process[i_proc].region_GID]["count_ct"]))
        L.append("\n")
        L.append("n_indegree_intra_regional_connection_types_per_process, ")
        L.append(str(process[i_proc].n_indegree_intra_regional_connection_types_per_process))
        L.append("\n")
        L.append("n_indegree_inter_regional_connection_types_per_process, ")
        L.append(str(rd_array[process[i_proc].region_GID]["n_indegree_inter_regional_connection_types_per_process"]))
        L.append("\n")
        L.append("n_indegree_dummy_connection_types_per_process, ")
        L.append(str(rd_array[process[i_proc].region_GID]["n_indegree_dummy_connection_types_per_process"]))
        
        L.append("\n")


        for bn in process[i_proc].bdl.keys():
            L.append("bundle_LID, ")
            L.append(str(process[i_proc].bdl[bn].bundle_LID))
            #print("B", process[i_proc].region_GID, process[i_proc].bdl[bn].bundle_name, process[i_proc].bdl[bn].bundle_LID)
            L.append(", bundle_name, ")
            L.append(str(process[i_proc].bdl[bn].bundle_name) )
            L.append(", n_indegree_inter_regional_connection_types_per_bundle, ")
            L.append(str(process[i_proc].bdl[bn].n_indegree_inter_regional_connection_types_per_bundle))
            L.append(", n_inter_regional_pre_process_per_bundle, ")
            L.append(str(len(process[i_proc].bdl[bn].inter_regional_pre_process_per_bundle) ))
            L.append(", inter_regional_pre_process_per_bundle")
            for pre_ID in process[i_proc].bdl[bn].inter_regional_pre_process_per_bundle:
                L.append(", ")
                L.append(str(pre_ID))
            L.append(", n_inter_regional_post_process_per_bundle, ")
            L.append(str(len(process[i_proc].bdl[bn].inter_regional_post_process_per_bundle) ))
            L.append(", inter_regional_post_process_per_bundle")
            for post_ID in process[i_proc].bdl[bn].inter_regional_post_process_per_bundle:
                L.append(", ")
                L.append(str(post_ID))
            L.append(", n_inter_regional_pre_post_address_per_bundle, ")
            L.append( str(len(process[i_proc].bdl[bn].inter_regional_pre_post_address_per_bundle) ))
            L.append(", inter_regional_pre_post_address_per_bundle")
            for address in process[i_proc].bdl[bn].inter_regional_pre_post_address_per_bundle:
                for element in address:
                    L.append(", ")
                    L.append(str(element))
            L.append(", n_pre_nt_at_bundle, " + str(process[i_proc].bdl[bn].n_pre_nt_at_bundle))
            L.append(", n_post_nt_per_pre_per_bundle, ")
            for pre_nt_LID_at_bundle in range(process[i_proc].bdl[bn].n_pre_nt_at_bundle):
                bundle_name = process[i_proc].bdl[bn].bundle_name
                pre_nt_SID_at_bundle = rd_array[process[i_proc].region_GID]["bdl"][bundle_name].pre_nt_SIDs_at_bundle[pre_nt_LID_at_bundle]
                L.append(str(len(rd_array[process[i_proc].region_GID]["bdl"][bundle_name].post_nt_per_pre_nt_at_bundle[ pre_nt_SID_at_bundle ] )))

            L.append("\n")

        # settings for indegree_connections_per_post_synaptic_neuron
        count_icppsn=0
        for i_icppsn, icppsn in enumerate(sd["indegree_connections_per_post_synaptic_neuron"]):
            if str(icppsn["region"]) == str(process[i_proc].region_name):
                L.append("indegree_connections_per_post_synaptic_neuron")
                L.append(", ")
                L.append(str(region_GID))
                L.append(", ")
                L.append(str(rd_array[region_GID]["structure_info"]["subregion_name_to_subregion_LID"][icppsn["subregion"]]))
                L.append(", ")
                L.append(str(rd_array[region_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][icppsn["subregion"]][icppsn["neuron_type"]]))
                L.append(", ")
                L.append(str(icppsn["n_recording_neurons"]))
                L.append(", ")
                L.append(str(icppsn["center_position"][0]))
                L.append(", ")
                L.append(str(icppsn["center_position"][1]))
                L.append(", ")
                L.append(str(icppsn["center_position"][2]))
                L.append(", ")
                L.append(str(icppsn["radius"]))
                L.append("\n")
                count_icppsn+=1
        L.append("n_indegree_connections_per_post_synaptic_neuron, " + str(count_icppsn) )

        # settings for outdegree_connections_per_pre_synaptic_neuron
        L.append("\n")
        count_ocppsn=0
        for i_ocppsn, ocppsn in enumerate(sd["outdegree_connections_per_pre_synaptic_neuron"]):
            if str(ocppsn["region"]) == str(process[i_proc].region_name):
                L.append("outdegree_connections_per_pre_synaptic_neuron")
                L.append(", ")
                L.append(str(region_GID))
                L.append(", ")
                L.append(str(rd_array[region_GID]["structure_info"]["subregion_name_to_subregion_LID"][ocppsn["subregion"]]))
                L.append(", ")
                L.append(str(rd_array[region_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][ocppsn["subregion"]][ocppsn["neuron_type"]]))
                L.append(", ")
                L.append(str(ocppsn["n_recording_neurons"]))
                L.append(", ")
                L.append(str(ocppsn["center_position"][0]))
                L.append(", ")
                L.append(str(ocppsn["center_position"][1]))
                L.append(", ")
                L.append(str(ocppsn["center_position"][2]))
                L.append(", ")
                L.append(str(ocppsn["radius"]))
                L.append("\n")
                count_ocppsn+=1
        L.append("n_outdegree_connections_per_pre_synaptic_neuron, " + str(count_ocppsn) )
        L.append("\n")

        count_rcw=0
        for i_rcw, rcw in enumerate(sd["record_connection_weights"]):
            if str(rcw["post"]["region"]) == str(process[i_proc].region_name):
                L.append("record_connection_weights")
                L.append(", ")
                pre_region_GID = sd["region_name_to_region_GID"][rcw["pre"]["region"]]
                L.append(str(pre_region_GID))
                L.append(", ")
                L.append(str(rd_array[pre_region_GID]["structure_info"]["subregion_name_to_subregion_LID"][rcw["pre"]["subregion"]]))
                L.append(", ")
                L.append(str(rd_array[pre_region_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][rcw["pre"]["subregion"]][rcw["pre"]["neuron_type"]]))
                L.append(", ")
                L.append(str(rcw["pre"]["center_position"][0]))
                L.append(", ")
                L.append(str(rcw["pre"]["center_position"][1]))
                L.append(", ")
                L.append(str(rcw["pre"]["center_position"][2]))
                L.append(", ")
                L.append(str(rcw["pre"]["radius"]))
                L.append(", ")
                post_region_GID = sd["region_name_to_region_GID"][rcw["pre"]["region"]]
                L.append(str(post_region_GID))
                L.append(", ")
                L.append(str(rd_array[post_region_GID]["structure_info"]["subregion_name_to_subregion_LID"][rcw["post"]["subregion"]]))
                L.append(", ")
                L.append(str(rd_array[post_region_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][rcw["post"]["subregion"]][rcw["post"]["neuron_type"]]))
                L.append(", ")                
                L.append(str(rcw["post"]["center_position"][0]))
                L.append(", ")
                L.append(str(rcw["post"]["center_position"][1]))
                L.append(", ")
                L.append(str(rcw["post"]["center_position"][2]))
                L.append(", ")
                L.append(str(rcw["post"]["radius"]))
                L.append(", ")
                L.append(str(rcw["n_record_post_synaptic_neurons"]))
                L.append(", ")
                L.append(str(rcw["n_record_connections"]))
                L.append(", ")

                for post_ct in rd_array[post_region_GID]["ct"]:
                    if post_ct.pre_rg_name == rcw["pre"]["region"] and post_ct.pre_sr_name == rcw["pre"]["subregion"] and post_ct.pre_nt_name == rcw["pre"]["neuron_type"] and post_ct.post_rg_name == rcw["post"]["region"] and post_ct.post_sr_name == rcw["post"]["subregion"] and post_ct.post_nt_name == rcw["post"]["neuron_type"]:
                        L.append(str(post_ct.connection_type_SerLID))
                        L.append("\n")
                count_rcw+=1
                
        L.append("n_record_connection_weights, " + str(count_rcw) )
        L.append("\n")

        #####################################################################
        # tentative DTI test
        if sd["use_DTI_data"] == True:
            L.append("n_voxel_sections_per_tile, ")
            L.append(str(len(process[i_proc].voxel_sections)))
            L.append("\n")
            
            L.append("length_on_a_side_of_voxel, ")
            L.append(str(50))
            L.append("\n")
            
            L.append("n_voxel_section_dimensions, ")
            L.append(str(len(process[i_proc].voxel_sections[0].position )))
            L.append("\n")

            # only S1 (region_GID == 1) due to only postsynaptic site having info of M1->S1 connections
            if process[i_proc].region_GID == 1:
                pre_region_GID = 0
                post_region_GID = 1
                L.append("initial_and_terminal_voxel_section_positions, ")
                for i_vs, vs in enumerate(process[pre_region_GID].voxel_sections):
                    for ip in vs.position:
                        for element in ip:
                            L.append(str(element))
                            L.append(", ")

                    for element in process[post_region_GID].voxel_sections[i_vs].terminals[0]:
                        L.append(str(element))
                        L.append(", ")
                L.append("\n")
        #######################################################################

        # settings for record_target_neuron_variables
        count_record_target_neuron_variables=0
        for i_rtnv, rtnv in enumerate(sd["record_target_neuron_variables"]):
            if str(rtnv["region"]) == str(process[i_proc].region_name):
                L.append("record_target_neuron_variables")
                L.append(", ")
                L.append(str(region_GID))
                L.append(", ")
                L.append(str(rd_array[region_GID]["structure_info"]["subregion_name_to_subregion_LID"][rtnv["subregion"]]))
                L.append(", ")
                L.append(str(rd_array[region_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][rtnv["subregion"]][rtnv["neuron_type"]]))
                L.append(", ")
                L.append(str(rtnv["n_recording_neurons"]))
                L.append(", ")
                L.append(str(rtnv["center_position"][0]))
                L.append(", ")
                L.append(str(rtnv["center_position"][1]))
                L.append(", ")
                L.append(str(rtnv["center_position"][2]))
                L.append(", ")
                L.append(str(rtnv["radius"]))
                L.append("\n")
                count_record_target_neuron_variables+=1

        L.append("n_record_target_neuron_variables, ")
        L.append(str(count_record_target_neuron_variables))
        L.append("\n")
        
        # switches
        L.extend(switch_strings)

        # finally add to output
        output += ''.join(L)
        
        with open("./process_"+str(process[i_proc].process_GID)+".dat", 'wt') as f:
            f.write(output)
            f.close()
            
def write_global_shared_info_file(sd, rd_array):
    output = ""
    L = []
    gsi_strings = make_global_shared_information_strings(sd, rd_array)
    L.extend(gsi_strings)
    switch_strings = make_switch_strings(sd)
    L.extend(switch_strings)
    hardware_info_strings = make_hardware_info_strings(sd)
    L.extend(hardware_info_strings)
    output += ''.join(L)

    file_name = "./global_shared_info.dat"
    with open(file_name, 'wt') as f:
        f.write(output)
        f.close()
        
def write_region_definition_file(sd, rd_array):
    
    for rd in rd_array:
        output = ""
        L = []
        L.extend(make_region_info_strings(rd, sd))
        L.extend(make_connection_type_strings(rd))
        output += ''.join(L)
        
        region_file_name = "./region_" + rd["structure_info"]["region_name"] + ".dat"
        with open(region_file_name, 'wt') as f:
            f.write(output)
            f.close()


def write_process_info_file(sd, process):
    output = ""
    L_temp = []

    # process and region ID
    for i_process in range(sd["n_process"]):
        L_temp.append("region_GID, " + str(i_process) + ", region_GID, " + str(process[i_process].region_GID) + ", region_name, " + str(process[i_process].region_name))

        # xy plane 4 verticies
        L_temp.append(", xy_plane_coordinates")
        for coordinate in process[i_process].xy_plane_4vertices:
            for element in coordinate:
                L_temp.append(", " + str(element))

        # xyz spatial extent for whloe tile
        L_temp.append(", xyz_cube_coordinates")
        for element in process[i_process].spatial_extent:
            L_temp.append(", " + str(element))
            
        # xyz spatial extents for subregions
        L_temp.append(", n_subregions, " + str(len(process[i_process].subregion_positions)))
        L_temp.append(", xyz_subregion_cube_coordinates" )

        for subregion_coordinates in process[i_process].subregion_positions:
            for element in subregion_coordinates:
                      L_temp.append(", " + str(element))

        L_temp.append("\n")

    output += ''.join(L_temp)
    
    region_file_name = "./process_info.dat"
    with open(region_file_name, 'wt') as f:
        f.write(output)
        f.close()

    """
    # process and region ID
    for i_process in range(sd["n_process"]):
        L_temp.append("region_GID, " + str(i_process) + ", region_GID, " + str(process[i_process].region_GID) + ", region_name, " + str(process[i_process].region_name))

        # xy plane 4 verticies
        L_temp.append(", xy_plane_coordinates")
        for coordinate in process[i_process].xy_plane_4vertices:
            for element in coordinate:
                L_temp.append(", " + str(element))

        # xyz spatial extent for whloe tile
        L_temp.append(", xyz_cube_coordinates")
        for element in process[i_process].spatial_extent:
            L_temp.append(", " + str(element))
            
        # xyz spatial extents for subregions
        L_temp.append(", n_subregions, " + str(len(process[i_process].subregion_positions)))
        L_temp.append(", xyz_subregion_cube_coordinates" )

        for subregion_coordinates in process[i_process].subregion_positions:
            for element in subregion_coordinates:
                      L_temp.append(", " + str(element))

        L_temp.append("\n")

    output += ''.join(L_temp)
    
    n_copy_files = int(sd["n_process"] / 1000) + 1
    for i_copy in range(n_copy_files):
        region_file_name = "./" + str(i_copy) + "/process_info" + "_" + str(i_copy) + ".dat"
        with open(region_file_name, 'wt') as f:
            f.write(output)
            f.close()
    """
    
def write_job_script_fugaku(sd, rd_array):

    # Fugaku settings
    small_upper_limit = 384
    large_upper_limit = 55296
    huge_upper_limit = 82944

    # resource group
    if sd["n_process"]/4 <= small_upper_limit:
        rscgrp = "small"
    elif small_upper_limit < sd["n_process"]/4 <= large_upper_limit:
        rscgrp = "large"
    elif large_upper_limit < sd["n_process"]/4 <= huge_upper_limit:
        rscgrp = "huge"
    else:
        print( "Error: n_process is wrong:", sd["n_process"]/4)
        sys.exit()
        
    absolute_path = os.getcwd()
    process_info_file_name = absolute_path + "/process_info.dat"
    gsi_file_name = absolute_path + "/global_shared_info.dat"
    for rd in rd_array:
        region_file_name = absolute_path + "/region_" + rd["structure_info"]["region_name"] + ".dat"
    
    output = ""
    L = []
    output += ''.join(L)
    L.append("#!/bin/bash\n")
    L.append("#PJM --rsc-list \"node=\"" +str(int(sd["n_process"]/4)) +"\n")
    L.append("#PJM --rsc-list \"rscunit=rscunit_ft01\"\n")
    L.append("#PJM --rsc-list \"rscgrp=" + rscgrp + "\"\n")
    L.append("#PJM --rsc-list \"elapse=0:20:00\"\n")
    L.append("#PJM --mpi \"max-proc-per-node=4\"\n")
    L.append("#PJM -S\n")
    L.append("\n")
    L.append("export PLE_MPI_STD_EMPTYFILE=off\n")
    
    L.append("llio_transfer " + absolute_path + "/monet_k\n")
    L.append("llio_transfer " + gsi_file_name + "\n")
    L.append("llio_transfer " + process_info_file_name + "\n")
    for rd in rd_array:
        region_file_name = absolute_path + "/region_" + rd["structure_info"]["region_name"] + ".dat"
        L.append("llio_transfer " + region_file_name + "\n")
    
    L.append("\n")
    L.append("module load lang\n")
    L.append("export OMP_NUM_THREADS=12\n")
    L.append("export PLE_MPI_STD_EMPTYFILE=off\n")
    L.append("\n")
    L.append("mpiexec -stdout-proc ./output.%j/%/1000r/stdout -stderr-proc ./output.%j/%/1000r/stderr ./monet_k " + absolute_path + "\n")
    L.append("\n")
    L.append("llio_transfer --purge " + absolute_path + "/monet_k\n")
    L.append("llio_transfer --purge " + gsi_file_name + " \n")
    L.append("llio_transfer --purge " + process_info_file_name + " \n")
    for rd in rd_array:
        region_file_name = absolute_path + "/region_" + rd["structure_info"]["region_name"] + ".dat"
        L.append("llio_transfer --purge " + region_file_name + "\n")

    output += ''.join(L)

    jobscript_name = "./run_monet_fugaku.sh"
    with open(jobscript_name, 'wt') as f:
        f.write(output)
        f.close()
        
def write_job_script_fugaku_detail_profile(sd, rd_array):

    # Fugaku settings
    small_upper_limit = 384
    large_upper_limit = 55296
    huge_upper_limit = 82944

    # simple profile 5 
    # standard profile 11 
    # detail profile 17

    n_sample_profile = 11
    
    sample_level = 2

    # resource group
    if sd["n_process"]/4 <= small_upper_limit:
        rscgrp = "small"
    elif small_upper_limit < sd["n_process"]/4 <= large_upper_limit:
        rscgrp = "large"
    elif large_upper_limit < sd["n_process"]/4 <= huge_upper_limit:
        rscgrp = "huge"
    else:
        print( "Error: n_process is wrong:", sd["n_process"]/4)
        sys.exit()
        
    absolute_path = os.getcwd()
    process_info_file_name = absolute_path + "/process_info.dat"
    gsi_file_name = absolute_path + "/global_shared_info.dat"
    for rd in rd_array:
        region_file_name = absolute_path + "/region_" + rd["structure_info"]["region_name"] + ".dat"
    
    output = ""
    L = []
    output += ''.join(L)
    L.append("#!/bin/bash\n")
    L.append("#PJM --rsc-list \"node=\"" +str(int(sd["n_process"]/4)) +"\n")
    L.append("#PJM --rsc-list \"rscunit=rscunit_ft01\"\n")
    L.append("#PJM --rsc-list \"rscgrp=" + rscgrp + "\"\n")
    L.append("#PJM --rsc-list \"elapse=3:00:00\"\n")
    L.append("#PJM --mpi \"max-proc-per-node=4\"\n")
    L.append("#PJM -S\n")
    L.append("\n")
    L.append("export PLE_MPI_STD_EMPTYFILE=off\n")
    
    L.append("llio_transfer " + absolute_path + "/monet_k\n")
    L.append("llio_transfer " + gsi_file_name + "\n")
    L.append("llio_transfer " + process_info_file_name + "\n")
    for rd in rd_array:
        region_file_name = absolute_path + "/region_" + rd["structure_info"]["region_name"] + ".dat"
        L.append("llio_transfer " + region_file_name + "\n")
    
    L.append("\n")
    L.append("module load lang\n")
    L.append("export OMP_NUM_THREADS=12\n")
    L.append("export PLE_MPI_STD_EMPTYFILE=off\n")
    L.append("\n")
    
    for i in range(1, n_sample_profile+1):
        L.append("fapp -C -Hevent=pa" + str(i) + " -d ./rep" + str(i) + " -L " + str(sample_level) + " mpiexec ./monet_k\n")
    
    L.append("\n")
    L.append("llio_transfer --purge " + absolute_path + "/monet_k\n")
    L.append("llio_transfer --purge " + gsi_file_name + " \n")
    L.append("llio_transfer --purge " + process_info_file_name + " \n")
    for rd in rd_array:
        region_file_name = absolute_path + "/region_" + rd["structure_info"]["region_name"] + ".dat"
        L.append("llio_transfer --purge " + region_file_name + "\n")

    output += ''.join(L)

    jobscript_name = "./run_monet_fugaku.sh"
    with open(jobscript_name, 'wt') as f:
        f.write(output)
        f.close()

        
def write_weight_distribution(rd, i):
    weight_distribution_settings = ""
    
    weight_distribution_settings += ", weight_distribution"

    if str(*rd["ct"][i].weight_distribution.keys()) == "uniform_dist":
        if isinstance(rd["ct"][i].weight_distribution[ str(*rd["ct"][i].weight_distribution.keys())]["uniform_value"], str):
            temp_coefficient_name = rd["ct"][i].weight_distribution[ str(*rd["ct"][i].weight_distribution.keys())]["uniform_value"]
            temp_coefficient = rd["intra_regional_connection"]["coefficients"][temp_coefficient_name]

            #print temp_coefficient_name, temp_coefficient 
            
            weight_distribution_settings += ", uniform_dist"
            weight_distribution_settings += ", uniform_value, " + str(temp_coefficient)
        else:
            weight_distribution_settings += ", uniform_dist"
            weight_distribution_settings += ", uniform_value, " + str(rd["ct"][i].weight_distribution[ str(*rd["ct"][i].weight_distribution.keys())]["uniform_value"])
        
    elif str(*rd["ct"][i].weight_distribution.keys()) == "gaussian_dist":
        weight_distribution_settings += ", gaussian_dist"
        weight_distribution_settings += ", mean, " + str(rd["ct"][i].weight_distribution[ str(*rd["ct"][i].weight_distribution.keys())]["mean"])
        weight_distribution_settings += ", sd, " + str(rd["ct"][i].weight_distribution[ str(*rd["ct"][i].weight_distribution.keys())]["sd"])
        
    elif str(*rd["ct"][i].weight_distribution.keys()) == "log_normal_dist":
        weight_distribution_settings += ", log_normal_dist"
        weight_distribution_settings += ", mean, " + str(rd["ct"][i].weight_distribution[ str(*rd["ct"][i].weight_distribution.keys())]["mean"])
        weight_distribution_settings += ", sd, " + str(rd["ct"][i].weight_distribution[ str(*rd["ct"][i].weight_distribution.keys())]["sd"])
    else:
        print( "Error: There is no matched name in weight_distribution: ", str(*rd["ct"][i].weight_distribution.keys()))
        sys.exit()

    return weight_distribution_settings


def write_synaptic_channel(rd, i):
    L_temp = []

    L_temp.append( ", synaptic_channel" )
    L_temp.append( ", colocalization" )
    L_temp.append( ", ")
    L_temp.append( str(rd["ct"][i].colocalization) )
    synaptic_channel_name = str(*rd["ct"][i].synaptic_channel.keys())

    L_temp.append( ", ")
    L_temp.append(str(synaptic_channel_name))
    if str(*rd["ct"][i].synaptic_channel[synaptic_channel_name]["dynamics"].keys()) == "alpha":
        L_temp.append( ", ")
        L_temp.append(str(*rd["ct"][i].synaptic_channel[synaptic_channel_name]["dynamics"].keys()))
        L_temp.append( ", tau, ")
        L_temp.append(str(rd["ct"][i].synaptic_channel[synaptic_channel_name]["dynamics"]["alpha"]["tau"]))
        L_temp.append( ", reversal_potential, " )
        L_temp.append( str(rd["ct"][i].synaptic_channel[synaptic_channel_name]["reversal_potential"]) )
        L_temp.append( write_weight_distribution(rd, i) )
        
    elif str(*rd["ct"][i].synaptic_channel[synaptic_channel_name]["dynamics"].keys()) == "beta":
        L_temp.append( ", ")
        L_temp.append(str(*rd["ct"][i].synaptic_channel[synaptic_channel_name]["dynamics"].keys()))
        L_temp.append( ", tau_r, ")
        L_temp.append( str(rd["ct"][i].synaptic_channel[synaptic_channel_name]["dynamics"]["beta"]["tau_r"]))
        L_temp.append( ", tau_d, ")
        L_temp.append( str(rd["ct"][i].synaptic_channel[synaptic_channel_name]["dynamics"]["beta"]["tau_d"]))
        L_temp.append( ", reversal_potential, " )
        L_temp.append( str(rd["ct"][i].synaptic_channel[synaptic_channel_name]["reversal_potential"]))
        L_temp.append( write_weight_distribution(rd, i))

    elif str(*rd["ct"][i].synaptic_channel[synaptic_channel_name]["dynamics"].keys()) == "exponential":
        L_temp.append( ", " )
        L_temp.append( str(*rd["ct"][i].synaptic_channel[synaptic_channel_name]["dynamics"].keys()))
        L_temp.append( ", tau, ")
        L_temp.append( str(rd["ct"][i].synaptic_channel[synaptic_channel_name]["dynamics"]["exponential"]["tau"]))
        L_temp.append( ", reversal_potential, ")
        L_temp.append( str(rd["ct"][i].synaptic_channel[synaptic_channel_name]["reversal_potential"]))
        L_temp.append( write_weight_distribution(rd, i))

    else:
        print( "Error: There is no matched name in  channel_dynamics: ", str(*rd["ct"][i].synaptic_channel[synaptic_channel_name]["dynamics"].keys()))
        sys.exit()
        
    return L_temp


def write_spatial_pattern(rd, i):
    spatial_pattern_settings = ""
    
    spatial_pattern_settings += ", spatial_pattern" 
    if str(*rd["ct"][i].spatial_pattern.keys()) == "two_dimensional_gaussian":

        if isinstance(rd["ct"][i].spatial_pattern["two_dimensional_gaussian"]["mu"], str):
            spatial_pattern_settings += ", " + str(*rd["ct"][i].spatial_pattern.keys())
            temp_coefficient_name = rd["ct"][i].spatial_pattern["two_dimensional_gaussian"]["mu"]
            #print temp_coefficient_name, rd["intra_regional_connection"]["coefficients"][temp_coefficient_name ]
            temp_coefficient = rd["intra_regional_connection"]["coefficients"][temp_coefficient_name ]
            
            spatial_pattern_settings += ", mu, " + str(rd["intra_regional_connection"]["coefficients"][temp_coefficient_name ])
            temp_coefficient_name = rd["ct"][i].spatial_pattern["two_dimensional_gaussian"]["sd"]
            spatial_pattern_settings += ", sd, " + str(rd["intra_regional_connection"]["coefficients"][temp_coefficient_name ])
            #print temp_coefficient_name, rd["intra_regional_connection"]["coefficients"][temp_coefficient_name ]

        else:
            spatial_pattern_settings += ", " + str(*rd["ct"][i].spatial_pattern.keys())
            spatial_pattern_settings += ", mu, " + str(rd["ct"][i].spatial_pattern["two_dimensional_gaussian"]["mu"])
            spatial_pattern_settings += ", sd, " + str(rd["ct"][i].spatial_pattern["two_dimensional_gaussian"]["sd"])
        
    elif str(*rd["ct"][i].spatial_pattern.keys()) == "orthogornal_cross":
        spatial_pattern_settings += ", " + str(*rd["ct"][i].spatial_pattern.keys())
        spatial_pattern_settings += ", a0, " + str(rd["ct"][i].spatial_pattern["orthogornal_cross"]["a0"])
        spatial_pattern_settings += ", a1, " + str(rd["ct"][i].spatial_pattern["orthogornal_cross"]["a1"])
        spatial_pattern_settings += ", pre_width, " + str(rd["ct"][i].spatial_pattern["orthogornal_cross"]["pre_width"])
        spatial_pattern_settings += ", post_width, " + str(rd["ct"][i].spatial_pattern["orthogornal_cross"]["post_width"])
        spatial_pattern_settings += ", probability, " + str(rd["ct"][i].spatial_pattern["orthogornal_cross"]["probability"])
        
    elif str(*rd["ct"][i].spatial_pattern.keys()) == "circular":
        spatial_pattern_settings += ", " + str(*rd["ct"][i].spatial_pattern.keys())
        spatial_pattern_settings += ", radius, " + str(rd["ct"][i].spatial_pattern["circular"]["radius"])
        spatial_pattern_settings += ", probability, " + str(rd["ct"][i].spatial_pattern["circular"]["probability"])
        
    elif str(*rd["ct"][i].spatial_pattern.keys()) == "square":
        spatial_pattern_settings += ", " + str(*rd["ct"][i].spatial_pattern.keys())
        spatial_pattern_settings += ", one_side, " + str(rd["ct"][i].spatial_pattern["square"]["one_side"])
        spatial_pattern_settings += ", probability, " + str(rd["ct"][i].spatial_pattern["square"]["probability"])
        
    elif str(*rd["ct"][i].spatial_pattern.keys()) == "ID_connect":
        spatial_pattern_settings += ", " + str(*rd["ct"][i].spatial_pattern.keys()) + ", "
        spatial_pattern_settings += "n_neuron_ID_pairs, " + str(len(rd["ct"][i].spatial_pattern["ID_connect"]["neuron_ID_pairs"]))
        for neuron_pairs in rd["ct"][i].spatial_pattern["ID_connect"]["neuron_ID_pairs"]:
            for np in neuron_pairs:
                for element in np:
                    spatial_pattern_settings += ", " + str(element)
    else:
        print( "Error: There is no matched name in spatial_pattern: ", str(*rd["ct"][i].spatial_pattern.keys()))
        sys.exit()
    
    return spatial_pattern_settings

def make_global_shared_information_strings(sd, rd_array):
    L_temp = []

    # system infomation
    L_temp.append("n_total_process, ")
    L_temp.append(str(sd["n_process"]))
    L_temp.append("\n")
    
    L_temp.append("n_OpenMP_threads, ")
    L_temp.append(str(sd["n_OpenMP_threads"]))
    L_temp.append("\n")
    
    L_temp.append("simulation_time, ")
    L_temp.append(str(sd["simulation_time"]))
    L_temp.append("\n")
    
    L_temp.append("dt, ")
    L_temp.append(str(sd["dt"]))
    L_temp.append("\n")

    L_temp.append("PRNG_seed, ")
    L_temp.append(str(sd["PRNG_seed"]))
    L_temp.append("\n")

    region_names = []
    subregion_names = []
    neuron_type_names = []

    n_regions = len(rd_array)
    n_subregions_per_region = []
    n_neuron_types_per_subregion = []
    n_neurons_per_neuron_type = []
    n_neurons_per_process = []

    subregion_GID_to_region_GID = []
    subregion_GID_to_subregion_LID = []
    neuron_type_GID_to_region_GID = []
    neuron_type_GID_to_subregion_LID = []
    neuron_type_GID_to_subregion_GID = []
    neuron_type_GID_to_neuron_type_LID = []
    region_subregion_neuron_type_names_to_neuron_type_GID = {}

    xy_lengths_per_tile = []

    count_sr_GID = 0    
    count_nt_GID = 0
    for i_rd, rd in enumerate(rd_array):
        rn = rd["structure_info"]["region_name"]
        n_subregions_per_region.append(len(rd["neuron_info"][rn]))
        region_names.append(rn)
        xy_lengths_per_tile.append(rd["structure_info"]["xy_length_per_tile"])
        region_subregion_neuron_type_names_to_neuron_type_GID[rn] = {}
        
        for i_sr, sr in enumerate(rd["neuron_info"][rn]):
            srn = str(*sr.keys())
            subregion_names.append(srn)
            n_neuron_types_per_subregion.append(len(rd["neuron_info"][rn][i_sr][srn]))
            subregion_GID_to_region_GID.append(i_rd)
            subregion_GID_to_subregion_LID.append(i_sr)
            region_subregion_neuron_type_names_to_neuron_type_GID[rn][srn]={}

            for i_nt, nt in enumerate(rd[ "neuron_info"][rn][i_sr][srn]):
                ntn = str(*nt.keys())
                neuron_type_names.append(ntn)
                neuron_type_GID_to_region_GID.append(i_rd)
                neuron_type_GID_to_subregion_LID.append(i_sr)
                neuron_type_GID_to_subregion_GID.append(count_sr_GID)
                neuron_type_GID_to_neuron_type_LID.append(i_nt)
                n_neurons_per_process.append(rd["structure_info"]["n_neurons_per_neuron_type"][i_nt])
                region_subregion_neuron_type_names_to_neuron_type_GID[rn][srn][ntn] = count_nt_GID

                count_nt_GID+=1
            count_sr_GID+=1

    #print region_subregion_neuron_type_names_to_neuron_type_GID


    
    L_temp.append("gsi_region_names")
    
    for rgn in region_names:
        L_temp.append(", ")
        L_temp.append(str(rgn))
    L_temp.append("\n")
    
    L_temp.append("gsi_subregion_names")
    for srn in subregion_names:
        L_temp.append(", ")
        L_temp.append(str(srn))
    L_temp.append("\n")

    L_temp.append("gsi_neuron_type_names")
    for ntn in neuron_type_names:
        L_temp.append(", ")
        L_temp.append(str(ntn))
    L_temp.append("\n")

    L_temp.append("gsi_n_regions, ")
    L_temp.append(str(n_regions))
    L_temp.append("\n")

    L_temp.append("gsi_n_subregions_per_region")
    for n_sr in n_subregions_per_region:
        L_temp.append(", ")
        L_temp.append(str(n_sr))
    L_temp.append("\n")

    L_temp.append("gsi_n_neuron_types_per_subregion")
    for n_nt in n_neuron_types_per_subregion:
        L_temp.append(", ")
        L_temp.append(str(n_nt))
    L_temp.append("\n")

    L_temp.append("gsi_n_all_subregions, ")
    L_temp.append(str(sum(n_subregions_per_region)))
    L_temp.append("\n")

    L_temp.append("gsi_n_all_neuron_types, ")
    L_temp.append(str(sum(n_neuron_types_per_subregion)))
    L_temp.append("\n")
    
    L_temp.append("gsi_subregion_GID_to_region_GID")
    for rg_GID in subregion_GID_to_region_GID:
        L_temp.append(", ")
        L_temp.append(str(rg_GID))
    L_temp.append("\n")
    
    L_temp.append("gsi_subregion_GID_to_subregion_LID")
    for sr_LID in subregion_GID_to_subregion_LID:
        L_temp.append(", ")
        L_temp.append(str(sr_LID))
    L_temp.append("\n")
    
    L_temp.append("gsi_neuron_type_GID_to_region_GID")
    for rg_GID in neuron_type_GID_to_region_GID:
        L_temp.append(", ")
        L_temp.append(str(rg_GID))
    L_temp.append("\n")
    
    L_temp.append("gsi_neuron_type_GID_to_subregion_LID")
    for sr_LID in neuron_type_GID_to_subregion_LID:
        L_temp.append(", ")
        L_temp.append(str(sr_LID))
    L_temp.append("\n")

    L_temp.append("gsi_neuron_type_GID_to_subregion_GID")
    for sr_GID in neuron_type_GID_to_subregion_GID:
        L_temp.append(", ")
        L_temp.append(str(sr_GID))
    L_temp.append("\n")
    
    L_temp.append("gsi_neuron_type_GID_to_neuron_type_LID")
    for nt_LID in neuron_type_GID_to_neuron_type_LID:
        L_temp.append(", ")
        L_temp.append(str(nt_LID))
    L_temp.append("\n")

    L_temp.append("gsi_n_neurons_per_process")
    for nnpp in n_neurons_per_process:
        L_temp.append(", ")
        L_temp.append(str(nnpp))
    L_temp.append("\n")
    
    L_temp.append("gsi_xy_lengths_per_tile")
    for xylpt in xy_lengths_per_tile:
        L_temp.append(", ")
        L_temp.append(str(xylpt))
    L_temp.append("\n")

    
    """
    L_temp.append("gsi_region_GID_to_region_name")
    L_temp.append("gsi_subregion_LID_to_subregion_name")
    L_temp.append("gsi_neuron_type_LID_to_neuron_type_name")
    """
    return L_temp

def write_neuron_parameter(sd, rd_array, region_GID, neuron_variable_name):

    region_name = str(*rd_array[region_GID]["neuron_info"].keys())
    L_temp = []
    for i in rd_array[region_GID]["neuron_info"][region_name]:
        subregion_name = str(*i.keys())
        for j in i[subregion_name]:
            neuron_type_name = str(*j.keys())
            L_temp.append(neuron_variable_name)

            L_temp.append(", ")
            L_temp.append(str(sd["region_name_to_region_GID"][region_name]))
            L_temp.append(", ")
            L_temp.append(str(rd_array[region_GID]["structure_info"]["subregion_name_to_subregion_LID"][subregion_name]))
            L_temp.append(", ")
            L_temp.append(str(rd_array[region_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][subregion_name][neuron_type_name]))
            L_temp.append(", ")
            L_temp.append(j[neuron_type_name][neuron_variable_name]["initial_setting_type"])

            if j[neuron_type_name][neuron_variable_name]["initial_setting_type"] == "homogeneous":
                L_temp.append(", value, ")
                L_temp.append(str(j[neuron_type_name][neuron_variable_name]["value"]))

            elif j[neuron_type_name][neuron_variable_name]["initial_setting_type"] == "gaussian":
                L_temp.append(", mean, ")
                L_temp.append(str(j[neuron_type_name][neuron_variable_name]["mean"]))
                L_temp.append(", sd, ")
                L_temp.append(str(j[neuron_type_name][neuron_variable_name]["sd"]))
                
            elif j[neuron_type_name][neuron_variable_name]["initial_setting_type"] == "random":
                L_temp.append(", upper_limit, ")
                L_temp.append(str(j[neuron_type_name][neuron_variable_name]["upper_limit"]))
                L_temp.append(", lower_limit, ")
                L_temp.append(str(j[neuron_type_name][neuron_variable_name]["lower_limit"]))
            
            else:
                "Error: There is no classification:", neuron_variable_name, j[neuron_type_name]["initial_setting_type"][neuron_variable_name]
                
            L_temp.append("\n")

    return L_temp

def write_neuron_parameter2(sd, rd, region_GID, neuron_variable_name):

    region_name = str(*rd["neuron_info"].keys())
    L_temp = []
    for i in rd["neuron_info"][region_name]:
        subregion_name = str(*i.keys())
        for j in i[subregion_name]:
            neuron_type_name = str(*j.keys())
            L_temp.append(neuron_variable_name)

            L_temp.append(", ")
            L_temp.append(str(sd["region_name_to_region_GID"][region_name]))
            L_temp.append(", ")
            L_temp.append(str(rd["structure_info"]["subregion_name_to_subregion_LID"][subregion_name]))
            L_temp.append(", ")
            L_temp.append(str(rd["structure_info"]["neuron_type_name_to_neuron_type_LID"][subregion_name][neuron_type_name]))
            L_temp.append(", ")
            L_temp.append(j[neuron_type_name][neuron_variable_name]["initial_setting_type"])

            if j[neuron_type_name][neuron_variable_name]["initial_setting_type"] == "homogeneous":
                L_temp.append(", value, ")
                L_temp.append(str(j[neuron_type_name][neuron_variable_name]["value"]))

            elif j[neuron_type_name][neuron_variable_name]["initial_setting_type"] == "gaussian":
                L_temp.append(", mean, ")
                L_temp.append(str(j[neuron_type_name][neuron_variable_name]["mean"]))
                L_temp.append(", sd, ")
                L_temp.append(str(j[neuron_type_name][neuron_variable_name]["sd"]))
                
            elif j[neuron_type_name][neuron_variable_name]["initial_setting_type"] == "random":
                L_temp.append(", upper_limit, ")
                L_temp.append(str(j[neuron_type_name][neuron_variable_name]["upper_limit"]))
                L_temp.append(", lower_limit, ")
                L_temp.append(str(j[neuron_type_name][neuron_variable_name]["lower_limit"]))
            
            else:
                "Error: There is no classification:", neuron_variable_name, j[neuron_type_name]["initial_setting_type"][neuron_variable_name]
                
            L_temp.append("\n")

    return L_temp

def make_switch_strings(sd):
    L_temp = []
    
    # switches
    L_temp.append("record_calculate_time, ")
    L_temp.append(str(sd["record_calculate_time"]["switch"]))
    if sd["record_calculate_time"]["switch"] == True:
       # master_zero, all, one_per_region, process_ID 
       L_temp.append(", process_rank, ")
       L_temp.append(str(sd["record_calculate_time"]["process_rank"]))
    L_temp.append("\n")
    
    L_temp.append("record_build_time, ")
    L_temp.append(str(sd["record_build_time"]["switch"]))
    if sd["record_build_time"]["switch"] == True:
       # master_zero, all, one_per_region, process_ID 
       L_temp.append(", process_rank, ")
       L_temp.append(str(sd["record_build_time"]["process_rank"]))
    L_temp.append("\n")
    
    L_temp.append("record_spike_times_per_process, ")
    L_temp.append(str(sd["record_spike_times_per_process"]))
    L_temp.append("\n")
    L_temp.append("record_neuron_positions_per_process, ")
    L_temp.append(str(sd["record_neuron_positions_per_process"]))
    L_temp.append("\n")
    L_temp.append("record_spike_times_per_neuron_type, ")
    L_temp.append(str(sd["record_spike_times_per_neuron_type"]))
    L_temp.append("\n")
    L_temp.append("record_neuron_positions_per_neuron_type, ")
    L_temp.append(str(sd["record_neuron_positions_per_neuron_type"]))
    L_temp.append("\n")        
    L_temp.append("record_neuron_state_varibles_per_neuron_type, ")
    L_temp.append(str(sd["record_neuron_state_varibles_per_neuron_type"]))
    L_temp.append("\n")
    L_temp.append("record_synaptic_conductance_buffers_per_neuron_type, ")
    L_temp.append(str(sd["record_synaptic_conductance_buffers_per_neuron_type"]))
    L_temp.append("\n")
    L_temp.append("record_spike_info, ")
    L_temp.append(str(sd["record_spike_info"]))
    L_temp.append("\n")
    L_temp.append("record_intra_regional_connection_positions, ")
    L_temp.append(str(sd["record_intra_regional_connection_positions"]))
    L_temp.append("\n")
    L_temp.append("record_intra_regional_connection_samples, ")
    L_temp.append(str(sd["record_intra_regional_connection_samples"]))
    L_temp.append("\n")
    L_temp.append("record_inter_regional_connection_positions, ")
    L_temp.append(str(sd["record_inter_regional_connection_positions"]))
    L_temp.append("\n")
    L_temp.append("record_inter_regional_connection_samples, ")
    L_temp.append(str(sd["record_inter_regional_connection_samples"]))
    L_temp.append("\n")
    L_temp.append("record_settings, ")
    L_temp.append(str(sd["record_settings"]))
    L_temp.append("\n")
    L_temp.append("record_spike_time_communications, ")
    L_temp.append(str(sd["record_spike_time_communications"]))
    L_temp.append("\n")
    L_temp.append("record_FR_subgrid_neuron_type_per_process, ")
    L_temp.append(str(sd["record_FR_subgrid_neuron_type_per_process"]))
    L_temp.append("\n")
    L_temp.append("record_FR_subgrid_per_process, ")
    L_temp.append(str(sd["record_FR_subgrid_per_process"]))
    L_temp.append("\n")
    L_temp.append("record_FR_process, ")
    L_temp.append(str(sd["record_FR_process"]))
    L_temp.append("\n")
    L_temp.append("terminal_output, ")
    L_temp.append(str(sd["terminal_output"]))
    L_temp.append("\n")
    
    for rfr in sd["record_firing_rate"].keys():
        L_temp.append("record_" + rfr + ", " )
        L_temp.append(str(sd["record_firing_rate"][rfr]))
        L_temp.append("\n")
        
    L_temp.append("cache_block_manual_setting, ")
    L_temp.append(str(sd["hardware_info"]["cache_block_manual_setting"]))
    L_temp.append("\n")

    return L_temp

def make_connection_type_strings(rd):
    L_temp = []
    
    for i in range(len(rd["ct"])):
        L_temp.append("connection_type_SerLID, ")
        L_temp.append(str(rd["ct"][i].connection_type_SerLID))
        L_temp.append(", pre_rg_GID, ")
        L_temp.append(str(rd["ct"][i].pre_rg_GID))
        L_temp.append(", pre_sr_LID, ")
        L_temp.append(str(rd["ct"][i].pre_sr_LID))
        L_temp.append(", pre_nt_LID, ")
        L_temp.append(str(rd["ct"][i].pre_nt_LID))
        L_temp.append(", pre_nt_SerLID, ")
        L_temp.append(str(rd["ct"][i].pre_nt_SerLID))
        L_temp.append(", post_rg_GID, ")
        L_temp.append(str(rd["ct"][i].post_rg_GID))
        L_temp.append(", post_sr_LID, ")
        L_temp.append(str(rd["ct"][i].post_sr_LID))
        L_temp.append(", post_nt_LID, ")
        L_temp.append(str(rd["ct"][i].post_nt_LID))
        L_temp.append(", post_nt_SerLID, ")
        L_temp.append(str(rd["ct"][i].post_nt_SerLID))
            
        L_temp.append(", pre_rg_name, ")
        L_temp.append(str(rd["ct"][i].pre_rg_name))
        L_temp.append(", pre_sr_name, ")
        L_temp.append(str(rd["ct"][i].pre_sr_name))
        L_temp.append(", pre_nt_name, ")
        L_temp.append(str(rd["ct"][i].pre_nt_name))
        L_temp.append(", post_rg_name, ")
        L_temp.append(str(rd["ct"][i].post_rg_name))
        L_temp.append(", post_sr_name, ")
        L_temp.append(str(rd["ct"][i].post_sr_name))
        L_temp.append(", post_nt_name, ")
        L_temp.append(str(rd["ct"][i].post_nt_name))
        L_temp.append(", connection_type_LID_in_nt, ")
        L_temp.append(str(rd["ct"][i].connection_type_LID_in_nt))
        L_temp.append(", inter_or_intra, ")
        L_temp.append(str(rd["ct"][i].inter_or_intra))
            
        if rd["ct"][i].inter_or_intra == "intra":
            L_temp.append(write_spatial_pattern(rd, i))
            L_temp.extend(write_synaptic_channel(rd, i))
            L_temp.append(", delay, ")
            L_temp.append(str(rd["ct"][i].delay))

            if rd["ct"][i].LTP == {}:
                L_temp.append(", LTP, false")
            # STDP
            elif rd["ct"][i].LTP["LTP_type"] == "STDP":
                L_temp.append(", LTP, ")
                L_temp.append(str(rd["ct"][i].LTP["LTP_type"]))
            
                # pre
                L_temp.append(", slave, region, ")
                L_temp.append(str(rd["ct"][i].LTP["pre"]["region"]))
                L_temp.append(", subregion, ")
                L_temp.append(str(rd["ct"][i].LTP["pre"]["subregion"]))
                L_temp.append(", neuron_type, ")
                L_temp.append(str(rd["ct"][i].LTP["pre"]["neuron_type"]))
                # post
                L_temp.append(", post, region, ")
                L_temp.append(str(rd["ct"][i].LTP["post"]["region"]))
                L_temp.append(", subregion, ")
                L_temp.append(str(rd["ct"][i].LTP["post"]["subregion"]))
                L_temp.append(", neuron_type, ")
                L_temp.append(str(rd["ct"][i].LTP["post"]["neuron_type"]))
                # STDP function (asymetric two exponential functions)
                L_temp.append(", LTP_parameters")
                L_temp.append(", time_window, ")
                L_temp.append(str(rd["ct"][i].LTP["LTP_parameters"]["time_window"]))
                L_temp.append(", right_coefficient, ")
                L_temp.append(str(rd["ct"][i].LTP["LTP_parameters"]["right_coefficient"]))
                L_temp.append(", right_time_constant, ")
                L_temp.append(str(rd["ct"][i].LTP["LTP_parameters"]["right_time_constant"]))
                L_temp.append(", letf_coefficient, ")
                L_temp.append(str(rd["ct"][i].LTP["LTP_parameters"]["letf_coefficient"]))
                L_temp.append(", left_time_constant, ")
                L_temp.append(str(rd["ct"][i].LTP["LTP_parameters"]["left_time_constant"]))
                    
            # PFCFP_MASTER
            elif rd["ct"][i].LTP["LTP_type"] == "PFCFP_master":
                L_temp.append(", LTP, ")
                L_temp.append(str(rd["ct"][i].LTP["LTP_type"]))
                # Climbing fiber
                L_temp.append(", master_pre_nt_SerLID, ")
                # convert GID, LID, and SerLID
                rg_GID = sd["region_name_to_region_GID"][rd["ct"][i].LTP["master_pre"]["region"]]
                sr_LID = rd_array[rg_GID]["structure_info"]["subregion_name_to_subregion_LID"][ rd["ct"][i].LTP["master_pre"]["subregion"] ]
                nt_LID = rd_array[rg_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][rd["ct"][i].LTP["master_pre"]["subregion"] ][rd["ct"][i].LTP["master_pre"]["neuron_type"] ]
                nt_SerLID = rd_array[rg_GID]["structure_info"]["LID_to_neuron_type_SerLID"][sr_LID][nt_LID ]
                L_temp.append(str(nt_SerLID))
                # Paralell fiber
                L_temp.append(", slave_pre_nt_SerLID, ")
                # convert GID, LID, and SerLID
                rg_GID = sd["region_name_to_region_GID"][rd["ct"][i].LTP["slave_pre"]["region"]]
                sr_LID = rd_array[rg_GID]["structure_info"]["subregion_name_to_subregion_LID"][ rd["ct"][i].LTP["slave_pre"]["subregion"] ]
                nt_LID = rd_array[rg_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][rd["ct"][i].LTP["slave_pre"]["subregion"] ][rd["ct"][i].LTP["slave_pre"]["neuron_type"] ]
                nt_SerLID = rd_array[rg_GID]["structure_info"]["LID_to_neuron_type_SerLID"][sr_LID][nt_LID ]
                L_temp.append(str(nt_SerLID))

            # PFCFP_SLAVE
            elif rd["ct"][i].LTP["LTP_type"] == "PFCFP_slave":
                L_temp.append(", LTP, ")
                L_temp.append(str(rd["ct"][i].LTP["LTP_type"]))
                # Climbing fiber
                L_temp.append(", master_pre_nt_SerLID, ")
                # convert GID, LID, and SerLID
                rg_GID = sd["region_name_to_region_GID"][rd["ct"][i].LTP["master_pre"]["region"]]
                sr_LID = rd_array[rg_GID]["structure_info"]["subregion_name_to_subregion_LID"][ rd["ct"][i].LTP["master_pre"]["subregion"] ]
                nt_LID = rd_array[rg_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][rd["ct"][i].LTP["master_pre"]["subregion"] ][rd["ct"][i].LTP["master_pre"]["neuron_type"] ]
                master_pre_nt_SerLID = rd_array[rg_GID]["structure_info"]["LID_to_neuron_type_SerLID"][sr_LID][nt_LID ]
                L_temp.append(str(master_pre_nt_SerLID))
                # Paralell fiber
                L_temp.append(", slave_pre_nt_SerLID, ")
                # convert GID, LID, and SerLID
                rg_GID = sd["region_name_to_region_GID"][rd["ct"][i].LTP["slave_pre"]["region"]]
                sr_LID = rd_array[rg_GID]["structure_info"]["subregion_name_to_subregion_LID"][ rd["ct"][i].LTP["slave_pre"]["subregion"] ]
                nt_LID = rd_array[rg_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][rd["ct"][i].LTP["slave_pre"]["subregion"] ][rd["ct"][i].LTP["slave_pre"]["neuron_type"] ]
                slave_pre_nt_SerLID = rd_array[rg_GID]["structure_info"]["LID_to_neuron_type_SerLID"][sr_LID][nt_LID ]
                L_temp.append(str(slave_pre_nt_SerLID))
                # Perkinjer cell
                L_temp.append(", post_nt_SerLID, ")
                rg_GID = sd["region_name_to_region_GID"][rd["ct"][i].LTP["post"]["region"]]
                sr_LID = rd_array[rg_GID]["structure_info"]["subregion_name_to_subregion_LID"][ rd["ct"][i].LTP["post"]["subregion"] ]
                nt_LID = rd_array[rg_GID]["structure_info"]["neuron_type_name_to_neuron_type_LID"][rd["ct"][i].LTP["post"]["subregion"] ][rd["ct"][i].LTP["post"]["neuron_type"] ]
                post_nt_SerLID = rd_array[rg_GID]["structure_info"]["LID_to_neuron_type_SerLID"][sr_LID][nt_LID ]
                L_temp.append(str(post_nt_SerLID))
                # master_connection_type_LID_in_nt
                for ct in rd["ct"]:
                    if ct.pre_nt_SerLID  == master_pre_nt_SerLID and ct.post_nt_SerLID == post_nt_SerLID:
                        if rd["ct"][i].LTP["master_pre"]["region"] == ct.pre_rg_name and rd["ct"][i].LTP["master_pre"]["subregion"] == ct.pre_sr_name and rd["ct"][i].LTP["master_pre"]["neuron_type"] == ct.pre_nt_name:
                            master_connection_type_LID_in_nt = ct.connection_type_LID_in_nt
                L_temp.append(", master_connection_type_LID_in_nt, ")
                L_temp.append(str(master_connection_type_LID_in_nt))
                
                # PFCFP function (rectangular shape)
                L_temp.append(", LTP_parameters")
                L_temp.append(", time_window, ")
                L_temp.append(str(rd["ct"][i].LTP["LTP_parameters"]["time_window"]))
                L_temp.append(", LTD_coefficient, ")
                L_temp.append(str(rd["ct"][i].LTP["LTP_parameters"]["LTD_coefficient"]))
                L_temp.append(", LTP_coefficient, ")
                L_temp.append(str(rd["ct"][i].LTP["LTP_parameters"]["LTP_coefficient"]))
                L_temp.append(", min, ")
                L_temp.append(str(rd["ct"][i].LTP["LTP_parameters"]["min"]))
                L_temp.append(", max, ")
                L_temp.append(str(rd["ct"][i].LTP["LTP_parameters"]["max"]))
  
        elif rd["ct"][i].inter_or_intra == "inter":
            L_temp.append(", bundle_name, ")
            L_temp.append(str(rd["ct"][i].bundle_name))
            L_temp.append(", bundle_LID, ")
            L_temp.append(str(rd["ct"][i].bundle_LID))
            #print("A", rd["structure_info"]["region_name"], rd["ct"][i].bundle_name, rd["ct"][i].bundle_LID)
            L_temp.append(", connection_type_LID_in_bundle, ")
            L_temp.append(str(rd["ct"][i].connection_type_LID_in_bundle))
            L_temp.append(", topology" )
            L_temp.append(", shift_x, ")
            L_temp.append(str(rd["bdl"][rd["ct"][i].bundle_name].topology["shift"][0]))
            L_temp.append(", shift_y, ")
            L_temp.append(str(rd["bdl"][rd["ct"][i].bundle_name].topology["shift"][1]))
            L_temp.append(", shift_z, ")
            L_temp.append(str(rd["bdl"][rd["ct"][i].bundle_name].topology["shift"][2]))
            L_temp.append(", scale, ")
            L_temp.append(str(rd["bdl"][rd["ct"][i].bundle_name].topology["scale"])) 
            L_temp.append(", rotate, ")
            L_temp.append(str(rd["bdl"][rd["ct"][i].bundle_name].topology["rotate"]))
            
            L_temp.append(write_spatial_pattern(rd, i))
            L_temp.extend(write_synaptic_channel(rd, i))

            L_temp.append(", delay, ")
            L_temp.append(str(rd["ct"][i].delay))
                
            if rd["ct"][i].LTP == {}:
                L_temp.append(", LTP, false") 
            else:
                L_temp.append(", LTP, ")
                L_temp.append(str(rd["ct"][i].LTP["LTP_type"]))
                
        elif rd["ct"][i].inter_or_intra == "dummy":
            L_temp.append(write_spatial_pattern(rd, i))
            L_temp.extend(write_synaptic_channel(rd, i))
            L_temp.append(", delay, ")
            L_temp.append(str(rd["ct"][i].delay))
            if rd["ct"][i].LTP == {}:
                L_temp.append(", LTP, false")

            
        else:
            print( "Error: There is no matched name in inter intra : ", rd["ct"][i].inter_or_intra)
            sys.exit()
            
        L_temp.append("\n")
        
    return L_temp

def make_region_info_strings(rd, sd):
    L = []
    
    region_name = rd["structure_info"]["region_name"]
    L.append("region_name, ")
    L.append(str(region_name))
    L.append("\n")
    
    L.append("region_GID, ")
    L.append(str(sd["region_name_to_region_GID"][region_name]))
    L.append("\n")
    
    L.append("tile_link_limit, ")
    L.append(str(rd["structure_info"]["tile_link_limit"]))
    L.append("\n")

    L.append("xy_length_per_tile, " +str(rd["structure_info"]["xy_length_per_tile"]))
    L.append("\n")

    L.append("n_same_region_process_GIDs, ")
    L.append(str(len(rd["process_ID"])))
    L.append("\n")
    
    L.append("same_region_process_GIDs")
    for i in rd["process_ID"]:
        L.append(", ")
        L.append(str(i))
    L.append("\n")

    # structure info
    L.append("n_subregions, ")
    L.append(str(rd["structure_info"]["n_subregions"]))
    L.append("\n")

    L.append("subregion_names")
    for i_sr, sr in enumerate(rd["neuron_info"][rd["structure_info"]["region_name"]]):
        L.append(", ")
        L.append(str(*sr.keys()))
    L.append("\n")

    L.append("n_neuron_types_per_subregion")
    for i in rd["structure_info"]["n_neuron_types_per_subregion"]:
        L.append(", ")
        L.append(str(i))
    L.append("\n")

    """
    L.append("subregion_positions")
    for i in process[i_proc].subregion_positions:
        for j in i:
            L.append(", ")
            L.append(str(j))
    L.append("\n")
    """

    L.append("neuron_type_names")
    for i in rd["structure_info"]["neuron_type_names"]:
        for j in i:
            L.append(", ")
            L.append(str(j))
    L.append("\n")

    L.append("n_neurons_per_neuron_type")
    for i in rd["structure_info"]["n_neurons_per_neuron_type"]:
        L.append(", ")
        L.append(str(i))
    L.append("\n")

    # neuron info
    L.append("neuron_model")
    for i in rd["structure_info"]["neuron_model_array"]:
        for j in i:
            L.append(", ")
            L.append(str(j))
    L.append("\n")
        
    L.append("E_or_I")
    for i in rd["structure_info"]["E_or_I_array"]:
        for j in i:
            L.append(", ")
            L.append(str(j))
    L.append("\n")
         
    L.append("membrane_time_constant")
    for i in rd["structure_info"]["membrane_time_constant_array"]:
        for j in i:
            L.append(", ")
            L.append(str(j))
    L.append("\n")
 
    L.append("spike_threshold")
    for i in rd["structure_info"]["spike_threshold_array"]:
        for j in i:
            L.append(", ")
            L.append(str(j))
    L.append("\n")
 
    L.append("reset_value")
    for i in rd["structure_info"]["reset_value_array"]:
        for j in i:
            L.append(", ")
            L.append(str(j))
    L.append("\n")
    
    L.append("E_rest")
    for i in rd["structure_info"]["E_rest_array"]:
        for j in i:
            L.append(", ")
            L.append(str(j))
    L.append("\n")

    
    L.append("absolute_refractory_period")
    for i in rd["structure_info"]["absolute_refractory_period_array"]:
        for j in i:
            L.append(", ")
            L.append(str(j))
    L.append("\n")

    
    for neuron_variable in ["I_ex", "membrane_potential"]:
        L += write_neuron_parameter2(sd, rd, sd["region_name_to_region_GID"][region_name], neuron_variable)
    
    
    # extended neuron position setting (this extended way arrowing hierachical description will be further extended to describe neuron information)
    for i in rd["neuron_info"][region_name]:
        subregion_name = str(*i.keys())
        for j in i[subregion_name]:
            neuron_type_name = str(*j.keys())
            L.append("position_type")
            L.append(", ")
            L.append(str(sd["region_name_to_region_GID"][region_name]))
            L.append(", ")
            L.append(str(rd["structure_info"]["subregion_name_to_subregion_LID"][subregion_name]))
            L.append(", ")
            L.append(str(rd["structure_info"]["neuron_type_name_to_neuron_type_LID"][subregion_name][neuron_type_name]))
            L.append(", ")
            L.append(j[neuron_type_name]["position_type"]["position_type_name"])

            # settings by different position types
            # random
            if j[neuron_type_name]["position_type"]["position_type_name"] == "random":
                pass
                    
                # gird square
            elif j[neuron_type_name]["position_type"]["position_type_name"] == "grid_square":
                # check the parameter
                n_neurons_from_grid_square = (j[neuron_type_name]["position_type"]["n_neurons_on_a_side"] **2) *j[neuron_type_name]["position_type"]["n_sheets"]
                n_neurons_from_density =  j[neuron_type_name]["n_cells_per_mm2"] * ((rd["structure_info"]["xy_length_per_tile"] / 1000.) ** 2) * rd["structure_info"]["n_neurons_factor"]
                    
                if (n_neurons_from_grid_square != n_neurons_from_density):
                    print("Error: mone_out:", region_name, subregion_name, neuron_type_name)
                    print( "n_neurons_from_grid_square", n_neurons_from_grid_square, "n_neurons_from_density", n_neurons_from_density)
                    sys.exit()
                    
                L.append(", n_sheets, ")
                L.append(str(j[neuron_type_name]["position_type"]["n_sheets"]))
                L.append(", n_neurons_on_a_side, ")
                L.append(str(int(j[neuron_type_name]["position_type"]["n_neurons_on_a_side"] * math.sqrt(rd["structure_info"]["n_neurons_factor"]))))

            # cluster (underdevelopment)
            elif j[neuron_type_name]["position_type"]["position_type_name"] == "cluster":
                L.append(", n_clusters, ")
                """
                L.append(str(len(process[i_proc].cluster_center)))
                
                for cc in process[i_proc].cluster_center:
                L.append(str(cc))
                """
                pass

                    
                # manual (underdevelopment)
            elif j[neuron_type_name]["position_type"]["position_type_name"] == "manual":
                L.append("")

            else:
                "Error: There is no classification for position type:", j[neuron_type_name]["position_type"]["position_type_name"]
                
            L.append("\n")

    # prepare counting of pre inputs per post
    intra_regional_pre_inputs_per_post = []
    n_indegree_intra_regional_connection_types_per_process = 0

    # count number of post neurons and initialize the arrays
    for i_post, row in enumerate(range(rd["structure_info"]["n_neuron_types_per_process"])):
        intra_regional_pre_inputs_per_post.append([])

    # count number of pre inputs per post
    i_connected_post_pre_neuron_type_pairs = 0
    connected_post_pre_neuron_type_pair_SerLIDs = [ [ [] for j in range(rd["structure_info"]["n_neuron_types_per_process"] ) ] for i in range(rd["structure_info"]["n_neuron_types_per_process"] ) ]
        
    for i_post, row in enumerate(rd["intra_regional_connection"]["connection_parameter_matrix"]):
        rank_synapse_in_post = 0
        for i_pre, cpm in enumerate(row):
            if len(cpm["synaptic_channel"]) > 0 and cpm["spatial_pattern"] != "None":
                
                # count number of synaptic channels per post
                n_indegree_intra_regional_connection_types_per_process += len(cpm["synaptic_channel"])

                # set running number of post pre of connections
                connected_post_pre_neuron_type_pair_SerLIDs[i_post][i_pre] = [ i_connected_post_pre_neuron_type_pairs, len(cpm["synaptic_channel"]), str(*cpm["synaptic_channel"][0].keys()), str(*cpm["synaptic_channel"][0][str(*cpm["synaptic_channel"][0].keys())]['dynamics'].keys()),  rank_synapse_in_post]

                i_connected_post_pre_neuron_type_pairs += 1
                rank_synapse_in_post += 1
                if len(cpm["synaptic_channel"]) == 1:
                    # set running number of pre according to number of channels
                    intra_regional_pre_inputs_per_post[i_post].append(i_pre)
                        
                elif len(cpm["synaptic_channel"]) == 2:
                    # set running number of pre according to number of channels
                    intra_regional_pre_inputs_per_post[i_post].append(i_pre)
                    intra_regional_pre_inputs_per_post[i_post].append(i_pre)

    n_connected_post_pre_neuron_type_pairs=i_connected_post_pre_neuron_type_pairs

    # wirte information of pre inputs per post
    L.append("n_connection_type_SerLIDs_per_post")
    #for ctspp in rd_array[process[i_proc].region_GID]["connection_type_SerLIDs_per_post"]:
    for nctppwp in rd["n_connection_types_per_post_with_padding"]:
        L.append(", ")
        L.append(str(nctppwp))
    L.append("\n")

    #print rd_array[process[i_proc].region_GID]["connection_type_SerLIDs_per_post"]
    for i_post, ct_SerLIDs in enumerate(rd["connection_type_SerLIDs_per_post"]):
        L.append("connection_type_SerLIDs_per_post, ")
        L.append(str(i_post))
        for ct_SerLID in ct_SerLIDs:
            L.append(", ")
            L.append(str(ct_SerLID[1]))
        L.append("\n")
                
    for i_pi, pre_inputs in enumerate(intra_regional_pre_inputs_per_post):
        L.append("intra_regional_pre_inputs_per_post, ")
        L.append(str(i_pi))
        for pre_input in pre_inputs:
            L.append(", ")
            L.append(str(pre_input))
        L.append("\n")

    L.append("n_connected_post_pre_neuron_type_pairs, ")
    L.append(str(n_connected_post_pre_neuron_type_pairs) )
    L.append("\n")

    
    intra_regional_pre_post_address = [ [] for i in rd["intra_regional_connection"]["connection_parameter_matrix"] ]
    for post_nt_SerLID, row in enumerate(rd["intra_regional_connection"]["connection_parameter_matrix"]):
        for pre_nt_SerLID, cpe in enumerate(row):
            intra_regional_pre_post_address[pre_nt_SerLID].append([
                cpe["pre_post_colocalization"][6],
                cpe["pre_post_colocalization"][7],
                cpe["pre_post_colocalization"][8],
                cpe["pre_post_colocalization"][0],
                cpe["pre_post_colocalization"][1],
                cpe["pre_post_colocalization"][2],
                cpe["pre_post_colocalization"][3],
                cpe["pre_post_colocalization"][4],
                cpe["pre_post_colocalization"][5],
            ])
            
    #print(intra_regional_pre_post_address)
    """
    for i_irppa, irppa in enumerate(intra_regional_pre_post_address):
        L.append("intra_regional_post_neuron_types_per_pre_neuron_type")
        for one_pair in irppa:
            L.append(", pre_nt_SerLID, " + str(i_irppa) )
            L.append(", " +str(one_pair[0]))
            L.append(", " +str(one_pair[1]))
            L.append(", " +str(one_pair[2]))
            L.append(", " +str(one_pair[3]))
            L.append(", " +str(one_pair[4]))
            L.append(", " +str(one_pair[5]))
            L.append(", " +str(one_pair[6]))
            L.append(", " +str(one_pair[7]))
            L.append(", " +str(one_pair[8]))
        L.append("\n")
    """     
    
    """
    # intra_regional_pre_neuron_types_in_pre_process_per_post_neuron_type
    for i_pre_process, pre_process in enumerate(process[i_proc].intra_regional_pre_post_address3):
        for pre_nt_SerLID, infos_per_pre_nt in enumerate(pre_process):
            L.append("intra_regional_post_neuron_types_per_pre_process_per_pre_neuron_type, i_pre_process, " + str(i_pre_process) )
            for post_nt_SerLID, info in enumerate(infos_per_pre_nt):
                L.append(", pre_nt_SerLID, " + str(pre_nt_SerLID) )
                L.append(", " + str(info[0]) + ", " + str(info[1]) + ", " + str(info[2]) + ", " + str(info[3]) + ", " + str(info[4]) + ", " + str(info[5]) + ", " + str(info[6]) + ", " + str(info[7]) + ", " + str(info[8]) + ", " +str(info[9]) + ", " + str(info[10]) )
            L.append("\n")
    """
    for pre_nt_SerLID, irpntppnt in enumerate(rd["intra_regional_post_neuron_types_per_pre_neuron_type"]):
        L.append("intra_regional_post_neuron_types_per_pre_neuron_type")

        for post_neuron_type in irpntppnt:
            L.append(", pre_nt_SerLID, " + str(pre_nt_SerLID) )

            for element in post_neuron_type:
                L.append(", " + str(element))
        L.append("\n")

    """
    # inter_regional_pre_neuron_types_in_pre_process_per_post_neuron_type
    for bundle_LID_at_process, bundle in enumerate(rd_array[process[i_proc].region_GID]["inter_regional_connection"]):
        bundle_name = "bundle_" + bundle["pre"]["region"] + "_to_" + bundle["post"]["region"]
            
        for pre_proc_LID_at_bundle, multiple_pre_addresses in enumerate(process[i_proc].bdl[bundle_name].inter_regional_post_neuron_types_per_pre_neuron_type_per_bundle):
            for pre_nt_LID_at_bundle, single_pre_addresses in enumerate(multiple_pre_addresses):
                L.append("inter_regional_post_neuron_types_per_pre_process_per_pre_neuron_type, bundle_LID_at_process, " + str(bundle_LID_at_process) + ", pre_proc_LID_at_bundle, " + str(pre_proc_LID_at_bundle) + ", pre_nt_LID_at_bundle, " + str(pre_nt_LID_at_bundle) )
                #print(rd_array[process[i_proc].region_GID]["bdl"][bundle_name].post_nt_per_pre_nt_at_bundle, pre_nt_LID_at_bundle)
                #print(process[i_proc].bdl[bundle_name].inter_regional_post_neuron_types_per_pre_neuron_type_per_bundle)
                #print(rd_array[process[i_proc].region_GID]["bdl"][bundle_name].pre_nt_SIDs_at_bundle[pre_nt_LID_at_bundle])
                pre_nt_SID_at_bundle = rd_array[process[i_proc].region_GID]["bdl"][bundle_name].pre_nt_SIDs_at_bundle[pre_nt_LID_at_bundle]
                #print(str(len(rd_array[process[i_proc].region_GID]["bdl"][bundle_name].post_nt_per_pre_nt_at_bundle[ pre_nt_SID_at_bundle ] )))
                L.append(", n_post_nt_per_pre_nt, " + str(len(rd_array[process[i_proc].region_GID]["bdl"][bundle_name].post_nt_per_pre_nt_at_bundle[ pre_nt_SID_at_bundle ] )) )
                for pre_post_pair_address in single_pre_addresses:
                    L.append(", " + str(pre_post_pair_address[0]) + ", " + str(pre_post_pair_address[1]) + ", " + str(pre_post_pair_address[2]) + ", " + str(pre_post_pair_address[3]) + ", " + str(pre_post_pair_address[4]) + ", " + str(pre_post_pair_address[5]) + ", " + str(pre_post_pair_address[6]) + ", " + str(pre_post_pair_address[7]) + ", " + str(pre_post_pair_address[8]) + ", " +str(pre_post_pair_address[9]) + ", " + str(pre_post_pair_address[10]) )
                L.append("\n")
    """
    
    for bundle_LID_at_region, bundle in enumerate(rd["inter_regional_connection"]):
        bundle_name = "bundle_" + bundle["pre"]["region"] + "_to_" + bundle["post"]["region"]

        for i_multiple, multiple in enumerate(rd["bdl"][bundle_name].new_inter_regional_post_neuron_types_per_pre_neuron_type_at_bundle):
            for i_single, single in enumerate(multiple):

                L.append("inter_regional_post_neuron_types_per_pre_neuron_type_at_bundle, bundle_LID_at_region, " + str(bundle_LID_at_region ) + ", pre_nt_LID_at_bundle, " + str(i_multiple) )

                for element in single:
                    L.append(", " + str(element))
                L.append("\n")

    
    # numbers of estimated_intra outdegree connections per one presynaptic neuron type
    L.append("n_estimated_intra_outdegree_connections_per_presynaptic_neuron_type")
    for neocppnt in rd["structure_info"]["n_estimated_intra_outdegree_connections_per_presynaptic_neuron_type"]:
        L.append(", " + str(int(neocppnt)))
    L.append("\n")

    # numbers of estimated_intra indegree connections per one postsynaptic neuron type
    L.append("n_estimated_intra_indegree_connections_per_postsynaptic_neuron_type")
    for neicppnt in rd["structure_info"]["n_estimated_intra_indegree_connections_per_postsynaptic_neuron_type"]:
        L.append(", " + str(int(neicppnt)))
    L.append("\n")
        
    # numbers of estimated_intra outdegree connections per one presynaptic neuron type
    L.append("n_estimated_intra_outdegree_connections_per_one_presynaptic_neuron")
    for neocppn in rd["structure_info"]["n_estimated_intra_outdegree_connections_per_one_presynaptic_neuron"]:
        L.append(", " + str(int(neocppn)))
    L.append("\n")

    # numbers of estimated_intra indegree connections per one postsynaptic neuron type
    L.append("n_estimated_intra_indegree_connections_per_one_postsynaptic_neuron")
    for neicppn in rd["structure_info"]["n_estimated_intra_indegree_connections_per_one_postsynaptic_neuron"]:
        L.append(", " + str(int(neicppn)))
    L.append("\n")

    # numbers of estimated_inter indegree connections per one postsynaptic neuron type
    L.append("n_estimated_inter_indegree_connections_per_postsynaptic_neuron_type")
    for neicppnt in rd["structure_info"]["n_estimated_inter_indegree_connections_per_postsynaptic_neuron_type"]:
        L.append(", " + str(int(neicppnt)))
    L.append("\n")

    # numbers of estimated_inter indegree connections per one postsynaptic neuron type
    L.append("n_estimated_inter_indegree_connections_per_one_postsynaptic_neuron")
    for neicppn in rd["structure_info"]["n_estimated_inter_indegree_connections_per_one_postsynaptic_neuron"]:
        L.append(", " + str(int(neicppn)))
    L.append("\n")

    # numbers of estimated_inter outdegree connections per one presynaptic neuron type
    #L.append("n_estimated_inter_outdegree_connections_per_presynaptic_neuron_type")
    for i_neiocppnt, neiocppnt in enumerate(rd["structure_info"]["n_estimated_inter_outdegree_connections_per_presynaptic_neuron_type"]):
        L.append("n_estimated_inter_outdegree_connections_per_presynaptic_neuron_type, bundle_LID, " + str(int(i_neiocppnt)))
        for i_neiocppntb, neiocppntb  in enumerate(neiocppnt):
            L.append(", pre_nt_LID_at_bundle, " + str(i_neiocppntb) + ", " + str(int(neiocppntb)))
            L.append("\n")
    
    # inter regional connection
    L.append("n_bundles, ")
    L.append(str(len(rd["bdl"])))
    L.append("\n")


    L.append("n_connection_type_per_process, ")
    L.append(str( rd["count_ct"]))
    L.append("\n")
    
    L.append("n_indegree_intra_regional_connection_types_per_process, ")
    L.append(str(rd["n_indegree_intra_regional_connection_types_per_process"]))
    L.append("\n")
    
    L.append("n_indegree_inter_regional_connection_types_per_process, ")
    L.append(str(rd["n_indegree_inter_regional_connection_types_per_process"]))
    L.append("\n")
    L.append("n_indegree_dummy_connection_types_per_process, ")
    L.append(str(rd["n_indegree_dummy_connection_types_per_process"]))
    L.append("\n")

    """
    for bn in process[i_proc].bdl.keys():
        L.append("bundle_LID, ")
        L.append(str(process[i_proc].bdl[bn].bundle_LID))
        #print("B", process[i_proc].region_GID, process[i_proc].bdl[bn].bundle_name, process[i_proc].bdl[bn].bundle_LID)
        L.append(", bundle_name, ")
        L.append(str(process[i_proc].bdl[bn].bundle_name) )
        L.append(", n_indegree_inter_regional_connection_types_per_bundle, ")
        L.append(str(process[i_proc].bdl[bn].n_indegree_inter_regional_connection_types_per_bundle))
        L.append(", n_inter_regional_pre_process_per_bundle, ")
        L.append(str(len(process[i_proc].bdl[bn].inter_regional_pre_process_per_bundle) ))
        L.append(", inter_regional_pre_process_per_bundle")
        for pre_ID in process[i_proc].bdl[bn].inter_regional_pre_process_per_bundle:
            L.append(", ")
            L.append(str(pre_ID))
        L.append(", n_inter_regional_post_process_per_bundle, ")
        L.append(str(len(process[i_proc].bdl[bn].inter_regional_post_process_per_bundle) ))
        L.append(", inter_regional_post_process_per_bundle")
        for post_ID in process[i_proc].bdl[bn].inter_regional_post_process_per_bundle:
            L.append(", ")
            L.append(str(post_ID))
        L.append(", n_inter_regional_pre_post_address_per_bundle, ")
        L.append( str(len(process[i_proc].bdl[bn].inter_regional_pre_post_address_per_bundle) ))
        L.append(", inter_regional_pre_post_address_per_bundle")
        for address in process[i_proc].bdl[bn].inter_regional_pre_post_address_per_bundle:
            for element in address:
                L.append(", ")
                L.append(str(element))
        L.append(", n_pre_nt_at_bundle, " + str(process[i_proc].bdl[bn].n_pre_nt_at_bundle))
        L.append(", n_post_nt_per_pre_per_bundle, ")
        for pre_nt_LID_at_bundle in range(process[i_proc].bdl[bn].n_pre_nt_at_bundle):
            bundle_name = process[i_proc].bdl[bn].bundle_name
            pre_nt_SID_at_bundle = rd_array[process[i_proc].region_GID]["bdl"][bundle_name].pre_nt_SIDs_at_bundle[pre_nt_LID_at_bundle]
            L.append(str(len(rd_array[process[i_proc].region_GID]["bdl"][bundle_name].post_nt_per_pre_nt_at_bundle[ pre_nt_SID_at_bundle ] )))
        
        L.append("\n")
    """
    
    for bn in rd["bdl"].keys():
        L.append("bundle_LID, ")

        L.append(str(rd["bdl"][bn].bundle_LID))
        #print("B", process[i_proc].region_GID, process[i_proc].bdl[bn].bundle_name, process[i_proc].bdl[bn].bundle_LID)
        L.append(", bundle_name, ")
        L.append(str(rd["bdl"][bn].bundle_name) )

        L.append(", n_indegree_inter_regional_connection_types_per_bundle, ")
        L.append(str(rd["bdl"][bn].n_indegree_inter_regional_connection_types_per_bundle))

        # pre region GID
        L.append(", pre_region_GID, ")
        L.append(str(rd["bdl"][bn].pre_region_GID))

        # post region GID
        L.append(", post_region_GID, ")
        L.append(str(rd["bdl"][bn].post_region_GID))

        # shift coordinate
        L.append(", shift, ")
        L.append(str(rd["bdl"][bn].shift[0]))
        L.append(", ")
        L.append(str(rd["bdl"][bn].shift[1]))
        L.append(", ")
        L.append(str(rd["bdl"][bn].shift[2]))
        
        # scale
        L.append(", scale, ")
        L.append(str(rd["bdl"][bn].scale))
        
        # tile_link_limit
        L.append(", tile_link_limit, ")
        L.append(str(rd["bdl"][bn].tile_link_limit))
        
        """
        L.append(", n_inter_regional_pre_process_per_bundle, ")
        L.append(str(len(rd["bdl"][bn].inter_regional_pre_process_per_bundle) ))
        L.append(", inter_regional_pre_process_per_bundle")
        for pre_ID in rd["bdl"][bn].inter_regional_pre_process_per_bundle:
            L.append(", ")
            L.append(str(pre_ID))
        L.append(", n_inter_regional_post_process_per_bundle, ")
        L.append(str(len(rd["bdl"][bn].inter_regional_post_process_per_bundle) ))
        L.append(", inter_regional_post_process_per_bundle")
        for post_ID in rd["bdl"][bn].inter_regional_post_process_per_bundle:
            L.append(", ")
            L.append(str(post_ID))

        L.append(", n_inter_regional_pre_post_address_per_bundle, ")
        L.append( str(len(rd["bdl"][bn].inter_regional_pre_post_address_per_bundle) ))
        L.append(", inter_regional_pre_post_address_per_bundle")
        for address in rd["bdl"][bn].inter_regional_pre_post_address_per_bundle:
            for element in address:
                L.append(", ")
                L.append(str(element))
        """
        
        L.append(", n_pre_nt_at_bundle, " + str(rd["bdl"][bn].n_pre_nt_at_bundle))
        L.append(", n_post_nt_per_pre_per_bundle, ")
        for pre_nt_LID_at_bundle in range(rd["bdl"][bn].n_pre_nt_at_bundle):
            bundle_name = rd["bdl"][bn].bundle_name
            pre_nt_SID_at_bundle = rd["bdl"][bundle_name].pre_nt_SIDs_at_bundle[pre_nt_LID_at_bundle]
            L.append(str(len(rd["bdl"][bundle_name].post_nt_per_pre_nt_at_bundle[ pre_nt_SID_at_bundle ] )))

        L.append("\n")

         
    
    # settings for indegree_connections_per_post_synaptic_neuron
    count_icppsn = 0
    for i_icppsn, icppsn in enumerate(sd["indegree_connections_per_post_synaptic_neuron"]):
        if str(icppsn["region"]) == region_name:
            L.append("indegree_connections_per_post_synaptic_neuron")
            L.append(", ")
            region_GID = sd["region_name_to_region_GID"][icppsn["region"]]
            L.append(str(region_GID))
            L.append(", ")
            L.append(str(rd["structure_info"]["subregion_name_to_subregion_LID"][icppsn["subregion"]]))
            L.append(", ")
            L.append(str(rd["structure_info"]["neuron_type_name_to_neuron_type_LID"][icppsn["subregion"]][icppsn["neuron_type"]]))
            L.append(", ")
            L.append(str(icppsn["n_recording_neurons"]))
            L.append(", ")
            L.append(str(icppsn["center_position"][0]))
            L.append(", ")
            L.append(str(icppsn["center_position"][1]))
            L.append(", ")
            L.append(str(icppsn["center_position"][2]))
            L.append(", ")
            L.append(str(icppsn["radius"]))
            L.append("\n")
            count_icppsn+=1
    L.append("n_indegree_connections_per_post_synaptic_neuron, " + str(count_icppsn) )

    # settings for outdegree_connections_per_pre_synaptic_neuron
    L.append("\n")
    count_ocppsn=0
    for i_ocppsn, ocppsn in enumerate(sd["outdegree_connections_per_pre_synaptic_neuron"]):
        if str(ocppsn["region"]) == region_name:
            L.append("outdegree_connections_per_pre_synaptic_neuron")
            L.append(", ")
            L.append(str(region_GID))
            L.append(", ")
            L.append(str(rd["structure_info"]["subregion_name_to_subregion_LID"][ocppsn["subregion"]]))
            L.append(", ")
            L.append(str(rd["structure_info"]["neuron_type_name_to_neuron_type_LID"][ocppsn["subregion"]][ocppsn["neuron_type"]]))
            L.append(", ")
            L.append(str(ocppsn["n_recording_neurons"]))
            L.append(", ")
            L.append(str(ocppsn["center_position"][0]))
            L.append(", ")
            L.append(str(ocppsn["center_position"][1]))
            L.append(", ")
            L.append(str(ocppsn["center_position"][2]))
            L.append(", ")
            L.append(str(ocppsn["radius"]))
            L.append("\n")
            count_ocppsn+=1
    L.append("n_outdegree_connections_per_pre_synaptic_neuron, " + str(count_ocppsn) )
    L.append("\n")

    count_rcw = 0
    for i_rcw, rcw in enumerate(sd["record_connection_weights"]):
        if str(rcw["post"]["region"]) == region_name:
            L.append("record_connection_weights")
            L.append(", ")
            pre_region_GID = sd["region_name_to_region_GID"][rcw["pre"]["region"]]
            L.append(str(pre_region_GID))
            L.append(", ")
            L.append(str(rd["structure_info"]["subregion_name_to_subregion_LID"][rcw["pre"]["subregion"]]))
            L.append(", ")
            L.append(str(rd["structure_info"]["neuron_type_name_to_neuron_type_LID"][rcw["pre"]["subregion"]][rcw["pre"]["neuron_type"]]))
            L.append(", ")
            L.append(str(rcw["pre"]["center_position"][0]))
            L.append(", ")
            L.append(str(rcw["pre"]["center_position"][1]))
            L.append(", ")
            L.append(str(rcw["pre"]["center_position"][2]))
            L.append(", ")
            L.append(str(rcw["pre"]["radius"]))
            L.append(", ")
            post_region_GID = sd["region_name_to_region_GID"][rcw["pre"]["region"]]
            L.append(str(post_region_GID))
            L.append(", ")
            L.append(str(rd["structure_info"]["subregion_name_to_subregion_LID"][rcw["post"]["subregion"]]))
            L.append(", ")
            L.append(str(rd["structure_info"]["neuron_type_name_to_neuron_type_LID"][rcw["post"]["subregion"]][rcw["post"]["neuron_type"]]))
            L.append(", ")                
            L.append(str(rcw["post"]["center_position"][0]))
            L.append(", ")
            L.append(str(rcw["post"]["center_position"][1]))
            L.append(", ")
            L.append(str(rcw["post"]["center_position"][2]))
            L.append(", ")
            L.append(str(rcw["post"]["radius"]))
            L.append(", ")
            L.append(str(rcw["n_record_post_synaptic_neurons"]))
            L.append(", ")
            L.append(str(rcw["n_record_connections"]))
            L.append(", ")

            for post_ct in rd["ct"]:
                if post_ct.pre_rg_name == rcw["pre"]["region"] and post_ct.pre_sr_name == rcw["pre"]["subregion"] and post_ct.pre_nt_name == rcw["pre"]["neuron_type"] and post_ct.post_rg_name == rcw["post"]["region"] and post_ct.post_sr_name == rcw["post"]["subregion"] and post_ct.post_nt_name == rcw["post"]["neuron_type"]:
                    L.append(str(post_ct.connection_type_SerLID))
                    L.append("\n")
            count_rcw+=1
                
    L.append("n_record_connection_weights, " + str(count_rcw) )
    L.append("\n")

    """
    #####################################################################
    # tentative DTI test
    if sd["use_DTI_data"] == True:
        L.append("n_voxel_sections_per_tile, ")
        L.append(str(len(process[i_proc].voxel_sections)))
        L.append("\n")
            
        L.append("length_on_a_side_of_voxel, ")
        L.append(str(50))
        L.append("\n")
            
        L.append("n_voxel_section_dimensions, ")
        L.append(str(len(process[i_proc].voxel_sections[0].position )))
        L.append("\n")

        # only S1 (region_GID == 1) due to only postsynaptic site having info of M1->S1 connections
        if process[i_proc].region_GID == 1:
            pre_region_GID = 0
            post_region_GID = 1
            L.append("initial_and_terminal_voxel_section_positions, ")
            for i_vs, vs in enumerate(process[pre_region_GID].voxel_sections):
                for ip in vs.position:
                    for element in ip:
                        L.append(str(element))
                        L.append(", ")

                for element in process[post_region_GID].voxel_sections[i_vs].terminals[0]:
                    L.append(str(element))
                    L.append(", ")
            L.append("\n")
        #######################################################################
    """
    
    # settings for record_target_neuron_variables
    count_record_target_neuron_variables = 0
    for i_rtnv, rtnv in enumerate(sd["record_target_neuron_variables"]):
        if str(rtnv["region"]) == region_name:
            L.append("record_target_neuron_variables")
            L.append(", ")
            region_GID = sd["region_name_to_region_GID"][region_name]
            L.append(str(region_GID))
            L.append(", ")
            L.append(str(rd["structure_info"]["subregion_name_to_subregion_LID"][rtnv["subregion"]]))
            L.append(", ")
            L.append(str(rd["structure_info"]["neuron_type_name_to_neuron_type_LID"][rtnv["subregion"]][rtnv["neuron_type"]]))
            L.append(", ")
            L.append(str(rtnv["n_recording_neurons"]))
            L.append(", ")
            L.append(str(rtnv["center_position"][0]))
            L.append(", ")
            L.append(str(rtnv["center_position"][1]))
            L.append(", ")
            L.append(str(rtnv["center_position"][2]))
            L.append(", ")
            L.append(str(rtnv["radius"]))
            L.append("\n")
            count_record_target_neuron_variables += 1

    L.append("n_record_target_neuron_variables, ")
    L.append(str(count_record_target_neuron_variables))
    L.append("\n")
    
    return L



def make_hardware_info_strings(sd):
    L_temp = []
    
    L_temp.append("CPU_name, ")
    L_temp.append(str(sd["hardware_info"]["CPU_name"]))
    L_temp.append("\n")
    
    L_temp.append("n_cache_blocks, ")
    L_temp.append(str(sd["hardware_info"]["n_cache_blocks"]))
    L_temp.append("\n")
        
    L_temp.append("L1_size(KB), ")
    L_temp.append(str(sd["hardware_info"]["L1_size(KB)"]))
    L_temp.append("\n")
        
    L_temp.append("L2_size(MB), ")
    L_temp.append(str(sd["hardware_info"]["L2_size(MB)"]))
    L_temp.append("\n")
    
    L_temp.append("DRAM_size_per_process(GB), ")
    L_temp.append(str(sd["hardware_info"]["DRAM_size_per_process(GB)"]))
    L_temp.append("\n")
    
    L_temp.append("GFLOPS_per_processor, ")
    L_temp.append(str(sd["hardware_info"]["GFLOPS_per_processor"]))
    L_temp.append("\n")

    L_temp.append("memory_bandwidth_per_processor(GB/s), ")
    L_temp.append(str(sd["hardware_info"]["memory_bandwidth_per_processor(GB/s)"]))
    L_temp.append("\n")
    
    return L_temp

          


def write_system_files(sd, rd_array, process, process_IDs):
    # preparetaion of gsi strings
    gsi_strings = make_global_shared_information_strings(sd, rd_array)

    # preparation of switch strings
    switch_strings = make_switche_strings(sd)

    # preparation of hardware info strings
    hardware_info_strings = make_hardware_info_strings(sd)
    
    output="""simulation_time, %(simulation_time)d
dt, %(dt)f
memory_per_node, %(memory_per_node[GB])d
""" % sd

    # shared global information that do not depends on process information
    L = []
    L.extend(gsi_strings)
    
    # system infomation
    L.append("n_total_process, ")
    L.append(str(sd["n_process"]))
    L.append("\n")

    L.append("n_OpenMP_threads, ")
    L.append(str(sd["n_OpenMP_threads"]))
    L.append("\n")

    
    
def make_region_info_strings_pre(sd, rd_array, process, process_IDs):
    L = []
    L.append("region_GID, ")
    L.append(str(process[i_proc].region_GID))
    L.append("\n")
    L.append("region_name, ")
    L.append(str(process[i_proc].region_name))
    L.append("\n")
    L.append("xy_length_per_tile, " +str(rd_array[region_GID]["structure_info"]["xy_length_per_tile"]))
    L.append("\n")
    
    L.append("n_same_region_process_GIDs, ")
    L.append(str(len(rd_array[process[i_proc].region_GID]["process_ID"])))
    L.append("\n")
    
    L.append("same_region_process_GIDs")
    for i in rd_array[process[i_proc].region_GID]["process_ID"]:
        L.append(", ")
        L.append(str(i))
    L.append("\n")

    # structure info
    L.append("n_subregions, ")
    L.append(str(rd_array[region_GID]["structure_info"]["n_subregions"]))
    L.append("\n")

    L.append("subregion_names")
    for i_sr, sr in enumerate(rd_array[region_GID]["neuron_info"][process[i_proc].region_name]):
        L.append(", ")
        L.append(str(*sr.keys()))
    L.append("\n")

    L.append("n_neuron_types_per_subregion")
    for i in rd_array[region_GID]["structure_info"]["n_neuron_types_per_subregion"]:
        L.append(", ")
        L.append(str(i))
    L.append("\n")
    
    L.append("n_neurons_per_neuron_type")
    for i in rd_array[region_GID]["structure_info"]["n_neurons_per_neuron_type"]:
        L.append(", ")
        L.append(str(i))
    L.append("\n")

    # neuron info
    L.append("neuron_model")
    for i in rd_array[region_GID]["structure_info"]["neuron_model_array"]:
        for j in i:
            L.append(", ")
            L.append(str(j))
    L.append("\n")
        
    L.append("E_or_I")
    for i in rd_array[region_GID]["structure_info"]["E_or_I_array"]:
        for j in i:
            L.append(", ")
            L.append(str(j))
    L.append("\n")
    
    L.append("membrane_time_constant")
    for i in rd_array[region_GID]["structure_info"]["membrane_time_constant_array"]:
        for j in i:
            L.append(", ")
            L.append(str(j))
    L.append("\n")
    
    L.append("spike_threshold")
    for i in rd_array[region_GID]["structure_info"]["spike_threshold_array"]:
        for j in i:
            L.append(", ")
            L.append(str(j))
    L.append("\n")
    
    L.append("reset_value")
    for i in rd_array[region_GID]["structure_info"]["reset_value_array"]:
        for j in i:
            L.append(", ")
            L.append(str(j))
    L.append("\n")
            
    L.append("E_rest")
    for i in rd_array[region_GID]["structure_info"]["E_rest_array"]:
        for j in i:
            L.append(", ")
            L.append(str(j))
    L.append("\n")
    
def write_process_files(sd, rd_array, process, process_IDs):
    L = []
    
    for i_proc in process_IDs:
        region_GID = process[i_proc].region_GID

        L.append("PRNG_seed, ")
        L.append(str(sd["PRNG_seed"]))
        L.append("\n")

        L.append("process_GID, ")
        L.append(str(process[i_proc].process_GID))
        L.append("\n")
        
        L.append("region_position")
        for i in process[i_proc].spatial_extent:
            L.append(", ")
            L.append(str(i))
        L.append("\n")
