# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 23:32:18 2023

@author: banbar
"""

import os
import sys
import geopandas as gpd
import shapely
from shapely.geometry import Point
from x2polygons import polygon_distance
from x2polygons import plot
from x2polygons import geometry
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def compute_distances(site_id, **kwargs):   
    # Header: results.txt
    # ID, OSM_ID, Ref_ID, Building_Type (from OSM) IoU, H, C, P, TF
    
    site_path = site_id + '/'
    script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
    
    
       
    data = dict()
    
    abs_osm_path = os.path.join(script_dir, site_path, 'matched_osm.shp')
    abs_reference_path = os.path.join(script_dir, site_path, 'matched_reference.shp')
    
    data['osm'] = gpd.read_file(abs_osm_path)
    data['reference'] = gpd.read_file(abs_reference_path)
    
    
    out_path = os.path.join(script_dir, site_path, "results.txt")
    
    c= 1
    source = 'osm'
    destination = 'reference'
    
    f_out = open(out_path, "w")
    header = "id " + "osm_gid " + "ref_gid " + "b_type " + " iou " + "hausdorff " + "chamfer " + "polis " + "tf \n" 
    f_out.write(header)
    
    for p_source in range(len(data[source])):
            for p_destination in range(len(data[destination])):
                if(data[source]['geometry'][p_source].intersects(data[destination]['geometry'][p_destination])):
                    p1 = data[source]['geometry'][p_source]
                    p2 = data[destination]['geometry'][p_destination]
                    
                    d_h = polygon_distance.hausdorff_distance(p1, p2, symmetrise = 'average')
                    d_c = polygon_distance.chamfer_distance(p1, p2, symmetrise = 'average')
                    d_p = polygon_distance.polis_distance(p1, p2, symmetrise = 'average')
                    d_tf = polygon_distance.turning_function_distance(p1, p2)
                    
    
                    
                    p_s_id = data[source]['gid'][p_source]
                    p_d_id = data[destination]['gid'][p_destination]
                    
                    
                    
                    intersect = p1.intersection(p2).area
                    union = p1.union(p2).area
                    iou = intersect / union
                    
                    out_case = str(c)+' '+str(p_s_id)+' '+str(p_d_id)+' '+data['osm']['building'][p_source]+' '+str(round(iou,2))+' '+str(round(d_h,2))+' '+str(round(d_c,2))+' '+str(round(d_p,2))+' '+str(round(d_tf,2))+'\n'
                    f_out.write(out_case)
                    
                    
                    if 'show_progress' in kwargs:
                        if(kwargs["show_progress"]):
                            print("Homomorph: ", c)
                    c+= 1
    f_out.close()


def obtain_distances(site_id):
    site_path = site_id + '/'
    script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
    input_path = os.path.join(script_dir, site_path, 'results.txt')
    
    f = open(input_path, "r")
    
    
    
    distance = dict()

    distance['hausdorff'] = []
    distance['polis'] = []
    distance['chamfer'] = []
    distance['turn_function'] = []
    
    header = f.readline()
    intersections = f.readlines()
    
    for i in intersections:
        record = i.rstrip('\n').split(' ')
        distance['hausdorff'].append(float(record[5]))
        distance['chamfer'].append(float(record[6]))
        distance['polis'].append(float(record[7]))
        distance['turn_function'].append(float(record[8]))
        
    
    f.close()
    return distance
    
    
def correlation_between_distance_measures(site_id, d1, d2, **kwargs):
    # Keyword argument: plot
    distance = obtain_distances(site_id)
    

    my_rho = np.corrcoef(distance[d1], distance[d2])
    
    
        
    if(kwargs["plot"] == True):
        print("Pearson correlation: ", my_rho)
        fig = plt.figure()
        ax = sns.regplot(x = distance[d1], y= distance[d2], marker='+',color='k')
    
        #ax.set(xlabel=, ylabel=distance_measure2.capitalize())
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        y_val_size = 14
        ax.tick_params(axis='both', which='major', labelsize = y_val_size)
        
        x_axis_label_size = 18
        y_axis_label_size = 18
        
        #yy = 'Turning Function'
        #plt.ylabel(yy, fontsize = y_axis_label_size)
        
        plt.xlabel(d1.capitalize(), fontsize = x_axis_label_size)
        plt.ylabel(d2.capitalize(), fontsize = y_axis_label_size)
    
        out_file_name = './' + site_id + '/' + 'correlation_' + d1 + '_' + d2 + '.svg'
        plt.tight_layout()
        plt.savefig(out_file_name, format='svg')
        
    return my_rho[0,1]

def correlation_between_distance_measures2(site_id, d1, d2, highlight_points=None, **kwargs):
    # Keyword argument: plot
    distance = obtain_distances(site_id)

    my_rho = np.corrcoef(distance[d1], distance[d2])
        
    if(kwargs.get("plot", False)):
        print("Pearson correlation: ", my_rho)
        fig = plt.figure()
        ax = sns.regplot(x=distance[d1], y=distance[d2], marker='+', color='k')
    
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        y_val_size = 14
        ax.tick_params(axis='both', which='major', labelsize=y_val_size)
        
        x_axis_label_size = 18
        y_axis_label_size = 18
        
        plt.xlabel(d1.capitalize(), fontsize=x_axis_label_size)
        plt.ylabel(d2.capitalize(), fontsize=y_axis_label_size)
        
        if highlight_points is not None:
            for point in highlight_points:
                x, y, label = point
                ax.scatter(x, y, color='red', zorder=5)
                ax.annotate(label, (x, y), textcoords="offset points", xytext=(10,10), ha='center', fontsize=10, color='red')
    
        out_file_name = f'./{site_id}/correlation_{d1}_{d2}.svg'
        plt.tight_layout()
        plt.savefig(out_file_name, format='svg')
        
    return my_rho[0,1]

def correlation_heatmap(site_id, corr):
    tick_labels = ['Hausdorff', 'Chamfer', 'PoLiS', 'Turn Function']
    d = {}
    fs = 14 #font size
    
    # Iterate over the rows of the correlation matrix
    for i, measure_name in enumerate(tick_labels):
        # Assign the correlation values from the i-th row to the respective distance measure name
        d[measure_name] = corr[i].tolist()
    
    df = pd.DataFrame(data=d)
    
    # Create the heatmap without annotations
    g = sns.heatmap(df, 
                    linewidths=1,
                    cmap='coolwarm'
                    )
    
    # Add annotations manually
    for i in range(len(df)):
        for j in range(len(df.columns)):
            text = "{:.2f}".format(df.iloc[i, j])  # Format the text to display two decimal places
            plt.text(j + 0.5, i + 0.5, text, ha='center', va='center', fontsize=fs)
    
    # Set tick labels and save the plot
    plt.xticks(ticks=[i + 0.5 for i in range(len(df.columns))], labels=df.columns, fontsize=fs)
    plt.yticks(ticks=[i + 0.5 for i in range(len(df))], labels=tick_labels, fontsize=fs-2)
    out_file_name = './' + site_id + '/' + 'correlation_heatmap.svg'
    plt.savefig(out_file_name)
    

def obtain_ranked_distances(site_id, distance_measures=None):
    # The last index of the results_ranked.txt is the max difference between ranks
    
    all_measures = ['hausdorff', 'chamfer', 'polis', 'turn_function']
    
    
    
    site_path = site_id + '/'
    script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
    input_path = os.path.join(script_dir, site_path, 'results.txt')
    
    out_path = os.path.join(script_dir, site_path, "results_ranked.txt")
    out_path_sorted = os.path.join(script_dir, site_path, "results_ranked_sorted.txt")
    
    f_out = open(out_path, "w")
    f_out_sorted = open(out_path_sorted, "w")
    
    f = open(input_path, "r")
    header = f.readline()
    header_out = "id\tosm_gid\tref_gid\tiou\tr_h\tr_c\tr_p\tr_tf\tr_difference"
    
    if distance_measures is None:
        header_out = header_out + '\n'
        f_out_sorted.write(header_out)
    else:
        header_out = header_out + '\t' + distance_measures[0][0:3] + '_' + distance_measures[1][0:3] + '\n'        
        f_out_sorted.write(header_out)
    
    distance = dict()

    distance["id"] = []

    distance["hausdorff"] = []
    distance["polis"] = []
    distance["chamfer"] = []
    distance["turn_function"] = []
    distance["iou"] = []

    intersections = f.readlines()
    
    for i in intersections:
        record = i.rstrip('\n').split(' ')
        distance["id"].append(int(record[0]))
        distance["hausdorff"].append(float(record[5]))
        distance["chamfer"].append(float(record[6]))
        distance["polis"].append(float(record[7]))
        distance["turn_function"].append(float(record[8]))
    
    # Sort and process distances - append the ranks for each distance measure
    
    distance["hausdorff"].sort()
    distance["chamfer"].sort()
    distance["polis"].sort()
    distance["turn_function"].sort()
    
    
    f.close()
    
    # Process the matches once more, and now include the ranks for each distance
    
    f = open(input_path, "r")
    header = f.readline()
     
    intersections = f.readlines()
    
    for i in intersections:
        record = i.rstrip('\n').split(' ')

        rh = distance["hausdorff"].index(float(record[5]))
        rc = distance["chamfer"].index(float(record[6]))
        rp = distance["polis"].index(float(record[7]))
        rtf = distance["turn_function"].index(float(record[8]))
        
        ranker = [rh, rc, rp, rtf]
        max_distance = 0

        for k1 in range(len(ranker)):
            for k2 in range(k1 + 1, len(ranker)):
                d = abs(ranker[k1] - ranker[k2])
                max_distance = max(max_distance, d)
        
        # If the user specified two distance measures
        if distance_measures is None:
            out_case = record[0] + ' ' + record[1] + ' ' + record[2] + ' ' + record[4] + ' ' + str(rh) + ' ' + str(rc) + ' ' + str(rp) + ' ' + str(rtf) + ' ' + str(max_distance) + '\n'
            f_out.write(out_case)
        else:
            first_rank = distance[distance_measures[0]].index(float(record[5 + all_measures.index(distance_measures[0])]))
            second_rank = distance[distance_measures[1]].index(float(record[5 + all_measures.index(distance_measures[1])]))
            chosen_difference = abs(first_rank - second_rank)
            out_case = record[0] + ' ' + record[1] + ' ' + record[2] + ' ' + record[4] + ' ' + str(rh) + ' ' + str(rc) + ' ' + str(rp) + ' ' + str(rtf) + ' ' + str(max_distance) + ' ' + str(chosen_difference) + '\n'
            f_out.write(out_case)
            
      
        
        
    f.close()
    f_out.close()
    
    # this time import as np array
    f_out = open(out_path, "r")
    
    if distance_measures is None:
        data = np.loadtxt(out_path, skiprows=0, dtype={'names': ('col1', 'col2', 'col3', 'col4', 'col5', 'col6', 'col7', 'col8', 'col9'),
                                                       'formats':(int, int, int, float, int, int, int, int, int)}) 
        #sorted_indices = np.argsort(data[:, 8])[::-1] #not working after the inclusion of specific data types / names
        sorted_indices = np.argsort(data, order='col9')[::-1]
        sd = data[sorted_indices] #sd: sorted_data
    else:
        data = np.loadtxt(out_path, skiprows=0, dtype={'names': ('col1', 'col2', 'col3', 'col4', 'col5', 'col6', 'col7', 'col8', 'col9', 'col10'),
                                                       'formats':(int, int, int, float, int, int, int, int, int, int)}) 
        #sorted_indices = np.argsort(data[:, 8])[::-1] #not working after the inclusion of specific data types / names
        sorted_indices = np.argsort(data, order='col10')[::-1]
        sd = data[sorted_indices] #sd: sorted_data
    
    c = 1
    for i in range(len(sd)):
        out_case = ''
        sd[i][0] = str(c)
        for j in range(len(sd[i])):
            out_case = out_case + str(sd[i][j]) + '\t' 
        out_case = out_case + '\n'
        c = c+1
        f_out_sorted.write(out_case)
    
    f_out_sorted.close()
    f_out.close()
    
    # Delete the results_ranked file - no further need
    os.remove(out_path)

def find_osm_ranked_row(site_id, gid):
    # Used in obtaining the centroids.shp
    # find_osm_ranked_row(site_id, 76) # finds the row with the OSM gid of 76, and retrieves its corresponding values
    site_path = site_id + '/'
    script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
    input_path = os.path.join(script_dir, site_path, 'results_ranked_sorted.txt')
    
    f = open(input_path, "r")
    header = f.readline()

    rows = f.readlines()
    for i in range(len(rows)):
        splitted = rows[i].rstrip('\n').split('\t')
        if (splitted[1] == str(gid)):
            f.close()
            return (splitted)
    
    print("The feature {} is not found!".format(str(gid)))
    f.close()
    

def heatmap_ranked_distances(site_id):
    site_path = site_id + '/'
    script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
    input_path = os.path.join(script_dir, site_path, 'results_ranked_sorted.txt')
    
    f = open(input_path, "r")
    header = f.readline()
    header = header.rstrip().split('\t')
    
    osm_path = os.path.join(script_dir, site_path, 'matched_osm.shp')
    data = dict()
    data['osm'] = gpd.read_file(osm_path)
    
    data['centroids'] = dict()
    data['centroids']['geom'] = []
    data['centroids']['building'] = []
    data['centroids'][header[1]] = [] #osm_gid
    data['centroids'][header[2]] = [] #ref_gid
    data['centroids'][header[3]] = [] #iou
    data['centroids'][header[4]] = [] #r_h
    data['centroids'][header[5]] = [] #r_c
    data['centroids'][header[6]] = [] #r_p
    data['centroids'][header[7]] = [] #r_tf
    data['centroids'][header[8]] = [] #r_max
    if(len(header) == 10):
        # the rank file also contains a selected pair of distance measures
        data['centroids'][header[9]] = [] #rank difference between chosen measures
    
    centroids_out = os.path.join(script_dir, site_path, 'matching_centroids.shp')

    
    for i in range(len(data['osm'])):
        # Obtain the centroid
        data['centroids']['geom'].append(data['osm']['geometry'][i].centroid)
        data['centroids']['building'].append(data['osm']['building'][i])    
        # Find the ranks
        data['centroids'][header[1]].append(data['osm']['gid'][i])
        matching_row = find_osm_ranked_row(site_id, data['osm']['gid'][i])
        data['centroids'][header[2]].append(matching_row[2])
        data['centroids'][header[3]].append(float(matching_row[3]))
        data['centroids'][header[4]].append(int(matching_row[4]))
        data['centroids'][header[5]].append(int(matching_row[5]))
        data['centroids'][header[6]].append(int(matching_row[6]))
        data['centroids'][header[7]].append(int(matching_row[7]))
        data['centroids'][header[8]].append(int(matching_row[8]))
        if(len(header) == 10): 
            data['centroids'][header[9]].append(int(matching_row[9]))
        
    
    # Retrieve the CRS from the OSM data
    osm_crs = data['osm']['geometry'].crs

    # Create a GeoDataFrame from the dictionary with the CRS from OSM data
    gdf = gpd.GeoDataFrame(data['centroids'], geometry='geom', crs=osm_crs)

    # Export
    gdf.to_file(centroids_out)  

    
    
    

def visualise_matching_polygons(site_id, id_pa, id_pb, **kwargs):
    site_path = site_id + '/'
    script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
    
    data = dict()
    
    file_name = str(id_pa) + '_' + str(id_pb) + '.emf'
    
    abs_osm_path = os.path.join(script_dir, site_path, 'matched_osm.shp')
    abs_reference_path = os.path.join(script_dir, site_path, 'matched_reference.shp')
    file_path = os.path.join(script_dir, site_path, file_name)
    
    data['osm'] = gpd.read_file(abs_osm_path)
    data['reference'] = gpd.read_file(abs_reference_path)
    
    for i in range(len(data["osm"])):
        if (data["osm"]["gid"][i] == id_pa):
                for j in range(len(data["reference"])):
                    if (data["reference"]["gid"][j] == id_pb):
                        poly_a = data["osm"]["geometry"][i]
                        poly_b = data["reference"]["geometry"][j]
                        
                        try:
                            if(kwargs["export_as_emf"]):
                                plot.plot_x2polygons(poly_a, poly_b, file_path = file_path)
                        except:
                            plot.plot_x2polygons(poly_a, poly_b)
                            
                
    
    
    
    
    


    
    
    
    

