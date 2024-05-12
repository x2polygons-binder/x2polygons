# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 23:32:18 2023

@author: banbar
"""

import os
import sys
import geopandas as gpd
import shapely

# Relative path of the OSM file

def read_data(site_id):   
    site_path = site_id + '/'
    
    script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
    
    data = dict()
    
    osm_path = os.path.join(script_dir, site_path, 'osm.shp')
    reference_path = os.path.join(script_dir, site_path, 'reference.shp')
    
    data['osm'] = gpd.read_file(osm_path)
    data['reference'] = gpd.read_file(reference_path)
    
    return data

def check_multipolygon(data):
    flag = 0
    # Check existence of multipolygons
    for p_source in data['osm'].geometry:
        # For OSM dataset
        if(isinstance(p_source, shapely.geometry.multipolygon.MultiPolygon)):
            flag = 1
            print("A multipolygon")
    
    for p_destination in data['reference'].geometry:
        # For reference dataset
        if(isinstance(p_destination, shapely.geometry.multipolygon.MultiPolygon)):
            flag = 1
            print("A multipolygon")
    
    if(flag==0):
        print("No multipolygons")
        
# check_multipolygon(data)



def find_intersecting_polygons(site_id, d1, d2, **kwargs):
    # Find the intesecting polygons
    # source, destination: source or destination
    #d1, d2: dataset1 and dataset2
    site_path = site_id + '/'
    
    script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
    
    out_file = script_dir + '/' + site_path + 'intersections.txt'
    out_path = os.path.join(script_dir, site_path, out_file)
    f = open(out_path, "w")
    if (d1 == 'osm'):
        # OSM has the additional feature - building type - include in the header
        header = "id " + str(d1) + " " + str(d2) + " " + "IoU " + "building_type \n" 
    else:
        header = "id " + str(d1) + " " + str(d2) + " " + "IoU " + " \n" 
    f.write(header)
    
    data = read_data(site_id)
    
    c= 1
    for p_d1 in range(len(data[d1])):
        for p_d2 in range(len(data[d2])):
            if(data[d1]['geometry'][p_d1].intersects(data[d2]['geometry'][p_d2])):
                if(kwargs["show_progress"]):
                    print("Intersection ", c)
                p1 = data[d1]['geometry'][p_d1]
                p2 = data[d2]['geometry'][p_d2]
                
                intersect = p1.intersection(p2).area
                union = p1.union(p2).area
                iou = intersect / union
                
                # Intersection ID, OSM gid, Reference ID, IoU
                if (d1 == 'osm'): # we have additional information that could be used - 'building' attribute
                    out_line = str(c) + ' ' + str(data[d1]['gid'][p_d1]) + ' ' + str(data[d2]['gid'][p_d2]) + ' ' + str(round(iou,2)) + ' ' + data[d1]['building'][p_d1] + '\n'
                else:
                    out_line = str(c) + ' ' + str(data[d1]['gid'][p_d1]) + ' ' + str(data[d2]['gid'][p_d2]) + ' ' + str(round(iou,2)) + '\n'
                f.write(out_line)
                
                c += 1
    
    f.close()


# Execute the iou threshold on the intersecting polygons
def filter_high_iou(input_file, iou_threshold):
    
    f_in = open(input_file, 'r')
    f_out = open('high_iou.txt', 'w')
    
    c = 1
    
    intersections = f_in.readlines()

    for i in intersections:
        iou = float(i.rstrip('\n').split(' ')[3])
        if(iou >= iou_threshold):
            new_line = str(c) + ' ' + i 
            f_out.write(new_line)
            c += 1
            print(i)
        
    f_in.close()
    f_out.close()


def identify_candidate_buildings(site_id):
    # Candidate buildings for a 1-1 match
    site_path = site_id + '/'
    
    script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
    
    input_path = script_dir + '/' + site_path + 'intersections.txt'
    f_in = open(input_path, 'r')
    header_in = f_in.readline()
    
    intersections = f_in.readlines()
    
    # Create an empty dictionary to store counts
    count_dict = dict()
    count_dict["dataset_1"] = {}
    count_dict["dataset_2"] = {}
    
    for i in intersections:
        
        building_id = i.rstrip('\n').split(' ')[1] # maybe <gid> would be a string? 
        if building_id in count_dict["dataset_1"]:
            count_dict["dataset_1"][building_id] += 1
        else:
            count_dict["dataset_1"][building_id] = 1
            
        building_id = i.rstrip('\n').split(' ')[2] # maybe <gid> would be a string? 
        if building_id in count_dict["dataset_2"]:
            count_dict["dataset_2"][building_id] += 1
        else:
            count_dict["dataset_2"][building_id] = 1
            
    # Convert the dictionary to a list of lists
    # result = [[item, count] for item, count in count_dict.items() if count==1]
    candidates = dict()
    
    candidates["dataset_1"] = [item for item, count in count_dict["dataset_1"].items() if count==1]
    candidates["dataset_2"] = [item for item, count in count_dict["dataset_2"].items() if count==1]
    
    f_in.close()
    return candidates
    

def identify_one_to_ones(site_id, candidates):    
    site_path = site_id + '/'    
    script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in    
    input_path = script_dir + '/' + site_path + "intersections.txt"
    f_in = open(input_path, "r")
    header_in = f_in.readline()
    
    out_path = script_dir + '/' + site_path + 'one_to_ones.txt'
    
    f_out = open(out_path, 'w')
    
    
    f_out.write(header_in)
    
    
    intersections = f_in.readlines()
    
    c = 0
    for i in intersections:
        splitted = i.rstrip('\n').split(' ')
        source_polygon = splitted[1]
        destination_polygon = splitted[2]
        if (source_polygon in candidates['dataset_1']) and (destination_polygon in candidates['dataset_2']):
            # save the line
            f_out.write(i)
            c+= 1
    
    print('Total 1-1s: ', c)
        
    f_in.close()
    f_out.close()
    

def export_one_to_ones(site_id, select_key):
    # Designed to work on QGIS
    # Assumption: select_key is common in both source & destination
    # does not work for geopackage
    
    site_path = site_id + '/'
    script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
    input_path =  script_dir + '/' + site_path + "one_to_ones.txt"
    f_in = open(input_path, "r")
    header_in = f_in.readline()
    

    data = dict()
    
    osm_path = os.path.join(script_dir, site_path, 'osm.shp')
    reference_path = os.path.join(script_dir, site_path, 'reference.shp')
    
    data['osm'] = gpd.read_file(osm_path)
    data['reference'] = gpd.read_file(reference_path)
    
    q_osm = ''
    q_reference = ''
    
    ones = f_in.readlines()
    for i in ones:
        splitted = i.rstrip('\n').split(' ')
        source_polygon = splitted[1]
        destination_polygon = splitted[2]
        
        q_osm = q_osm + "(data['osm']['gid'] == " + str(source_polygon) + ") | "
        q_reference = q_reference + "(data['reference']['gid'] == " + str(destination_polygon) + ") | "
        
    
    #delete the last | 
    q_osm = q_osm[:-2]
    q_reference = q_reference[:-2]
    
    osm_out = os.path.join(script_dir, site_path, 'matched_osm.shp')
    reference_out = os.path.join(script_dir, site_path, 'matched_reference.shp')
    
    data['osm'][eval(q_osm)].to_file(osm_out)
    data['reference'][eval(q_reference)].to_file(reference_out)
    
    
    f_in.close()
    

    
