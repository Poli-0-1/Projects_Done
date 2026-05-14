

#### LSTM-only code for the thesis project on Model developing for  
#### classification of differentiating and non-differentiating hiPSCs cells in cardiomyocytes


### Andrea Policano


# final update in 05/12/2026




##########################################################################################################################################

######################### Importing multiple libraries



from random import choices
import re
import time
from ast import List, arg
from doctest import debug
from email import message
from math import isnan
from shutil import which
import time
from tkinter import font
from tracemalloc import start
from turtle import width

from anyio import open_file
import numpy as np
import matplotlib.pyplot as plt
import os
from skimage import io

from scipy.interpolate import make_interp_spline

os.environ["TF_ENABLE_ONEDNN_OPTS"] = "0"  # Disable oneDNN optimizations to avoid issues with certain operations


import PIL
import skimage
from sklearn.semi_supervised import LabelSpreading
import tensorflow as tf
import pathlib
import cv2
import tifffile as tiff
import seaborn as sns
import mahotas
import tqdm
import InquirerPy
import yaspin
import rich
import pandas as pd


from mahotas import median_filter
from scipy.ndimage import binary_fill_holes, gaussian_filter
from termcolor import colored
from rich.progress import SpinnerColumn, Progress, TextColumn
from rich.spinner import Spinner, SPINNERS
from yaspin import yaspin
from tqdm import tqdm
from InquirerPy import inquirer, prompt
from skimage.segmentation import watershed
from scipy import ndimage  
from skimage import morphology, measure, segmentation
from gc import callbacks
from genericpath import isdir
from tabnanny import verbose
from matplotlib import units
from skimage.measure import label, regionprops, shannon_entropy
from skimage.filters import threshold_otsu
from os import listdir
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, Embedding, GRU, Bidirectional
from tensorflow.keras.preprocessing import sequence
from tensorflow.keras.metrics import AUC
from sklearn.model_selection import permutation_test_score, train_test_split
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import f1_score
from tensorflow.keras.losses import BinaryCrossentropy
from tensorflow.keras.callbacks import EarlyStopping
from torch import dropout, gru
from sklearn.metrics import confusion_matrix, roc_curve
from skimage.morphology import local_maxima











##########################################################################################################################################
##########################################################################################################################################

    ######################### Explain to the user the requirements before using the algorithm 
    
start_time = time.time()



# Print the introduction to the algorithm and the requirements to use it
print("\n####################################################################################")
print("########################################################\n")


print(colored("\nIs this your first time using this tool?", 'green', attrs = ['bold']))
first_time = inquirer.select(
    message = "\n",
    choices = ["Yes", "No"],
    default = "Yes"
).execute()


if first_time == "Yes":
    print("\nGood Day to you.    \U0001f44b\U0001f62f\n\nThe following program is a tool with the main objective of detecting and assert successfull differentiation events of hiPSCs-CMs taking place in 2D, 2.5D and 3D cultures.\n\nThe algorithm exploits time-series of gray-scale fluorescence images in order to derive data and information from each single image, later implementing the resulting info by \x1B[3mfeeding\x1B[0m them to a Long Short-Term Memory Neural Network.\n")
    print("\nOnce the model is trained, it will be able to classify the samples in two main categories:\n\n \t \u2022 Differentiating Samples  \u26ab  \u27f6  \u2714\n\n \t \u2022 Non-Differentiating Samples  \u26ab  \u27f6  \u274c\n")
    print("\nIn the end the algorithm will be able to provide a final assertion regarding the succes of the whole differentiation procedure, based on the classification of the samples and their respective time-series of images.\n")
    print("\nThe following code is a work in progress, so please be patient and report any bugs or issues you may encounter to the author of the code: \n \t \u27a1  \x1B[3m(Andrea Policano).\x1B[0m\n\nThank you for your collaboration, I wish you a good analysis and a good day!\n")
    print("\nPlease, before the tool is put to work check the following constraints:\n- Presence of marker name in the .TIFF format images \n- Equal amount of images in each subfolder\n- When asked to give the marker order and names please write them exactly as present in the image filename\n\n\n####################################################################################\n")


print("\n\n\n")




#User_quest = inquirer.confirm(
#    message = "\nAre you ready to start the analysis considering the previous advices?\n>", default = True).execute()     # Ask the user if he is ready to start the analysis, if not then the program will stop and ask to check the requirements

print(colored("\nAre you ready to start the analysis ?", 'green', attrs = ['blink']))

User_quest = inquirer.select(
    message = "\n",
    choices = ["Yes", "No"],
    default = "Yes"
).execute()  # Ask the user if he is ready to start the analysis, if not then the program will stop and ask to check the requirements





if User_quest == "Yes":  # If the user is ready to start the analysis then we can go on:


    ######################### Specify the markers present in the analysis and their order 
    
    markers_std = input(colored("\n\nPlease specify the markers and order (comma separated):", 'green', attrs = ['bold']) + "\n>") 
    markers_order = [x.strip().lower() for x in markers_std.split(",")]
    
    print("\n\n\n")
    
    
    
    print(colored("\nDo you have markers that can be considered as aliases of others?", 'green', attrs = ['bold']))
    add_alias = inquirer.select(
        message = "\n",
        choices = ["Yes", "No"],
        default = "Yes"
        ).execute()       # Check for markers that have tha same purpose and are substituted to others troughout te samples 
    
    
    marker_aliases = {}
    
    
    if add_alias == "Yes":
        input_aliases = input(colored("\nPlease, provide the markers (--> Separated, original --> alias):\t   (Use Comma in case of multiple markers)", 'green', attrs = ['bold']) + "\n>")      # Specify the markers in question to later consider them as the "same" (at least at a purpose level)
        input_aliases = input_aliases.strip().lower().split(",")                                ###################
                                                                                                                  #  Split the multiple aliases for each marker and check if there are more
        for aliases in input_aliases:                                                                             #
            sections = [s.strip() for s in aliases.split("-->") if s.strip()]                                     #
                                                                                                                  #  Put the associated markers to the same list, to later use in the parsing loop 
            if len(sections) >= 2:                                                                                #
                marker = sections[0]                                                                      #
                single_alias = sections[1:]                           #
                marker_aliases[marker] = single_alias                                           ###################

    else:
        print("\nPerfect, analysis will continue")        # If no markers for the same purpose are present then we can go on 
    
    
    print("\n\n\n")
    
    
    nuclear_input = input(colored("\nPlease identify the nuclear environment: \t(specify the marker associated to the nuclear environment", 'green', attrs = ['bold']) + "\n>")       
    
    
    
    single_marker = input(colored("\nAre there specific markers you would like to solely work on?", 'green', attrs = ['bold']) + "\n>")
    
    
    if single_marker in ["yes", "Yes", "Y", "y"]:
        print("\nPlease specify the marker :")
        sole_markers = inquirer.checkbox(
            message = "\n",
            choices = markers_order,
            validate = lambda result: len(result) > 0,
            invalid_message = "Please select at least one marker."
        ).execute()
        
        markers_order = sole_markers

    
    
    ### Input data import and specification of the usable extensions
    data_folder = input(colored("\n\nInsert Path for the folder containing the data:", 'green', attrs = ['bold']) + "\n>")
    print("\n\n\n")
    
    
    overlap_check = input(colored("\nActivating overlap", 'green', attrs = ['bold']) + "\n>")       # If activated it will produce the overlap between the original image and the segmented bdies found 
    
    
    


    
    
    
    
    
########################################################################################################################################################################
########################################################################################################################################################################
    
    
    
    ######################### Importing the time-series of images
    
    
    time_series = []            # List that will be filled with the arrays of the images in their respective temporal order 
    time_series_meta = []       # For metadata purposes and histogram build, just to correctly keep track in plots visualization

    

    ###Reading in an iterative way using the path as reading index
    for sub_folder in sorted(os.listdir(data_folder), key = lambda x: int(x[1:].split("_")[0])):     # Sort by considering the format of the file name, consider "D" and later on check for the int value of the timepoint, i.e. D7 --> 7, D18 --> 18, etc., this is done to have a correct sorting of the folders and images, otherwise t would cause problems at a lexicographic level
        subfold_pathway = os.path.join(data_folder, sub_folder)

        if not os.path.isdir(subfold_pathway):
            continue 
        
        
        for a in sorted(os.listdir(subfold_pathway)):       # Sort aplhabetically the images, quite important to do
            sample_path = os.path.join(subfold_pathway, a)
            
            if not os.path.isdir(sample_path):              # Iterate through folders and subfolders until we supposedly reach the images 
                continue
            
            sample_folder = [None] * len(markers_order)     # Produce temporary folders to store variable and later work on them
            sample_meta = [None] * len(markers_order)
            
            
            for images in sorted(os.listdir(sample_path)):      # Divide name and extensions, in order to use them later and check if everything seems to be in order and ready for extraction 
                
                _, extensions = os.path.splitext(images)
                
                if extensions.lower() in [".tif", ".tiff"]:
                    img_path = os.path.join(sample_path, images)
                    img = tiff.imread(img_path)                     # "Read" the image using a tiff-specific package
                    img_name = images.lower()
                    
                    
                    
                    
                    

                    for index, marker in enumerate(markers_order):      # Check for the images names and control whether or not they are in the marker aliases
                        if marker in img_name:
                            sample_folder[index] = img
                            sample_meta[index] = (sub_folder, marker)
                            break
                        
                        
                                                        # Checck if the markers with the same function are present and add them anyway without overlapping errors
                        
                        elif marker in marker_aliases:
                            for alias in marker_aliases[marker]:
                                if alias in img_name:
                                    sample_folder[index] = img
                                    sample_meta[index] = (sub_folder, marker)
                                    break
                                    
                 
                        
            if all(i is not None for i in sample_folder):       # Finally append the images's arrays and "metadata" to their respective list
                time_series.append(sample_folder)    
                time_series_meta.append(sample_meta)
                
                print(f"\n\nLoading Sequence: {len(time_series)} --> {sub_folder}")      # Print the timepoint and the marker associated to the image, so as to keep track of the images loaded
                print("\n\tMarker and Shape of the image:")
                
                for m_img, (_, id_marker) in zip(sample_folder, sample_meta):

                    img_shape = img.shape                      # Check images shape in order to see if the dimensions are correct
                    print(f"\t\t>{id_marker}: {img_shape}")
                

    print(f"\n###########################################################################\n\n\nLoaded TIFF Images\t {len(time_series)} sequences (time_points x markers)\n\n")     # Print the time_series variable's structure and its length (which is the sum of Sample across all subfolders)
    print(time_series)
    print("\n\n\n\n")










########################################################################################################################################################################
########################################################################################################################################################################



    ######################### Plot Additional Statistics related to image quality (Luminosity, Blur, Contrast)
    
    
    print("\n\n\n###########################################################################")
    
    print(colored("\nWould you like to see all the statistics related to the image quality? (if no only partial infos will be shown)", 'green', attrs = ['bold']))
    
    IQA_input = inquirer.select(
        message = "\n",
        choices = ["Yes", "No"],
        default = "Yes").execute()      # Ask the user whether or not he wants to see the statistics related to the image quality
    
    
    if IQA_input in ["Y", "y", "YES", "yes", "Yes"]:
        
        
        List_Blur  = []  # List to store the blur values for each image in the sequence
        List_Contrast = []  # List to store the contrast values for each image in the sequence
        List_Luminosity = []  # List to store the luminosity values for each image in the sequence
        List_Edges = []  # List to store the edge density values for each image in the sequence
        List_Tenegrad = []  # List to store the Tenengrad values for each image in the sequence
        
        
        for index, (single_sequence, related_metadata) in enumerate(tqdm(zip(time_series, time_series_meta),
                                                                         desc = "Importing Images",
                                                                         total = len(time_series),
                                                                         colour = "#900000",
                                                                         dynamic_ncols = True)):            # Loop through all the images and save single element and their indexes
            for image_index, (array, (time_point, marker)) in enumerate(zip(single_sequence, related_metadata)):
                
                img_array = array
                
                Blur = cv2.Laplacian(img_array, cv2.CV_8U).var() # Calculate the Laplacian to measure blur
                List_Blur.append(Blur)
                
                
                Contrast = shannon_entropy(img_array)  # Calculate Shannon Entropy for contrast
                List_Contrast.append(Contrast)
                
                        
                Luminosity = np.mean(img_array)  # Calculate luminosity 
                List_Luminosity.append(Luminosity)
                
                
                Edges = cv2.Canny(img_array.astype(np.uint8), 50, 150)  # Apply Canny edge detection
                N_Edges = np.sum(Edges > 0)             # Count the number of edge pixels(those with value > 0) obtained from the previous function
                Edges_Ratio = N_Edges / Edges.size   # Calculate the ratio of edge pixels to the total number of pixels in the image
                List_Edges.append(Edges_Ratio)
                
                
                G_x = cv2.Sobel(img_array, cv2.CV_64F, 1, 0, ksize = 3)   # Sobel operator in the X direction
                G_y = cv2.Sobel(img_array, cv2.CV_64F, 0, 1, ksize = 3)   # Sobel operator in the Y direction
                Grad_Magnitude = np.sqrt(G_x**2 + G_y**2)          # Gradient magnitude
                Tenegrad_val = np.mean(Grad_Magnitude)
                List_Tenegrad.append(Tenegrad_val)
                
                
                
                
        min_blur = np.min(List_Blur) if List_Blur else 0
        max_blur = np.max(List_Blur) if List_Blur else 0
        
        min_luminosity = np.min(List_Luminosity) if List_Luminosity else 0
        max_luminosity = np.max(List_Luminosity) if List_Luminosity else 0
        
        min_contrast = np.min(List_Contrast) if List_Contrast else 0
        max_contrast = np.max(List_Contrast) if List_Contrast else 0
        
        min_edges = np.min(List_Edges) if List_Edges else 0
        max_edges = np.max(List_Edges) if List_Edges else 0
        
        min_tenegrad = np.min(List_Tenegrad) if List_Tenegrad else 0
        max_tenegrad = np.max(List_Tenegrad) if List_Tenegrad else 0
        
        

        # Verify mapping of metrics to images 
        
        rows = []
        metric_idx = 0

        for seq_idx, (img_list, meta_list) in enumerate(zip(time_series, time_series_meta), start=1):
            for img_idx, (img, (timepoint, marker)) in enumerate(zip(img_list, meta_list), start=1):
                shape = img.shape
                blur = List_Blur[metric_idx] if metric_idx < len(List_Blur) else None
                lumin = List_Luminosity[metric_idx] if metric_idx < len(List_Luminosity) else None
                contrast = List_Contrast[metric_idx] if metric_idx < len(List_Contrast) else None
                edges = List_Edges[metric_idx] if metric_idx < len(List_Edges) else None
                tenegrad = List_Tenegrad[metric_idx] if metric_idx < len(List_Tenegrad) else None
                
                rows.append((seq_idx, timepoint, marker, img_idx, shape, blur, lumin, contrast, edges, tenegrad))
                metric_idx += 1



        print("\n\n\n")
        print("Seq | Timepoint                | Marker | Img# | Shape       | Blur     | Lumin    | Contrast    | Edges_Density    | Sharpness")
        print("----+--------------------------+--------+------+-------------+----------+----------+-------------+------------------+----------")
        for r in rows:
            print(f"{r[0]:>3} | {r[1]:24} | {r[2]:6} | {r[3]:>3} | {r[4]} | {r[5]:8.3f} | {r[6]:8.3f} | {r[7]:8.3f} | {r[8]:8.3f} | {r[9]:8.3f}")


        
        print(f"\nImage Quality Analysis Results:\n\n\u2022 List of Blur Values: \n\n{List_Blur}\n\n\u2022 List of Luminosity Values: \n\n{List_Luminosity}\n\n\u2022 List of Contrast Values: \n\n{List_Contrast}\n\n\u2022 List of Edge Density Values: \n\n{List_Edges}\n\n\u2022 List of Sharpness Values: \n\n{List_Tenegrad}\n") 
        print(f"\n \t \u2022 Minimum Blur Value: {min_blur}\n \t \u2022 Maximum Blur Value: {max_blur}\n")
        print(f"\n \t \u2022 Minimum Luminosity Value: {min_luminosity}\n \t \u2022 Maximum Luminosity Value: {max_luminosity}\n")        
        print(f"\n \t \u2022 Minimum Contrast Value: {min_contrast}\n \t \u2022 Maximum Contrast Value: {max_contrast}\n")
        print(f"\n \t \u2022 Minimum Edge Density Value: {min_edges}\n \t \u2022 Maximum Edge Density Value: {max_edges}\n")
        print(f"\n \t \u2022 Minimum Sharpness Value: {min_tenegrad}\n \t \u2022 Maximum Sharpness Value: {max_tenegrad}\n")



    
        
        
        
        fig, axs = plt.subplots(5, 1, figsize = (10, 15), sharex = True)             # Create subplots for each metric, sharing the x-axis for better comparison
        fig.subplots_adjust(hspace = 0.4, bottom = 0.15)
        
        
        axs[0].plot(List_Blur, label = "Blur", marker = "p", color = "blue")          # Plot the blur values
        axs[0].set_ylabel("Blur", fontsize = 10)
        axs[0].grid(True, color = "#EEEEEE", linestyle = "--", linewidth = 0.5)

        
        axs[1].plot(List_Luminosity, label = "Luminosity", marker = "p", color = "orange")  # Plot the luminosity values  
        axs[1].set_ylabel("Luminosity", fontsize = 10)
        axs[1].grid(True, color = "#EEEEEE", linestyle = "--", linewidth = 0.5)
        
        axs[2].plot(List_Contrast, label = 'Contrast', color = 'green', marker = "p")  # Plot the contrast values
        axs[2].set_ylabel("Contrast", fontsize = 10)
        axs[2].grid(True, color = "#EEEEEE", linestyle = "--", linewidth = 0.5)
        
        axs[3].plot(List_Edges, label = 'Edges', marker = 'p', color = 'purple') # Plot the edge density values
        axs[3].set_ylabel("Edge Density", fontsize = 10)
        axs[3].grid(True, color = "#EEEEEE", linestyle = "--", linewidth = 0.5)
        
        
        axs[4].plot(List_Tenegrad, label = 'Sharpness', marker = 'p', color = 'red') # Plot the Tenengrad values
        axs[4].set_ylabel("Sharpness", fontsize = 10)
        axs[4].grid(True, color = "#EEEEEE", linestyle = "--", linewidth = 0.5)
        axs[4].set_xlabel("Image Index", fontsize = 10)
        
        plt.tight_layout(rect = [0, 0.03, 1, 0.95])  # Adjust layout to make room for the title
        plt.show()
        
        
        
        

        ## Plot the IQA metrics with better visualization and grids 
        
        plt.plot(List_Blur, label = "Blur", marker = "p", color = "blue")          # Plot the blur values
        plt.xlabel("Image Index", fontsize = 10)
        plt.xticks(np.arange(len(List_Blur)), [f"{i+1}" for i in range(len(List_Blur))])  # Set x-ticks to be the index of the images
        plt.ylabel("Blur Values", fontsize = 10)
        plt.grid(True, color = "#EEEEEE", linestyle = "--", linewidth = 0.5)  # Add grid for better readability
        plt.title("IQA - Images Blurring", fontsize = 12)
        plt.show()
        
        
        
        plt.plot(List_Luminosity, label = "Luminosity", marker = "p", color = "orange")  # Plot the luminosity values
        plt.xlabel("Image Index", fontsize = 10)
        plt.ylabel("Luminosity Values", fontsize = 10)
        plt.xticks(np.arange(len(List_Luminosity)), [f"{i+1}" for i in range(len(List_Luminosity))])
        plt.grid(True, color = "#EEEEEE", linestyle = "--", linewidth = 0.5)  
        plt.title("IQA - Images Luminosity", fontsize = 12)
        plt.show()
        
        
        
        plt.plot(List_Contrast, label = 'Contrast', color = 'green', marker = "p")  # Plot the contrast values
        plt.xlabel("Image Index", fontsize = 10)
        plt.ylabel("Contrast Values", fontsize = 10)
        plt.xticks(np.arange(len(List_Contrast)), [f"{i + 1}" for i in range(len(List_Contrast))])
        plt.grid(True, color = "#EEEEEE", linestyle = "--", linewidth = 0.5)
        plt.title("IQA - Images Contrast", fontsize = 12)
        plt.show()
        
        
        
        
        
        plt.plot(List_Edges, label = 'Edges', marker = 'p', color = 'purple') # Plot the edge density values
        plt.xlabel("Image Index", fontsize = 10)
        plt.ylabel("Edge Density Values", fontsize = 10)
        plt.xticks(np.arange(len(List_Edges)), [f"{i + 1}" for i in range(len(List_Edges))])
        plt.grid(True, color = "#EEEEEE", linestyle = "--", linewidth = 0.5)
        plt.title("IQA - Images Edge Density", fontsize = 12)
        plt.show()
        
        
        
        plt.plot(List_Tenegrad, label = 'Sharpness', marker = 'p', color = 'red') # Plot the Tenengrad values
        plt.xlabel("Image Index", fontsize = 10)
        plt.ylabel("Sharpness Values", fontsize = 10)
        plt.xticks(np.arange(len(List_Tenegrad)), [f"{i + 1}" for i in range(len(List_Tenegrad))])
        plt.grid(True, color = "#EEEEEE", linestyle = "--", linewidth = 0.5)
        plt.title("IQA - Images Sharpness", fontsize = 12)
        plt.show()
        
        
        

        # Upper and Lower bounds for the Metrics related to Image Quality Analysis (IQA)
        # The values obtained are derived from the following calculations applied on the control images
        # The values were then taken and inserted as "representatives" of the acceptable ranges for our expectations
        
        '''
        min_bound_Blur = min_blur * 0.9
        max_bound_Blur = max_blur * 1.1
        
        min_bound_Lum = min_luminosity * 0.9
        max_bound_Lum = max_luminosity * 1.1
        
        min_bound_Contr = min_contrast * 0.9
        max_bound_Contr = max_contrast * 1.1
        
        min_bound_Edg = min_edges * 0.9
        max_bound_Edg = max_edges * 1.1
        
        
        min_bound_Teneg = min_tenegrad * 0.9
        max_bound_Teneg = max_tenegrad * 1.1
        
        '''
        
        
        min_bound_Blur = 1.54576
        max_bound_Blur = 279.84862
        
        min_bound_Lum = 4.05424
        max_bound_Lum = 200.0
        
        min_bound_Contr = 1.10639
        max_bound_Contr = 7.86189
        
        min_bound_Edg = 1.5419407894736843e-05
        max_bound_Edg = 0.13861932754516604
        
        
        min_bound_Sharp = 2.0002292961445325
        max_bound_Sharp = 50.69124668147343
        
        
        print("\n\n[----------------------------------------------------------------------------------]")
        print("[----------------------------------------------------------------------------------]\n\n")
        
        
        
        print(f"\n \t \u2022 Minimum Bound for Blur: {min_bound_Blur}\n \t \u2022 Maximum Bound for Blur: {max_bound_Blur}\n")
        print(f"\n \t \u2022 Minimum Bound for Luminosity: {min_bound_Lum}\n \t \u2022 Maximum Bound for Luminosity: {max_bound_Lum}\n")
        print(f"\n \t \u2022 Minimum Bound for Contrast: {min_bound_Contr}\n \t \u2022 Maximum Bound for Contrast: {max_bound_Contr}\n")
        print(f"\n \t \u2022 Minimum Bound for Edges: {min_bound_Edg}\n \t \u2022 Maximum Bound for Edges: {max_bound_Edg}\n")
        print(f"\n \t \u2022 Minimum Bound for Sharpness: {min_bound_Sharp}\n \t \u2022 Maximum Bound for Sharpness: {max_bound_Sharp}\n")
        
        
        
        import csv
        
        filename = "IQA_Values.csv"
        
        rows = zip(List_Blur, List_Luminosity, List_Contrast, List_Edges, List_Tenegrad)
        
        with open(filename, 'w', newline = '') as f:
            
            writer = csv.writer(f)
            
            writer.writerow(["Blur", "Luminosity", "Contrast", "Edges", "Sharpness"])
            writer.writerows(rows)
            
            print(f"\n\nValues saved to {filename}\n in {data_folder}")
        
        
        
        def check_bounds(x, y, z, g, h):
            
            flag_blur = []
            flag_contrast = []
            flag_luminosity = []
            flag_edges = []
            flag_sharp = []
        
            
            
        def check_bounds(List_Blur, List_Luminosity, List_Contrast, List_Edges, List_Tenegrad):
            
            flag_blur = []
            flag_contrast = []
            flag_luminosity = []
            flag_edges = []
            flag_sharp = []

            # Find indices where values are out of bounds
            for idx, value in enumerate(List_Contrast):
                if value < min_bound_Contr or value > max_bound_Contr:
                    flag_contrast.append(idx)

            for idx, value in enumerate(List_Luminosity):
                if value < min_bound_Lum or value > max_bound_Lum:
                    flag_luminosity.append(idx)

            for idx, value in enumerate(List_Blur):
                if value < min_bound_Blur or value > max_bound_Blur:
                    flag_blur.append(idx)

            for idx, value in enumerate(List_Edges):
                if value < min_bound_Edg or value > max_bound_Edg:
                    flag_edges.append(idx)

            for idx, value in enumerate(List_Tenegrad):
                if value < min_bound_Sharp or value > max_bound_Sharp:
                    flag_sharp.append(idx)

            # Check if any flags are not empty
            if flag_contrast or flag_luminosity or flag_blur or flag_edges or flag_sharp:
                print("\n \x1B[3mWarning\x1B[0m: Some values are out of the expected bounds, please check the images and their quality before proceeding with the analysis")
                print(f"\n \t \u2022 Contrast Values out of bounds: {flag_contrast}")
                print(f"\t \u2022 Luminosity Values out of bounds: {flag_luminosity}")
                print(f"\t \u2022 Blur Values out of bounds: {flag_blur}")
                print(f"\t \u2022 Edges Values out of bounds: {flag_edges}")
                print(f"\t \u2022 Sharpness Values out of bounds: {flag_sharp}\n")

            return flag_contrast, flag_luminosity, flag_blur, flag_edges, flag_sharp


        flag_contrast, flag_luminosity, flag_blur, flag_edges, flag_sharp = check_bounds(List_Blur, List_Luminosity, List_Contrast, List_Edges, List_Tenegrad)  # Check if the values are out of bounds
        
     

    
    
    else:
                
        List_Blur  = []  # List to store the blur values for each image in the sequence
        List_Contrast = []  # List to store the contrast values for each image in the sequence
        List_Luminosity = []  # List to store the luminosity values for each image in the sequence
        List_Edges = []  # List to store the edge density values for each image in the sequence
        List_Tenegrad = []  # List to store the Tenengrad values for each image in the sequence        
        
        
        for index, (single_sequence, related_metadata) in enumerate(tqdm(zip(time_series, time_series_meta),
                                                                         desc = "Importing Images",
                                                                         colour = "#900000",
                                                                         total = len(time_series),
                                                                         dynamic_ncols = True)):            # Loop through all the images and save single element and their indexes
            for image_index, (array, (time_point, marker)) in enumerate(zip(single_sequence, related_metadata)):
                
                img_array = array
                
                Blur = cv2.Laplacian(img_array, cv2.CV_8U).var() # Calculate the Laplacian to measure blur
                List_Blur.append(Blur)
                
                Contrast = shannon_entropy(img_array)  # Calculate Shannon Entropy for contrast
                List_Contrast.append(Contrast)
                        
                Luminosity = np.mean(img_array)  # Calculate luminosity 
                List_Luminosity.append(Luminosity)
                
                Edges = cv2.Canny(img_array.astype(np.uint8), 50, 150)  # Apply Canny edge detection
                N_Edges = np.sum(Edges > 0)             # Count the number of edge pixels(those with value > 0) obtained from the previous function
                Edges_Ratio = N_Edges / Edges.size   # Calculate the ratio of edge pixels to the total number of pixels in the image
                List_Edges.append(Edges_Ratio)
                
                
                G_x = cv2.Sobel(img_array, cv2.CV_64F, 1, 0, ksize = 3)   # Sobel operator in the X direction
                G_y = cv2.Sobel(img_array, cv2.CV_64F, 0, 1, ksize = 3)   # Sobel operator in the Y direction
                Grad_Magnitude = np.sqrt(G_x**2 + G_y**2)          # Gradient magnitude
                Tenegrad_val = np.mean(Grad_Magnitude)
                List_Tenegrad.append(Tenegrad_val)
                
                
                
        min_blur = np.min(List_Blur) if List_Blur else 0
        max_blur = np.max(List_Blur) if List_Blur else 0
        
        min_luminosity = np.min(List_Luminosity) if List_Luminosity else 0
        max_luminosity = np.max(List_Luminosity) if List_Luminosity else 0
        
        min_contrast = np.min(List_Contrast) if List_Contrast else 0
        max_contrast = np.max(List_Contrast) if List_Contrast else 0
        
        min_edges = np.min(List_Edges) if List_Edges else 0
        max_edges = np.max(List_Edges) if List_Edges else 0
        
        min_tenegrad = np.min(List_Tenegrad) if List_Tenegrad else 0
        max_tenegrad = np.max(List_Tenegrad) if List_Tenegrad else 0
        
        
        
        # === Verify mapping of metrics to images ===
        rows = []
        metric_idx = 0

        for seq_idx, (img_list, meta_list) in enumerate(zip(time_series, time_series_meta), start=1):
            for img_idx, (img, (timepoint, marker)) in enumerate(zip(img_list, meta_list), start=1):
                shape = img.shape
                blur = List_Blur[metric_idx] if metric_idx < len(List_Blur) else None
                lumin = List_Luminosity[metric_idx] if metric_idx < len(List_Luminosity) else None
                contrast = List_Contrast[metric_idx] if metric_idx < len(List_Contrast) else None
                edges = List_Edges[metric_idx] if metric_idx < len(List_Edges) else None
                tenegrad = List_Tenegrad[metric_idx] if metric_idx < len(List_Tenegrad) else None
                                
                rows.append((seq_idx, timepoint, marker, img_idx, shape, blur, lumin, contrast, edges, tenegrad))
                metric_idx += 1



        print("\n\n\n")
        print("Seq | Timepoint                | Marker | Img# | Shape       | Blur     | Lumin    | Contrast    | Edges_Density    | Sharpness")
        print("----+--------------------------+--------+------+-------------+----------+----------+-------------+------------------+----------")
        for r in rows:
            print(f"{r[0]:>3} | {r[1]:24} | {r[2]:6} | {r[3]:>3} | {r[4]} | {r[5]:8.3f} | {r[6]:8.3f} | {r[7]:8.3f} | {r[8]:8.3f} | {r[9]:8.3f}")

        
        


        min_bound_Blur = 1.54576
        max_bound_Blur = 279.84862
        
        min_bound_Lum = 4.05424
        max_bound_Lum = 200.0
        
        min_bound_Contr = 1.10639
        max_bound_Contr = 7.86189
        
        min_bound_Edg = 1.5419407894736843e-05
        max_bound_Edg = 0.13861932754516604
        
        
        min_bound_Sharp = 2.0002292961445325
        max_bound_Sharp = 50.69124668147343
        





        def check_bounds(List_Blur, List_Luminosity, List_Contrast, List_Edges, List_Tenegrad):
            
            flag_blur = []
            flag_contrast = []
            flag_luminosity = []
            flag_edges = []
            flag_sharp = []

            # Find indices where values are out of bounds
            for idx, value in enumerate(List_Contrast):
                if value < min_bound_Contr or value > max_bound_Contr:
                    flag_contrast.append(idx)

            for idx, value in enumerate(List_Luminosity):
                if value < min_bound_Lum or value > max_bound_Lum:
                    flag_luminosity.append(idx)

            for idx, value in enumerate(List_Blur):
                if value < min_bound_Blur or value > max_bound_Blur:
                    flag_blur.append(idx)

            for idx, value in enumerate(List_Edges):
                if value < min_bound_Edg or value > max_bound_Edg:
                    flag_edges.append(idx)

            for idx, value in enumerate(List_Tenegrad):
                if value < min_bound_Sharp or value > max_bound_Sharp:
                    flag_sharp.append(idx)

            # Check if any flags are not empty
            if flag_contrast or flag_luminosity or flag_blur or flag_edges or flag_sharp:
                print("\n \x1B[3mWarning\x1B[0m: Some values are out of the expected bounds, please check the images and their quality before proceeding with the analysis")
                print(f"\n \t \u2022 Contrast Values out of bounds: {flag_contrast}")
                print(f"\t \u2022 Luminosity Values out of bounds: {flag_luminosity}")
                print(f"\t \u2022 Blur Values out of bounds: {flag_blur}")
                print(f"\t \u2022 Edges Values out of bounds: {flag_edges}")
                print(f"\t \u2022 Sharpness Values out of bounds: {flag_sharp}\n")

            return flag_contrast, flag_luminosity, flag_blur, flag_edges, flag_sharp


        flag_contrast, flag_luminosity, flag_blur, flag_edges, flag_sharp = check_bounds(List_Blur, List_Luminosity, List_Contrast, List_Edges, List_Tenegrad)  # Check if the values are out of bounds
        
     











########################################################################################################################################################################
########################################################################################################################################################################




    ######################### Plotting of Pixels's Intensity Histograms (Whole Sequence / Single Element)
    
    print("\n\n\n\n\n###########################################################################")


    print(colored("\nWould you like to see the histograms related to the pixels intensity values of each sample or whole sequences of timepoints?", 'green', attrs = ['bold']))

    hist_input = inquirer.select(
        message = "\n\n\u2022N --> No, keep going with the normal analysis\n\u20221 --> Plot Histograms of complete sequences of samples for single timepoint\n\u20222--> PLot Histograms for each single sample\n",
        choices = ["N", "1", "2"],
        default = "N").execute()


    if hist_input in ["2"]:             # Check the input answer for 1 HIstograms for each single image marker

        for index, (single_sequence, related_metadata) in enumerate(zip(time_series, time_series_meta)):            # Loop through all the images and save single element and their indexes
            for image_index, (array, (time_point, marker)) in enumerate(zip(single_sequence, related_metadata)):
                fig, ax = plt.subplots(figsize = (10, 4))

                pixels = array.flatten()            # Flatten the array to ease the plotting
                
                    
                n, bins, patches = ax.hist(pixels, bins = range(256), edgecolor = 'none')                       # Plot the Histograms considering the pixels intensity range 
                ax.set_title(f"Histogram Pixel Intensities, Sequence --> {time_point}, Marker --> {marker}")       # Consider also the associated Time Point and Marker      (TOOK A REALLY LONG WHILE TO MAKE IT OUT)
                ax.set_xlabel("Pixel Intensities")
                ax.set_ylabel("Pixel Frequencies")
                #ax.set_xlim(0,255)
                plt.tight_layout()
                plt.show()




    elif hist_input in ["1"]:           # Check the input for 1 Histogram for a complete sequence per sample
        
            for index, single_sequence in enumerate(time_series):          # Loop trough the samples
                fig, ax = plt.subplots( figsize = (10, 4))      # Produce the sequence of markers plot
                pixels = []
                for image in single_sequence:                       # Loop through the single sequence of markers and append the flattened arrays to a single list
                    pixels.extend(image.flatten())
                
                
                n, bins, patches = ax.hist(pixels, bins = range(256), edgecolor = 'none') 
                ax.set_title(f"Histogram Pixel Intensities, Sequence --> {index}")
                ax.set_xlabel("Pixel Intensities")
                ax.set_ylabel("Pixel Frequencies")
                #ax.set_xlim(0,255)
                plt.show()
    
      

    
    
    
    
    
##########################################################################################################################################
###########################################################################################################################################



    ######################### Preprocessing of the images 

    ### Develop a funtion to apply the preprocessing to all sequences in the time-series list directly


    ## Applied some changes to resizing, beneath there is the original code line

    ##    def preprocess_single_seq(image_stack, target_size = (256,256), debug = False):

    def preprocess_single_seq(image_stack, metadata_stack, seq_idx, global_starting_index = 0, debug = False):

        extracted_features = []
        
        raw_rp = {}     # Dictionary to later store th features and use them later in the STATISTICS section
        
        for i, array in enumerate(image_stack):   #iterate

            global_index = global_starting_index + i  # To keep track of the global index across multiple sequences

            print(f"\tFiltering, Cleaning, Threhsolding and Features Extraction of array {i+1} of {len(image_stack)}:")
            
            ### Filter application to remove possible noise, while keeping intact, as much as reasonably possible, the biological relevan data
            transformed_layer = cv2.normalize(array, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)  #change format to allow filtering operation
            
            
            # Small Section for creating histograms related to the behaviour of pixel intensities in the original array and filtered cases (Bilateral, Median) 
            
            '''
            plt.figure(figsize=(15, 4))

            # Original behaviour
            
            plt.subplot(1,3,1)
            plt.hist(array.flatten(), bins = 256, color = "red")
            plt.title("Original Pixel Intensity Behaviour")
            plt.xlabel("Pizel Intensity")
            plt.ylabel("Frequency")

            # Bilateral filtered image histogram
            bilateral_img = cv2.bilateralFilter(array.astype(np.uint8), d=8, sigmaColor=75, sigmaSpace=75)
            plt.subplot(1, 3, 2)
            plt.hist(bilateral_img.flatten(), bins=256, color = "blue")
            plt.title("Bilateral Filtered Histogram")
            plt.xlabel("Pixel Intensity")
            plt.ylabel("Frequency")

            # Median filtered image histogram
            median_img = median_filter(array, size=3)
            plt.subplot(1, 3, 3)
            plt.hist(median_img.flatten(), bins=256, color = "green")
            plt.title("Median Filtered Histogram")
            plt.xlabel("Pixel Intensity")
            plt.ylabel("Frequency")

            plt.tight_layout()
            plt.show()
            '''



            
            # Create some parameters for adaptive changes according to the images
            layer_std = np.std(transformed_layer)
            
            base_sigma_color = 75
            
            sig_col = max(20, min(150, base_sigma_color * (layer_std/30)))
            
            
            # Apply the filtering, BILATERAL for edge preservation and noise reduction
            filtered_layer = cv2.bilateralFilter(transformed_layer, d = 8, sigmaColor = sig_col, sigmaSpace = 35)
            rescaled_layer = filtered_layer.astype(np.float32) / 255.0 # rechange format to the original so as to not lose biological relevant data
            
            
            
            
            ### Resizing of images
            #resizing_img = cv2.resize(rescaled_layer, target_size)
            
            
            
            normalized_img = rescaled_layer.astype(np.float32) / rescaled_layer.max()

            
            
            
            ## Apply local thresholding if image is over specific threshold values (considering Contrast, Luminosity and Edge Sharpness)
            ## This part is born from trying to "retrieve" biological information from images that can be considered as "low quality"
            ## or difficult to segment/deal with. These challenging instances, at least in the cardiac context, are due to asynchronous beating of the cardiomyocytes
            ## once they reach a certain differentiation phase (Early CMs phase, changes according to the culture's protocol but it should generally be around day 10-to-15).
            ## The control values are based on the "SMALL" datasets that the code was based upon, therefore they are highly subject to change/alteration for optimal control     
         
         
         
            if (List_Contrast[global_index] > 7.07) and (List_Luminosity[global_index] > 65.0) and (List_Edges[global_index] < 0.005):        #Derived "control" values from the IQA analysis
                
                ### Detect and Specify the image which is over the threshold values 
                print(f"\t\t  Low quality detected! (Image Global Index: {global_index + 1}) \n\t\tApplying adaptive thresholding...")
                
                
                
                
                ### Threshold image using Adaptive Mean thresholding
                thresholded_img = cv2.adaptiveThreshold((normalized_img * 255).astype(np.uint8), 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY, 31, 0)
                
                # Here multiple values for the blockSize and C were tested, the ones above seemed to be the most effective
                # If curious, these are the other values that were tested that were found to be "of interest"
                # for the local thresholding(Considering the positive results in both Adaptive_Mean and Adaptive_Gaussian):
                            # (15, 0), (21, -1), (21, 0), (31, -1), (31, 0),
                            # (35, -1), (41, -1), (45, -1), (51, -1), (51, -3)
                
                
                # Convert to binary mask, just to be consistent 
                binary = (thresholded_img > 0)
                
                
                
                
            # In case the value control are not exceeded, the code will try to apply a more "standard" segmentation approach    
            
            else: 
            ### Apply segmentation using Kittler-Illingworth's-Otsu fallback mechanism
                try:
                    img_8 = (normalized_img * 255).astype(np.uint8)
                    thresholded_img = mahotas.thresholding.kittler_illingworth(img_8)
                    binary = normalized_img > (thresholded_img / 255)
                    
                except:
                    otsu_threshold = skimage.filters.threshold_otsu(normalized_img)
                    binary = normalized_img > otsu_threshold
                
            

            
            ### Morphological cleaning of binary mask
            cleaned_img = morphology.remove_small_objects(binary, min_size = 128) 
            cleaned_img = morphology.remove_small_holes(cleaned_img, area_threshold = 96) 
            

    
            
            ### Distance Transform 
            dist_tran_img = ndimage.distance_transform_edt(cleaned_img)
            dist_smooth = gaussian_filter(dist_tran_img, sigma = 1)
            
            
            ### Local Maximas to use for the Watershed map
            local_peaks = local_maxima(dist_smooth)
            peaks = label(local_peaks)


            ### Watershed to split touching object
            limits = watershed(-dist_smooth, peaks, mask = cleaned_img)
        
            
            
            ### Feature extraction section 
            layer_props = regionprops(limits, intensity_image = cleaned_img)
            
            
            
            
            
            ### Section specific for Region of Interest (ROI) investigation, used to check segmentation and local bodies separation
            '''
            
            input_ROI = input("Specific ROI?")
            
             
            ### Info of specific ROI
            if input_ROI in ["y", "Y", "yes", "Yes", "YES"]:
                
                coordinates = input("Dimensions of ROI? --> (Format: x_start, x_end, y_start, y_end)")
                
                x1 = int(coordinates.split(",")[0])
                x2 = int(coordinates.split(",")[1])
                y1 = int(coordinates.split(",")[2])
                y2 = int(coordinates.split(",")[3])
                
                
                info_ROI = limits[y1:y2, x1:x2]
                
                num_bodies_ROI = len(np.unique(info_ROI)) - 1  # Subtract 1 to exclude the background label (0)
                
                print(f"\n\t\t  --> Number of bodies detected in the specified ROI: {num_bodies_ROI}\n")
                
            
                
                fig, axs = plt.subplots(1, 2, figsize = (12, 6))
                axs[0].imshow(normalized_img[y1:y2, x1:x2], cmap = "gray")
                axs[0].set_title("Normalized_Img - zoom")
                axs[1].imshow(info_ROI, cmap = "prism")
                axs[1].set_title("Labeled_Mask - zoom")
                axs[0].axis('off')
                axs[1].axis('off')
                plt.tight_layout()
                plt.show()
                
            '''

            
            
            
            
            if overlap_check in ["y", "Y", "yes", "Yes", "YES"]:


                fig, axes = plt.subplots(1, 3, figsize = (18, 6))


                axes[0].imshow(array, cmap = 'gray')
                axes[0].set_title("Original Fluorescence Image")
                axes[0].axis('off')


                # Mask only
                axes[1].imshow(limits, cmap = 'nipy_spectral')
                axes[1].set_title("Segmentation Mask")
                axes[1].axis('off')


                # Overlay: Original image with mask on top (semi-transparent)
                axes[2].imshow(array, cmap = 'gray')
                axes[2].imshow(limits, cmap = 'nipy_spectral', alpha = 0.3)  # alpha controls transparency
                axes[2].set_title("Overlay (Mask on Fluorescence)")
                axes[2].axis('off')

                plt.tight_layout()
                plt.show()

                
                height, width = array.shape
                dpi = 100  # Dots per inch
                
                fig_overlap = plt.figure(figsize = (width / dpi, height / dpi), dpi = dpi)
                
                ax_overlap = fig_overlap.add_axes([0, 0, 1, 1])
                ax_overlap.imshow(array, cmap = 'gray')
                ax_overlap.imshow(limits, cmap = 'nipy_spectral', alpha = 0.3)
                ax_overlap.axis('off')
                
                
                #plt.savefig(f"C:/Users/User/Desktop/Tesi_Magistrale/Code/Overlap/singles_base/Overlap_Img_{seq_idx + 1}_Array_{i + 1}.png", dpi = dpi, bbox_inches = 'tight', pad_inches = 0)
                

            
            ### Store the raw version of regionprops for later optional use in the Statistics section 
            if seq_idx not in raw_rp:
                raw_rp[seq_idx] = {}
            raw_rp[seq_idx][i] = layer_props
                        
            
            
            # Extracted features from regionprops function
            # Each feature is calculated as the mean value across all detected regions in the layer
            # If no regions are detected, default values are set to 0 or defaults
            
            # Each feature is explained in the README file
            Num_bodies = len(layer_props)
            
            print(f"\t{Num_bodies} bodies detected\n")
            
            
            # Area features
            if layer_props:  # Check if there are any detected regions
                areas = [p.area for p in layer_props]
                Area = np.mean(areas)
                Area_std = np.std(areas)
                Area_Min = np.min(areas)
                Area_Max = np.max(areas)
            else: 
                Area = Area_std = Area_Min = Area_Max = 0  
                
                
                
            # Intensity features
            if layer_props:
                Intensity = np.mean([p.mean_intensity for p in layer_props])
                Intensity_std = np.std([p.mean_intensity for p in layer_props])
                Intensity_Min = np.min([p.mean_intensity for p in layer_props])
                Intensity_Max = np.max([p.mean_intensity for p in layer_props])
                
            else:
                Intensity = 0
                Intensity_std = 0
                Intensity_Min = 0
                Intensity_Max = 0
                
       
            
            # Solidity features
            if layer_props:
                solidities = [p.solidity for p in layer_props]
                Solidity = np.mean(solidities)
                Solidity_std = np.std(solidities)
            else:
                Solidity = Solidity_std = 0
            
            
            
            # Major Axis Length to get the diameter of the cells
            if layer_props:
                Major_axes = [p.major_axis_length for p in layer_props]
                
                Average_axes = np.mean(Major_axes)
                Major_axes_min = np.min(Major_axes)
                Major_axes_max = np.max(Major_axes)
                Major_axes_std = np.std(Major_axes)
                
                quan_25 = np.percentile(Major_axes, 25)
                quan_75 = np.percentile(Major_axes, 75)
                
                IQR = quan_75 - quan_25
                
                lower, upper = quan_25 - (1.5*IQR), quan_75 + (1.5*IQR)
            
            else:
                Average_axes = Major_axes_min = Major_axes_max = Major_axes_std = 0
                        
            
            
            # Perimeter features
            if layer_props:    
                perimeters = [p.perimeter for p in layer_props] 
                Perimeter = np.mean(perimeters)
                Perimeter_std = np.std(perimeters)
            else:
                Perimeter, Perimeter_std = 0
                
                
                
                
            # Diameter features
            if layer_props:
                diameters = [p.equivalent_diameter_area for p in layer_props]
                Diameter = np.mean(diameters)
                Diameter_std = np.std(diameters)
            else:
                Diameter, Diameter_std = 0 
            
            
            
            # Orientation features
            if layer_props:
                orientations = [p.orientation for p in layer_props]
                Orientation = np.mean(orientations)
                Orientation_std = np.std(orientations)
            else:
                Orientation, Orientation_std = 0    
             
            
            
            # Eccentricity features
            if layer_props:
                eccentricities = [p.eccentricity for p in layer_props]
                Eccentricity = np.mean(eccentricities)
                Eccentricity_std = np.std(eccentricities)
            else:
                Eccentricity, Eccentricity_std = 0
            
            
            # Centroids = np.mean([p.centroid for p in layer_props], axis=0) if layer_props else np.zeros(2)  # Centroid coordinates
            # Moments = np.mean([p.moments for p in layer_props], axis=0) if layer_props else np.zeros((3, 3))  # Raw moments
            # Moments_Central = np.mean([p.moments_central for p in layer_props], axis=0) if layer_props else np.zeros((3, 3))  # Central moments
            
            # Moments features 
            if layer_props:
                HU_list = []
                for giggino in layer_props:
                    hu_moments = giggino.moments_hu
                    log_hu_moments = -np.sign(hu_moments) * np.log(np.abs(hu_moments) + 1e-10)  # Avoid log(0) by adding a small constant 
                    HU_list.append(log_hu_moments)
                    Moments_Hu = np.mean(HU_list, axis=0)
            else:
                Moments_Hu = np.zeros(7)
        

            extracted_features.append([Num_bodies, Area, Area_Min, Area_Max, Area_std, Solidity, Solidity_std, Diameter, Diameter_std, Perimeter, Perimeter_std,
                                       Orientation, Orientation_std, Eccentricity, Eccentricity_std, Average_axes, Major_axes_max,
                                       Major_axes_min, Major_axes_std,*Moments_Hu])        #list of feature later transformed into an array
            
            
            if debug:
                fig, axs = plt.subplots(1, 4, figsize = (20, 6))
                axs[0].imshow(normalized_img, cmap = "gray")
                axs[0].set_title("Normalized_Img")
                
                axs[1].imshow(binary, cmap = "gray")
                axs[1].set_title("Img_Mask")
                
                axs[2].imshow(dist_tran_img, cmap = "magma")
                axs[2].set_title("Distance_Transform")
                
                axs[3].imshow(limits, cmap = "nipy_spectral")
                axs[3].set_title("Watershed_Img")
                
                for g in axs:
                    g.axis("off")
                    
                plt.show()
             
                
                
            
                plt.imshow(limits, cmap = "nipy_spectral")
                plt.axis("off")
                
                plt.tight_layout()
                plt.show()





        return np.array(extracted_features), raw_rp  # Return the extracted features as a numpy array



    ### Function to apply preprocessing to all images in the folder 
    # remember to put "target_size = (256,256)" if you want to resize the images to a specific size, after "time series"
    
    
    def all_seq_preprocessing(time_series):
        
        preprocessed_elements = []  # Final list of vectors extracted from all sequences
        global_idx = 0  # Initialize a global index to keep track of the image index across sequences
        global_rp = {}  # Global dict for the regionprops of all sequences, needed to conduct the Statistics section 
        
        
        
        
        print("\n\n\n###########################################################################")    
        
        print(colored("\nWould you like to visualize the Plots related to the different Samples throughout the Preprocessing section?", 'green', attrs = ['bold']))
            
        debug_inp = inquirer.select(
            message = "\n",
            choices = ["Yes", "No"],
            default = "Yes").execute()  # Ask the user if he wants to see the plots of the image processing steps
        
        
        debug = debug_inp in ["Y", "y", "YES", "yes", "Yes"]
        
        
        
        
        SPINNERS["custom_spinner"] = {
            
            "interval" : 100,  
            "frames" : ["▱▱▱▱▱▱▱","▰▱▱▱▱▱▱","▰▰▱▱▱▱▱","▰▰▰▱▱▱▱",
                        "▰▰▰▰▱▱▱","▰▰▰▰▰▱▱","▰▰▰▰▰▰▱","▰▰▰▰▰▰▰"]
            
        } 
        
        
        
        
        with Progress(
            
            SpinnerColumn(spinner_name = "custom_spinner", style = 'dark_red'),
            TextColumn("[progress.description]{task.description}"),
        ) as progress:  # Use the Progress context manager to handle the progress bar
            
            
            
            task = progress.add_task("[dark_red]Preprocessing Images...", total = len(time_series))  # Add a task to the progress bar
        
        
            
            # Functional loop
            for i, (seq, seq_meta) in enumerate(zip(time_series, time_series_meta)):  # Loop through each sequence in the time series
                print(colored("\n\n\nWorking with Sequence: ", 'red') + f"{i+1}" + colored(" of ", 'red') + f"{len(time_series)}")
                
                
                
                #remember to put "target_size = (256,256)" if you want to resize the images to a specific size, after "seq"
                seq_features, raw_rp_seq = preprocess_single_seq(seq, seq_meta, seq_idx = i, global_starting_index = global_idx,debug = debug)
                preprocessed_elements.append(seq_features)          # Append the resulting features in an list, while being associated to the respective array
                global_rp.update(raw_rp_seq)                        # Update the global regionprops dictionary with the current sequence's data
                



                global_idx += len(seq)  # Update the global index by adding the number of images in the current sequence



            progress.update(task, description = "[dark_red italic]Completed")  # Update the progress bar after processin each sequence
        time.sleep(1)  # Add a small delay to allow the progress bar to update smoothly


        return np.array(preprocessed_elements), global_rp     # Final resutl should be a 3-dimensional array  (num_samples, time_step, feature_vector)  + regionprops dictionary









    ready_to_be_fed, raw_rp = all_seq_preprocessing(time_series)

    
    
    if ready_to_be_fed.ndim == 3:
        print(f"\n\n\nThe shape of the final array coming from the Preprocessing phase is:\t{ready_to_be_fed.shape}")
        print(f"\n\n\nVariable Type:\t{type(ready_to_be_fed)}")
        print("\nThe data seems \x1B[3mready to be fed\x1B[0m to the LSTM neural network, analysis will continue. \n")
        print(ready_to_be_fed)
    
    else:
        print(f"\n\n\nIt seems there is a problem with the shape of the data obtained from the Preprocessing phase\n >{ready_to_be_fed.shape}")
        print("\nPlease check again the data or the processing of the images")



  
else:
    print("\n\n\nIt seems you are not ready to start the analysis, please check the requirements and try again later")









############################################################################################################################################################
############################################################################################################################################################

    ######################### Statistics of the Extracted features from input data 




print("\n\n\n" + "#" * 75 + "\n\n\n")



## Dictionary preparation and fixing/formatting

def DF_prep(raw_rp, time_series_meta):

    raw_props = {}

    for seq_index, marker_dict in raw_rp.items():

        raw_props[seq_index] = {}
        timepoint = time_series_meta[seq_index][0][0]  # Get the timepoint from the first entry

        for marker, props in marker_dict.items():

            marker_name = time_series_meta[seq_index][marker][1]  # Get the marker name from metadata

            raw_props[seq_index][marker_name] = {
                
                "Timepoint": timepoint,
                "Area": [p.area for p in props],
                "Solidity": [p.solidity for p in props],
                "Major_Axes_Lengths": [p.major_axis_length for p in props],
                "Perimeter": [p.perimeter for p in props],
                "Diameter": [p.equivalent_diameter_area for p in props],
                "Orientation": [p.orientation for p in props],
                "Eccentricity": [p.eccentricity for p in props],
                "Num_bodies": [len(props)]
            }
    
    
    
    return raw_props




## Statistics section designed to produce boxplots over temporal courses for each feature and marker
## The boxplots are designed to show the distribution of the features across different timepoints, allowing to visualize the temporal evolution of the features for each marker.
## Temporal evolution plots are also produced, focusing on the mean value of one specific feature and its evolution over time
## The Kruskal-Wallis and Wilcoxon signed-rank tests were implemented for features intra- and inter-model comparisons.

def Statistics(raw_props, time_series_meta, feature_names):
        

    def Boxplot_feat(raw_props, time_series_meta, feature_names):

        # Collect all markers from the dictionary
        markers = sorted({marks for seq in raw_props.values() for marks in seq.keys()})

        for marker in markers:
            for feature in feature_names:

                data = []
                labels = []

                for seq_idx in sorted(raw_props.keys()):

                    # Skip if marker missing, came up with this due to presence of other datasets lacking certain timepoints of interest
                    if marker not in raw_props[seq_idx]:
                        continue


                    # Extract day from dict
                    timepoint = time_series_meta[seq_idx][0][0]
                    match = re.search(r'[dD](\d+)', timepoint)
                    tp_day = int(match.group(1)) if match else None


                    # Get asociated values from the extracted day
                    vals = raw_props[seq_idx][marker][feature]
                    if len(vals) == 0:
                        continue


                    data.append(vals)
                    labels.append(f"D{tp_day}")


                if len(data) == 0:
                    continue



                # Plotting the boxplots for eahc feature and marker across timepoints 
                plt.figure(figsize=(10, 6))
                bp = plt.boxplot(data, labels=labels,
                                patch_artist=True, showmeans=True, meanline=True)




                colors = sns.color_palette("Set2", len(data))
                
                for patch, color in zip(bp['boxes'], colors):
                    patch.set_facecolor(color)
                    patch.set_alpha(0.65)


                plt.title(f"{marker.capitalize()} — {feature} Across Timepoints",
                        fontsize = 15, weight = 'bold')
                
                plt.xlabel("Timepoint (Days)", fontsize = 12)
                plt.ylabel(feature, fontsize = 12)
                plt.grid(True, axis = 'y', alpha = 0.3)


                plt.tight_layout()
                plt.show(block = False)
                plt.pause(0.1)
                plt.close()








            

    def Temporal_Evolution(raw_props, time_series_meta, feature):
        
 
        
        
        
        ############################################### FAVOURITE ONE
        
        
        # Collect unique markers
        all_markers = sorted(set(marker for seq in raw_props.values() for marker in seq.keys()))
        
        fig, ax = plt.subplots(figsize = (14, 6))
        colors = sns.color_palette("husl", len(all_markers))
        
        # Build list of all sequence labels in order
        sequence_labels = []
        x_positions = []
        
        for seq_idx in sorted(raw_props.keys()):
            timepoint = time_series_meta[seq_idx][0][0]
            
            # Extract day number
            match = re.search(r'[dD](\d+)', timepoint)
            if not match:
                continue
            
            tp_day = int(match.group(1))
            
            # Count previous occurrences of this day to detect replicates
            prev_count = 0
            for prev_idx in range(seq_idx):
                prev_tp = time_series_meta[prev_idx][0][0]
                prev_match = re.search(r'[dD](\d+)', prev_tp)
                if prev_match and int(prev_match.group(1)) == tp_day:
                    prev_count += 1
            
            # Create label
            if prev_count > 0:
                label = f"D{tp_day}_S{prev_count + 1}"
            else:
                label = f"D{tp_day}"
            
            sequence_labels.append(label)
            x_positions.append(seq_idx)
            
            
        
        # Plot each marker
        for m_idx, marker in enumerate(all_markers):
            
            means = []
            mins = []
            maxs = []
            x_vals = []
            labels_for_marker = []
            
            for seq_idx in sorted(raw_props.keys()):
                
                if marker not in raw_props[seq_idx]:
                    continue
                
                # Get raw values (individual cell measurements)
                raw_vals = raw_props[seq_idx][marker][feature]
                
                if len(raw_vals) == 0:
                    continue
                
                means.append(np.mean(raw_vals))
                mins.append(np.min(raw_vals))
                maxs.append(np.max(raw_vals))
                x_vals.append(seq_idx)
                labels_for_marker.append(sequence_labels[seq_idx])
            
            if len(means) == 0:
                continue
            
            x_vals = np.array(x_vals)
            means = np.array(means)
            mins = np.array(mins)
            maxs = np.array(maxs)
            
            # Plot mean line
            ax.plot(x_vals, means, marker='o', markersize = 6, linewidth = 2,
                label = marker.capitalize(), color = colors[m_idx], zorder = 5)
            
            # Plot min/max range as subtle shaded area
            ax.fill_between(x_vals, mins, maxs, color = colors[m_idx],
                        alpha = 0.25, zorder = 1)
            
            # Add small vertical error bars for better visibility, comment the line to remove it,  I personally found it slightly useful
            for x, mean, mn, mx in zip(x_vals, means, mins, maxs):
                ax.plot([x, x], [mn, mx], color = colors[m_idx],
                    linewidth = 1.5, alpha = 0.4, zorder = 2)
        
        
        # Set x-axis days, adn slightly modify formatting for better visibility
        ax.set_xticks(x_positions)
        ax.set_xticklabels(sequence_labels, rotation = 20, ha = 'right')
        
        
        
        # Add vertical lines to separate different days, took an hell of a lot of time to do
        current_day = None
        for seq_idx in sorted(raw_props.keys()):
            timepoint = time_series_meta[seq_idx][0][0]
            match = re.search(r'[dD](\d+)', timepoint)
            if match:
                day = int(match.group(1))
                if current_day is not None and day != current_day:
                    ax.axvline(x = seq_idx - 0.5, color = 'gray', linestyle='--',
                            alpha = 0.25, linewidth = 1, zorder = 0)
                current_day = day
        
        
        # Further styling and formatting of the plot
        ax.set_title(f"Temporal Evolution: {feature}\n(Shaded area = min/max range from individual cells)",
                    fontsize = 14, weight = 'bold', pad = 15)
        
        ax.set_xlabel("Sequence (Timepoint)", fontsize = 10, weight = 'bold')
        
        ax.set_ylabel(f"{feature} Value", fontsize = 10, weight = 'bold')
        
        ax.legend(title = "Marker", fontsize = 11, title_fontsize = 12,
                loc = 'best', frameon = True, shadow = True)
        
        ax.grid(True, alpha = 0.25, axis = 'y')
        
        plt.tight_layout()
        plt.show()

        









    Boxplot_feat(raw_props, time_series_meta, feature_names)


    for feats in feature_names:
        Temporal_Evolution(raw_props, time_series_meta, feats)










print("\n\n\n" + "#" * 75 + "\n\n\n")

print("Execute Statistics?\n")

stat_inp = inquirer.select(
    message = "\n",
    choices = ["Yes" , "No"],
    default = "No"
).execute()


if stat_inp in ["Yes"]:
    
    
    feature_names = ["Num_bodies", "Area", "Solidity", "Diameter", 
                     "Perimeter", "Orientation", "Eccentricity", "Major_Axes_Lengths"]
    
    print("\n\n\nExecuting Statistics Section...\n\n")
    
    with Progress(
        SpinnerColumn(spinner_name = "custom_spinner", style = 'dark_red'),
        TextColumn("[progress.description]{task.description}"),
    ) as progress:
        
        task = progress.add_task("[dark_red]Converting to DataFrame...", total = 2)
        
        raw_props = DF_prep(raw_rp, time_series_meta)
        
        progress.update(task, advance = 1, description = "[dark_red italic]Features Plotting...")
        
        Statistics(raw_props = raw_props, time_series_meta = time_series_meta, feature_names = feature_names)
        
        progress.update(task, description = "[dark_red italic]Completed")
        
    time.sleep(1)

else:
    print("\n\n\nStatistics section skipped as per user request.\n\n")
    
    
    

############################################################################################################################################################



### Part related to the saving of data for later comparison with the CELLPOSE VERSION of the tool on the google colab environment
'''

import pickle

output_dir = "data_to_plot"
os.makedirs(output_dir, exist_ok = True)

method_oi = "Model_1"

data_to_save = {
    
    "raw_props" : raw_props,
    "feature_names" : feature_names,
    "time_series_meta" : time_series_meta,
    "ready_to_be_fed" : ready_to_be_fed,
    "method" : method_oi
}



all_markers = []

for seq_idx in data_to_save["raw_props"]:
  all_markers.extend(data_to_save["raw_props"][seq_idx].keys())


uniques_marks = sorted(set(all_markers))


data_to_save["markers"] = uniques_marks




filename = f"{output_dir}/{method_oi}_data.pkl"

with open(filename, 'wb') as file:
    pickle.dump(data_to_save, file)





    
end_time = time.time() - start_time

print(colored("\n\nTime Spent for Analysis", 'green', attrs = ['bold']) + f"\n\t> {end_time/60:.2f} minutes\n\n")
    
'''



############################################################################################################################################################
print("\n\n\n" + "#" * 75 + "\n\n\n")






######################### GRU model development, or at least the first try, due to limited dataset the development halted

########### First create the datasets related to training and testing 
#just a try, will have to modify later with the real inputs
num_sequences = ready_to_be_fed.shape[0]

y = np.random.randint(0, 2, size = (num_sequences,))

x_train, x_test, y_train, y_test = train_test_split(                # Training and testing groups creation
    ready_to_be_fed, y, test_size = 0.2, random_state = 42, stratify = y  
)


########### GRU Architecture

#identify the interesting parameters of the array
time_steps = ready_to_be_fed.shape[1]
num_features = ready_to_be_fed.shape[2]


### Model architecture creation
GRU_model = Sequential()    #Sequential so that it can interact with additional layers

#First two layers are GRU, both with 64 HUs
GRU_model.add(GRU(units = 64, return_sequences = True, input_shape = (time_steps, num_features), dropout = 0.2))
    
GRU_model.add(GRU(units = 64, return_sequences = False, dropout = 0.2))

#Third layer is fully connected and has 64 HUs, plus it works with the ReLUs
GRU_model.add(Dense(32, activation = 'relu'))

#Addition of "Dropout" in order to reduce the overfitting
GRU_model.add(Dropout(0.2))

#Final layer for binary classification
GRU_model.add(Dense(1, activation = 'sigmoid'))



### Model Compiling 
GRU_optimizer = tf.keras.optimizers.Adam(
    learning_rate = 0.001,
    clipnorm = 1.0)

GRU_model.compile(optimizer = GRU_optimizer, loss = 'binary_crossentropy',        # Model Compiling
                metrics = ['accuracy', AUC(name = 'auc')])



######################################################################################################################################


######################### Training phase of the model


#Create an early stop checkpoint in case something goes south
stopping_point = EarlyStopping(monitor = 'val_loss', patience = 10, reset_best_weights = True)

#fit the model
model_fitting = GRU_model.fit(x_train, y_train, validation_data = (x_test, y_test),
                            epochs = 100, batch = 8, verbose = 1, callbacks = [stopping_point])



############## Evaluation of training schedule of the model

#create a function to evaluate the reults from the fitting of the input data
def train_plotting_graphs(model_fitting):
    
    #define range
    epochs_range = range(1, len(model_fitting.history['loss']) + 1)
    
    fig, axis = plt.subplots(1,3, figsize = (20, 6)) 
    
    ## Loss plotting
    axis[0].plot(epochs_range, model_fitting.history['loss'], label = 'Training_Loss')
    axis[0].plot(epochs_range, model_fitting.history['val_loss'], label = 'Validation_Loss')
    axis[0].set_title = 'Learning Curvature'
    axis[0].set_xlabel = 'Epochs'
    axis[0].set_ylabel = 'Loss'
    axis[0].legend()
    axis[0].grid(True)
    
    
    ## Accuracy plotting
    axis[1].plot(epochs_range, model_fitting.history['accuracy'], label = 'Training_Accuracy')
    axis[1].plot(epochs_range, model_fitting.history['val_accuracy'], label = 'Validation_Accuracy')
    axis[1].set_title = 'Accuracy throughout the Epoches'
    axis[1].set_xlabel = 'Epochs'
    axis[1].set_ylabel = 'Accuracy'
    axis[1].legend()
    axis[1].grid(True)
    
    
    ## AUC plotting
    axis[2].plot(epochs_range, model_fitting.history['auc'], label = 'AUC')
    axis[2].plot(epochs_range, model_fitting.history['val_auc'], label = 'Validation AUC')
    axis[2].set_title = 'Area Under the Curve'
    axis[2].set_xlabel = 'Epochs'
    axis[2].set_ylabel = 'AUC Score'
    axis[2].legend()
    axis[2].grid(True)
    
    
    
    #show the otained plots
    plt.tight_layout()      #useful to give us a better visualization 
    plt.show()


#apply function to plot
train_plotting_graphs(model_fitting)

######################################################################################################################################

######################### Model Evaluation

model_results = GRU_model.evaluate(x_test, y_test, verbose = 1)

print('Test Loss results:', model_results[0])
print('Test Accuracy results:', model_results[1]) 

y_pred_prob = GRU_model.predict(x_test)
y_pred_class = (y_pred_prob > 0.5).astype(int)

f1_model_score = f1_score(y_test, y_pred_class)
print('F1 Score results:', f1_model_score)


## Confusion Matrix
cm = confusion_matrix(y_test, y_pred_class)
figure, ax = plt.subplot()
figure.set_size_inches(12, 8)
sns.heatmap(cm, annot = True, fmt = 'd', cmap = plt.cm.Blues, cbar = False,
            xticklabels = ['Not Differentiated', 'Differentiated'], yticklabels = ['Not Differentiated', 'Differentiated'])

plt.xlabel('Predicted')
plt.ylabel('Effective')
plt.title('Differentiation Confusion Matrix')
plt.tight_layout()
plt.show()


######################################################################################################################################

######################### Predict New Inputs

#define a function to predict new unseen data
def predict_new_inputs(sequence_batch):
    new_seq = scaler.transform(sequence_batch.reshape(-1, sequence_batch.shape[-1])).reshape(sequence_batch.shape)
    return GRU_model.predict(new_seq)


######################################################################################################################################


