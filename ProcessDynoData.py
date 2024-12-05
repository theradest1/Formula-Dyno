# +
"""This program takes in a folder of log files of formats .mdf, .csv, and .xls, 
and outputs a heatmap of all the cells in the VE map visited in those logs. Log files can be added to the folder and the
program can be rerun to update the heatmap. The heatmap is cumulative.

LIMITATIONS:
1. This program can only handle these log file formats: .mdf, .csv, and .xls
2. At least the following channels MUST be selected when the log is exported: 't[s]', 'Plenum_MAP[kPa]', 'RPM[RPM]', 
"CLVEAdapting[bool]", "Eng_SteadyState[bool]
3. If the names of those channels change, the code may need to be modified.
4. The parent folder of the logs must be specified if it is different form the value in "logs_folder"

NOTES:
1. Progress over time is saved in the folder "Processed_Data". If you want to recalculate everything, delete this folder.
2. Other file file types may be saved in the logs folder. This program will ignore them. If they are of a readable, type, 
however, the program will attempt to process them. 
3. The final output of the program is the .png file in Processed_Data.
4. If only error messages are displayed before displaying output, you should still trust the output.
5. If multiple files of the same name with different extensions are available, the mdf, then xls, then csv will be tried."""

#import required modules. if missing any, do "pip install <module name>" in your terminal.
import mdfreader
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

"""Recursively searches through the parent of all log folders to be looked through. 
The log folder must be located in the same directory as this script. The log folder can be infinitely deep."""
def get_relative_paths(folder_path):
    relative_paths = []

    for root, dirs, files in os.walk(folder_path):
        for file in files:
            file_path = os.path.join(root, file)
            relative_path = os.path.relpath(file_path, folder_path)
            
            # only selects files of supported types (.xls, .csv, and .mdf)
            if relative_path.upper().endswith(".MDF"):
                relative_paths.append(os.path.join(folder_path, relative_path))

    return relative_paths

#paths
logs_folder = "Dyno"  # NOTE: may need to change depending on location of log files.
completed_stuff_folder = 'Processed_Data' # folder that files from this program will be stored to
completed_logs_path = 'completed_logs.txt' # list of files that have been parsed successfully and are represented in the heatmap

counts_path = "counts.csv" # current heatmap data as a dataframe
counts_heatmap_path = 'counts_heatmap_output.png' # current heatmap as a png

cylPhiDiff_path = "cylPhiDiff.csv" # current heatmap data as a dataframe
cylPhiDiff_heatmap_path = 'cylPhiDiff_heatmap.png' # current heatmap as a png

phiOffset_path = "phiOffset.csv" # current heatmap data as a dataframe
phiOffset_heatmap_path = 'phiOffset_heatmap.png' # current heatmap as a png

count = 0 # counter for number of cells visited in a file

filepaths = get_relative_paths(logs_folder)

# bins from the VE map. add a dummy first bin to make pandas happy.
loadBreaks = [0.00, 7.69, 15.38, 23.08, 30.77, 38.46, 46.15, 53.85, 61.54, 69.23, 76.92, 84.62, 92.31, 100.00, 105.13, 110.27, 115.40]
rpmBreaks = [1750.00, 2359.38, 2968.75, 3578.13, 4187.50, 4796.88, 5406.25, 6015.63, 6625.00, 7234.38, 7843.75, 8453.13, 9062.50, 9671.88, 10281.25, 10890.63, 11500.00]
mapBreaks = [18.00, 23.06, 28.13, 33.19, 38.25, 43.31, 48.38, 53.44, 58.50, 63.56, 68.63, 73.69, 78.75, 83.81, 88.88, 93.94, 99.00]

if not os.path.exists(completed_stuff_folder):  # first execution of program
    os.mkdir(completed_stuff_folder)
    completed_logs = []
    with open(os.path.join(completed_stuff_folder, completed_logs_path), 'w'):  # create a blank file
        pass
else:
    try:
        with open(os.path.join(completed_stuff_folder,completed_logs_path), 'r') as file:
            completed_logs = [line.strip() for line in file]
    except Exception as e:
        print("""Check that the VEVisits_data folder exists in this directory and that the contents include 
        completed_logs.txt, counts.csv, and heatmap_output.png. If not, delete VEVisits_data and rerun.""")

if not os.path.exists(os.path.join(completed_stuff_folder,counts_path)):  # first execution of program
    # Create an empty dataframe with RPM and MAP buckets as columns
    df_counts = pd.DataFrame(index=rpmBreaks, columns=mapBreaks, dtype=float)
    df_counts = df_counts.fillna(0)  # Initialize all counts to 0
    df_counts.to_csv(os.path.join(completed_stuff_folder,counts_path), index=True, header=True)
    
    df_cylPhiDiff = pd.DataFrame(index=rpmBreaks, columns=loadBreaks, dtype=float)
    df_cylPhiDiff = df_cylPhiDiff.fillna(0)  # Initialize all counts to 0
    df_cylPhiDiff.to_csv(os.path.join(completed_stuff_folder,cylPhiDiff_path), index=True, header=True)
    
    df_phiOffset = pd.DataFrame(index=rpmBreaks, columns=mapBreaks, dtype=float)
    df_phiOffset = df_phiOffset.fillna(0)  # Initialize all counts to 0
    df_phiOffset.to_csv(os.path.join(completed_stuff_folder,phiOffset_path), index=True, header=True)
else:
    try:
        df_counts = pd.read_csv(os.path.join(completed_stuff_folder, counts_path), index_col=0)
        df_counts.index = df_counts.index.astype(float)
        df_counts.columns = df_counts.columns.astype(float)
        
        df_cylPhiDiff = pd.read_csv(os.path.join(completed_stuff_folder, cylPhiDiff_path), index_col=0)
        df_cylPhiDiff.index = df_cylPhiDiff.index.astype(float)
        df_cylPhiDiff.columns = df_cylPhiDiff.columns.astype(float)
        
        df_phiOffset = pd.read_csv(os.path.join(completed_stuff_folder, phiOffset_path), index_col=0)
        df_phiOffset.index = df_phiOffset.index.astype(float)
        df_phiOffset.columns = df_phiOffset.columns.astype(float)
    except Exception as e:
        print("""Check that the VEVisits_data folder exists in this directory and that the contents include 
        completed_logs.txt, counts.csv, and heatmap_output.png. If not, delete VEVisits_data and rerun.""")

"""For a given row of a log, find out which cell in the VE Map it corresponds to."""
def get_bucket(row, buckets, channel):
    val = row[channel] # will be either RPM or MAP
    my_buckets = buckets.copy() #remove dummy bin 
    # find bin above and below val. return the bin that val is closer to.
    for pot_val_idx in range(len(buckets)):
        try:
            lower_bound = buckets[pot_val_idx]
            upper_bound = buckets[pot_val_idx+1]
            if val >= lower_bound and val <= upper_bound:
                rng = upper_bound - lower_bound
                crit = rng/2 + lower_bound
                if val > crit:
                    return upper_bound
                else:
                    return lower_bound
        except IndexError as e: # val was either bigger or smaller than the range of buckets.
            if val >= buckets[len(buckets)-1]:
                return buckets[len(buckets)-1]
            else:
                return buckets[0]
            
"""Go through the rows in a file and figure out what cell to modify and if adaptation was possible"""
def update_dataframes(row):
    global count
    
    #get bucket
    rpm_bucket = get_bucket(row, rpmBreaks, 'RPM[RPM]')
    map_bucket = get_bucket(row, mapBreaks, 'Plenum_MAP[kPa]')
    load_bucket = get_bucket(row, loadBreaks, 'Load[%]')

    #if in steady state and engine is running
    if row["Steady_State_Op[bool]"] == 1.0 and row["RPM[RPM]"] >= rpmBreaks[0]:
        count += 1
        
        pastBucketCount = df_counts.loc[rpm_bucket, map_bucket]
        df_counts.loc[rpm_bucket, map_bucket] = pastBucketCount + 1
        
        cylPhiDiff = row["Cyl1_Phi[EqRatio]"] / row["Cyl2_Phi[EqRatio]"]
        newCylPhiDiffAvg = (df_cylPhiDiff.loc[rpm_bucket, load_bucket]*pastBucketCount + cylPhiDiff)/(pastBucketCount + 1)
        df_cylPhiDiff.loc[rpm_bucket, load_bucket] = newCylPhiDiffAvg
        
        phiOffset = (row["Cyl1_Phi[EqRatio]"] + row["Cyl2_Phi[EqRatio]"])/2 - row["Target_Phi[EqRatio]"]
        newPhiOffsetAvg = (df_phiOffset.loc[rpm_bucket, map_bucket]*pastBucketCount + phiOffset)/(pastBucketCount + 1)
        df_phiOffset.loc[rpm_bucket, map_bucket] = newPhiOffsetAvg

"""Takes in a path to an mdf file and returns a formatted dataframe to be processed in parse_data()"""
def deal_with_mdf(input_path):
    # Read mdf file
    mdf_file = mdfreader.Mdf(input_path)

    # convert to csv temporarily
    mdf_file.export_to_csv('temp.csv')
    
    # try to read the csv with one set of channels. use xls nomenclature
    df = pd.read_csv('temp.csv', low_memory=False, usecols=["t", "RPM", "Plenum_MAP", "Steady_State_Op", "Cyl1_Phi", "Cyl2_Phi", "Target_Phi", "Load"])[1:]
    
    df.rename(columns={'t': 't[s]',
                        'RPM': 'RPM[RPM]',
                        'Plenum_MAP' : 'Plenum_MAP[kPa]',
                        'Steady_State_Op': 'Steady_State_Op[bool]',
                        'Cyl1_Phi': 'Cyl1_Phi[EqRatio]', 
                        'Cyl2_Phi': 'Cyl2_Phi[EqRatio]', 
                        'Target_Phi': 'Target_Phi[EqRatio]',
                        'Load': 'Load[%]'
                        }, inplace=True)
    
    # get rid of the temporary file
    os.remove('temp.csv')
    
    # make all values numbers
    df = df.apply(pd.to_numeric)
    
    return df

"""Checks if a filepath has already been successfully read"""
def already_done(filepath):
    for pot_filepath in completed_logs:
        if filepath == pot_filepath:
            return True
    return False

"""Parses a given file path by filetype. Updates the heatmap dataframe accordingly."""
def parse_data(input_path):
    sensor_data = deal_with_mdf(input_path)
    print("Done parsing.")
    
    # start counting how many cells in this file
    count_before = count
    
    # calculate the visits and add to heatmap dataframe
    sensor_data.apply(update_dataframes, axis=1)
    count_after = count
    
    print(f"{count_after-count_before} datapoints were added for a total of {count_after}\n")

"""Makes and saves the heatmap."""
def make_heatmap(df, title, savedImg_path):
    fig, ax = plt.subplots(figsize=(25, 15))  # set the figure size (adjust as needed)
    ax.xaxis.tick_top()  # move x-axis ticks to the top
    heatmap = sns.heatmap(df, annot=True, cmap='coolwarm', fmt='g', cbar = False, ax=ax)
    
    # adjust axis tick label orientation and size
    heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=0, fontsize=13)
    heatmap.set_yticklabels(heatmap.get_yticklabels(), rotation=0, fontsize=13)
    
    plt.suptitle(title, fontsize=24, fontweight='bold', y=0.95)
    
    # Save the heatmap to a file
    plt.savefig(os.path.join(completed_stuff_folder, savedImg_path), bbox_inches='tight')
    plt.show()


"""Runs all code in this file."""
def main():
    # parse all filepaths in the logs folder
    
    addedFiles = 0
    for filepath in filepaths:
        if not already_done(filepath):
            try:
                print("Now parsing: " + filepath)
                parse_data(filepath)
                completed_logs.append(filepath)
                addedFiles += 1
            except ValueError as e:
                print(e)
                print(filepath + " could not be parsed. The required channels may not be available in this log.\n")
                continue
    
    if addedFiles == 0:
        print("No unparsed mdf files found.\nDelete the processed data folder if you updated a file.")

    # store the seen logs
    with open(os.path.join(completed_stuff_folder,completed_logs_path), 'w') as file:
        for filepath in completed_logs:
            file.write(f"{filepath}\n")

    # store the current heatmap df
    df_counts.to_csv(os.path.join(completed_stuff_folder,counts_path), index=True, header=True)
    df_cylPhiDiff.to_csv(os.path.join(completed_stuff_folder, cylPhiDiff_path), index=True, header=True)
    df_phiOffset.to_csv(os.path.join(completed_stuff_folder, phiOffset_path), index=True, header=True)

    # display the heatmap
    make_heatmap(df_counts, 'Visits Per Cell (RPM vs MAP)', counts_heatmap_path)
    make_heatmap(df_cylPhiDiff, 'Cylinder Phi difference (RPM vs Load)', cylPhiDiff_heatmap_path)
    make_heatmap(df_phiOffset, 'Phi Offset (RPM vs MAP)', phiOffset_heatmap_path)

main()


