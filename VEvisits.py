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
1. Progress over time is saved in the folder "VEVisits_data". If you want to recalculate everything, delete this folder.
2. Other file file types may be saved in the logs folder. This program will ignore them. If they are of a readable, type, 
however, the program will attempt to process them. 
3. The final output of the program is the .png file in VEVisits_data.
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

logs_folder = "Dyno"  # NOTE: may need to change depending on location of log files.
completed_stuff_folder = 'VEVisits_data' # folder that files from this program will be stored to
completed_logs_path = 'completed_logs.txt' # list of files that have been parsed successfully and are represented in the heatmap
current_counts_path = "current_counts.csv" # current heatmap data as a dataframe
current_heatmap_path = 'heatmap_output.png' # current heatmap as a png

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
        completed_logs.txt, current_counts.csv, and heatmap_output.png. If not, delete VEVisits_data and rerun.""")

if not os.path.exists(os.path.join(completed_stuff_folder,current_counts_path)):  # first execution of program
    # Create an empty dataframe with RPM and MAP buckets as columns
    df_counts = pd.DataFrame(index=rpmBreaks[1:], columns=mapBreaks[1:])
    print(df_counts)
    df_counts = df_counts.fillna(0)  # Initialize all counts to 0
    df_counts.to_csv(os.path.join(completed_stuff_folder,current_counts_path), index=True, header=True)
else:
    try:
        df_counts = pd.read_csv(os.path.join(completed_stuff_folder,current_counts_path), index_col=0)
        # Convert indices to integers
        df_counts.index = df_counts.index.astype(int)
        # Convert all data to floats
        df_counts.columns = df_counts.columns.astype(float)
    except Exception as e:
        print("""Check that the VEVisits_data folder exists in this directory and that the contents include 
        completed_logs.txt, current_counts.csv, and heatmap_output.png. If not, delete VEVisits_data and rerun.""")

"""For a given row of a log, find out which cell in the VE Map it corresponds to."""
def get_bucket(row, buckets, channel):
    val = row[channel] # will be either RPM or MAP
    my_buckets = buckets.copy()[1:] #remove dummy bin 
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
def update_counts(row):
    global count
    rpm_bucket = get_bucket(row, rpmBreaks, 'RPM[RPM]')
    map_bucket = get_bucket(row, mapBreaks, 'Plenum_MAP[kPa]')

    if row["Steady_State_Op[bool]"] == 1.0: # adaptation possible
        # increment the corresponding cell in the heatmap dataframe
        df_counts.loc[rpm_bucket, map_bucket] += 1
        count += 1

"""Takes in a path to an mdf file and returns a formatted dataframe to be processed in parse_data()"""
def deal_with_mdf(input_path):
    # Read mdf file
    mdf_file = mdfreader.Mdf(input_path)

    # convert to csv temporarily
    mdf_file.export_to_csv('temp.csv')
    
    # try to read the csv with one set of channels. use xls nomenclature
    df = pd.read_csv('temp.csv', low_memory=False, usecols=["t", "RPM", "Plenum_MAP", "Steady_State_Op"])[1:]
    df.rename(columns={'t': 't[s]',
                        'RPM': 'RPM[RPM]',
                        "Plenum_MAP" : "Plenum_MAP[kPa]",
                        'Steady_State_Op': 'Steady_State_Op[bool]'}, inplace=True)
    #finally:
    #    os.remove('temp.csv') # get rid of this temporary file
        
    # if these reading steps don't work, this file will be aborted
    df = df.apply(pd.to_numeric) # make all values numbers
    #df['CLVEAdapting[bool]'] = np.where(df['CLVEAdapting[bool]'] > 0.5, 1.0, 0.0) #if booleans are not boolean, interpolate
    
    return df

"""Checks if a filepath has already been successfully read, just in another format."""
def already_done(filepath):
    #for pot_filepath in completed_logs:
    #    if filepath[:-4] in pot_filepath[:-4]:
    #        return True
    return False

"""Parses a given file path by filetype. Updates the heatmap dataframe accordingly."""
def parse_data(input_path):
    global count
    sensor_data = deal_with_mdf(input_path)
    print("Done parsing.")
    
    # start counting how many cells in this file
    count_before = count
    
    # calculate the visits and add to heatmap dataframe
    sensor_data.apply(update_counts, axis=1)
    count_after = count
    
    print("Done adding to heatmap. " + str(count_after-count_before) + " cells were visited in this file.\n")

"""Makes and saves the heatmap."""
def make_heatmap():
    print("Generating heatmap.")
    
    title = 'Visits Per Cell in VE Map (RPM vs MAP)'
    
    fig, ax = plt.subplots(figsize=(25, 15))  # set the figure size (adjust as needed)
    ax.xaxis.tick_top()  # move x-axis ticks to the top
    heatmap = sns.heatmap(df_counts, annot=True, cmap='coolwarm', fmt='g', cbar = False, ax=ax)
    # adjust axis tick label orientation and size
    heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=0, fontsize=13)
    heatmap.set_yticklabels(heatmap.get_yticklabels(), rotation=0, fontsize=13)
    
    plt.suptitle(title, fontsize=24, fontweight='bold', y=0.95)
    
    # Save the heatmap to a file
    plt.savefig(os.path.join(completed_stuff_folder,current_heatmap_path), bbox_inches='tight')
    plt.show()


# -

"""Runs all code in this file."""
def main():
    # parse all filepaths in the logs folder
    for filepath in filepaths:
        if not already_done(filepath):
            try:
                print("Now parsing: " + filepath)
                parse_data(filepath)
                completed_logs.append(filepath)
            except ValueError as e:
                print(e)
                print(filepath + " could not be parsed. The required channels (Time, RPM, Plenum MAP, Steady State, and CLVEAdapting) may not be available in this log.\n")
                continue
    # temp file can stick around sometimes
    #if os.path.exists('temp.csv'):
        #os.remove('temp.csv')

    # store the seen logs
    with open(os.path.join(completed_stuff_folder,completed_logs_path), 'w') as file:
        for filepath in completed_logs:
            file.write(f"{filepath}\n")

    # store the current heatmap df
    df_counts.to_csv(os.path.join(completed_stuff_folder,current_counts_path), index=True, header=True)

    # display the heatmap
    make_heatmap()

main()



