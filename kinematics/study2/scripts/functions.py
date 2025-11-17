# loading in some basic packages
import numpy as np                           # basic data operations
import pandas as pd                          # data wrangling
import os                                    # for foldering       
import glob                                  # for file handling    
from tqdm import tqdm                        # for progress bars
from scipy.ndimage import gaussian_filter1d  # for smoothing
from scipy.spatial import distance           # for calculating distances
from dtw_functions import get_keypoints


### Functions to use when preprocessing of the time series data ###
def merge_ts(file, body_cols_to_keep): 
    # load in the body and hand time series
    try:
        body_ts = pd.read_csv(file + "_body.csv")
        hand_ts = pd.read_csv(file + "_hands.csv")
    except:
        # use engine='python' to avoid parser error
        # print(f"Using engine='python' to read the csv files: {file}")
        body_ts = pd.read_csv(file + "_body.csv", engine='python')
        hand_ts = pd.read_csv(file + "_hands.csv", engine='python')
        
    # keep only the columns we need
    body_ts = body_ts[body_cols_to_keep]
    # merge body and hand time series
    merged_ts = pd.merge(body_ts, hand_ts, on="time", how="left")

    # check if the number of rows in the merged time series is the same as the body and hand time series
    if len(merged_ts) != len(body_ts) or len(merged_ts) != len(hand_ts):
        print(f"Warning: The number of rows in the merged time series is different from the body and hand time series for file: {file}")

    return merged_ts

def fill_missing_values(df):
    for column in df.columns:
        df[column] = df[column].interpolate(method='linear', x=df['time'])
        # df[column] = df[column].interpolate(method='spline', order=3) # cubic spline interpolation
    return df

def adjust_aspect_ratio(df, width=16, height=9):
    # calculate scale factors
    scale_x = width / height
    
    # multiply x-coordinates by scale factor
    keypoint_names = [col for col in list(df.columns) if col.startswith('X_')]
    for keypoint in keypoint_names:
        df[keypoint] = df[keypoint] * scale_x

    return df

def flip_y_axis(df):
    y_cols = [col for col in df.columns if col.startswith('Y_')]
    for col in y_cols:
        df[col] = df[col] * -1
    return df

def normalize_size(df, target_height=1):
    # Get shoulder and hip coordinates
    left_shoulder = np.array([np.median(df['X_LEFT_SHOULDER']), 
                              np.median(df['Y_LEFT_SHOULDER']), 
                              np.median(df['Z_LEFT_SHOULDER'])])
    right_shoulder = np.array([np.median(df['X_RIGHT_SHOULDER']), 
                               np.median(df['Y_RIGHT_SHOULDER']), 
                               np.median(df['Z_RIGHT_SHOULDER'])])
    left_hip = np.array([np.median(df['X_LEFT_HIP']), 
                         np.median(df['Y_LEFT_HIP']), 
                         np.median(df['Z_LEFT_HIP'])])
    right_hip = np.array([np.median(df['X_RIGHT_HIP']), 
                          np.median(df['Y_RIGHT_HIP']), 
                          np.median(df['Z_RIGHT_HIP'])])
    
    # Calculate midpoints
    shoulder_midpoint = (left_shoulder + right_shoulder) / 2
    hip_midpoint = (left_hip + right_hip) / 2
    
    # Calculate current height (torso height)
    current_height = np.linalg.norm(shoulder_midpoint - hip_midpoint)
    
    # Calculate scaling factor
    scale_factor = target_height / current_height if current_height > 0 else 1.0
    df['scale_factor'] = scale_factor

    keypoint_names = [col for col in list(df.columns) if not col.startswith('visibility') and col!='time' and col!='scale_factor']
    # Apply scaling to all keypoints
    for keypoint in keypoint_names:
        df[keypoint] = df[keypoint] * scale_factor
        
    return df

def normalize_position(df, ref_point = "mid-torso"):
    keypoint_names = [col for col in list(df.columns) if not col.startswith('visibility') and col!='time']

    for idx, row in df.iterrows():
        if ref_point == 'mid-shoulder':
            ref_x = np.mean([row['X_LEFT_SHOULDER'], row['X_RIGHT_SHOULDER']])
            ref_y = np.mean([row['Y_LEFT_SHOULDER'], row['Y_RIGHT_SHOULDER']])
            ref_z = np.mean([row['Z_LEFT_SHOULDER'], row['Z_RIGHT_SHOULDER']])

        elif ref_point == 'mid-torso':
            midshoulder_x = np.mean([row['X_LEFT_SHOULDER'], row['X_RIGHT_SHOULDER']])
            midshoulder_y = np.mean([row['Y_LEFT_SHOULDER'], row['Y_RIGHT_SHOULDER']])
            midshoulder_z = np.mean([row['Z_LEFT_SHOULDER'], row['Z_RIGHT_SHOULDER']])
            midhip_x = np.mean([row['X_LEFT_HIP'], row['X_RIGHT_HIP']])
            midhip_y = np.mean([row['Y_LEFT_HIP'], row['Y_RIGHT_HIP']])
            midhip_z = np.mean([row['Z_LEFT_HIP'], row['Z_RIGHT_HIP']])
            ref_x = np.mean([midshoulder_x, midhip_x])
            ref_y = np.mean([midshoulder_y, midhip_y])
            ref_z = np.mean([midshoulder_z, midhip_z])
        
        for keypoint in keypoint_names:
            if keypoint.startswith('X_'):
                df.at[idx, keypoint] = row[keypoint] - ref_x
            elif keypoint.startswith('Y_'):
                df.at[idx, keypoint] = row[keypoint] - ref_y
            elif keypoint.startswith('Z_'):
                df.at[idx, keypoint] = row[keypoint] - ref_z
    return df
    

def preprocess_ts(file, body_cols_to_keep, merged_folder, interpolated_folder, norm_smooth_folder, flipped = False):
    filename = file.split("/")[-1]

    # Merge body and hands time series
    merged_df = merge_ts(file, body_cols_to_keep)
    merged_file = merged_folder + filename + "_merged.csv" if not flipped else merged_folder + "flipped/" + filename + "_merged.csv"
    merged_df.to_csv(merged_file, index=False)

    # Interpolate missing values
    interpolated_df = fill_missing_values(merged_df)
    interpolated_file = interpolated_folder + filename + "_interpolated.csv" if not flipped else interpolated_folder + "flipped/" + filename + "_interpolated.csv"
    # interpolated_df.to_csv(interpolated_file, index=False)

    # Apply smoothing
    xyz_cols = [col for col in interpolated_df.columns if col.startswith("X_") or col.startswith("Y_") or col.startswith("Z_")]
    smoothed_df = interpolated_df[xyz_cols].apply(lambda x: gaussian_filter1d(x, sigma=2))
    smoothed_df.insert(0, 'time', interpolated_df["time"])

    # Normalize time series
    normalized_df = normalize_size(smoothed_df.copy())
    normalized_df = normalize_position(normalized_df)

    # Adjust for aspect ratio
    normalized_df = adjust_aspect_ratio(normalized_df)

    # Flip y-axis
    normalized_df = flip_y_axis(normalized_df)

    # calculate relative finger positions
    normalized_df = calculate_relative_position(normalized_df, get_keypoints())
    normalized_df = change_wrist_col_name(normalized_df)

    # drop unnecessary columns
    xyz_cols = [col for col in normalized_df.columns if col.startswith("x_") or col.startswith("y_") or col.startswith("z_")]
    normalized_df = normalized_df[xyz_cols]
    normalized_df.insert(0, 'time', interpolated_df["time"])

    # Save normalized and smoothed time series
    norm_smooth_file = norm_smooth_folder + filename + "_ns.csv" if not flipped else norm_smooth_folder + "flipped/" + filename + "_ns.csv"
    normalized_df.to_csv(norm_smooth_file, index=False)


### Functions to use when merging the elan annotation with the time series data ###
def calculate_relative_position(df, keypoints):
    relative_keypoints = [keypoint for keypoint in keypoints if "WRIST" not in keypoint]
    for keypoint in relative_keypoints:
        df[keypoint.lower() + "_relative"] = df[keypoint] - df["_".join(keypoint.split("_")[:2]) + "_WRIST"]
    return df

def change_wrist_col_name(df):
    wrist_cols = [col for col in df.columns if "WRIST" in col]
    new_wrist_cols = [col.split("_y")[0].lower() for col in wrist_cols]
    df.rename(columns=dict(zip(wrist_cols, new_wrist_cols)), inplace=True)
    return df

def center_wrist(df):
    wrist_cols = [col for col in df.columns if "wrist" in col if "_change" not in col]
    for col in wrist_cols:
        change_col = col + "_centered"
        df[change_col] = df[col] - df[col].mean()
    return df


# this function loads in annotations and the original time of the timeseries dataframe, and returns annotations for the time series dataframe
def load_in_event(time_original, adata, col, speaker, i, adj_dur):
    output = np.full(len(time_original), np.nan, dtype=object)  # Initialize output array with NaN values
    speaker_1 = adata.loc[i, 'speaker_1']
    speaker_2 = adata.loc[i, 'speaker_2']
    if adj_dur == True:
        begin_1_col = 'begin_time_1_adj'
        end_1_col = 'end_time_1_adj'
        begin_2_col = 'begin_time_2_adj'
        end_2_col = 'end_time_2_adj'
    else:
        begin_1_col = 'begin_time_1'
        end_1_col = 'end_time_1'
        begin_2_col = 'begin_time_2'
        end_2_col = 'end_time_2'

    if speaker == 'A':
        if speaker_1 == 'A':
            output[(time_original >= adata.loc[i, begin_1_col]) & (time_original <= adata.loc[i, end_1_col])] = adata.iloc[i, adata.columns.get_loc(col)]
        elif speaker_2 == 'A':
            output[(time_original >= adata.loc[i, begin_2_col]) & (time_original <= adata.loc[i, end_2_col])] = adata.iloc[i, adata.columns.get_loc(col)]
    elif speaker == 'B':
        if speaker_1 == 'B':
            output[(time_original >= adata.loc[i, begin_1_col]) & (time_original <= adata.loc[i, end_1_col])] = adata.iloc[i, adata.columns.get_loc(col)]
        elif speaker_2 == 'B':
            output[(time_original >= adata.loc[i, begin_2_col]) & (time_original <= adata.loc[i, end_2_col])] = adata.iloc[i, adata.columns.get_loc(col)]
    
    return output


def export_merge_annot(MT_files, anno, ts_annot_folder, adj_dur=False):
    # merged_data_list = []
    keypoints = get_keypoints()
    merged_data = pd.DataFrame()  # Initialize an empty DataFrame

    n_comparisons = len(anno)
    print("Number of comparisons: " + str(n_comparisons))

    skip_count = 0
    skip_files = []

    for mt_file in tqdm(MT_files):
        fname = os.path.basename(mt_file).split('_ns')[0]
        # check if the processed file already exists
        if glob.glob(ts_annot_folder + fname + "*.csv"):
            skip_count += 1
            skip_files.append(fname)
            continue # Skip the file if it already exists. 

        speaker = fname.split('_')[1].upper() # Extract the speaker from the file name
        pair_n = int(fname.split('_')[0]) # Extract the pair number from the file name and convert it to an integer
        # anno["pairnr"] = anno["pair"].str[-2:].astype(int)
        # if pair_n not in anno['pairnr'].values:
        if pair_n not in anno['pair'].values:
            print("Pair number " + str(pair_n) + " not found in the annotation file.")
            continue

        # print("Now processing: " + mt_file + " for speaker " + speaker + "...")
        try:
            mdata = pd.read_csv(mt_file)
        except:
            # use engine='python' to avoid parser error
            mdata = pd.read_csv(mt_file, engine='python')
        adata = anno[anno['pair'] == pair_n].reset_index(drop=True)
        # mdata = calculate_relative_position(mdata, keypoints)
        # change_wrist_col_name(mdata)
        merged_data = mdata.copy()

        cols = ["comparison_id"]

        for i in range(len(adata)):
            comparison_id = adata.loc[i, 'comparison_id']
            hands_dtw = adata.loc[i, 'hands_dtw']
            # Initialize a pandas dataframe and a dictionary to hold new columns (this approach is faster than appending to a DataFrame in a loop)
            new_cols_df = pd.DataFrame()
            new_cols = {}
            # Add the new column to the dictionary
            new_cols['File'] = [fname] * len(merged_data)
            new_cols['Speaker'] = [speaker] * len(merged_data)

            # Apply the function to each column and store the result in the dictionary
            for col in cols:
                new_cols[col] = load_in_event(merged_data['time'], adata, col, speaker, i, adj_dur)

            # Create a new DataFrame from the dictionary
            new_cols_df = pd.DataFrame(new_cols)
            # Concatenate the original DataFrame with the new DataFrame
            final_merged_data = pd.concat([merged_data, new_cols_df], axis=1)
            final_merged_data = final_merged_data[final_merged_data['comparison_id'].notna()]
            # drop unnecesary columns and export the merged data to a csv file
            cols_to_drop = [col for col in final_merged_data.columns if "X" in col or "Y" in col or "Z" in col]
            final_merged_data.drop(columns=cols_to_drop, inplace=True)
            final_merged_data.to_csv(ts_annot_folder + fname + "_" + str(comparison_id) + "_" + hands_dtw + ".csv", index=False)
            
    if skip_count > 0:
        print(f"{skip_count} files skipped because they already exist")
        print(f"Files skipped: {skip_files}")



def export_merge_annot_size(MT_files, anno, ts_annot_folder):
    # merged_data_list = []
    keypoints = get_keypoints()
    merged_data = pd.DataFrame()  # Initialize an empty DataFrame

    n_gestures = len(anno)
    print("Number of gestures: " + str(n_gestures))

    skip_count = 0
    skip_files = []

    for mt_file in tqdm(MT_files):
        fname = os.path.basename(mt_file).split('_ns')[0]
        # check if the processed file already exists
        if glob.glob(ts_annot_folder + fname + "*.csv"):
            skip_count += 1
            skip_files.append(fname)
            continue # Skip the file if it already exists. 

        speaker = fname.split('_')[1].upper() # Extract the speaker from the file name
        pair_n = int(fname.split('_')[0]) # Extract the pair number from the file name and convert it to an integer
        # anno["pairnr"] = anno["pair"].str[-2:].astype(int)
        # if pair_n not in anno['pairnr'].values:
        if pair_n not in anno['pair'].values:
            print("Pair number " + str(pair_n) + " not found in the annotation file.")
            continue

        # print("Now processing: " + mt_file + " for speaker " + speaker + "...")
        try:
            mdata = pd.read_csv(mt_file)
        except:
            # use engine='python' to avoid parser error
            mdata = pd.read_csv(mt_file, engine='python')
        adata = anno[anno['pair'] == pair_n].reset_index(drop=True)
        mdata = calculate_relative_position(mdata, keypoints)
        change_wrist_col_name(mdata)
        merged_data = mdata.copy()

        cols = ["comparison_id"]

        for i in range(len(adata)):
            gesturer = adata.loc[i, 'gesturer']
            if speaker != gesturer:
                continue

            comparison_id = adata.loc[i, 'comparison_id']
            hands = adata.loc[i, 'A_hands'] if gesturer == 'A' else adata.loc[i, 'B_hands']
            # Initialize a pandas dataframe and a dictionary to hold new columns (this approach is faster than appending to a DataFrame in a loop)
            new_cols_df = pd.DataFrame()
            new_cols = {}
            # Add the new column to the dictionary
            new_cols['File'] = [fname] * len(merged_data)
            new_cols['Speaker'] = [speaker] * len(merged_data)

            # Apply the function to each column and store the result in the dictionary
            for col in cols:
                time_original = merged_data['time']
                output = np.full(len(time_original), np.nan, dtype=object)  # Initialize output array with NaN values
                output[(time_original >= adata.loc[i, 'begin_time']) & (time_original <= adata.loc[i, 'end_time'])] = adata.iloc[i, adata.columns.get_loc(col)]
                new_cols[col] = output

            # Create a new DataFrame from the dictionary
            new_cols_df = pd.DataFrame(new_cols)
            # Concatenate the original DataFrame with the new DataFrame
            final_merged_data = pd.concat([merged_data, new_cols_df], axis=1)
            final_merged_data = final_merged_data[final_merged_data['comparison_id'].notna()]
            # drop unnecesary columns and export the merged data to a csv file
            cols_to_drop = [col for col in final_merged_data.columns if "X" in col or "Y" in col or "Z" in col]
            final_merged_data.drop(columns=cols_to_drop, inplace=True)
            final_merged_data.to_csv(ts_annot_folder + fname + "_" + str(comparison_id) + "_" + hands + ".csv", index=False)
            
    if skip_count > 0:
        print(f"{skip_count} files skipped because they already exist")
        print(f"Files skipped: {skip_files}")


def make_export_size_df(output_folder, ts_annot_folder, keypoints, aligned=False):
    output_filename = "gesture_size.csv" if not aligned else "aligned_gesture_size.csv"
    # make an empty dataframe to store the results
    df_size = pd.DataFrame(columns=["pair", "comparison_id", "hands",
                                    "average_size", "size_left", "size_right"])

    # specify columns we want to keep in the timeseries dataframe (before merging with annotations)
    cols_to_keep = ["File", "Speaker", "comparison_id", "time"]
    cols_to_keep.extend(keypoints)

    error_count = 0
    error_files = []

    if os.path.exists(output_folder + output_filename):
        print(f"The file {output_folder + output_filename} already exists. Skipping the export.")
        
    else:
        ts_annot_folder_files = [file for file in os.listdir(ts_annot_folder) if file.endswith(".csv")]
        for filename in tqdm(ts_annot_folder_files):
            pair = filename.split("_")[0]
            comparison_id = filename.split("_")[2]
            hands = filename.split("_")[3].split(".")[0]
            size_array = np.array([pair, comparison_id, hands])

            MT = pd.read_csv(ts_annot_folder + filename)
                
            try:
                ### retrieve the min and max for the wrists
                min_left = (MT["x_left_wrist"].min(), MT["y_left_wrist"].min(), MT["z_left_wrist"].min())
                max_left = (MT["x_left_wrist"].max(), MT["y_left_wrist"].max(), MT["z_left_wrist"].max())
                min_right = (MT["x_right_wrist"].min(), MT["y_right_wrist"].min(), MT["z_right_wrist"].min())
                max_right = (MT["x_right_wrist"].max(), MT["y_right_wrist"].max(), MT["z_right_wrist"].max())
                
                ### calculate the Euclidean distance for each hand
                left_size = distance.euclidean(min_left, max_left)
                right_size = distance.euclidean(min_right, max_right)
                average_size = (left_size + right_size) / 2 if hands == "both" else left_size if hands == "left" else right_size

                size_array = [pair, comparison_id, hands, average_size, left_size, right_size]
                
                ### append the size_array to the df_size dataframe
                df_size = pd.concat([df_size, pd.DataFrame([size_array], columns=df_size.columns)])

            except:
                error_count += 1
                error_files.append(filename)
                pass # do nothing and continue to the next line


        # sort the dataframe by comparison_id
        df_size["comparison_id"] = df_size["comparison_id"].astype(int)
        df_size = df_size.sort_values(by=["comparison_id"])

        # save the dataframe to a csv file
        df_size.to_csv(output_folder + output_filename, index=False)

        # check the shape of the dataframe
        print(f"The follwing {error_count} files were skipped. The files might contain missing values for the keypoints or have too few datapoints.")
        print(error_files)
        print(df_size.shape)