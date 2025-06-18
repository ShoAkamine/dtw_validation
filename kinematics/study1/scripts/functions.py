# loading in some basic packages
import numpy as np                           # basic data operations
import pandas as pd                          # data wrangling
import os                                    # for foldering       
import glob                                  # for file handling    
from tqdm import tqdm                        # for progress bars
from dtw_functions import retrieve_keypoints

### Functions to use when preprocessing of the time series data ###
def merge_ts(file, body_cols_to_keep): 
    # load in the body and hand time series
    try:
        body_ts = pd.read_csv(file + "_body.csv")
        hand_ts = pd.read_csv(file + "_hands.csv")
    except:
        # use engine='python' to avoid parser error
        print(f"Using engine='python' to read the csv files: {file}")
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


### Functions to use when merging the elan annotation with the time series data ###
def calculate_relative_position(df, keypoints):
    relative_keypoints = [keypoint for keypoint in keypoints if "WRIST" not in keypoint]
    for keypoint in relative_keypoints:
        df[keypoint.lower() + "_relative"] = df[keypoint] - df["_".join(keypoint.split("_")[:2]) + "_WRIST"]
    return df

def change_wrist_col_name(df):
    wrist_cols = [col for col in df.columns if "WRIST" in col]
    new_wrist_cols = [col.lower() for col in wrist_cols]
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
    keypoints = retrieve_keypoints()
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

        speaker = fname.split('pp')[1] # Extract the speaker from the file name
        pair_n = int(fname.split('pair')[1].split('_')[0]) # Extract the pair number from the file name and convert it to an integer
        anno["pairnr"] = anno["pair"].str[-2:].astype(int)
        if pair_n not in anno['pairnr'].values:
            # print("Pair number " + str(pair_n) + " not found in the annotation file.")
            continue

        # print("Now processing: " + mt_file + " for speaker " + speaker + "...")
        mdata = pd.read_csv(mt_file)
        adata = anno[anno['pairnr'] == pair_n].reset_index(drop=True)
        mdata = flip_sign_for_y(mdata)
        mdata = calculate_relative_position(mdata, keypoints)
        change_wrist_col_name(mdata)
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

