# loading in some basic packages
import numpy as np                           # basic data operations
import pandas as pd                          # data wrangling
import os                                    # for foldering  
import itertools                             # for list combinations  
from tqdm import tqdm                        # for progress bars

# loading packages for dynamic time warpingfrom dtw import dtw
from dtw import dtw


def retrieve_keypoints():
    # identify which column indices in a dataframe correspond to which body parts
    keypoint_thumb_right = ["X_RIGHT_THUMB_TIP", "Y_RIGHT_THUMB_TIP", "Z_RIGHT_THUMB_TIP"]
    keypoint_thumb_left = ["X_LEFT_THUMB_TIP", "Y_LEFT_THUMB_TIP", "Z_LEFT_THUMB_TIP"]
    keypoint_index_right =  ["X_RIGHT_INDEX_FINGER_TIP", "Y_RIGHT_INDEX_FINGER_TIP", "Z_RIGHT_INDEX_FINGER_TIP"]
    keypoint_index_left = ["X_LEFT_INDEX_FINGER_TIP", "Y_LEFT_INDEX_FINGER_TIP", "Z_LEFT_INDEX_FINGER_TIP"]
    keypoint_middle_right = ["X_RIGHT_MIDDLE_FINGER_TIP", "Y_RIGHT_MIDDLE_FINGER_TIP", "Z_RIGHT_MIDDLE_FINGER_TIP"]
    keypoint_middle_left = ["X_LEFT_MIDDLE_FINGER_TIP", "Y_LEFT_MIDDLE_FINGER_TIP", "Z_LEFT_MIDDLE_FINGER_TIP"]
    keypoint_ring_right = ["X_RIGHT_RING_FINGER_TIP", "Y_RIGHT_RING_FINGER_TIP", "Z_RIGHT_RING_FINGER_TIP"]
    keypoint_ring_left = ["X_LEFT_RING_FINGER_TIP", "Y_LEFT_RING_FINGER_TIP", "Z_LEFT_RING_FINGER_TIP"]
    keypoint_pinky_right =  ["X_RIGHT_PINKY_FINGER_TIP", "Y_RIGHT_PINKY_FINGER_TIP", "Z_RIGHT_PINKY_FINGER_TIP"]
    keypoint_pinky_left = ["X_LEFT_PINKY_FINGER_TIP", "Y_LEFT_PINKY_FINGER_TIP", "Z_LEFT_PINKY_FINGER_TIP"]
    keypoint_wrist_right = ["X_RIGHT_WRIST_y", "Y_RIGHT_WRIST_y", "Z_RIGHT_WRIST_y"] # the y is added to the wrist keypoint because the wrist keypoint is present in both body and hands time series
    keypoint_wrist_left = ["X_LEFT_WRIST_y", "Y_LEFT_WRIST_y", "Z_LEFT_WRIST_y"]

    # make one 1 list with all keypoints
    keypoints = list(itertools.chain(keypoint_thumb_right, keypoint_thumb_left, 
                                    keypoint_index_right, keypoint_index_left, 
                                    keypoint_middle_right, keypoint_middle_left, 
                                    keypoint_ring_right, keypoint_ring_left, 
                                    keypoint_pinky_right, keypoint_pinky_left, 
                                    keypoint_wrist_right, keypoint_wrist_left))
    return keypoints


#make a function that takes in two timeseries and produces a normalized dtw distance
def dtw_distance_normalized_asymmetric(ts1, ts2):
    ts1 = np.array(ts1)
    ts2 = np.array(ts2)
    
    window_size = min(len(ts1), len(ts2)) // 10 # the window size is 10% of the length of the shortest timeseries
    if window_size < 3:
        window_size = 3 # the window size should be at least 3 frames

    res = dtw(ts1, ts2, 
              keep_internals = True,
              step_pattern = "asymmetric",
              open_begin = True, 
              open_end = True,
              window_args= {"window_type": "sakoechiba", "window_size": window_size})

    normalized_distance = res.normalizedDistance
    
    return normalized_distance


def dtw_distance_normalized(ts1, ts2):
    ts1 = np.array(ts1)
    ts2 = np.array(ts2)

    window_size = min(len(ts1), len(ts2)) // 10 # the window size is 10% of the length of the shortest timeseries
    if window_size < 3:
        window_size = 3 # the window size should be at least 3 frames

    res = dtw(ts1, ts2, 
                keep_internals = True,
                step_pattern = "symmetric2",
                window_args= {"window_type": "sakoechiba", "window_size": window_size})
    normalized_distance = res.normalizedDistance
    
    return normalized_distance


### Initiate the dtw distance function

# identify which column indices in a dataframe correspond to which body parts
keypoint_thumb_right = ["x_right_thumb_tip_relative", "y_right_thumb_tip_relative"]
keypoint_thumb_left = ["x_left_thumb_tip_relative", "y_left_thumb_tip_relative"]
keypoint_index_right =  ["x_right_index_finger_tip_relative", "y_right_index_finger_tip_relative"]
keypoint_index_left = ["x_left_index_finger_tip_relative", "y_left_index_finger_tip_relative"]
keypoint_middle_right = ["x_right_middle_finger_tip_relative", "y_right_middle_finger_tip_relative"]
keypoint_middle_left = ["x_left_middle_finger_tip_relative", "y_left_middle_finger_tip_relative"]
keypoint_ring_right = ["x_right_ring_finger_tip_relative", "y_right_ring_finger_tip_relative"]
keypoint_ring_left = ["x_left_ring_finger_tip_relative", "y_left_ring_finger_tip_relative"]
keypoint_pinky_right =  ["x_right_pinky_finger_tip_relative", "y_right_pinky_finger_tip_relative"]
keypoint_pinky_left = ["x_left_pinky_finger_tip_relative", "y_left_pinky_finger_tip_relative"]
keypoint_wrist_right = ["x_right_wrist", "y_right_wrist"] # the y is added to the wrist keypoint because the wrist keypoint is present in both body and hands time series
keypoint_wrist_left = ["x_left_wrist", "y_left_wrist"]

# make one 1 list with all keypoints
keypoints_relative = list(itertools.chain(keypoint_thumb_right, keypoint_thumb_left, 
                                 keypoint_index_right, keypoint_index_left, 
                                 keypoint_middle_right, keypoint_middle_left, 
                                 keypoint_ring_right, keypoint_ring_left, 
                                 keypoint_pinky_right, keypoint_pinky_left, 
                                 keypoint_wrist_right, keypoint_wrist_left))


def retrieve_relative_keypoints():
    return keypoints_relative


# make a dependent dtw such that each keypoint dtw distance is added up and divided by the number of keypoints
def dtw_distance_dependent(MT1, MT2, distance_array, hands_dtw):
    dis_left_thumb = dtw_distance_normalized(MT1[keypoint_thumb_left], MT2[keypoint_thumb_left])
    dis_right_thumb = dtw_distance_normalized(MT1[keypoint_thumb_right], MT2[keypoint_thumb_right])
    dis_left_index = dtw_distance_normalized(MT1[keypoint_index_left], MT2[keypoint_index_left])
    dis_right_index = dtw_distance_normalized(MT1[keypoint_index_right], MT2[keypoint_index_right])
    dis_left_middle = dtw_distance_normalized(MT1[keypoint_middle_left], MT2[keypoint_middle_left])
    dis_right_middle = dtw_distance_normalized(MT1[keypoint_middle_right], MT2[keypoint_middle_right])
    dis_left_ring = dtw_distance_normalized(MT1[keypoint_ring_left], MT2[keypoint_ring_left])
    dis_right_ring = dtw_distance_normalized(MT1[keypoint_ring_right], MT2[keypoint_ring_right])
    dis_left_pinky = dtw_distance_normalized(MT1[keypoint_pinky_left], MT2[keypoint_pinky_left])
    dis_right_pinky = dtw_distance_normalized(MT1[keypoint_pinky_right], MT2[keypoint_pinky_right])
    dis_left_wrist = dtw_distance_normalized(MT1[keypoint_wrist_left], MT2[keypoint_wrist_left])
    dis_right_wrist = dtw_distance_normalized(MT1[keypoint_wrist_right], MT2[keypoint_wrist_right])

    ## calculate the mean distance for the hands
    distance_left = np.mean([dis_left_thumb, dis_left_index, dis_left_middle, 
                             dis_left_ring, dis_left_pinky, dis_left_wrist])
    distance_right = np.mean([dis_right_thumb, dis_right_index, dis_right_middle, 
                              dis_right_ring, dis_right_pinky, dis_right_wrist])
    
    ### calculate the distance for opposite hands (left thumb with right thumb, left index with right index, etc.)
    if "_" in hands_dtw: # if the hands_dtw is 'left_right' or 'right_left'
        dis_left_right_thumb = dtw_distance_normalized(MT1[keypoint_thumb_left], MT2[keypoint_thumb_right])
        dis_left_right_index = dtw_distance_normalized(MT1[keypoint_index_left], MT2[keypoint_index_right])
        dis_left_right_middle = dtw_distance_normalized(MT1[keypoint_middle_left], MT2[keypoint_middle_right])
        dis_left_right_ring = dtw_distance_normalized(MT1[keypoint_ring_left], MT2[keypoint_ring_right])
        dis_left_right_pinky = dtw_distance_normalized(MT1[keypoint_pinky_left], MT2[keypoint_pinky_right])
        dis_left_right_wrist = dtw_distance_normalized(MT1[keypoint_wrist_left], MT2[keypoint_wrist_right])
        dis_right_left_thumb = dtw_distance_normalized(MT1[keypoint_thumb_right], MT2[keypoint_thumb_left])
        dis_right_left_index = dtw_distance_normalized(MT1[keypoint_index_right], MT2[keypoint_index_left])
        dis_right_left_middle = dtw_distance_normalized(MT1[keypoint_middle_right], MT2[keypoint_middle_left])
        dis_right_left_ring = dtw_distance_normalized(MT1[keypoint_ring_right], MT2[keypoint_ring_left])
        dis_right_left_pinky = dtw_distance_normalized(MT1[keypoint_pinky_right], MT2[keypoint_pinky_left])
        dis_right_left_wrist = dtw_distance_normalized(MT1[keypoint_wrist_right], MT2[keypoint_wrist_left])

        ## calculate the mean distance for the opposite hands
        distance_left_right = np.mean([dis_left_right_thumb, dis_left_right_index, dis_left_right_middle,
                                dis_left_right_ring, dis_left_right_pinky, dis_left_right_wrist])
        distance_right_left = np.mean([dis_right_left_thumb, dis_right_left_index, dis_right_left_middle,
                                    dis_right_left_ring, dis_right_left_pinky, dis_right_left_wrist])
    else:
        dis_left_right_thumb = np.nan
        dis_left_right_index = np.nan
        dis_left_right_middle = np.nan
        dis_left_right_ring = np.nan
        dis_left_right_pinky = np.nan
        dis_left_right_wrist = np.nan
        dis_right_left_thumb = np.nan
        dis_right_left_index = np.nan
        dis_right_left_middle = np.nan
        dis_right_left_ring = np.nan
        dis_right_left_pinky = np.nan
        dis_right_left_wrist = np.nan


    if hands_dtw == "both":
        distance = np.mean([distance_left, distance_right])
    elif hands_dtw == "left":
        distance = distance_left
    elif hands_dtw == "right":
        distance = distance_right
    elif hands_dtw == "left_right":
        distance = distance_left_right
    elif hands_dtw == "right_left":
        distance = distance_right_left

    #make a np array and append the each distance to the dataframe
    distance_array = np.append(distance_array, 
                               [distance,
                                dis_left_thumb, dis_right_thumb, dis_left_index, dis_right_index, 
                                dis_left_middle, dis_right_middle, dis_left_ring, dis_right_ring, 
                                dis_left_pinky, dis_right_pinky, dis_left_wrist, dis_right_wrist,
                                dis_left_right_thumb, dis_left_right_index, dis_left_right_middle,
                                dis_left_right_ring, dis_left_right_pinky, dis_left_right_wrist,
                                dis_right_left_thumb, dis_right_left_index, dis_right_left_middle,
                                dis_right_left_ring, dis_right_left_pinky, dis_right_left_wrist])

    return distance_array


def make_export_dtw_df(dtw_folder, ts_annot_folder, keypoints, anno):
    # make an empty dataframe to store the results
    df_distance = pd.DataFrame(columns=["pair", "comparison_id", "average_distance", 
                                        "left_thumb", "right_thumb", "left_index", "right_index", 
                                        "left_middle", "right_middle", "left_ring", "right_ring", 
                                        "left_pinky", "right_pinky", "left_wrist", "right_wrist",
                                        "left_right_thumb", "left_right_index", "left_right_middle",
                                        "left_right_ring", "left_right_pinky", "left_right_wrist",
                                        "right_left_thumb", "right_left_index", "right_left_middle",
                                        "right_left_ring", "right_left_pinky", "right_left_wrist"])

    # specify columns we want to keep in the timeseries dataframe (before merging with annotations)
    cols_to_keep = ["File", "Speaker", "comparison_id", "time"]
    cols_to_keep.extend(keypoints)

    error_count = 0
    error_files = []

    # using dtw to compare distances and show a warping line
    if os.path.exists(dtw_folder + "dtw_distance.csv"):
        print("DTW distance file already exists.")
    else:
        ts_annot_folder_files = [file for file in os.listdir(ts_annot_folder) if file.endswith(".csv")]
        for filename in tqdm(ts_annot_folder_files):
            pair = filename.split("_")[0]
            speaker = filename.split("_")[2]
            comparison_id = filename.split("_")[3]
            hands_dtw = anno[anno["comparison_id"] == int(comparison_id)]["hands_dtw"].values[0]
            distance_array = np.array([pair, comparison_id])

            if speaker == "ppB":
                continue
            else:
                ### make the earlier gesture as the referent gesture (makes a difference for asymmetric DTW)
                MT_A = pd.read_csv(ts_annot_folder + filename)
                MT_B = pd.read_csv(ts_annot_folder + filename.replace("ppA", "ppB"))
                if MT_A.loc[0, "time"] < MT_B.loc[0, "time"]:
                    MT1 = MT_A[cols_to_keep]
                    MT2 = MT_B[cols_to_keep]
                else:
                    MT1 = MT_B[cols_to_keep]
                    MT2 = MT_A[cols_to_keep]

                ### compute DTW distance
                try:
                    distance_array = dtw_distance_dependent(MT1, MT2, distance_array, hands_dtw)
                    df_distance = pd.concat([df_distance, pd.DataFrame([distance_array], columns=df_distance.columns)])
                except:
                    error_count += 1
                    error_files.append(filename)
                    # print(f"[DTW]Error in: {filename}. The file might contain missing values for the keypoints or have too few datapoints. Skipping this file...")
                    pass # do nothing and continue to the next line


        # sort the dataframe by comparison_id
        df_distance["comparison_id"] = df_distance["comparison_id"].astype(int)
        df_distance = df_distance.sort_values(by=["comparison_id"])

        df_distance = pd.merge(df_distance, anno, on="comparison_id", how="left")

        # save the dataframe to a csv file
        df_distance.to_csv(dtw_folder + "dtw_distance.csv", index=False)

        # check the shape of the dataframe
        print(f"The follwing {error_count} files were skipped. The files might contain missing values for the keypoints or have too few datapoints.")
        print(error_files)
        print(df_distance.shape)
    



def make_export_fake_dtw_df(dtw_folder, ts_annot_folder, keypoints, anno):
    # make an empty dataframe to store the results
    df_distance = pd.DataFrame(columns=["pair", "comparison_id", "average_distance",
                                        "left_thumb", "right_thumb", "left_index", "right_index", 
                                        "left_middle", "right_middle", "left_ring", "right_ring", 
                                        "left_pinky", "right_pinky", "left_wrist", "right_wrist"])

    # specify columns we want to keep in the timeseries dataframe (before merging with annotations)
    cols_to_keep = ["File", "Speaker", "comparison_id", "time"]
    cols_to_keep.extend(keypoints)

    error_count = 0
    error_files = []

    # using dtw to compare distances and show a warping line
    if os.path.exists(dtw_folder + "dtw_distance_fake.csv"):
        print("DTW distance file already exists.")
    else:
        ts_annot_folder_files = [file for file in os.listdir(ts_annot_folder) if file.endswith(".csv")]
        for filename in tqdm(ts_annot_folder_files):
            pair = filename.split("_")[0]
            speaker = filename.split("_")[2]
            comparison_id = filename.split("_")[3].split(".")[0]
            hands_dtw = anno[anno["comparison_id"] == int(comparison_id)]["hands_dtw"].values[0]
            fake_speaker = "ppB" if speaker == "ppA" else "ppA"
            try:
                fake_filename_list = [filename for filename in ts_annot_folder_files if (filename.split("_")[2] == fake_speaker) and (filename.split("_")[3] == comparison_id)]
                fake_filename = fake_filename_list[0]
            except:
                print(f"Could not find a fake file for {filename}. Skipping this file...")
                continue
            distance_array = np.array([pair, comparison_id])

            if speaker == "ppB":
                continue
            else:
                MT1 = pd.read_csv(ts_annot_folder + filename)
                MT2 = pd.read_csv(ts_annot_folder + fake_filename)
                # select only the columns that we want to keep
                MT1 = MT1[cols_to_keep]
                MT2 = MT2[cols_to_keep]
                try:
                    distance_array = dtw_distance_dependent(MT1, MT2, distance_array, hands_dtw)
                    df_distance = pd.concat([df_distance, pd.DataFrame([distance_array], columns=df_distance.columns)])
                except:
                    error_count += 1
                    error_files.append(filename)
                    # print(f"[DTW]Error in: {filename}. The file might contain missing values for the keypoints or have too few datapoints. Skipping this file...")
                    pass # do nothing and continue to the next line


        # sort the dataframe by comparison_id
        df_distance["comparison_id"] = df_distance["comparison_id"].astype(int)
        df_distance = df_distance.sort_values(by=["comparison_id"])

        df_distance = pd.merge(df_distance, anno, on="comparison_id", how="left")

        # save the dataframe to a csv file
        df_distance.to_csv(dtw_folder + "dtw_distance_fake.csv", index=False)

        print(f"The follwing {error_count} files were skipped. The files might contain missing values for the keypoints or have too few datapoints.")
        print(error_files)
        print(df_distance.shape)


def make_export_non_aligned_dtw_df(dtw_folder, ts_annot_folder, keypoints, anno):
    # make an empty dataframe to store the results
    df_distance = pd.DataFrame(columns=["pair", "comparison_id", "average_distance", 
                                        "left_thumb", "right_thumb", "left_index", "right_index", 
                                        "left_middle", "right_middle", "left_ring", "right_ring", 
                                        "left_pinky", "right_pinky", "left_wrist", "right_wrist",
                                        "left_right_thumb", "left_right_index", "left_right_middle",
                                        "left_right_ring", "left_right_pinky", "left_right_wrist",
                                        "right_left_thumb", "right_left_index", "right_left_middle",
                                        "right_left_ring", "right_left_pinky", "right_left_wrist"])

    # specify columns we want to keep in the timeseries dataframe (before merging with annotations)
    cols_to_keep = ["File", "Speaker", "comparison_id", "time"]
    cols_to_keep.extend(keypoints)

    error_count = 0
    error_files = []

    # using dtw to compare distances and show a warping line
    if os.path.exists(dtw_folder + "dtw_distance_non_aligned_gestures.csv"):
        print("DTW distance file already exists.")
    else:
        ts_annot_folder_files = [file for file in os.listdir(ts_annot_folder) if file.endswith(".csv")]
        for filename in tqdm(ts_annot_folder_files):
            pair = filename.split("_")[0]
            speaker = filename.split("_")[2]
            comparison_id = filename.split("_")[3]
            hands_dtw = anno[anno["comparison_id"] == int(comparison_id)]["hands_dtw"].values[0]
            distance_array = np.array([pair, comparison_id])

            if speaker == "ppB":
                continue
            else:
                ### make the earlier gesture as the referent gesture (makes a difference for asymmetric DTW)
                MT_A = pd.read_csv(ts_annot_folder + filename)
                MT_B = pd.read_csv(ts_annot_folder + filename.replace("ppA", "ppB"))
                if MT_A.loc[0, "time"] < MT_B.loc[0, "time"]:
                    MT1 = MT_A[cols_to_keep]
                    MT2 = MT_B[cols_to_keep]
                else:
                    MT1 = MT_B[cols_to_keep]
                    MT2 = MT_A[cols_to_keep]
                try:
                    distance_array = dtw_distance_dependent(MT1, MT2, distance_array, hands_dtw)
                    df_distance = pd.concat([df_distance, pd.DataFrame([distance_array], columns=df_distance.columns)])
                except:
                    error_count += 1
                    error_files.append(filename)
                    pass # do nothing and continue to the next line


        # sort the dataframe by comparison_id
        df_distance["comparison_id"] = df_distance["comparison_id"].astype(int)
        df_distance = df_distance.sort_values(by=["comparison_id"])

        df_distance = pd.merge(df_distance, anno, on="comparison_id", how="left")

        # save the dataframe to a csv file
        df_distance.to_csv(dtw_folder + "dtw_distance_non_aligned_gestures.csv", index=False)

        print(f"The follwing {error_count} files were skipped. The files might contain missing values for the keypoints or have too few datapoints.")
        print(error_files)
        print(df_distance.shape)