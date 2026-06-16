# -*- coding: utf-8 -*-
"""
Created on Thu May  1 11:43:13 2025

@author: Mathieu_Thabault
"""
import numpy as np
import pandas as pd

# Define parameters 
csv_path = r"#paste path here"
time_column = 't'  # Event time column (in seconds)
max_lag = 500  # ms
bin_size = 5  # ms
lags = np.arange(-bin_size / 2, max_lag + bin_size, bin_size)
center_bins = (lags[:-1] + lags[1:]) / 2

# Keep only center_bins >= 0
pos_mask = center_bins >= 0
center_bins = center_bins[pos_mask]


df = pd.read_csv(csv_path)

# Add recording ID
df['recording_id'] = df['mouse'].astype(str) + '_' + df['cell_number'].astype(str)

# Define groups
groups = {
    "10w": {
        "WT": (df['genotype'] == 'wt') & (df['age'] == '10w'),
        "KO": (df['genotype'] == 'ko') & (df['age'] == '10w')
    },
    "20w": {
        "WT": (df['genotype'] == 'wt') & (df['age'] == '20w'),
        "KO": (df['genotype'] == 'ko') & (df['age'] == '20w')
    }
}

# Autocorrelogram function
def compute_autocorrelogram(event_times, max_lag=5000, bin_size=100):
    """
    Compute autocorrelogram using fast counting (searchsorted).
    event_times: sorted array of event times in ms.
    """
    diffs_all = []

    for i in range(len(event_times)):
        low = event_times[i] - max_lag
        high = event_times[i] + max_lag

        # Get indices of events in window
        start_idx = np.searchsorted(event_times, low, side='left')
        end_idx = np.searchsorted(event_times, high, side='right')

        # Get diffs (exclude self)
        diffs = event_times[start_idx:end_idx] - event_times[i]
        diffs = diffs[diffs != 0]  # remove zero

        diffs_all.append(diffs)

    if len(diffs_all) == 0:
        return None

    # Concatenate all diffs and make histogram
    diffs_all = np.concatenate(diffs_all)
    hist, _ = np.histogram(diffs_all, bins=lags)

    # Keep only positive lags
    hist = hist[pos_mask]
    norm = np.sum(hist)

    if norm > 0:
        return hist / norm
    else:
        return None

# Calculate autocorrelograms
results = {}

for age in groups.keys():
    results[age] = {}
    for genotype in groups[age].keys():
        subset = df[groups[age][genotype]]
        autocorr_all = []

        for rec_id, group_data in subset.groupby('recording_id'):
            event_times = np.sort(group_data[time_column].values) * 1000  # Convert to ms
            if len(event_times) < 2:
                continue

            autocorr = compute_autocorrelogram(event_times, max_lag=max_lag, bin_size=bin_size)
            if autocorr is not None:
                autocorr_all.append(autocorr)

        if len(autocorr_all) > 0:
            stack = np.array(autocorr_all)
            results[age][genotype] = {
                "mean": np.mean(stack, axis=0),
                "sem": np.std(stack, axis=0) / np.sqrt(stack.shape[0]),
                "n": stack.shape[0]
            }
        else:
            print(f"No valid recordings for {genotype} at {age}")

# Prep dataset for export
export_rows = []

for age in groups.keys():
    for genotype in groups[age].keys():
        subset = df[groups[age][genotype]]
        autocorr_all = []

        for rec_id, group_data in subset.groupby('recording_id'):
            event_times = np.sort(group_data[time_column].values) * 1000  # Convert to ms
            if len(event_times) < 2:
                continue

            autocorr = compute_autocorrelogram(event_times, max_lag=max_lag, bin_size=bin_size)
            if autocorr is not None:
                for lag, value in zip(center_bins, autocorr):
                    export_rows.append({
                        "recording_id": rec_id,
                        "age": age.lower(),
                        "genotype": genotype.lower(),
                        "lag_ms": lag,
                        "autocorr_value": value
                    })

# Create dataframe
df_export = pd.DataFrame(export_rows)

# Export to CSV
output_path = r"paste path here"
df_export.to_csv(output_path, index=False)

print(f"Saved individual autocorrelogram data to {output_path}")
