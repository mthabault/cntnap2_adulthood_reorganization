import pandas as pd

# Define paths
input_path = "paste path here"
output_path = "paste path here"

# Load data
df = pd.read_csv(input_path)
df = df.sort_values(["video", "mouse", "start_time"]).reset_index(drop=True)

# Classify behaviours
exploratory = {"digging", "supported_rearing", "unsupported_rearing"}
rsb = {"scratching", "grooming_long", "grooming_short", "head_body_twitch"}
empty = "empty"

def get_class(beh):
    if beh in exploratory:
        return "exploratory"
    if beh in rsb:
        return "rsb"
    if beh == empty:
        return "empty"
    return "other"

# Define transition classifier
def classify_transition(b1_row, next_rows):
    t1_end = b1_row["end_time"]
    b1 = b1_row["behavior"]
    class_b1 = get_class(b1)

    for idx, row in next_rows.iterrows():
        beh2 = row["behavior"]
        start2 = row["start_time"]
        delta = start2 - t1_end

        # empty segment
        if beh2 == empty:
            if row["duration"] >= 3:
                return "isolated", 0, 0
            else:
                continue

        # real behaviour
        if delta >= 3:
            return "isolated", 0, 0

        class_b2 = get_class(beh2)

        # repeat flag
        repeat_flag = 1 if beh2 == b1 else 0

        # rsb destination flag
        rsb_dest_flag = 1 if beh2 in rsb else 0

        # same vs diff class
        if class_b1 == class_b2 and class_b1 not in {"empty", "other"}:
            class_cat = "same"
        else:
            class_cat = "diff"

        return class_cat, repeat_flag, rsb_dest_flag

    return "isolated", 0, 0

# Compute transitions per mouse
results = []

for (video, mouse, genotype, sex, age), subdf in df.groupby(["video", "mouse", "genotype", "sex", "age"]):
    real_events = subdf[subdf["behavior"] != empty]
    total_events = len(real_events)
    if total_events == 0:
        continue

    counts_class = {"same": 0, "diff": 0, "isolated": 0}
    count_repeat = 0
    count_rsb_dest = 0

    for idx, row in real_events.iterrows():
        next_rows = subdf[subdf.index > idx]
        class_cat, repeat_flag, rsb_dest_flag = classify_transition(row, next_rows)

        counts_class[class_cat] += 1
        count_repeat += repeat_flag
        count_rsb_dest += rsb_dest_flag

    p_same = counts_class["same"] / total_events
    p_diff = counts_class["diff"] / total_events
    p_isolated = counts_class["isolated"] / total_events
    p_repeat = count_repeat / total_events
    p_rsb_dest = count_rsb_dest / total_events

    results.append({
        "video": video,
        "mouse": mouse,
        "genotype": genotype,
        "sex": sex,
        "age": age,
        "p_repeat": p_repeat,
        "p_same": p_same,
        "p_diff": p_diff,
        "p_isolated": p_isolated,
        "p_rsb_dest": p_rsb_dest,
        "p_sum_same_diff_isolated": p_same + p_diff + p_isolated,
        "n_events": total_events
    })

# Save dataset as CSV
results_df.to_csv(output_path, index=False)

print("Saved to:", output_path)
