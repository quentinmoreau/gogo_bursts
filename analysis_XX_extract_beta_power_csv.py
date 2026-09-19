import os
import numpy as np
import pandas as pd
from pathlib import Path

data_dir  = Path("/home/qmoreau/gogo_bursts/stats_new")
deriv_dir = Path("/home/common/bonaiuto/gogo_bursts/derivatives_v2/processed")
demo_dir  = Path("/home/common/bonaiuto/gogo_bursts/data_v2")

BASE_WINDOW = (-1.5, -1.0)
EPOCH_WINDOWS = {
    "STIM": (-1.5, 1.5),
    "RESP": (-1.5, 1.5),
}

demo = pd.concat([
    pd.read_csv(demo_dir / "GOGO_Demographics_2025_COMO.csv"),
    pd.read_csv(demo_dir / "GOGO_Demographics_2025_Driving.csv"),
], ignore_index=True)
demo = demo[demo["Status"] == "TD"]


def extract_epoch(epoch, window):
    burst_csv = data_dir / f"PC_12_{epoch}_trial_burst_counts.csv"
    out_csv   = data_dir / f"beta_power_{epoch}_contra_trials.csv"

    if out_csv.exists():
        print(f"Already exists: {out_csv}")
        return

    rt_df = (
        pd.read_csv(burst_csv)
        .query("tertile == 1")[["subject_id", "condition", "trial_idx", "response_time"]]
        .assign(subject_id=lambda d: d["subject_id"].astype(str))
    )

    all_rows = []

    for _, row in demo.iterrows():
        subject_id = str(row["ParticipantID"])
        fname      = deriv_dir / subject_id / f"{subject_id}_beta_power.npz"

        if not fname.exists():
            continue

        print(f"  [{epoch}] Loading {subject_id}")
        data       = np.load(fname, allow_pickle=True)
        subj_pow   = data["all_beta_pow"].item()
        power_time = data["time"]

        base_idx  = np.where((power_time >= BASE_WINDOW[0]) & (power_time <= BASE_WINDOW[1]))[0]
        epoch_idx = np.where((power_time >= window[0]) & (power_time <= window[1]))[0]
        epoch_time = power_time[epoch_idx]

        for condition in ["SHORT", "LONG"]:
            pow_mat = subj_pow["contra"][epoch][condition]  # trials x time
            if pow_mat is None or len(pow_mat) == 0:
                continue

            for trial_i in range(pow_mat.shape[0]):
                trial_pow = pow_mat[trial_i, :]
                base_mean = trial_pow[base_idx].mean()
                if base_mean == 0 or np.isnan(base_mean):
                    continue

                epoch_pow_bc = (trial_pow[epoch_idx] - base_mean) / base_mean * 100

                row_dict = {
                    "subject_id":    subject_id,
                    "condition":     condition,
                    "trial_idx":     trial_i,
                    "response_time": np.nan,
                }
                for t, v in zip(epoch_time, epoch_pow_bc):
                    row_dict[f"time_{round(t, 3)}"] = v

                all_rows.append(row_dict)

    df_out = pd.DataFrame(all_rows)

    df_out["subject_id"] = df_out["subject_id"].astype(str)
    df_out["trial_idx"]  = df_out["trial_idx"].astype(int)
    rt_df["trial_idx"]   = rt_df["trial_idx"].astype(int)

    df_out = df_out.drop(columns=["response_time"]).merge(
        rt_df, on=["subject_id", "condition", "trial_idx"], how="left"
    )

    df_out.to_csv(out_csv, index=False)
    print(f"Saved -> {out_csv}  ({len(df_out)} rows)")


for epoch, window in EPOCH_WINDOWS.items():
    extract_epoch(epoch, window)