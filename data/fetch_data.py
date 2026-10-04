"""Download Experiment 4 of Heffner et al. (2021) and rebuild the analysis files.

Source: https://github.com/jpheffne/epe (data/study4_mdd/). The upstream repository
carries no licence, so the raw data are not redistributed here.

Steps (ported from the original preprocessing notebook):
  1. Exclude participants whose rating of "neutral" on the emotion circumplex falls
     outside the preregistered 100 x 100 px square around the centre (508 -> 351).
  2. Divide each signal by its SD over the 7,020 retained trials (no centring,
     population SD), as sklearn.preprocessing.scale(with_mean=False) does.
  3. Split by CES-D group (cutoff 16) into study4_HC.csv (205) and study4_dep.csv (146).

The output is checked against SHA-256 hashes of the files the paper was fitted on.
Usage:  python data/fetch_data.py
"""
import hashlib
import io
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

COMMIT = "3031323d4140067b102d6db4fc5e7312038fa745"
BASE = f"https://raw.githubusercontent.com/jpheffne/epe/{COMMIT}/data/study4_mdd/"
OUT = Path(__file__).resolve().parent
SHA256 = {
    "study4_HC.csv": "290b2ab9a3559702334159a21118ddf4ac8e69a44d80c615eacf8aa6f9f86b0c",
    "study4_dep.csv": "d51d78d0a9b0890e5843a93d802f1b3e66f090a4b60f2d43f71e7759545fb44f",
}
SCALE = ["EPE_X1", "EPE_Y1", "RPE", "EP_X1", "EP_Y1", "RP",
         "feedback_X1", "feedback_Y1", "unfairness"]


def get(name):
    with urllib.request.urlopen(BASE + name) as r:
        return pd.read_csv(io.BytesIO(r.read()))


def main():
    beh, emo = get("study4_ug.csv"), get("study4_ec.csv")

    check = emo[emo["emotion"] == "neutral"].copy()
    inside = check["X1"].between(-50, 50) & check["Y1"].between(-50, 50)
    check["exclude"] = np.where(inside, "keep", "remove")
    beh = pd.merge(beh, check[["sub", "exclude"]], on="sub").query('exclude == "keep"')

    for col in SCALE:
        new = "unfairness_norm" if col == "unfairness" else f"{col}_scale"
        x = beh[col].to_numpy(dtype=float)
        beh[new] = x / x.std(ddof=0)

    files = {
        "study4_HC.csv": beh[beh["cesDepress"] == "Healthy_Controls"],
        "study4_dep.csv": beh[beh["cesDepress"] == "Clinically_Depressed"],
    }
    ok = True
    for name, df in files.items():
        df.to_csv(OUT / name)
        h = hashlib.sha256((OUT / name).read_bytes()).hexdigest()
        match = SHA256[name] is None or h == SHA256[name]
        ok &= match
        print(f"{name}: {df['sub'].nunique()} participants, {len(df)} trials, "
              f"sha256 {h[:16]}... {'OK' if match else 'MISMATCH'}")
    if not ok:
        sys.exit("Checksum mismatch: output differs from the files used in the paper.")


if __name__ == "__main__":
    main()
