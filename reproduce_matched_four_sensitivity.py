#!/usr/bin/env python3
"""Recompute the V8.9.1 post hoc matched-four descriptive sensitivity from frozen pairwise CSVs.
No model fitting or simulation is performed.
"""
from pathlib import Path
import pandas as pd, re
ROOT=Path(__file__).resolve().parent
COMMON={1001,2002,3003,5005}

def seed(x):
    m=re.search(r"seed(\d+)",str(x)); return int(m.group(1)) if m else None
rows=[]
def add(fn, metrics, family=None):
    df=pd.read_csv(ROOT/fn)
    if "seed_i" not in df:
        df["seed_i"]=df["identity_i"].map(seed); df["seed_j"]=df["identity_j"].map(seed)
    groups=df[family].drop_duplicates().tolist() if family else [None]
    for g in groups:
        sub=df if g is None else df[df[family]==g]
        p1=sub[sub.arm=="P1"]; p4=sub[sub.arm=="P4"]
        p1m=p1[p1.seed_i.isin(COMMON)&p1.seed_j.isin(COMMON)]
        p4m=p4[p4.seed_i.isin(COMMON)&p4.seed_j.isin(COMMON)]
        for metric in metrics:
            a=float(p1m[metric].mean()); b_all=float(p4[metric].mean()); b=float(p4m[metric].mean())
            rows.append(dict(source_file=fn,family=g or "aggregate",metric=metric,p1_common4_mean=a,p4_all5_mean=b_all,p4_common4_mean=b,pct_change_p4_all5_vs_p1=(b_all/a-1)*100,pct_change_p4_common4_vs_p1=(b/a-1)*100,p1_pair_count_common4=len(p1m),p4_pair_count_all5=len(p4),p4_pair_count_common4=len(p4m)))
add("G2_LOCATION_PAIRWISE.csv",["location_l2_distance"])
add("G2_PROJECTION_DISTANCE.csv",["projection_frobenius_distance"],"family")
add("G3_SIGMA_PAIRWISE.csv",["relative_frobenius"])
add("G3_CORRELATION_CONSEQUENCE.csv",["correlation_rmse"])
add("G3_SPECTRUM_CONSEQUENCE.csv",["spectrum_l2_distance"])
add("G3_COMPONENT_PAIRWISE.csv",["relative_frobenius_B","relative_frobenius_L"])
out=ROOT/"POSTHOC_MATCHED_FOUR_SENSITIVITY_V8_9_1.csv"
pd.DataFrame(rows).to_csv(out,index=False)
print(out)
