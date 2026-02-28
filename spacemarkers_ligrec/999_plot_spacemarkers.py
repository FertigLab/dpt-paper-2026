import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append('..')
from xsample.xsample import plot_hist


sm = pd.read_csv("FIBROBLASTS_near_PDAC_sig_LR_scores.csv", index_col=0)
df = sm[sm.columns[sm.columns.str.startswith('HC')]]
plot_hist(df, title="SpaceMarkers Significant Fibroblasts near PDAC")


sm = pd.read_csv("PDAC_near_FIBROBLASTS_sig_LR_scores.csv", index_col=0)
df = sm[sm.columns[sm.columns.str.startswith('HC')]]
plot_hist(df, title="SpaceMarkers Significant PDAC near FIBROBLASTS")