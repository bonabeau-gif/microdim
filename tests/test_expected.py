from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
def test_hcom2_thresholds():
 n=pd.read_csv(ROOT/'expected/hcom2_nmf_threshold.csv'); assert n.loc[n['rank']==52,'seed_pass_fraction'].iloc[0]==.5
 s=pd.read_csv(ROOT/'expected/hcom2_support_strength.csv'); assert s.loc[s.latent_dimension==6,'seed_pass_fraction'].iloc[0]==.6
def test_ho():
 h=pd.read_csv(ROOT/'expected/ho_support_controls.csv'); assert h.loc[h.support_model=='Species-label shuffle','primary_error'].iloc[0]>h.loc[h.support_model=='Published support','primary_error'].iloc[0]
