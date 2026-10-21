#!/usr/bin/env python3
from pathlib import Path
import argparse, hashlib, json, os, subprocess, sys, zipfile, io
import requests
ROOT=Path(__file__).resolve().parents[1]; LOCK=json.loads((ROOT/'sources.lock.json').read_text()); DATA=ROOT/'data'/'third_party'; DATA.mkdir(parents=True,exist_ok=True)
COHORTS=json.loads((ROOT/'config/cohorts.json').read_text())
def dl(url,path,md5=None):
    path.parent.mkdir(parents=True,exist_ok=True); r=requests.get(url,timeout=120); r.raise_for_status(); path.write_bytes(r.content)
    if md5 and hashlib.md5(r.content).hexdigest()!=md5: raise RuntimeError(f'MD5 mismatch {path.name}')
def fetch_borenstein():
    s=LOCK['borenstein_curated']; base=f"https://raw.githubusercontent.com/{s['repo']}/{s['commit']}/data/processed_data"
    for c in COHORTS:
        for fn in ['genera.tsv','mtb.tsv','metadata.tsv']:
            dl(f'{base}/{c}/{fn}',DATA/'borenstein'/c/fn)
def fetch_github_archive(key):
    s=LOCK[key]; url=f"https://codeload.github.com/{s['repo']}/zip/{s['commit']}"; r=requests.get(url,timeout=180); r.raise_for_status(); z=zipfile.ZipFile(io.BytesIO(r.content)); out=DATA/key; out.mkdir(exist_ok=True); z.extractall(out)
def fetch_ho():
    s=LOCK['ho_resource_competition']; meta=requests.get(f"https://zenodo.org/api/records/{s['record']}",timeout=120).json(); out=DATA/'ho'; out.mkdir(exist_ok=True)
    known=s.get('known_md5',{})
    for f in meta['files']:
        name=f['key']; url=f['links']['self']; md5=known.get(name) or (f.get('checksum','').split(':',1)[1] if f.get('checksum','').startswith('md5:') else None); dl(url,out/name,md5)
    (out/'zenodo_record.json').write_text(json.dumps(meta,indent=2))
if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--all',action='store_true'); ap.add_argument('--cohorts',action='store_true'); ap.add_argument('--hcom2',action='store_true'); ap.add_argument('--ho',action='store_true'); ap.add_argument('--han',action='store_true'); ap.add_argument('--jin',action='store_true'); a=ap.parse_args()
    if a.all or a.cohorts: fetch_borenstein()
    if a.all or a.hcom2: fetch_github_archive('mahajan_hcom2')
    if a.all or a.han: fetch_github_archive('han_metabolomics')
    if a.all or a.jin: fetch_github_archive('jin_hcom2')
    if a.all or a.ho: fetch_ho()
