#!/usr/bin/env python3
"""Audit the public Ho data files and record the primary evaluation-set IDs.

Assembly/control prediction is executed with the published model program.
"""
exec(open(__file__.replace('run_ho_validation.py','_bootstrap.py')).read())
from pathlib import Path
import hashlib,json,pandas as pd, argparse
parser=argparse.ArgumentParser()
parser.add_argument('--source-dir', type=Path, default=ROOT/'data/third_party/ho', help='Unpacked Ho et al. model archive')
a=parser.parse_args()
p=a.source_dir
if not p.exists(): raise SystemExit('Run scripts/fetch_sources.py --ho first')
lock=json.loads((ROOT/'sources.lock.json').read_text())['ho_resource_competition']
required=['BHI_assems.mat','BHI_niches.mat','PredictAssemblies_BHI.m']
missing=[name for name in required if not (p/name).is_file()]
if missing: raise SystemExit('Source directory does not contain the required Ho archive files: '+', '.join(missing))
(ROOT/'results').mkdir(exist_ok=True)
rows=[]
for f in sorted(p.glob('*')):
 if f.is_file():
  data=f.read_bytes()
  expected=lock.get('known_md5',{}).get(f.name)
  actual=hashlib.md5(data).hexdigest()
  if expected and actual!=expected: raise SystemExit('Checksum mismatch: '+f.name)
  rows.append({'file':f.name,'bytes':f.stat().st_size,'sha256':hashlib.sha256(data).hexdigest(),'md5':actual,'expected_md5':expected or ''})
pd.DataFrame(rows).to_csv(ROOT/'results/ho_source_inventory.csv',index=False)
json.dump({'primary_random_3_to_13_assembly_ids':list(range(106,170)),'n':64},open(ROOT/'results/ho_primary_set.json','w'),indent=2)
print('Ho source integrity inventory and primary 64-assembly identifiers written to results/.')
