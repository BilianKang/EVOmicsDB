#!/usr/bin/env python3
"""Rebuild one entry using the release analysis and entry-specific configuration."""
import argparse, hashlib, importlib.util, json, subprocess, sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()
def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--id',type=int,required=True);p.add_argument('--data-root',type=Path,required=True)
 p.add_argument('--out',type=Path,required=True);p.add_argument('--rscript',default='Rscript')
 p.add_argument('--grouping-profile',choices=['original','minimal'],default='original');p.add_argument('--dry-run',action='store_true')
 a=p.parse_args();out=a.out.resolve();data=a.data_root.resolve()
 if out.exists():p.error('Use a new output directory; existing results are never overwritten.')
 rows=json.loads((ROOT/'config/differential_analysis.json').read_text())['datasets']
 row=next((dict(r) for r in rows if r['dataset_id']==a.id),None)
 if row is None:p.error('Entry ID must be in the 110-entry registry.')
 assets={r['path']:r for r in json.loads((ROOT/'data/input_assets.json').read_text())['files']}
 checked=[]
 def resolve(rel):
  item=assets[rel];f=data/rel
  if not f.is_file():raise FileNotFoundError(f)
  actual=digest(f)
  if actual!=item['sha256']:raise ValueError(f'Input SHA-256 mismatch: {rel}')
  checked.append({'path':rel,'sha256':actual});return str(f)
 counts=resolve(row['data_file'])
 if a.grouping_profile=='minimal':
  p.error('This repository does not include separate grouping exports. Use the original group file from the external inputs.')
 elif row.get('group_file'):row['group_file']=resolve(row['group_file'])
 if row.get('probe_map_file'):row['probe_map_file']=resolve(row['probe_map_file'])
 f=ROOT/'scripts/differential_config.py'
 spec=importlib.util.spec_from_file_location('evomics_release_config',f);mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)
 cfg=mod.DifferentialDatasetConfig.from_dict(row);params=cfg.r_arguments()
 params.update({'--counts_file':counts,'--output_rda':str(out/f'ID{a.id}.rda'),'--output_all_csv':str(out/f'ID{a.id}_all.csv'),'--output_csv':str(out/f'ID{a.id}_significant.csv'),'--is_output_csv':True})
 script=ROOT/'R/differential_analysis.r';command=[a.rscript,str(script)]
 for k,v in params.items():
  if v is not None:command.extend([k,('TRUE' if v else 'FALSE') if isinstance(v,bool) else str(v)])
 record={'dataset_id':a.id,'grouping_profile':a.grouping_profile,'inputs':checked,'script_sha256':digest(script),'effective_configuration':row,'command':command}
 if a.dry_run:print(json.dumps(record,indent=2));return 0
 out.mkdir(parents=True)
 with (out/'run.log').open('w') as log:r=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT)
 record['returncode']=r.returncode;record['outputs']={f.name:digest(f) for f in out.iterdir() if f.is_file()}
 (out/'run.json').write_text(json.dumps(record,indent=2)+'\n')
 print(f'ID{a.id}: exit {r.returncode}; see {out / "run.log"}')
 return r.returncode
if __name__=='__main__':sys.exit(main())
