#!/usr/bin/env python3
"""Recount first-column identifiers within each molecular category.
Run from the archive root with --data-root pointing to externally retrieved inputs.
No expression values, statistical models or differential results are changed.
"""
import argparse,csv,hashlib,json
from pathlib import Path
from collections import defaultdict

def main():
 p=argparse.ArgumentParser();p.add_argument('--archive-root',type=Path,default=Path(__file__).resolve().parents[1]);p.add_argument('--data-root',required=True,type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 reg=list(csv.DictReader((a.archive_root/'data/entries.tsv').open(),delimiter='\t'))
 assets=json.loads((a.archive_root/'data/input_assets.json').read_text())['files']
 hashes={x['path']:x['sha256'] for x in assets}
 all_raw=defaultdict(set);all_folded=defaultdict(lambda:defaultdict(set));entry_rows=[]
 labels={'Protein':'Proteome','Phosphoprotein':'Phosphoproteome','Metabolite':'Metabolites'}
 for row in reg:
  rel=row['Data File'];path=a.data_root/rel
  digest=hashlib.sha256(path.read_bytes()).hexdigest()
  if hashes.get(rel)!=digest:raise ValueError(f"Checksum mismatch: entry {row['dataset_id']}, {rel}")
  category=labels.get(row['Molecule Type'],row['Molecule Type'].removeprefix('cfrna-'))
  names=set()
  with path.open(encoding='utf-8',newline='') as f:
   next(f,None)
   for line in f:
    sep='\t' if '\t' in line else ','
    name=next(csv.reader([line],delimiter=sep),[''])[0].strip()
    if name and name.upper()!='ID':names.add(name)
  for name in names:all_folded[category][name.casefold()].add(name)
  all_raw[category].update(names)
  entry_rows.append({'dataset_id':int(row['dataset_id']),'molecular_category':category,'registered_molecule_type':row['Molecule Type'],'input_path':rel,'input_sha256':digest,'case_sensitive_identifiers':len(names),'case_insensitive_identifiers':len({n.casefold() for n in names})})
 a.out.mkdir(parents=True,exist_ok=False)
 def write(name,rows):
  with (a.out/name).open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
 order=['mRNA','lncRNA','miRNA','Proteome','Phosphoproteome','Metabolites']
 summary=[{'molecular_category':c,'distinct_identifiers':len(all_folded[c]),'matrix_entries':sum(r['molecular_category']==c for r in entry_rows),'case_sensitive_identifiers':len(all_raw[c]),'case_only_duplicates_removed':len(all_raw[c])-len(all_folded[c])} for c in order]
 write('counts.tsv',summary);write('entry_audit.tsv',entry_rows)
 unique=[{'molecular_category':c,'casefold_identifier':k,'original_spellings':' | '.join(sorted(v))} for c in order for k,v in sorted(all_folded[c].items())]
 write('identifier_union.tsv',unique)
 note={'registry_entries_checked':len(reg),'all_input_sha256_verified':True,'rule':'Read first column after header; trim surrounding whitespace; discard empty identifiers and ID sentinel; Unicode casefold; set union within molecular category. cfRNA entries contribute to their corresponding RNA categories. No cross-category merging, gene alias mapping, isoform collapsing, precursor/arm collapsing or biological entity deduplication. circRNA is audited per entry but excluded from the six-category plot; no harmonized junction catalog is asserted.','counts':{r['molecular_category']:r['distinct_identifiers'] for r in summary}}
 (a.out/'counting_method.json').write_text(json.dumps(note,indent=2)+'\n')
 print(json.dumps(note,indent=2))
if __name__=='__main__':main()
