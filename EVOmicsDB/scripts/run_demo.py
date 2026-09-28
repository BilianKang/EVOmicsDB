#!/usr/bin/env python3
"""Run the current differential, volcano and ROC scripts on artificial inputs."""
import argparse
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--rscript', default='Rscript')
    args = parser.parse_args()
    out = args.out.resolve()
    if out.exists():
        parser.error('Use a new output directory; existing files are not overwritten.')
    out.mkdir(parents=True)
    scripts = ROOT / 'R'
    commands = [
        [args.rscript, str(scripts / 'differential_analysis.r'),
         '--counts_file', str(ROOT / 'data/example/synthetic_log2_expression.tsv'),
         '--group_file', str(ROOT / 'data/example/synthetic_groups.tsv'),
         '--sample_column', 'sample', '--group_column', 'group',
         '--case_code', 'C', '--control_code', 'N',
         '--analysis_profile', 'transcriptomics_normalized',
         '--input_scale', 'log2', '--normalization_todo', 'none',
         '--transform_todo', 'none', '--feature_filter_mode', 'none',
         '--model_imputation', 'none', '--visual_imputation', 'none',
         '--output_rda', str(out / 'demo.rda'),
         '--output_all_csv', str(out / 'demo_all.csv'),
         '--output_csv', str(out / 'demo_significant.csv')],
        [args.rscript, str(scripts / 'volcano_plot.r'),
         '--input_rda', str(out / 'demo.rda'),
         '--output_png', str(out / 'demo_volcano.png')],
        [args.rscript, str(scripts / 'roc.r'),
         '--utils_r', str(scripts / 'roc_utils.R'),
         '--input_rda', str(out / 'demo.rda'),
         '--group_mapping_case_label', 'case',
         '--group_mapping_control_label', 'control',
         '--group_mapping_source', 'synthetic_groups.tsv_C_case_N_control',
         '--custom_gene', 'DEMO_001', 'DEMO_011',
         '--output_png', str(out / 'demo_roc.png'),
         '--output_csv', str(out / 'demo_roc_coordinates.csv'),
         '--output_summary_csv', str(out / 'demo_roc_summary.csv')],
    ]
    records = []
    for i, command in enumerate(commands, 1):
        with (out / f'step_{i}.log').open('w') as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
        records.append({'command': command, 'returncode': result.returncode})
        (out / 'run.json').write_text(json.dumps(records, indent=2) + '\n')
        if result.returncode:
            raise SystemExit(f'Step {i} failed; inspect {out / f"step_{i}.log"}')
    print(f'PASS: synthetic differential, volcano and ROC runs; outputs: {out}')


if __name__ == '__main__':
    main()
