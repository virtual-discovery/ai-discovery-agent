"""Write design.json for fresh-li-bcc-soft-shear (decision points, graded fields, calibration from calibration.json).

    python3 -B benchmarks/construction/fresh-li-bcc-soft-shear/write_design.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/fresh-li-bcc-soft-shear'
CAL = json.loads((HERE/'calibration.json').read_text())
NOTES = json.loads((HERE/'design_notes.json').read_text())


def calibration():
    out = {}
    for field, row in CAL.items():
        if 'robustness' not in row:
            out[field] = row
            continue
        nearest = sorted(row['ablations_failing_this_field'].items(), key=lambda kv: abs(kv[1]['error']))[:3]
        out[field] = {'reference': row['reference'], 'tolerance': row['tolerance'],
                      'robustness': {k: round(v['value'], 5) for k, v in row['robustness'].items()},
                      'max_robustness_error': round(row['max_robust_error'], 5),
                      'nearest_failing_ablations': {k: round(v['value'], 5) for k, v in nearest}}
    return out


def main():
    design = {
        'task_id': 'fresh-li-bcc-soft-shear', 'catalog_id': 'fresh-li-bcc-soft-shear', 'start': 'fresh',
        'source_group': 'fresh_li_bcc_elastic',
        'title': 'bcc Li: soft tetragonal shear, Zener anisotropy and compression (PBE)',
        'decision_points': NOTES['decision_points'],
        'graded_fields': NOTES['graded_fields'],
        'length_estimate': NOTES['length_estimate'],
        'runtime_dependencies': [
            'Quantum ESPRESSO 6.7 pw.x, OpenMPI, numpy (all in the runtime image; no extra layer)',
            'One pseudopotential file vendored in environment/lab/ (0.38 MB)'],
        'data_provenance': {
            'archive_url': 'https://archive.materialscloud.org/records/rcyfm-68h65',
            'doi': '10.24435/materialscloud:f3-ym',
            'files_used': ['SSSP_1.3.0_PBE_precision.tar.gz', 'SSSP_1.3.0_PBE_precision.json', 'LICENSE.txt'],
            'md5': {'SSSP_1.3.0_PBE_precision.tar.gz': 'fde94756886f32ada7bf597547557eb5',
                    'SSSP_1.3.0_PBE_precision.json': '1692c5c9ce89e1c7c783f8f0eee0cbfa',
                    'LICENSE.txt': 'c0f9df05fb4cd4921f4caa25606815ef'},
            'vendored_sha256': {'li_pbe_v1.4.uspp.F.UPF': '02cc4b3810e28a43b570277840d9631ca24f634d922928200c047b72c0a86bf5'},
            'license': 'Archive CC BY 4.0; the GBRV 1.4 Li UPF file GPL v3 (per the archive LICENSE.txt, vendored as SSSP_LICENSE.txt)',
            'hidden_truth': 'Construction calculations through task_lab.run (reference_jobs.py, run_batch.sh, '
                            'reference_results.json), stored only in environment/lab/private/'},
        'exceptions': {},
        'notes': NOTES['notes'],
        'calibration': calibration(),
    }
    (TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')
    print('design.json written:', len(design['decision_points']), 'decision points')


if __name__ == '__main__':
    main()
