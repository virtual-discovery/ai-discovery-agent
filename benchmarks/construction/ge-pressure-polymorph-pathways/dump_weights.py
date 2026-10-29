"""Extract the DeePMD network weights from Fantasy-frz.pb (read from the archive zip) into the lab's
potential/dp_weights.npz + dp_meta.json, using the minimal protobuf reader pbparse.py (no TensorFlow).

    python3 -B dump_weights.py [zip]
"""
import json
import sys
import zipfile
from pathlib import Path
import numpy as np
from pbparse import graph_from_bytes

HERE = Path(__file__).resolve().parent
ZIP = Path(sys.argv[1] if len(sys.argv) > 1 else
           Path.home()/'.cache/ai-theorist/ge-pressure-polymorph-pathways/allfiles.zip')
OUT = HERE.parents[1]/'tasks/ge-pressure-polymorph-pathways/environment/lab/potential'


def main():
    raw = zipfile.ZipFile(ZIP).read('allfiles/Fantasy-frz.pb')
    g = graph_from_bytes(raw)
    script = json.loads(g['train_attr/training_script'][2][0].decode())
    d0, d1 = script['model']['descriptor']['list']
    assert d0['type'] == 'se_e3' and d1['type'] == 'se_e2_a' and script['model']['fitting_net']['atom_ener'] == [-4.52]
    W = {}
    for k, (op, inp, c, a) in g.items():
        if op == 'Const' and c is not None and hasattr(c, 'dtype') and c.dtype == np.float64 and \
                (k.startswith(('filter_type', 'layer_', 'final_layer')) or k.endswith(('t_std', 't_avg'))):
            W[k] = np.asarray(c, float)
    assert not np.any(W['descrpt_attr_0/t_avg']) and not np.any(W['descrpt_attr_1/t_avg'])
    meta = {'d0': {'rcut': d0['rcut'], 'rcut_smth': d0['rcut_smth'], 'sel': d0['sel'][0], 'neuron': d0['neuron']},
            'd1': {'rcut': d1['rcut'], 'rcut_smth': d1['rcut_smth'], 'sel': d1['sel'][0], 'neuron': d1['neuron'],
                   'axis': d1['axis_neuron']},
            'bias_atom_e': float(g['fitting_attr/t_bias_atom_e'][2][0]),
            'source': 'allfiles/Fantasy-frz.pb of Materials Cloud archive 238ad-1bg44 (CC BY 4.0)',
            'model_version': g['model_attr/model_version'][2][0].decode()}
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT/'dp_weights.npz', **W)
    (OUT/'dp_meta.json').write_text(json.dumps(meta, indent=1) + '\n')
    print(len(W), 'arrays;', meta)


if __name__ == '__main__':
    main()
