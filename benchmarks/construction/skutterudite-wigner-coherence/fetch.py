"""Refetch the pinned members of the Materials Cloud archive 5jsfj-cdg70 (no full download).

    python3 -B fetch.py [DEST]        # default DEST: ~/.cache/ai-theorist/skutterudite-wigner-coherence/raw

Archive: materials_cloud_archive_submission.zip (1,213,559,335 bytes, md5 9dc2bd3b8659022e5703862f82adaa40),
https://archive.materialscloud.org/records/5jsfj-cdg70 (CC BY 4.0). Each member is read with HTTP Range
requests (remotezip.py), checked against the zip CRC-32 and size, and then against members.json
(sha256), which this script writes on the first run and verifies afterwards. About 0.61 GB.
"""
import hashlib
import json
import sys
from pathlib import Path
import remotezip

HERE = Path(__file__).resolve().parent
DEST = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/skutterudite-wigner-coherence/raw'
P = 'materials_cloud_archive_submission/'
RUNS = {
    'CoSb3/tc_calc_no_isotope': 'Co_based/CoSb3/5_phono3py/tc_calc_no_isotope/',
    'CoSb3/tc_calc_with_isotope': 'Co_based/CoSb3/5_phono3py/tc_calc_with_isotope/',
    'CoSb3/v_operator_mean_free_path_isotope': 'Co_based/CoSb3/5_phono3py/v_operator_mean_free_path_isotope/',
    'YbCo4Sb12/hiphive': 'Co_based/YbCo4Sb12/3_hiphive_FCS/5_phono3py/',
    'FeSb3/hiphive': 'Fe_based/FeSb3/4_hiphive_FCS/5_phono3py/',
    'YbFe4Sb12-NM/hiphive_fc2hip_fc3hip': 'Fe_based/YbFe4Sb12/NM/6_hiphive/5_phono3py/fc2hip_fc3hip/',
    'YbFe4Sb12-NM/thirdorder_no_isotope': 'Fe_based/YbFe4Sb12/NM/3_FCS_thirdorder/5_phono3py/no_isotope/',
    'YbFe4Sb12-NM/thirdorder_isotope': 'Fe_based/YbFe4Sb12/NM/3_FCS_thirdorder/5_phono3py/isotope/',
    'IrSb3/hiphive': 'Ir_based/IrSb3/3_hiphive_FCS/5_phono3py/',
    'YbIr4Sb12/hiphive': 'Ir_based/YbIr4Sb12/3_hiphive_FCS/5_phono3py/',
}
RUN_FILES = {
    'CoSb3/tc_calc_no_isotope': ['kappa-m151515.hdf5', 'phono3py.yaml', '4_cmd_launch_tc.sh', 'tc_15x15x15_CoSb3_test.out'],
    'CoSb3/tc_calc_with_isotope': ['kappa-m151515.hdf5', 'phono3py.yaml', '1_cmd_launch_tc_workstation.sh', 'tc_CoSb3.out'],
    'CoSb3/v_operator_mean_free_path_isotope': ['kappa-m151515.hdf5', 'phono3py.yaml', '1_compute_velocity_operator.sh',
                                                'nohup.out', 'tc_data.dat', '2_plot_MFP.py'],
    'YbCo4Sb12/hiphive': ['kappa-m999.hdf5', 'phono3py.yaml', '1_cmd_launch_tc_workstation.sh', 'tc_YbCo4Sb12.out'],
    'FeSb3/hiphive': ['kappa-m999.hdf5', 'phono3py.yaml', '1_cmd_launch_tc_workstation.sh', 'tc_FeSb3.out'],
    'YbFe4Sb12-NM/hiphive_fc2hip_fc3hip': ['kappa-m999.hdf5', 'phono3py.yaml', '1_cmd_launch_tc_workstation.sh', 'tc_YbFe4Sb12.out'],
    'YbFe4Sb12-NM/thirdorder_no_isotope': ['kappa-m999.hdf5', 'phono3py.yaml', '3_cmd_launch_tc_workstation.sh', 'tc_YbFe4Sb12.out'],
    'YbFe4Sb12-NM/thirdorder_isotope': ['kappa-m999.hdf5', 'phono3py.yaml', '3_cmd_launch_tc_workstation.sh', 'tc_YbFe4Sb12.out',
                                        'tc_YbFe4Sb12_with_isotope.out'],
    'IrSb3/hiphive': ['kappa-m999.hdf5', 'phono3py.yaml', '1_cmd_launch_tc_workstation.sh', 'tc_IrSb3.out'],
    'YbIr4Sb12/hiphive': ['kappa-m999.hdf5', 'phono3py.yaml', '1_cmd_launch_tc_workstation.sh', 'tc_YbIr4Sb12.out'],
}
HARMONIC = {
    'CoSb3': ('Co_based/CoSb3/7_phonon_PDOS/', 'Co_based/CoSb3/6_MSD/'),
    'YbCo4Sb12': ('Co_based/YbCo4Sb12/5_phonon_PDOS/', 'Co_based/YbCo4Sb12/4_MSD/'),
    'FeSb3': ('Fe_based/FeSb3/6_phonon_PDOS/', 'Fe_based/FeSb3/5_MSD/'),
    'YbFe4Sb12-NM': ('Fe_based/YbFe4Sb12/NM/7_phonon_PDOS/', 'Fe_based/YbFe4Sb12/NM/4_MSD/'),
    'IrSb3': ('Ir_based/IrSb3/5_phonon_PDOS/', 'Ir_based/IrSb3/4_MSD/'),
    'YbIr4Sb12': ('Ir_based/YbIr4Sb12/5_phonon_PDOS/', 'Ir_based/YbIr4Sb12/4_MSD/'),
}
EXTRA = ['README.txt',
         'Fe_based/YbFe4Sb12/AFM/2_phonons/post_processing_matdyn_crystal/matdyn__crystal.freq',
         'Fe_based/YbFe4Sb12/NM/2_phonons/post_processing_matdyn_crystal/matdyn__crystal.freq',
         'Fe_based/YbFe4Sb12/AFM/1_vc-relax/vc_relax.in']


def member_list():
    names = [RUNS[r] + f for r, files in RUN_FILES.items() for f in files]
    for comp, (pdos, msd) in HARMONIC.items():
        names += [pdos + 'projected_dos.dat', pdos + 'phonopy.yaml', msd + 'POSCAR',
                  msd + 'thermal_displacements_calcs/mesh_50/thermal_displacements.yaml',
                  msd + 'thermal_displacements_calcs/mesh_50/1_cmd_MSD.sh']
        if comp != 'YbFe4Sb12-NM':          # this PDOS folder has no command files
            names += [pdos + 'cmd.sh', pdos + 'settings.conf']
    return [P + n for n in names + EXTRA]


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def main():
    names = member_list()
    record = HERE/'members.json'
    known = json.loads(record.read_text()) if record.exists() else {}
    directory = remotezip.central_directory()
    out = {}
    for name in names:
        info = directory[name]
        path = DEST/name[len(P):]
        if path.exists() and path.stat().st_size == info['usize'] and (name not in known or sha256(path.read_bytes()) == known[name]['sha256']):
            data = path.read_bytes()
        else:
            data = remotezip.read_member(info)          # checks size and CRC-32
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
            print('fetched', name, info['usize'])
        digest = sha256(data)
        if name in known and known[name]['sha256'] != digest:
            raise SystemExit(f'checksum mismatch for {name}')
        out[name] = {'bytes': info['usize'], 'crc32': '%08x' % info['crc'], 'sha256': digest}
    if set(out) != set(known):          # first run, or members added to the list: record the new checksums
        record.write_text(json.dumps(out, indent=1) + '\n')
    print(f'{len(out)} members verified in {DEST}')


if __name__ == '__main__':
    main()
