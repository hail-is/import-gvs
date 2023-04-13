import logging
import os
from import_gvs import import_gvs
import hail as hl
import sys


dirname = os.path.dirname(__file__)
test_resources_dir = os.path.join(dirname, "../build")


def run_import(out_path, tmpdir):
    logging.info(f"running import and writing VDS to {out_path}")
    def file_to_list(path):
        path = os.path.join(test_resources_dir, 'avro', path)
        a = []
        with hl.hadoop_open(path, 'r') as f:
            for line in f:
                a.append(os.path.join(test_resources_dir, line.strip()))
        return a

    refs = file_to_list('ref_files.txt')
    sample_mapping = file_to_list('sample_mapping.txt')
    site_filtering = file_to_list('site_filtering.txt')
    vet = file_to_list('vet_files.txt')
    vqsr_filt = file_to_list('vqsr_filtering.txt')
    vqsr_tranche = file_to_list('vqsr_tranche.txt')

    import_gvs(
        refs=[refs],
        vets=[vet],
        sample_mapping=sample_mapping,
        site_filtering_data=site_filtering,
        vqsr_filtering_data=vqsr_filt,
        vqsr_tranche_data=vqsr_tranche,
        final_path=out_path,
        tmp_dir=tmpdir,
        partitions_per_sample=4,
        unphase=True)

[_, vds_path, tmp_dir] = sys.argv

hl.init(tmp_dir=tmp_dir, log=os.path.join(tmp_dir, 'hail.log'))

run_import(vds_path, tmp_dir)