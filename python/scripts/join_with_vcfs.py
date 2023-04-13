import logging
import os
import sys

import hail as hl

dirname = os.path.dirname(__file__)
test_resources_dir = os.path.join(dirname, '../test/resources')


def join_with_vcfs(dense_path, out_path):
    logging.info(f"densifying GVS VCFs and writing to {out_path}")
    vcf_mt = hl.import_vcf(os.path.join(test_resources_dir, 'gvs-vcfs', '*.vcf.gz'),
                           force_bgz=True, reference_genome='GRCh38').key_rows_by('locus')

    vds_mt = hl.read_matrix_table(dense_path).key_rows_by('locus')

    joined = hl.experimental.full_outer_join_mt(vcf_mt, vds_mt)

    joined = joined.rename(
        {
            'left_col': 'vcf_col',
            'right_col': 'vds_col',
            'left_row': 'vcf_row',
            'right_row': 'vds_row',
            'left_entry': 'vcf_entry',
            'right_entry': 'vds_entry',
        })

    joined.write(out_path, overwrite=True)


[_, dense_path, out_path, tmp_dir] = sys.argv

hl.init(tmp_dir=tmp_dir, log=os.path.join(tmp_dir, 'hail.log'))

join_with_vcfs(dense_path, out_path)