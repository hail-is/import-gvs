import logging
import os
import sys

import hail as hl

dirname = os.path.dirname(__file__)

[_, vds_path, dense_path, tmp_dir] = sys.argv

logging.info(f"densifying VDS and writing to {dense_path}")
hl.init(tmp_dir=tmp_dir, log=os.path.join(tmp_dir, 'hail.log'))

vds = hl.vds.read_vds(vds_path)
vds.variant_data = vds.variant_data.drop('GT')
mt = hl.vds.to_dense_mt(vds)

mt.write(dense_path, overwrite=True)
