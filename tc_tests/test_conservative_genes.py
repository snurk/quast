import os
from common import *

name = os.path.basename(__file__)[5:-3]

run_quast(name, contigs=[contigs_1k_2], params='--find-conserved-genes')
assert_metric(name, 'Complete conserved orthologs (%)', ['0.00'])