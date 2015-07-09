#!/usr/bin/python

import os
from common import *

name = os.path.basename(__file__)[5:-3]
only_ns_in_contigs = 'only_Ns_in_sequence.fasta'


run_quast(name, contigs=[contigs_1k_1, only_ns_in_contigs, contigs_1k_2], params='-R ' + reference_1k)
check_report_files(name)
assert_report_header(name, [contigs_1k_1, only_ns_in_contigs, contigs_1k_2])
assert_metric(name, 'Largest contig', ['1000', '616', '760'])
assert_metric(name, 'Largest alignment', ['1000', '-', '760'])