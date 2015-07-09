#!/usr/bin/python

import os
from common import *

name = os.path.basename(__file__)[5:-3]
incorrect_chars_contigs = 'incorrect_chars_in_sequence.fasta'

run_quast(name, contigs=[incorrect_chars_contigs])
assert_report_header(name, [incorrect_chars_contigs])  # without ref - no check for non-ACGTN characters
assert_metric(name, 'Largest contig', ['3980'])