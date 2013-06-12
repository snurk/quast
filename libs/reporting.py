############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

from libs import qconfig

references = []

# Here you can modify content and order of metrics in QUAST reports and names of metrcis as well


class Order(object):
    def __init__(self, ref_free_fields, ref_dependant_fields, is_groupped=False):
        self.is_groupped = is_groupped
        self.ref_free_fields = ref_free_fields
        self.ref_dependant_fields = ref_dependant_fields


class Fields:
####################################################################################
###########################  CONFIGURABLE PARAMETERS  ##############################
####################################################################################
    ### for indent before submetrics
    TAB = '    '

    ### List of available fields for reports. Values (strings) should be unique! ###

    # Header
    NAME = 'Assembly'

    # Basic statistics
    NUMCONTIGS = '# contigs'
    CONTIGS = ('# contigs (>= %d bp)', tuple(qconfig.contig_thresholds))
    LARGCONTIG = 'Largest contig'
    TOTALLEN = 'Total length'
    TOTALLENS = ('Total length (>= %d bp)', tuple(qconfig.contig_thresholds))
    N50 = 'N50'
    N75 = 'N75'
    L50 = 'L50'
    L75 = 'L75'
    GC = 'GC (%)'

    # Misassemblies
    MISASSEMBL = '# misassemblies'
    MISCONTIGS = '# misassembled contigs'
    MISCONTIGSBASES = 'Misassembled contigs length'
    MISINTERNALOVERLAP = 'Misassemblies inter-contig overlap'
    ### additional list of metrics for detailed misassemblies report
    MIS_ALL_EXTENSIVE = '# misassemblies'
    MIS_RELOCATION = TAB + '# relocations'
    MIS_TRANSLOCATION = TAB + '# translocations'
    MIS_INVERTION = TAB + '# inversions'
    MIS_EXTENSIVE_CONTIGS = '# misassembled contigs'
    MIS_EXTENSIVE_BASES = 'Misassembled contigs length'
    MIS_LOCAL = '# local misassemblies'

    # Unaligned
    UNALIGNED = '# unaligned contigs'
    UNALIGNEDBASES = 'Unaligned contigs length'
    AMBIGUOUS = '# ambiguously mapped contigs'
    AMBIGUOUSEXTRABASES = 'Extra bases in ambiguously mapped contigs'
    MISLOCAL = '# local misassemblies'
    ### additional list of metrics for detailed unaligned report
    UNALIGNED_FULL_CNTGS = '# fully unaligned contigs'
    UNALIGNED_FULL_LENGTH = 'Fully unaligned length'
    UNALIGNED_PART_CNTGS = '# partially unaligned contigs'
    UNALIGNED_PART_WITH_MISASSEMBLY = TAB + '# with misassembly'
    UNALIGNED_PART_SIGNIFICANT_PARTS = TAB + '# both parts are significant'
    UNALIGNED_PART_LENGTH = 'Partially unaligned length'

    # Indels and mismatches
    MISMATCHES = '# mismatches'
    INDELS = '# indels'
    INDELSBASES = 'Indels length'
    SUBSERROR = '# mismatches per 100 kbp'
    INDELSERROR = '# indels per 100 kbp'
    MIS_SHORT_INDELS = TAB + '# short indels (<= %d bp)' % qconfig.SHORT_INDEL_THRESHOLD
    MIS_LONG_INDELS = TAB + '# long indels (> %d bp)' % qconfig.SHORT_INDEL_THRESHOLD
    UNCALLED = "# N's"
    UNCALLED_PERCENT = "# N's per 100 kbp"

    # Genome statistics
    MAPPEDGENOME = 'Genome fraction (%)'
    DUPLICATION_RATIO = 'Duplication ratio'
    GENES = '# genes'
    OPERONS = '# operons'
    AVGIDY = 'Average %IDY'  # Deprecated
    LARGALIGN = 'Largest alignment'
    NG50 = 'NG50'
    NA50 = 'NA50'
    NGA50 = 'NGA50'
    LG50 = 'LG50'
    LA50 = 'LA50'
    LGA50 = 'LGA50'
    NG75 = 'NG75'
    NA75 = 'NA75'
    NGA75 = 'NGA75'
    LG75 = 'LG75'
    LA75 = 'LA75'
    LGA75 = 'LGA75'

    # Predicted genes
    PREDICTED_GENES_UNIQUE = '# predicted genes (unique)'
    PREDICTED_GENES = ('# predicted genes (>= %d bp)', tuple(qconfig.genes_lengths))

    # Reference statistics
    REFLEN = 'Reference length'
    ESTREFLEN = 'Estimated reference length'
    REFGC = 'Reference GC (%)'
    REF_GENES = 'Reference genes'
    REF_OPERONS = 'Reference operons'

    # TODO: GC for aligned

    ### content and order of metrics in MAIN REPORT (<quast_output_dir>/report.txt, .tex, .tsv).
    # The first list is for reference-free metrics, the second is for each reference
    order = Order([NAME, NUMCONTIGS, CONTIGS, LARGCONTIG, TOTALLEN, TOTALLENS,
                   N50, N75, L50, L75, GC, ESTREFLEN, PREDICTED_GENES, PREDICTED_GENES_UNIQUE],

                  [REFLEN, REFGC, NG50, NG75, MISASSEMBL, MISCONTIGSBASES, MISLOCAL, UNALIGNED, UNALIGNEDBASES, MAPPEDGENOME,
                   DUPLICATION_RATIO, UNCALLED_PERCENT, SUBSERROR, INDELSERROR, GENES, OPERONS,
                   LARGALIGN, NA50, NGA50, NA75, NGA75])
    # order = [NAME, CONTIGS, TOTALLENS, NUMCONTIGS, LARGCONTIG, TOTALLEN, REFLEN, ESTREFLEN, GC, REFGC,
    #          N50, NG50, N75, NG75, MISASSEMBL, MISLOCAL, UNALIGNED, UNALIGNEDBASES, MAPPEDGENOME, DUPLICATION_RATIO,
    #          UNCALLED_PERCENT, SUBSERROR, INDELSERROR, GENES, OPERONS, PREDICTED_GENES_UNIQUE, PREDICTED_GENES,
    #          LARGALIGN, NA50, NGA50, NA75, NGA75]

    # content and order of metrics in DETAILED MISASSEMBLIES REPORT (<quast_output_dir>/contigs_reports/misassemblies_report.txt, .tex, .tsv)
    misassemblies_order = Order([NAME],
                                [MIS_ALL_EXTENSIVE, MIS_RELOCATION, MIS_TRANSLOCATION, MIS_INVERTION,
                                 MIS_EXTENSIVE_CONTIGS, MIS_EXTENSIVE_BASES, MIS_LOCAL, MISMATCHES,
                                 INDELS, MIS_SHORT_INDELS, MIS_LONG_INDELS, INDELSBASES])

    # content and order of metrics in DETAILED UNALIGNED REPORT (<quast_output_dir>/contigs_reports/unaligned_report.txt, .tex, .tsv)
    unaligned_order = Order([NAME],
                            [UNALIGNED_FULL_CNTGS, UNALIGNED_FULL_LENGTH, UNALIGNED_PART_CNTGS,
                             UNALIGNED_PART_WITH_MISASSEMBLY, UNALIGNED_PART_SIGNIFICANT_PARTS,
                             UNALIGNED_PART_LENGTH, UNCALLED])

    ### list of GAGE metrics (--gage option)
    GAGE_NUMCONTIGS = 'Contigs #'
    GAGE_MINCONTIG = 'Min contig'
    GAGE_MAXCONTIG = 'Max contig'
    GAGE_N50 = 'N50'
    GAGE_GENOMESIZE = 'Genome size'
    GAGE_ASSEMBLY_SIZE = 'Assembly size'
    GAGE_CHAFFBASES = 'Chaff bases'
    GAGE_MISSINGREFBASES = 'Missing reference bases'
    GAGE_MISSINGASMBLYBASES = 'Missing assembly bases'
    GAGE_MISSINGASMBLYCONTIGS = 'Missing assembly contigs'
    GAGE_DUPREFBASES = 'Duplicated reference bases'
    GAGE_COMPRESSEDREFBASES = 'Compressed reference bases'
    GAGE_BADTRIM = 'Bad trim'
    GAGE_AVGIDY = 'Avg idy'
    GAGE_SNPS = 'SNPs'
    GAGE_SHORTINDELS = 'Indels < 5bp'
    GAGE_LONGINDELS = 'Indels >= 5'
    GAGE_INVERSIONS = 'Inversions'
    GAGE_RELOCATION = 'Relocation'
    GAGE_TRANSLOCATION = 'Translocation'
    GAGE_NUMCORCONTIGS = 'Corrected contig #'
    GAGE_CORASMBLYSIZE = 'Corrected assembly size'
    GAGE_MINCORCONTIG = 'Min correct contig'
    GAGE_MAXCORCOTING = 'Max correct contig'
    GAGE_CORN50 = 'Corrected N50'

    # content and order of metrics in GAGE report (<quast_output_dir>/gage_report.txt, .tex, .tsv)
    gage_order = Order([NAME, GAGE_NUMCONTIGS, GAGE_MINCONTIG, GAGE_MAXCONTIG],
                       [GAGE_N50, GAGE_GENOMESIZE, GAGE_ASSEMBLY_SIZE,
                        GAGE_CHAFFBASES, GAGE_MISSINGREFBASES, GAGE_MISSINGASMBLYBASES,
                        GAGE_MISSINGASMBLYCONTIGS, GAGE_DUPREFBASES, GAGE_COMPRESSEDREFBASES, GAGE_BADTRIM,
                        GAGE_AVGIDY, GAGE_SNPS, GAGE_SHORTINDELS, GAGE_LONGINDELS, GAGE_INVERSIONS,
                        GAGE_RELOCATION, GAGE_TRANSLOCATION, GAGE_NUMCORCONTIGS, GAGE_CORASMBLYSIZE,
                        GAGE_MINCORCONTIG, GAGE_MAXCORCOTING, GAGE_CORN50])

    ### Grouping of metrics and set of main metrics for HTML version of main report
    grouped_order = Order([
            ('Basic statistics', [NUMCONTIGS, CONTIGS, LARGCONTIG, TOTALLEN, TOTALLENS,
                                  N50, N75, L50, L75, GC,]),
            ('Predicted genes', [PREDICTED_GENES_UNIQUE, PREDICTED_GENES,]),
        ], [
            ('Misassemblies', [MIS_ALL_EXTENSIVE, MIS_RELOCATION, MIS_TRANSLOCATION, MIS_INVERTION,
                               MIS_EXTENSIVE_CONTIGS, MIS_EXTENSIVE_BASES, MIS_LOCAL]),
            ('Unaligned', [UNALIGNED_FULL_CNTGS, UNALIGNED_FULL_LENGTH, UNALIGNED_PART_CNTGS,
                           UNALIGNED_PART_WITH_MISASSEMBLY, UNALIGNED_PART_SIGNIFICANT_PARTS,
                           UNALIGNED_PART_LENGTH,]),
            ('Mismatches', [MISMATCHES, INDELS, INDELSBASES, SUBSERROR, INDELSERROR,
                            MIS_SHORT_INDELS, MIS_LONG_INDELS, UNCALLED, UNCALLED_PERCENT,]),
            ('Genome statistics', [MAPPEDGENOME, DUPLICATION_RATIO, GENES, OPERONS, LARGALIGN,
                                   NG50, NG75, NA50, NA75, NGA50, NGA75, LG50, LG75, LA50, LA75, LGA50, LGA75,]),
            ('Reference statistics', [REFLEN, ESTREFLEN, REFGC, REF_GENES, REF_OPERONS,])
        ],
        is_groupped=True
    )

    # for a "short" version of HTML report
    main_metrics = [NUMCONTIGS, LARGCONTIG, TOTALLEN, N50,
                    MIS_ALL_EXTENSIVE, MIS_EXTENSIVE_BASES,
                    SUBSERROR, INDELSERROR, UNCALLED_PERCENT,
                    MAPPEDGENOME, DUPLICATION_RATIO, GENES, OPERONS, NGA50,
                    PREDICTED_GENES_UNIQUE, PREDICTED_GENES,]

####################################################################################
########################  END OF CONFIGURABLE PARAMETERS  ##########################
####################################################################################


    class Quality:
        MORE_IS_BETTER = 'More is better'
        LESS_IS_BETTER = 'Less is better'
        EQUAL = 'Equal'

    quality_dict = {
        Quality.MORE_IS_BETTER:
            [LARGCONTIG, TOTALLEN, TOTALLENS, N50, NG50, N75, NG75, NA50, NGA50, NA75, NGA75, LARGALIGN,
             MAPPEDGENOME, GENES, OPERONS, PREDICTED_GENES_UNIQUE, PREDICTED_GENES, AVGIDY],
        Quality.LESS_IS_BETTER:
            [NUMCONTIGS, CONTIGS, L50, LG50, L75, LG75,
             MISLOCAL, MISASSEMBL, MISCONTIGS, MISCONTIGSBASES, MISINTERNALOVERLAP,
             UNALIGNED, UNALIGNEDBASES, AMBIGUOUS, AMBIGUOUSEXTRABASES,
             UNCALLED, UNCALLED_PERCENT,
             LA50, LGA50, LA75, LGA75, DUPLICATION_RATIO, INDELS, INDELSERROR, MISMATCHES, SUBSERROR,
             MIS_SHORT_INDELS, MIS_LONG_INDELS, INDELSBASES],
        Quality.EQUAL:
            [REFLEN, ESTREFLEN, GC, REFGC],
        }

    for name, metrics in filter(lambda (name, metrics): name in ['Misassemblies', 'Unaligned'], grouped_order.ref_dependant_fields):
        quality_dict['Less is better'].extend(metrics)


#################################################

import os
import logging
from libs.qutils import print_timestamp, warning

####################################################################################
# Reporting module (singleton) for QUAST
#
# See class Fields to available fields for report.
# Usage from QUAST modules:
#  from libs import reporting
#  report = reporting.get(fasta_filename)
#  report.add_field(reporting.Field.N50, n50)
#
# Import this module only after final changes in qconfig!
#
####################################################################################

reports = {}  # basefilename -> Report
assemblies_order = []  # for printing in appropriate order

#################################################


def get_main_metrics():
    lists = map(take_tuple_metric_apart, Fields.main_metrics)
    m_metrics = []
    for l in lists:
        for m in l:
            m_metrics.append(m)
    return m_metrics


def take_tuple_metric_apart(field):
    metrics = []

    if isinstance(field, tuple): # TODO: rewrite it nicer
        thresholds = map(int, ''.join(field[1]).split(','))
        for i, feature in enumerate(thresholds):
            metrics.append(field[0] % feature)
    else:
        metrics = [field]

    return metrics


def get_quality(metric):
    for quality, metrics in Fields.quality_dict.iteritems():
        if metric in Fields.quality_dict[quality]:
            return quality
    return Fields.Quality.EQUAL


# Report for one assembly. Structure:
# common_fields = {
#   field1: value1,
#   field2: value2
# }
# reference_fields = {
#   reference1: {
#     field3: value3,
#     field4: value4,
#   }
#   reference2: {
#     field3: value5,
#     field4: value6,
#   }
# }
class Report(object):
    def __init__(self, name):
        self.basic_fields = {}
        self.reference_fields = dict([(ref, {}) for ref in references])

        self.add_field(Fields.NAME, name)

    def add_field(self, field, value, reference=None):
        assert field in Fields.__dict__.itervalues(), 'Unknown field: %s' % field
        if reference:
            self.reference_fields[reference][field] = value
        else:
            self.basic_fields[field] = value

    def append_field(self, field, value, reference=None):
        assert field in Fields.__dict__.itervalues(), 'Unknown field: %s' % field
        if reference:
            self.reference_fields[reference].setdefault(field, []).append(value)
        else:
            self.basic_fields.setdefault(field, []).append(value)

    def get_field(self, field, reference=None):
        assert field in Fields.__dict__.itervalues(), 'Unknown field: %s' % field
        if reference:
            return self.reference_fields[reference].get(field, None)
        else:
            return self.basic_fields.get(field, None)


def get(filename):
    filename = os.path.basename(filename)
    if filename not in assemblies_order:
        assemblies_order.append(filename)
    return reports.setdefault(filename, Report(filename))


def delete(filename, reference=None):
    filename = os.path.basename(filename)
    if filename in assemblies_order:
        assemblies_order.remove(filename)
    if filename in reports.keys():
        reports.pop(filename)


# ATTENTION! Contents numeric values, needed to be converted into strings
def get_table(order=Fields.order):
    # if not isinstance(order[0], tuple):  # is not a groupped metrics order
    #     order = [('', order)]

    table = []

    def append_line(reference, rows, field, pattern=None, feature=None, i=None):
        quality = get_quality(field)
        values = []

        for assembly_name in assemblies_order:
            report = get(assembly_name)
            value = report.get_field(field, reference)

            if feature is None or value is None:
                values.append(value)
            else:
                values.append(value[i] if i < len(value) else None)

        if filter(lambda v: value is not None, values):
            metric_name = field if (feature is None) else pattern % feature
            # ATTENTION! Contents numeric values, needed to be converted to strings.
            rows.append({
                'metricName': metric_name,
                'quality': quality,
                'values': values,
                'isMain': field in Fields.main_metrics,
            })

    for reference in [None] + references:
        groups = []
        table.append((reference.name if reference else None, groups))

        if reference:
            suborder = order.ref_dependant_fields
        else:
            suborder = order.ref_free_fields

        if not order.is_groupped:
            suborder = [('', suborder[:])]

        for group_name, fields in suborder:
            group_rows = []
            groups.append((group_name, group_rows))

            for field in fields:
                if isinstance(field, tuple):
                    for i, feature in enumerate(field[1]):
                        append_line(reference, group_rows, field, field[0], feature, i)
                else:
                    append_line(reference, group_rows, field)

    # if not isinstance(order[0], tuple):  # is not a groupped metrics order
    #     group_name, group_rows = table[0]
    #     return group_rows
    # else:
    return table


# def is_groupped_table(table):
#     return isinstance(table[0], tuple)
#
#
# def get_all_rows_out_of_table(table, is_transposed=False):
#     all_rows = []
#     if is_groupped_table(table):
#         for group_name, rows in table:
#             all_rows += rows
#     else:
#         all_rows = table
#
#     return all_rows


def val_to_str(val):
    if val is None:
        return '-'
    else:
        return str(val)


def save_txt(filename, table, is_transposed=False):
    txt_file = open(filename, 'w')
    if qconfig.min_contig:
        print >> txt_file, 'All statistics are based on contigs of size >= %d bp, unless otherwise noted ' % qconfig.min_contig + \
                           '(e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).'
        print >> txt_file

    all_rows = []
    for ref_name, groups in table:
        for group_name, group_rows in groups:
            all_rows.extend(group_rows)

    # determine width of columns for nice spaces
    colwidths = [0] * (len(all_rows[0]['values']) + 1)
    for row in all_rows:
        for i, cell in enumerate([row['metricName']] + map(val_to_str, row['values'])):
            colwidths[i] = max(colwidths[i], len(cell))

    for ref_name, groups in table:
        all_ref_rows = []
        for group_name, group_rows in groups:
            all_ref_rows.extend(group_rows)

        if ref_name and len(references) > 1:
            print >> txt_file, ref_name

        for row in all_ref_rows:
            print >> txt_file, '  '.join('%-*s' % (colwidth, cell)
                                         for colwidth, cell
                                         in zip(colwidths,
                                                [row['metricName']] +
                                                map(val_to_str, row['values'])))
        if len(references) > 1 and len(all_ref_rows) > 1:  # TODO: make (... > 0 or ... > 1 and metricName is NAME)
            print >> txt_file

    txt_file.close()


def save_tsv(filename, table, is_transposed=False):
    tsv_file = open(filename, 'w')

    for ref_name, ref_groups in table:
        all_rows = []
        for group_name, group_rows in ref_groups:
            all_rows.extend(group_rows)

        if ref_name and len(references) > 1:
            print >> tsv_file, ref_name

        for row in all_rows:
            print >> tsv_file, '\t'.join([row['metricName']] + map(val_to_str, row['values']))

        if len(references) > 1:
            print >> tsv_file

    tsv_file.close()


def parse_number(val):
    # Float?
    try:
        num = int(val)
    except ValueError:
        # Int?
        try:
            num = float(val)
        except ValueError:
            num = None

    return num


def get_num_from_table_value(val):
    if isinstance(val, int) or isinstance(val, float):
        num = val

    elif isinstance(val, basestring) and len(val.split()) > 0:
                                                                      # 'x + y part' format?
        tokens = val.split()                                          # tokens = [x, +, y, part]
        if len(tokens) >= 3:                                          # Yes, 'y + x part' format
            x, y = parse_number(tokens[0]), parse_number(tokens[2])
            if x is None or y is None:
                num = None
            else:
                num = (x, y)                                          # Tuple value. Can be compared lexicographically.
        else:
            num = parse_number(tokens[0])
    else:
        num = val

    return num


def __prettify_latex_row(row):
    # escape characters
    for esc_char in "\\ % $ # _ { } ~ ^".split():
        row = row.replace(esc_char, '\\' + esc_char)

    # pretty '>=' and '<=', '>'
    row = row.replace('>=', '$\\geq$')
    row = row.replace('<=', '$\\leq$')
    row = row.replace('>', '$>$')

    # pretty indent
    if row.startswith(Fields.TAB):
        row = "\hspace{5mm}" + row.lstrip()

    # pretty highlight
    row = row.replace('HIGHLIGHTEDSTART', '{\\bf ')
    row = row.replace('HIGHLIGHTEDEND', '}')

    return row


def save_tex(filename, table, is_transposed=False):
    tex_file = open(filename, 'w')
    # Header
    print >>tex_file, '\\documentclass[12pt,a4paper]{article}'
    print >>tex_file, '\\begin{document}'
    print >>tex_file, '\\begin{table}[ht]'
    print >>tex_file, '\\begin{left}'
    print >>tex_file, '\\caption{All statistics are based on contigs of size $\geq$ %d bp, unless otherwise noted ' % qconfig.min_contig + \
                      '(e.g., "\# contigs ($\geq$ 0 bp)" and "Total length ($\geq$ 0 bp)" include all contigs).}'
    print >>tex_file, '\\vspace{10pt}'

    all_rows = []
    for ref_name, ref_groups in table:
        for group_name, group_rows in ref_groups:
            all_rows.extend(group_rows)

    assemblies_n = len(all_rows[0]['values'])
    print >>tex_file, '\\begin{tabular}{|l*{' + str(assemblies_n) + '}{|r}|}'

    # Table header
    assembly_names = all_rows[0]['values']
    header_row = ' & '.join([''] + assembly_names)
    print >>tex_file, '\\hline'
    print >>tex_file, __prettify_latex_row(header_row), '\\\\'

    for i, (ref_name, ref_groups) in enumerate(table):
        ref_rows = []

        for group_name, group_rows in ref_groups:
            if i == 0:  # Skipping header (Assembly: contigs_1, contigs_2...) - it's being printed above
                group_rows = group_rows[1:]
            ref_rows.extend(group_rows)

        if not ref_rows:
            continue

        # Reference name
        if ref_name and len(references) > 1:
            row = ref_name
            print >>tex_file, '\\multicolumn{' + str(assemblies_n) + '}{l}{} \\\\'
            print >>tex_file, '\\multicolumn{' + str(assemblies_n) + '}{l}' + \
                              '{' + __prettify_latex_row(row) + '}', '\\\\'

        # Table body
        print >>tex_file, '\\hline'
        for row in ref_rows:
            values = row['values']
            quality = row['quality'] if ('quality' in row) else Fields.Quality.EQUAL

            if is_transposed or quality not in [Fields.Quality.MORE_IS_BETTER, Fields.Quality.LESS_IS_BETTER]:
                cells = map(val_to_str, values)
            else:
                # Checking the first value, assuming the others are the same type and format
                num = get_num_from_table_value(values[0])
                if num is None:  # Not a number
                    cells = map(val_to_str, values)
                else:
                    nums = map(get_num_from_table_value, values)
                    best = None
                    if quality == Fields.Quality.MORE_IS_BETTER:
                        best = max(nums)
                    if quality == Fields.Quality.LESS_IS_BETTER:
                        best = min(nums)

                    if len([num for num in nums if num != best]) == 0:
                        cells = map(val_to_str, values)
                    else:
                        cells = ['HIGHLIGHTEDSTART' + val_to_str(v) + 'HIGHLIGHTEDEND'
                                 if get_num_from_table_value(v) == best
                                 else val_to_str(v)
                                 for v in values]

            row = ' & '.join([row['metricName']] + cells)
            print >>tex_file, __prettify_latex_row(row), '\\\\'

        print >>tex_file, '\\hline'

    # Footer
    print >>tex_file, '\\end{tabular}'
    print >>tex_file, '\\end{left}'
    print >>tex_file, '\\end{table}'
    print >>tex_file, '\\end{document}'
    tex_file.close()


def save(output_dirpath, report_name, transposed_report_name, order):
    # Where total report will be saved
    table = get_table(order)

    log = logging.getLogger('quast')
    log.info('  Creating total report...')
    report_txt_filename = os.path.join(output_dirpath, report_name) + '.txt'
    report_tsv_filename = os.path.join(output_dirpath, report_name) + '.tsv'
    report_tex_filename = os.path.join(output_dirpath, report_name) + '.tex'
    save_txt(report_txt_filename, table)
    save_tsv(report_tsv_filename, table)
    save_tex(report_tex_filename, table)
    log.info('    saved to ' + report_txt_filename + ', ' + os.path.basename(report_tsv_filename) + \
             ', and ' + os.path.basename(report_tex_filename))

    if transposed_report_name:  # TODO: make for multiple references
        log.info('  Transposed version of total report...')

        all_rows = []
        for ref_name, ref_groups in table:
            for group_name, group_rows in ref_groups:
                all_rows.extend(group_rows)

        if all_rows[0]['metricName'] != Fields.NAME:
            warning('transposed version can\'t be created! First column have to be assemblies names')
        else:
            # Transposing table
            # initializing with a header, that contains metric names
            transposed_table = \
                [('', [('', [{'metricName': all_rows[0]['metricName'],
                              'values': [all_rows[i]['metricName'] for i in xrange(1, len(all_rows))]
                             }]
                       )]
                 )]
            # filling up with values
            for i in range(len(all_rows[0]['values'])):
                values = []
                for j in range(1, len(all_rows)):
                    values.append(all_rows[j]['values'][i])
                transposed_table[0][1][0][1].append(
                    {'metricName': all_rows[0]['values'][i],  # name of assembly, assuming the first line is assemblies names
                     'values': values}
                )

            report_txt_filename = os.path.join(output_dirpath, transposed_report_name) + '.txt'
            report_tsv_filename = os.path.join(output_dirpath, transposed_report_name) + '.tsv'
            report_tex_filename = os.path.join(output_dirpath, transposed_report_name) + '.tex'
            save_txt(report_txt_filename, transposed_table, is_transposed=True)
            save_tsv(report_tsv_filename, transposed_table, is_transposed=True)
            save_tex(report_tex_filename, transposed_table, is_transposed=True)
            log.info('    saved to ' + report_txt_filename + ', ' + os.path.basename(report_tsv_filename) + \
                     ', and ' + os.path.basename(report_tex_filename))


def save_gage(output_dirpath):
    save(output_dirpath, qconfig.gage_report_prefix + qconfig.report_prefix,
         qconfig.gage_report_prefix + qconfig.transposed_report_prefix, Fields.gage_order)


def save_total(output_dirpath):
    print_timestamp()
    log = logging.getLogger('quast')
    log.info('Summarizing...')
    save(output_dirpath, qconfig.report_prefix, qconfig.transposed_report_prefix, Fields.order)


def save_misassemblies(output_dirpath):
    save(output_dirpath, "misassemblies_report", "", Fields.misassemblies_order)


def save_unaligned(output_dirpath):
    save(output_dirpath, "unaligned_report", "", Fields.unaligned_order)
