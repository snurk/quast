############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import datetime
import os
from libs.qutils import warning


output_dirpath = None


simplejson_error = False
try:
    import json
except:
    try:
        import simplejson as json
    except:
        warning('Can\'t build html report - please install python-simplejson')
        simplejson_error = True

total_report_fn       = '/report.json'
contigs_lengths_fn    = '/contigs_lengths.json'
ref_length_fn         = '/ref_length.json'
aligned_contigs_fn    = '/aligned_contigs_lengths.json'
assemblies_lengths_fn = '/assemblies_lengths.json'
in_contigs_suffix_fn  = '_in_contigs.json'
gc_fn                 = '/gc.json'

prefix_fn             = '/'
suffix_fn             = '.json'


def save(fpath, what):
    if simplejson_error:
        return None

    if os.path.exists(fpath):
        os.remove(fpath)

    json_file = open(fpath, 'w')
    json.dump(what, json_file, separators=(',', ':'))
    json_file.close()
    return fpath


def save_total_report(min_contig, output_dirpath=output_dirpath):
    if output_dirpath:
        from libs import reporting
        assemblies_names = reporting.assemblies_order
        report = reporting.get_table(reporting.Fields.grouped_order)
        t = datetime.datetime.now()

        return save(output_dirpath + total_report_fn, {
            'date': t.strftime('%d %B %Y, %A, %H:%M:%S'),
            'assembliesNames': assemblies_names,
            'report': report,
            'minContig': min_contig,
        })
    else:
        return None

#def save_old_total_report(output_dir, min_contig):
#    from libs import reporting
#    table = reporting.table()
#
#    def try_convert_back_to_number(str):
#        try:
#            val = int(str)
#        except ValueError:
#            try:
#                val = float(str)
#            except ValueError:
#                val = str
#
#        return val
#
#
#    table = [[try_convert_back_to_number(table[i][j]) for i in xrange(len(table))] for j in xrange(len(table[0]))]
#
#    # TODO: check correctness, not sure that header and result are correct:
#    header = table[0]
#    results = table[1:]
#
#    t = datetime.datetime.now()
#
#    return save(output_dir + total_report_fn, {
#        'date' : t.strftime('%d %B %Y, %A, %H:%M:%S'),
#        'header' : header,
#        'results' : results,
#        'min_contig' : min_contig,
#        })

def save_contigs_lengths(filenames, lists_of_lengths, output_dirpath=output_dirpath):
    if output_dirpath:
        lists_of_lengths = [sorted(list, reverse=True) for list in lists_of_lengths]

        return save(output_dirpath + contigs_lengths_fn, {
            'filenames': map(os.path.basename, filenames),
            'lists_of_lengths': lists_of_lengths
        })
    else:
        return None


def save_reference_length(reference_length, output_dirpath=output_dirpath):
    if output_dirpath:
        return save(output_dirpath + ref_length_fn, {'reflen': reference_length })
    else:
        return None


def save_aligned_contigs_lengths(filenames, lists_of_lengths, output_dirpath=output_dirpath):
    if output_dirpath:
        lists_of_lengths = [sorted(list, reverse=True) for list in lists_of_lengths]

        return save(output_dirpath + aligned_contigs_fn, {
            'filenames': map(os.path.basename, filenames),
            'lists_of_lengths': lists_of_lengths
        })
    else:
        return None


def save_assembly_lengths(filenames, assemblies_lengths, output_dirpath=output_dirpath):
    if output_dirpath:
        return save(output_dirpath + assemblies_lengths_fn, {
            'filenames': map(os.path.basename, filenames),
            'assemblies_lengths': assemblies_lengths
        })
    else:
        return None


def save_features_in_contigs(filenames, feature_name, features_in_contigs, ref_features_num, output_dirpath=output_dirpath):
    if output_dirpath:
        return save(output_dirpath + prefix_fn + feature_name + in_contigs_suffix_fn, {
            'filenames': map(os.path.basename, filenames),
            feature_name + '_in_contigs': dict((os.path.basename(fn), feature_amounts) for (fn, feature_amounts) in features_in_contigs.items()),
            'ref_' + feature_name + '_number': ref_features_num,
            })
    else:
        return None


def save_GC_info( filenames, list_of_GC_distributions, output_dirpath=output_dirpath):
    if output_dirpath:
        return save(output_dirpath + gc_fn, {
            'filenames': map(os.path.basename, filenames),
            'list_of_GC_distributions': list_of_GC_distributions,
            'lists_of_gc_info': None,
            })
    else:
        return None


















