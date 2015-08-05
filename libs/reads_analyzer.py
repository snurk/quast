import pipes
import shutil
import stat
import subprocess
from libs import reporting, qconfig, qutils

from libs.log import get_logger

logger = get_logger(qconfig.LOGGER_DEFAULT_NAME)
import shlex
import os

bwa_dirpath = os.path.join(qconfig.LIBS_LOCATION, 'bwa-master')
samtools_dirpath = os.path.join(qconfig.LIBS_LOCATION, 'samtools-1.2')
bedtools_dirpath = os.path.join(qconfig.LIBS_LOCATION, 'bedtools-2.17.0')
manta_dirpath = os.path.join(qconfig.LIBS_LOCATION, 'manta')

def process_single_file(index, ref_fpath, bwa_threads, reads_fpath, output_dirpath, res_path, log_path, err_path):
    assembly_name = qutils.name_from_fpath(ref_fpath)
    res_fpath = os.path.join(res_path, assembly_name + '.res')
    logger.info('  ' + qutils.index_to_str(index) + assembly_name)

    if os.path.isfile(res_fpath) and os.path.getsize(res_fpath) > 0:
        logger.info('  ' + qutils.index_to_str(index) + 'Using existing BWA alignments...')
        return res_fpath
    sam_fpath = os.path.join(output_dirpath, assembly_name + '.sam')
    bam_fpath = os.path.join(output_dirpath, assembly_name + '.bam')
    logger.info('  ' + qutils.index_to_str(index) + 'Running BWA...')
    qutils.call_subprocess([bin_fpath('bwa'), 'index', ref_fpath], stdout=open(log_path, 'a'),
                           stderr=open(err_path, 'a'))
    cmd = bin_fpath('bwa') + ' mem -t' + str(bwa_threads) + ' ' + ref_fpath + ' ' + ' '.join(reads_fpath)
    qutils.call_subprocess(shlex.split(cmd), stdout=open(sam_fpath, 'w'), stderr=open(err_path, 'a'))
    qutils.call_subprocess([samtools_fpath('samtools'), 'faidx', ref_fpath], stderr=open(err_path, 'a'))
    qutils.call_subprocess([samtools_fpath('samtools'), 'view', '-@', str(bwa_threads), '-bS', sam_fpath], stdout=open(bam_fpath, 'w'),
                           stderr=open(err_path, 'a'))
    qutils.assert_file_exists(bam_fpath, 'bam file')
    qutils.call_subprocess([samtools_fpath('samtools'), 'flagstat', bam_fpath], stdout=open(res_fpath, 'w'),
                           stderr=open(err_path, 'a'))
    logger.info('  ' + qutils.index_to_str(index) + 'Analysis is finished.')
    return res_fpath


def get_coverage(output_dirpath, ref_name, bam_fpath, err_path, cov_fpath):
    if not os.path.isfile(cov_fpath):
        bamsorted_fpath = os.path.join(output_dirpath, ref_name + '.sorted')
        if not os.path.isfile(bamsorted_fpath):
            qutils.call_subprocess([samtools_fpath('samtools'), 'sort',  '-@', str(qconfig.max_threads), bam_fpath, bamsorted_fpath], stdout=open(err_path, 'w'),
                               stderr=open(err_path, 'a'))
        qutils.call_subprocess([samtools_fpath('samtools'), 'depth', bamsorted_fpath + '.bam'], stdout=open(cov_fpath, 'w'),
                               stderr=open(err_path, 'a'))
        qutils.assert_file_exists(cov_fpath, 'coverage file')
    return cov_fpath


def get_insert_size(ref_fpath, output_dirpath, err_path, res_fpath):
    insert_size = 0
    if os.path.isfile(res_fpath) and os.path.getsize(res_fpath) > 0:
        first_line = open(res_fpath).readline().rstrip()
        if 'insert' in first_line:
            insert_size = int(first_line.split(':')[1])
        return insert_size
    ref_name = qutils.name_from_fpath(ref_fpath)
    bam_fpath = os.path.join(output_dirpath, ref_name + '.bam')
    ls = subprocess.Popen([samtools_fpath('samtools'), 'view', '-@', str(qconfig.max_threads), '-f66', bam_fpath],
                          stdout=subprocess.PIPE,
                          stderr=open(err_path, 'a'))
    ls = subprocess.Popen(['cut', '-f', '9'], stdin=ls.stdout, stdout=subprocess.PIPE,
                          stderr=open(err_path, 'a'))
    ls = subprocess.Popen(['head', '-n', '10000'], stdin=ls.stdout, stdout=subprocess.PIPE,
                          stderr=open(err_path, 'a'))
    ls = subprocess.Popen(['sed', "s/^-//"], stdin=ls.stdout, stdout=subprocess.PIPE,
                          stderr=open(err_path, 'a'))
    output = ls.communicate()[0]
    if output:
        insert_sizes = output.split('\n')
        insert_sizes = [int(i) for i in insert_sizes if i != '']
        if len(insert_sizes) > 0:
            insert_size = int(sum(insert_sizes) / len(insert_sizes))
    with open(res_fpath, 'w') as res_file:
        res_file.writelines('Mean insert size:' + str(insert_size) + '\n')
    return insert_size


def run_manta(ref_fpath, output_dirpath, res_path, err_path):
    ref_name = qutils.name_from_fpath(ref_fpath)

    bam_fpath = os.path.join(output_dirpath, ref_name + '.bam')
    bamsorted_fpath = os.path.join(output_dirpath, ref_name + '.sorted')
    bed_fpath = os.path.join(res_path, ref_name + '.bed')
    cov_fpath = os.path.join(res_path, ref_name + '.cov')
    vcfoutput_dirpath = os.path.join(output_dirpath, 'manta_output')
    logger.info('  ' + 'Pre-processing for searching structural variations...' % err_path)
    logger.info('  ' + 'Logging to %s...' % err_path)
    if os.path.isfile(bed_fpath):
        logger.info('  Using existing BED-file for %s...' % ref_name)
        if not os.path.isfile(cov_fpath):
            cov_fpath = get_coverage(output_dirpath, ref_name, bam_fpath, err_path, cov_fpath)
        return bed_fpath, cov_fpath

    qutils.call_subprocess([samtools_fpath('samtools'), 'faidx', ref_fpath], stderr=open(err_path, 'a'))
    qutils.call_subprocess([samtools_fpath('samtools'), 'sort', '-@', str(qconfig.max_threads), bam_fpath, bamsorted_fpath],
                           stderr=open(err_path, 'a'))
    qutils.call_subprocess([samtools_fpath('samtools'), 'index', bamsorted_fpath + '.bam'], stderr=open(err_path, 'a'))
    cov_fpath = get_coverage(output_dirpath, ref_name, bam_fpath, err_path, cov_fpath)

    logger.info('  Running Manta for %s...' % ref_name)
    qutils.call_subprocess([os.path.join(manta_dirpath, 'bin/configManta.py'), '--normalBam', bamsorted_fpath + '.bam',
                            '--referenceFasta', ref_fpath, '--runDir', vcfoutput_dirpath],
                           stderr=open(err_path, 'a'))
    qutils.call_subprocess([os.path.join(vcfoutput_dirpath, 'runWorkflow.py'), '-m', 'local',
                            '-j', '1', '-g', str(qconfig.max_threads)],
                           stderr=open(err_path, 'a'))
    temp_fpath = os.path.join(vcfoutput_dirpath, 'results/variants/diploidSV.vcf.gz')
    unpacked_fpath = temp_fpath + '.unpacked'
    cmd = 'gunzip -c %s' % temp_fpath
    qutils.call_subprocess(shlex.split(cmd), stdout=open(unpacked_fpath, 'w'), stderr=open(err_path, 'a'), logger=logger)
    from manta import vcfToBedpe
    vcfToBedpe.vcfToBedpe(open(unpacked_fpath), open(bed_fpath, 'w'))
    if os.path.exists(bed_fpath):
        logger.info('  Structural variations saved to ' + bed_fpath)
        return bed_fpath, cov_fpath
    else:
        logger.info('  Failed searching structural variations.')
        return None, cov_fpath

def bin_fpath(fname):
    return os.path.join(bwa_dirpath, fname)


def samtools_fpath(fname):
    return os.path.join(samtools_dirpath, fname)


def bed_fpath(fname):
    return os.path.join(bedtools_dirpath, 'bin', fname)


def all_required_binaries_exist(bin_dirpath, binary):
    if not os.path.isfile(os.path.join(bin_dirpath, binary)):
        return False
    return True


def do(ref_fpath, contigs_fpaths, reads_fpaths, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    from libs import search_references_meta
    if search_references_meta.is_quast_first_run:
        output_dir = os.path.join(output_dir, 'aux')
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)

    logger.print_timestamp()
    logger.info('Running Reads analyzer...')

    if not all_required_binaries_exist(bwa_dirpath, 'bwa'):
        # making
        logger.info('Compiling BWA (details are in ' + os.path.join(bwa_dirpath, 'make.log') + ' and make.err)')
        return_code = qutils.call_subprocess(
            ['make', '-C', bwa_dirpath],
            stdout=open(os.path.join(bwa_dirpath, 'make.log'), 'w'),
            stderr=open(os.path.join(bwa_dirpath, 'make.err'), 'w'), )

        if return_code != 0 or not all_required_binaries_exist(bwa_dirpath, 'bwa'):
            logger.error('Failed to compile BWA (' + bwa_dirpath + ')! '
                                                                   'Try to compile it manually. ' + (
                             'You can restart QUAST with the --debug flag '
                             'to see the command line.' if not qconfig.debug else ''))
            logger.info('Failed aligning the reads.')
            return None, None

    if not all_required_binaries_exist(samtools_dirpath, 'samtools'):
        # making
        logger.info(
            'Compiling SAMtools (details are in ' + os.path.join(samtools_dirpath, 'make.log') + ' and make.err)')
        return_code = qutils.call_subprocess(
            ['make', '-C', samtools_dirpath],
            stdout=open(os.path.join(samtools_dirpath, 'make.log'), 'w'),
            stderr=open(os.path.join(samtools_dirpath, 'make.err'), 'w'), )

        if return_code != 0 or not all_required_binaries_exist(samtools_dirpath, 'samtools'):
            logger.error('Failed to compile SAMtools (' + samtools_dirpath + ')! '
                                                                             'Try to compile it manually. ' + (
                             'You can restart QUAST with the --debug flag '
                             'to see the command line.' if not qconfig.debug else ''))
            logger.info('Failed aligning the reads.')
            return None, None

    temp_output_dir = os.path.join(output_dir, 'temp_output')
    final_output_dir = os.path.join(output_dir, 'output')
    vcfoutput_dirpath = os.path.join(temp_output_dir, 'manta_output')

    if not os.path.isdir(temp_output_dir):
        os.mkdir(temp_output_dir)
    if not os.path.isdir(final_output_dir):
        os.mkdir(final_output_dir)
    if not os.path.isdir(vcfoutput_dirpath):
        os.mkdir(vcfoutput_dirpath)

    proc_files = contigs_fpaths[:]
    if ref_fpath:
        chromosomes = 0
        ref_file = open(ref_fpath)
        for line in ref_file:
            if line[0] == '>':
                chromosomes += 1
        ref_file.close()
        proc_files.append(ref_fpath)
    n_jobs = min(qconfig.max_threads, len(proc_files))
    bwa_threads = qconfig.max_threads // n_jobs
    log_path = os.path.join(temp_output_dir, 'align_reads.log')
    err_path = os.path.join(temp_output_dir, 'align_reads.err')
    logger.info('  ' + 'Logging to files %s and %s...' % (log_path, err_path))
    from joblib import Parallel, delayed

    res_fpaths = Parallel(n_jobs=n_jobs)(delayed(process_single_file)(index, fpath, bwa_threads, reads_fpaths,
                        temp_output_dir, final_output_dir, log_path, err_path) for (index, fpath) in enumerate(proc_files))
    insert_size_fpath = os.path.join(final_output_dir, 'ins_size.txt')
    insert_size = get_insert_size(ref_fpath, temp_output_dir, err_path, insert_size_fpath)
    if insert_size > 1000:
        qconfig.extensive_misassembly_threshold = insert_size
    if ref_fpath:
        bed_fpath, cov_fpath = run_manta(ref_fpath, temp_output_dir, final_output_dir, err_path)
    else:
        bed_fpath, cov_fpath = None, None
    if ref_fpath:
        assembly_name = qutils.name_from_fpath(ref_fpath)
        ref_respath = os.path.join(final_output_dir, assembly_name + '.res')
        ref_results = open(ref_respath)
        for line in ref_results:
            if 'properly paired' in line:
                ref_paired_reads = line.split()[0]
            elif 'mapped' in line and '%' in line:
                ref_mapped_reads = line.split()[0]
            elif 'singletons' in line:
                ref_singletons = line.split()[0]
            elif 'different' in line and 'mapQ' not in line and chromosomes > 1:
                ref_diffchrom = line.split()[0]

    # process all contigs files
    for index, contigs_fpath in enumerate(contigs_fpaths):
        report = reporting.get(contigs_fpath)
        results = open(res_fpaths[index])
        assembly_label = qutils.label_from_fpath(contigs_fpath)
        for line in results:
            if 'total' in line:
                report.add_field(reporting.Fields.TOTALREADS, line.split()[0])
                report.add_field(reporting.Fields.REF_READS, line.split()[0])
            elif 'read1' in line:
                report.add_field(reporting.Fields.LEFT_READS, line.split()[0])
            elif 'read2' in line:
                report.add_field(reporting.Fields.RIGHT_READS, line.split()[0])
            elif 'properly paired' in line:
                report.add_field(reporting.Fields.PROPERLYPAIR_READS, line.split()[0])
            elif 'mapped' in line and '%' in line:
                report.add_field(reporting.Fields.MAPPED_READS, line.split()[0])
                if line.split()[0] == '0':
                    logger.info('  ' + qutils.index_to_str(index) + 'BWA: nothing aligned for ' + '\'' + assembly_label + '\'.')
            elif 'singletons' in line:
                report.add_field(reporting.Fields.SINGLETONS, line.split()[0])
            elif 'different' in line and 'mapQ' not in line and len(contigs_fpaths) > 1:
                report.add_field(reporting.Fields.READS_DIFFCHROM, line.split()[0])

        if ref_fpath:
            report.add_field(reporting.Fields.REFPROPERLYPAIR_READS, ref_paired_reads)
            if chromosomes > 1:
                report.add_field(reporting.Fields.REFREADS_DIFFCHROM, ref_diffchrom)
            report.add_field(reporting.Fields.REFSINGLETONS, ref_singletons)
            report.add_field(reporting.Fields.REFMAPPED_READS, ref_mapped_reads)
            report.add_field(reporting.Fields.REFCHROMOSOMES, chromosomes)

    if not qconfig.debug:
        shutil.rmtree(temp_output_dir, ignore_errors=True)

    reporting.save_reads(output_dir)
    logger.info('Done.')
    return bed_fpath, cov_fpath
