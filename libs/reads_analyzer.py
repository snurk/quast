import shutil
from libs import reporting, qconfig, qutils

from libs.log import get_logger

logger = get_logger(qconfig.LOGGER_DEFAULT_NAME)
import shlex
import os

bwa_dirpath = os.path.join(qconfig.LIBS_LOCATION, 'bwa-master')
samtools_dirpath = os.path.join(qconfig.LIBS_LOCATION, 'samtools-1.2')
bedtools_dirpath = os.path.join(qconfig.LIBS_LOCATION, 'bedtools-2.17.0')
lumpytools_dirpath = os.path.join(qconfig.LIBS_LOCATION, 'lumpy')


def process_single_file(index, ref_fpath, bwa_threads, reads_fpath, output_dirpath, res_path, log_path, err_path):
    assembly_name = qutils.name_from_fpath(ref_fpath)
    res_fpath = os.path.join(res_path, assembly_name + '.res')
    logger.info('  ' + qutils.index_to_str(index) + assembly_name)
    cov_fpath = os.path.join(res_path, assembly_name + '.cov')

    if os.path.isfile(res_fpath):
        logger.info('  ' + qutils.index_to_str(index) + 'Using existing BWA alignments...')
        return res_fpath, cov_fpath
    sam_fpath = os.path.join(output_dirpath, assembly_name + '.sam')
    bam_fpath = os.path.join(output_dirpath, assembly_name + '.bam')
    bamsorted_fpath = os.path.join(output_dirpath, assembly_name + '.sorted.bam')
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
    qutils.call_subprocess([samtools_fpath('samtools'), 'sort', bam_fpath, bamsorted_fpath], stdout=open(cov_fpath, 'w'),
                           stderr=open(err_path, 'a'))
    qutils.call_subprocess([samtools_fpath('samtools'), 'depth', bamsorted_fpath], stdout=open(cov_fpath, 'w'),
                           stderr=open(err_path, 'a'))
    return res_fpath, cov_fpath

def run_lumpy(ref_fpath, output_dirpath, res_path, err_path):

    ref_name = qutils.name_from_fpath(ref_fpath)
    logger.info('  Running LUMPY for %s...' % ref_name)
    bed_fpath = os.path.join(res_path, ref_name + '.bed')
    logger.info('  ' + 'Logging to %s...' % (err_path))
    if os.path.isfile(bed_fpath):
        logger.info('  Using existing bed-file for %s...' % ref_name)
        return bed_fpath
    bam_fpath = os.path.join(output_dirpath, ref_name + '.bam')
    ##preprocessing for lumpy
    vcfoutput_dirpath = os.path.join(output_dirpath, 'lumpy_output')

    bamdiscordants_fpath = os.path.join(vcfoutput_dirpath, ref_name + '.discordants.bam')
    splitreads_fpath = os.path.join(vcfoutput_dirpath, ref_name + '.split')
    bamsplitter_fpath = os.path.join(vcfoutput_dirpath, ref_name + '.splitters.bam')
    qutils.call_subprocess([samtools_fpath('samtools'), 'faidx', ref_fpath], stderr=open(err_path, 'a'))
    discordsorted_fpath= os.path.join(vcfoutput_dirpath, ref_name + '.discordants.sorted')
    splitsorted_fpath = os.path.join(vcfoutput_dirpath, ref_name + '.splitters.sorted')
    readgroup_fpath = os.path.join(vcfoutput_dirpath, ref_name + '.readgroup')
    histo_fpath = os.path.join(vcfoutput_dirpath, ref_name + '.histo')

    temp_fpath = os.path.join(vcfoutput_dirpath, ref_name + '.temp')

    cmd = samtools_fpath('samtools') + ' view -@ %s -b -F 1294 ' % qconfig.max_threads + bam_fpath
    qutils.call_subprocess(shlex.split(cmd), stdout=open(bamdiscordants_fpath, 'w'),
                           stderr=open(err_path, 'a'))
    qutils.call_subprocess([samtools_fpath('samtools'), 'view', '-@', str(qconfig.max_threads), '-h', bam_fpath], stdout=open(splitreads_fpath, 'w'),
                           stderr=open(err_path, 'a'))

    qutils.call_subprocess([os.path.join(lumpytools_dirpath, 'scripts/extractSplitReads_BwaMem'), '-i', splitreads_fpath],
                           stdout=open(temp_fpath, 'w'),
                           stderr=open(err_path, 'a'))
    qutils.call_subprocess([samtools_fpath('samtools'), 'view', '-@', str(qconfig.max_threads), '-Sb', temp_fpath], stdout=open(bamsplitter_fpath, 'w'),
                           stderr=open(err_path, 'a'))
    qutils.call_subprocess([samtools_fpath('samtools'), 'sort', '-@', str(qconfig.max_threads), bamsplitter_fpath, splitsorted_fpath,],
                           stderr=open(err_path, 'a'))
    qutils.call_subprocess([samtools_fpath('samtools'), 'sort', '-@', str(qconfig.max_threads), bamdiscordants_fpath, discordsorted_fpath],
                           stderr=open(err_path, 'a'))
    qutils.call_subprocess([samtools_fpath('samtools'), 'view', '-@', str(qconfig.max_threads), '-r', 'readgroup1', bam_fpath],
                           stdout=open(temp_fpath, 'w'),
                           stderr=open(err_path, 'a'))
    qutils.call_subprocess(['tail', temp_fpath, '-n', '100000'], stdout=open(readgroup_fpath, 'w'),
                           stderr=open(err_path, 'a'))

    bam_info = open(temp_fpath).readline().split()
    if len(bam_info) < 10:
        logger.info('  Failed searching structural variations.')
        return
    read_length = str(len(bam_info[9]))
    qutils.call_subprocess(
        [os.path.join(lumpytools_dirpath, 'scripts/pairend_distro.py'), '-r', read_length, '-X', '4', '-N', '10000',
         '-o', histo_fpath],
        stdin=open(readgroup_fpath), stdout=open(temp_fpath, 'w'), stderr=open(err_path, 'a'))

    bam_statistics = open(temp_fpath).readline()
    if len(bam_statistics) < 4:
        logger.info('  Failed searching structural variations.')
        return
    mean = bam_statistics.split()[1]
    stdev = bam_statistics.split()[3]
    pe_files = 'id:%s,bam_file:%s.bam,histo_file:%s,mean:%s,stdev:%s,read_length:%s,min_non_overlap:%s,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:2' %\
               (ref_name, discordsorted_fpath, histo_fpath, mean, stdev, read_length, read_length)
    sr_files = 'id:%s,bam_file:%s.bam,back_distance:10,weight:1,min_mapping_threshold:20' % (ref_name, splitsorted_fpath)
    qutils.call_subprocess([os.path.join(lumpytools_dirpath, 'bin/lumpy'), '-mw', '4', '-tt', '0', '-pe', pe_files, '-sr', sr_files],
                           stdout=open(temp_fpath, 'w'),
                           stderr=open(err_path, 'a'))
    qutils.call_subprocess([os.path.join(lumpytools_dirpath, 'scripts/vcfToBedpe'), '-i', temp_fpath, '-o', bed_fpath],
        stderr=open(err_path, 'a'))
    logger.info('  Structural variations saved to ' + bed_fpath)
    return bed_fpath

def bin_fpath(fname):
    return os.path.join(bwa_dirpath, fname)


def samtools_fpath(fname):
    return os.path.join(samtools_dirpath, fname)


def bed_fpath(fname):
    return os.path.join(bedtools_dirpath, 'bin', fname)


def lumpy_fpath(fname):
    return os.path.join(lumpytools_dirpath, 'bin', fname)


def all_required_binaries_exist(bin_dirpath, binary):
    if not os.path.isfile(os.path.join(bin_dirpath, binary)):
        return False
    return True


def do(ref_fpath, contigs_fpaths, reads_fpaths, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    logger.print_timestamp()
    logger.info('Running reads aligner...')

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
                             'You can restart Quast with the --debug flag '
                             'to see the command line.' if not qconfig.debug else ''))
            logger.info('Failed aligning the reads.')
            return

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
                             'You can restart Quast with the --debug flag '
                             'to see the command line.' if not qconfig.debug else ''))
            logger.info('Failed aligning the reads.')
            return

    if not all_required_binaries_exist(lumpytools_dirpath, 'bin/lumpy'):
        # making
        logger.info(
            'Compiling LUMPY (details are in ' + os.path.join(lumpytools_dirpath, 'make.log') + ' and make.err)')
        return_code = qutils.call_subprocess(
            ['make', '-C', lumpytools_dirpath],
            stdout=open(os.path.join(lumpytools_dirpath, 'make.log'), 'w'),
            stderr=open(os.path.join(lumpytools_dirpath, 'make.err'), 'w'), )

        if return_code != 0 or not all_required_binaries_exist(lumpytools_dirpath, 'bin/lumpy'):
            logger.error('Failed to compile LUMPY (' + lumpytools_dirpath + ')! '
                                                                               'Try to compile it manually. ' + (
                             'You can restart Quast with the --debug flag '
                             'to see the command line.' if not qconfig.debug else ''))
            logger.info('Failed searching structural variations.')
            return

    temp_output_dir = os.path.join(output_dir, 'temp_output')
    final_output_dir = os.path.join(output_dir, 'output')
    vcfoutput_dirpath = os.path.join(temp_output_dir, 'lumpy_output')

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
    res_fpaths = Parallel(n_jobs=n_jobs)(delayed(process_single_file)(index, fpath, bwa_threads,
                                                                        reads_fpaths, temp_output_dir, final_output_dir, log_path,
                                                                        err_path) for (index, fpath) in enumerate(proc_files))
    cov_fpaths = [res_fpaths[i][1] for i in range(len(res_fpaths))]

    if ref_fpath:
        bed_fpath = run_lumpy(ref_fpath, temp_output_dir, final_output_dir, err_path)
    else:
        bed_fpath = None
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
        results = open(res_fpaths[0][index])
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
    return bed_fpath, cov_fpaths
