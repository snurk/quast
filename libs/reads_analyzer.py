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


def process_single_file(ref_fpath, bwa_threads, reads_fpath, output_dirpath, res_path, log_path, err_path, is_reference=False):
    assembly_name = qutils.name_from_fpath(ref_fpath)
    res_dirpath = os.path.join(res_path, assembly_name + '.res')

    if os.path.isfile(res_dirpath):
        logger.info('Using existing bwa alignments for ' + assembly_name)
        return res_dirpath
    sam_dirpath = os.path.join(output_dirpath, assembly_name + '.sam')
    bam_dirpath = os.path.join(output_dirpath, assembly_name + '.bam')
    bams_dirpath = os.path.join(output_dirpath, assembly_name + '.sorted')
    qutils.call_subprocess([bin_fpath('bwa'), 'index', ref_fpath], stdout=open(log_path, 'a'),
                           stderr=open(err_path, 'a'))
    cmd = bin_fpath('bwa') + ' mem -t' + str(bwa_threads) + ' ' + ref_fpath + ' ' + ' '.join(reads_fpath)
    qutils.call_subprocess(shlex.split(cmd), stdout=open(sam_dirpath, 'w'), stderr=open(err_path, 'a'))
    qutils.call_subprocess([sam_fpath('samtools'), 'faidx', ref_fpath], stderr=open(err_path, 'a'))
    qutils.call_subprocess([sam_fpath('samtools'), 'view', '-bS', sam_dirpath], stdout=open(bam_dirpath, 'w'),
                           stderr=open(err_path, 'a'))

    qutils.assert_file_exists(bam_dirpath, 'bam file')
    qutils.call_subprocess([sam_fpath('samtools'), 'sort', bam_dirpath, bams_dirpath],
                           stderr=open(err_path, 'a'))

    qutils.call_subprocess([sam_fpath('samtools'), 'flagstat', bam_dirpath], stdout=open(res_dirpath, 'w'),
                           stderr=open(err_path, 'a'))


    return res_dirpath

def create_vcf(ref_fpath, output_dirpath, res_path, log_path, err_path):

    ref_name = qutils.name_from_fpath(ref_fpath)
    bed_dirpath = os.path.join(res_path, ref_name + '.bed')
    if os.path.isfile(bed_dirpath):
        logger.info('Using existing bed-file for ' + ref_name)
        return bed_dirpath
    bam_dirpath = os.path.join(output_dirpath, ref_name + '.bam')
    ##preprocessing for lumpy
    vcfoutput_dirpath = os.path.join(output_dirpath, 'lumpy_output')

    bamdiscordants_dirpath = os.path.join(vcfoutput_dirpath, ref_name + '.discordants.bam')
    splitreads_dirpath = os.path.join(vcfoutput_dirpath, ref_name + '.split')
    bamsplitter_dirpath = os.path.join(vcfoutput_dirpath, ref_name + '.splitters.bam')
    qutils.call_subprocess([sam_fpath('samtools'), 'faidx', ref_fpath], stderr=open(err_path, 'a'))
    discordsorted_dirpath = os.path.join(vcfoutput_dirpath, ref_name + '.discordants.sorted')
    splitsorted_dirpath = os.path.join(vcfoutput_dirpath, ref_name + '.splitters.sorted')
    readgroup_dirpath = os.path.join(vcfoutput_dirpath, ref_name + '.readgroup')
    histo_dirpath = os.path.join(vcfoutput_dirpath, ref_name + '.histo')

    temp_file = os.path.join(vcfoutput_dirpath, ref_name + '.temp')

    cmd = sam_fpath('samtools') + ' view -b -F 1294 ' + bam_dirpath
    qutils.call_subprocess(shlex.split(cmd), stdout=open(bamdiscordants_dirpath, 'w'),
                           stderr=open(err_path, 'a'))
    qutils.call_subprocess([sam_fpath('samtools'), 'view', '-h', bam_dirpath], stdout=open(splitreads_dirpath, 'w'),
                           stderr=open(err_path, 'a'))

    qutils.call_subprocess([os.path.join(lumpytools_dirpath, 'scripts/extractSplitReads_BwaMem'), '-i', splitreads_dirpath],
                           stdout=open(temp_file, 'w'),
                           stderr=open(err_path, 'a'))
    qutils.call_subprocess([sam_fpath('samtools'), 'view', '-Sb', temp_file], stdout=open(bamsplitter_dirpath, 'w'),
                           stderr=open(err_path, 'a'))
    qutils.call_subprocess([sam_fpath('samtools'), 'sort', bamsplitter_dirpath, splitsorted_dirpath],
                           stderr=open(err_path, 'a'))
    qutils.call_subprocess([sam_fpath('samtools'), 'sort', bamdiscordants_dirpath, discordsorted_dirpath],
                           stderr=open(err_path, 'a'))
    qutils.call_subprocess([sam_fpath('samtools'), 'view', '-r', 'readgroup1', bam_dirpath],
                           stdout=open(temp_file, 'w'),
                           stderr=open(err_path, 'a'))
    qutils.call_subprocess(['tail', temp_file, '-n', '+100000'], stdout=open(readgroup_dirpath, 'w'),
                           stderr=open(err_path, 'a'))

    bam_info = open(temp_file).readline()
    if len(bam_info) < 4:
        return
    read_length = bam_info.split()[4]
    qutils.call_subprocess(
        [os.path.join(lumpytools_dirpath, 'scripts/pairend_distro.py'), '-r', read_length, '-X', '4', '-N', '10000',
         '-o', histo_dirpath],
        stdin=open(readgroup_dirpath), stdout=open(temp_file, 'w'), stderr=open(err_path, 'a'))
    bam_statistics = open(temp_file).readline()
    if len(bam_statistics) < 4:
        return
    mean = bam_statistics.split()[1]
    stdev = bam_statistics.split()[3]
    pe_files = 'id:%s,bam_file:%s.bam,histo_file:%s,mean:%s,stdev:%s,read_length:%s,min_non_overlap:%s,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:2' %\
               (ref_name, discordsorted_dirpath, histo_dirpath, mean, stdev, read_length, read_length)
    sr_files = 'id:%s,bam_file:%s.bam,back_distance:10,weight:1,min_mapping_threshold:20' % (ref_name, splitsorted_dirpath)
    qutils.call_subprocess([os.path.join(lumpytools_dirpath, 'bin/lumpy'), '-mw', '4', '-tt', '0', '-pe', pe_files, '-sr', sr_files],
                           stdout=open(temp_file, 'w'),
                           stderr=open(err_path, 'a'))
    qutils.call_subprocess(
        [os.path.join(lumpytools_dirpath, 'scripts/vcfToBedpe'), '-i', temp_file, '-o', bed_dirpath],
        stderr=open(err_path, 'a'))

    return bed_dirpath

def bin_fpath(fname):
    return os.path.join(bwa_dirpath, fname)


def sam_fpath(fname):
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
        logger.info('Compiling bwa (details are in ' + os.path.join(bwa_dirpath, 'make.log') + ' and make.err)')
        return_code = qutils.call_subprocess(
            ['make', '-C', bwa_dirpath],
            stdout=open(os.path.join(bwa_dirpath, 'make.log'), 'w'),
            stderr=open(os.path.join(bwa_dirpath, 'make.err'), 'w'), )

        if return_code != 0 or not all_required_binaries_exist(bwa_dirpath, 'bwa'):
            logger.error('Failed to compile bwa (' + bwa_dirpath + ')! '
                                                                   'Try to compile it manually. ' + (
                             'You can restart Quast with the --debug flag '
                             'to see the command line.' if not qconfig.debug else ''))
            logger.info('Failed aligning the reads.')
            return

    if not all_required_binaries_exist(samtools_dirpath, 'samtools'):
        # making
        logger.info(
            'Compiling samtools (details are in ' + os.path.join(samtools_dirpath, 'make.log') + ' and make.err)')
        return_code = qutils.call_subprocess(
            ['make', '-C', samtools_dirpath],
            stdout=open(os.path.join(samtools_dirpath, 'make.log'), 'w'),
            stderr=open(os.path.join(samtools_dirpath, 'make.err'), 'w'), )

        if return_code != 0 or not all_required_binaries_exist(samtools_dirpath, 'samtools'):
            logger.error('Failed to compile samtools (' + samtools_dirpath + ')! '
                                                                             'Try to compile it manually. ' + (
                             'You can restart Quast with the --debug flag '
                             'to see the command line.' if not qconfig.debug else ''))
            logger.info('Failed aligning the reads.')
            return

    if not all_required_binaries_exist(lumpytools_dirpath, 'bin/lumpy'):
        # making
        logger.info(
            'Compiling lumpytools (details are in ' + os.path.join(lumpytools_dirpath, 'make.log') + ' and make.err)')
        return_code = qutils.call_subprocess(
            ['make', '-C', lumpytools_dirpath],
            stdout=open(os.path.join(lumpytools_dirpath, 'make.log'), 'w'),
            stderr=open(os.path.join(lumpytools_dirpath, 'make.err'), 'w'), )

        if return_code != 0 or not all_required_binaries_exist(lumpytools_dirpath, 'bin/lumpy'):
            logger.error('Failed to compile lumpytools (' + lumpytools_dirpath + ')! '
                                                                               'Try to compile it manually. ' + (
                             'You can restart Quast with the --debug flag '
                             'to see the command line.' if not qconfig.debug else ''))
            logger.info('Failed aligning the reads.')
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
    log_path = os.path.join(temp_output_dir, 'bwa.log')
    err_path = os.path.join(temp_output_dir, 'bwa.err')

    from joblib import Parallel, delayed
    res_dirpaths = Parallel(n_jobs=n_jobs)(delayed(process_single_file)(fpath, bwa_threads,
                                                                        reads_fpaths, temp_output_dir, final_output_dir, log_path,
                                                                        err_path) for fpath in proc_files)
    create_vcf(ref_fpath, temp_output_dir, final_output_dir, log_path, err_path)

    if ref_fpath:
        assembly_name = qutils.name_from_fpath(ref_fpath)
        ref_dirpath = os.path.join(final_output_dir, assembly_name + '.res')
        ref_results = open(ref_dirpath)
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
        results = open(res_dirpaths[index])
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

    #if not qconfig.debug:
        #shutil.rmtree(temp_output_dir, ignore_errors=True)

    reporting.save_reads(output_dir)
    logger.info('Done.')
