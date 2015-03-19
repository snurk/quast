from libs import reporting, qconfig, qutils

from libs.log import get_logger
import subprocess

logger = get_logger(qconfig.LOGGER_DEFAULT_NAME)
import shlex
import os

bwa_dirpath = os.path.join(qconfig.LIBS_LOCATION, 'bwa-master')
samtools_dirpath = os.path.join(qconfig.LIBS_LOCATION, 'samtools-1.2')


def process_single_file(ref_fpath, reads_fpath, output_dirpath):
    cmd = sam_fpath('samtools') + ' faidx ' + ref_fpath
    subprocess.call(shlex.split(cmd), stdout=open(os.path.join(output_dirpath, 'bwa.log'), 'a+'),
                    stderr=open(os.path.join(output_dirpath, 'bwa.err'), 'a+'))
    bam_dirpath = output_dirpath + '/out.bam'
    sam_dirpath = output_dirpath + '/out.sam'

    sam = open(sam_dirpath, 'w+')
    cmd = bin_fpath('bwa') + ' mem -t' + str(qconfig.max_threads) + ' ' + ref_fpath + ' ' + ' '.join(reads_fpath)
    subprocess.call(shlex.split(cmd), stdout=sam, stderr=open(os.path.join(output_dirpath, 'bwa.err'), 'a+'))
    sam.close()

    cmd = sam_fpath('samtools') + ' import ' + ref_fpath + ' ' + sam_dirpath + ' ' + bam_dirpath
    subprocess.call(shlex.split(cmd), stdout=open(os.path.join(output_dirpath, 'bwa.log'), 'a+'),
                    stderr=open(os.path.join(output_dirpath, 'bwa.err'), 'a+'))
    cmd = sam_fpath('samtools') + ' flagstat ' + ' ' + bam_dirpath
    results = subprocess.check_output(shlex.split(cmd))

    return results


def bin_fpath(fname):
    return os.path.join(bwa_dirpath, fname)


def sam_fpath(fname):
    return os.path.join(samtools_dirpath, fname)


def all_required_binaries_exist(bin_dirpath, binary):
    if not os.path.isfile(os.path.join(bin_dirpath, binary)):
        return False
    return True


def do(ref_fpath, reads_fpaths, output_dir):
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

    bwa_output_dir = os.path.join(output_dir, 'bwa_output')

    if not os.path.isdir(bwa_output_dir):
        os.mkdir(bwa_output_dir)

    cmd = bin_fpath('bwa') + ' index ' + ref_fpath
    subprocess.call(shlex.split(cmd), stdout=open(os.path.join(bwa_output_dir, 'bwa.log'), 'a+'),
                    stderr=open(os.path.join(bwa_output_dir, 'bwa.err'), 'a+'))

    # process all reads files
    res = process_single_file(ref_fpath, reads_fpaths, bwa_output_dir).split('\n')
    report = reporting.get_reads()

    for line in res:
        if 'total' in line:
            report.add_field(reporting.Fields.TOTALREADS, line.split()[0])
        elif 'properly paired' in line:
            report.add_field(reporting.Fields.PROPERLYPAIR_READS, line.split()[0])
        elif 'read1' in line:
            report.add_field(reporting.Fields.LEFT_READS, line.split()[0])
        elif 'read2' in line:
            report.add_field(reporting.Fields.RIGHT_READS, line.split()[0])
        elif 'mapped' in line and '%' in line:
            report.add_field(reporting.Fields.MAPPED_READS, line.split()[0])

    reporting.save_reads(output_dir)
    logger.info('Done.')
