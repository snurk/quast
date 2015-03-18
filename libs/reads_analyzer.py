from libs import reporting, qconfig, qutils

from libs.log import get_logger
import subprocess
logger = get_logger(qconfig.LOGGER_DEFAULT_NAME)
import shlex
import os

bwa_dirpath = os.path.join(qconfig.LIBS_LOCATION, 'bwa-0.7.12')
samtools_dirpath = os.path.join(qconfig.LIBS_LOCATION, 'samtools-1.2')
def process_single_file(ref_fpath, reads_fpath,output_dirpath):

    cmd=sam_fpath('samtools') + ' faidx ' + ref_fpath
    subprocess.call(shlex.split(cmd))
    bam_dirpath = output_dirpath + "/out.bam"
    sam_dirpath = output_dirpath + "/out.sam"

    sam=open(sam_dirpath, 'w+')
    cmd=bin_fpath('bwa') + " mem -t" + str(qconfig.max_threads) + " " + ref_fpath + " " + " ".join(reads_fpath)
    subprocess.call(shlex.split(cmd),stdout = sam)
    sam.close()

    cmd=sam_fpath('samtools') + ' import '+ ref_fpath + " " + sam_dirpath + " " + bam_dirpath
    subprocess.call(shlex.split(cmd))
    cmd=sam_fpath('samtools') + ' flagstat ' + " " + bam_dirpath
    results = subprocess.check_output(shlex.split(cmd))

    return results

def bin_fpath(fname):
    return os.path.join(bwa_dirpath, fname)

def sam_fpath(fname):
    return os.path.join(samtools_dirpath, fname)

def do(ref_fpath, reads_fpaths, output_dirpath):

    logger.print_timestamp()
    logger.info('Running reads aligner...')

    cmd=bin_fpath('bwa') + " index " + ref_fpath
    subprocess.call(shlex.split(cmd))

    # process all reads files
    res=process_single_file(ref_fpath, reads_fpaths, output_dirpath).split('\n')
    report = reporting.getReads()

    for line in res:
        if 'total' in line:
            report.add_field(reporting.Fields.TOTALREADS,line[0])
        elif 'properly paired' in line:
            report.add_field(reporting.Fields.PROPERLYPAIR_READS,line[0])
        elif 'read1' in line:
            report.add_field(reporting.Fields.LEFT_READS,line[0])
        elif 'read2' in line:
            report.add_field(reporting.Fields.RIGHT_READS,line[0])
        elif 'mapped' in line and '%' in line:
            report.add_field(reporting.Fields.MAPPED_READS,line[0])


    reporting.save_reads(output_dirpath)
    logger.info('Done.')
