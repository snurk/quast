from libs import reporting, qconfig, qutils

from libs.log import get_logger
import subprocess
logger = get_logger(qconfig.LOGGER_DEFAULT_NAME)
import shlex

def process_single_file(ref_fpath, reads_fpath,output_dirpath):

    cmd='samtools faidx ' + ref_fpath
    subprocess.call(shlex.split(cmd))
    bam_dirpath = output_dirpath + "/out.bam"
    sam_dirpath = output_dirpath + "/out.sam"

    sam=open(sam_dirpath, 'w+')
    cmd="./bwa mem -t" + str(qconfig.max_threads) + " " + ref_fpath + " " + " ".join(reads_fpath)
    subprocess.call(shlex.split(cmd),stdout = sam)
    sam.close()

    cmd='samtools import '+ ref_fpath + " " + sam_dirpath + " " + bam_dirpath
    subprocess.call(shlex.split(cmd))
    cmd='samtools flagstat ' + " " + bam_dirpath
    results = subprocess.check_output(shlex.split(cmd))

    return results

def do(ref_fpath, reads_fpaths, output_dirpath):

    logger.print_timestamp()
    logger.info('Running reads aligner...')

    cmd="./bwa index " + ref_fpath
    subprocess.call(shlex.split(cmd))

    # process all contig files
    res=process_single_file(ref_fpath, reads_fpaths, output_dirpath).split('\n')

    report = reporting.get(reads_fpaths[0])

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



    logger.info('Done.')
