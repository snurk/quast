import os
import platform
from Busco import busco_v1
from libs import reporting, qconfig, qutils
from libs.log import get_logger
logger = get_logger(qconfig.LOGGER_DEFAULT_NAME)

busco_dirpath = os.path.join(qconfig.LIBS_LOCATION, 'Busco')
blast_dirpath = os.path.join(busco_dirpath, 'ncbi-blast-2.2.28+')
augustus_dirpath = os.path.join(busco_dirpath, 'augustus-3.0.3')

if platform.system() == 'Darwin':
    hmmer_dirpath = os.path.join(busco_dirpath, 'hmmer-3.1b2-mac')
else:
    hmmer_dirpath = os.path.join(busco_dirpath, 'hmmer-3.1b2')

def all_required_binaries_exist(bin_dirpath, binary):
    if not os.path.isfile(os.path.join(bin_dirpath, binary)):
        return False
    return True

def do(contigs_fpaths, out_dirpath):
    logger.print_timestamp()
    logger.info('Running BUSCO...')

    if not all_required_binaries_exist(augustus_dirpath, 'bin/augustus'):
        # making
        logger.info('Compiling augustus (details are in ' + os.path.join(augustus_dirpath, 'make.log') + ' and make.err)')
        return_code = qutils.call_subprocess(
            ['make', augustus_dirpath],
            stdout=open(os.path.join(augustus_dirpath, 'make.log'), 'w'),
            stderr=open(os.path.join(augustus_dirpath, 'make.err'), 'w'), )

        if return_code != 0 or not all_required_binaries_exist(augustus_dirpath, 'bin/augustus'):
            logger.error('Failed to compile augustus (' + augustus_dirpath + ')! '
                                                                   'Try to compile it manually. ' + (
                             'You can restart Quast with the --debug flag '
                             'to see the command line.' if not qconfig.debug else ''))
            logger.info('Failed aligning the reads.')
        scripts_path = os.path.join(augustus_dirpath, 'scripts')
        return_code = qutils.call_subprocess(
            ['make', scripts_path],
            stdout=open(os.path.join(augustus_dirpath, 'make.log'), 'w'),
            stderr=open(os.path.join(augustus_dirpath, 'make.err'), 'w'), )

        if return_code != 0 or not all_required_binaries_exist(augustus_dirpath, 'bin/augustus'):
            logger.error('Failed to compile genewise (' + augustus_dirpath + ')! '
                                                                   'Try to compile it manually. ' + (
                             'You can restart Quast with the --debug flag '
                             'to see the command line.' if not qconfig.debug else ''))
            logger.info('Failed aligning the reads.')
    if not all_required_binaries_exist(hmmer_dirpath, 'src/hmmsearch'):
        os.chdir(hmmer_dirpath)
        # making
        logger.info('Compiling hmmer (details are in ' + os.path.join(hmmer_dirpath, 'make.log') + ' and make.err)')
        qutils.call_subprocess(['./configure'],stdout=open(os.path.join(hmmer_dirpath, 'make.log'), 'w'),
            stderr=open(os.path.join(hmmer_dirpath, 'make.err'), 'w'), )

        return_code = qutils.call_subprocess(
            ['make', '-C', hmmer_dirpath],
            stdout=open(os.path.join(hmmer_dirpath, 'make.log'), 'w'),
            stderr=open(os.path.join(hmmer_dirpath, 'make.err'), 'w'), )
        if return_code != 0 or not all_required_binaries_exist(hmmer_dirpath, 'src/hmmsearch'):
            logger.error('Failed to compile hmmer (' + hmmer_dirpath + ')! '
                                                                   'Try to compile it manually. ' + (
                             'You can restart Quast with the --debug flag '
                             'to see the command line.' if not qconfig.debug else ''))
            logger.info('Failed aligning the reads.')
            return

    if not os.path.isdir(out_dirpath):
        os.makedirs(out_dirpath)

    n_jobs = min(len(contigs_fpaths), qconfig.max_threads)
    busco_threads = qconfig.max_threads//n_jobs

    from joblib import Parallel, delayed
    results = Parallel(n_jobs=n_jobs)(delayed(busco_v1.do)(['-in', contigs_fpath,'-o', qutils.name_from_fpath(contigs_fpath), '-l', 'E', '-m', 'genome', '-f', '-c', str(busco_threads)],
                                                           out_dirpath) for index, contigs_fpath in enumerate(contigs_fpaths))

    # saving results
    for i, contigs_fpath in enumerate(contigs_fpaths):
        report = reporting.get(contigs_fpath)
        total, complete, part = results[i]
        report.add_field(reporting.Fields.CORE_COMPLETE, ('%.2f' % (float(complete)*100.0/total)))
        report.add_field(reporting.Fields.CORE_PART, ('%.2f' % (float(part)*100.0/total)))
    logger.print_timestamp()
    logger.info('Done.')

