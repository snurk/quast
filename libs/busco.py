import os
import platform
from Busco import busco_v1
from libs import reporting, qconfig, qutils
from libs.log import get_logger
logger = get_logger(qconfig.LOGGER_DEFAULT_NAME)

busco_dirpath = os.path.join(qconfig.LIBS_LOCATION, 'Busco')
if platform.system() == 'Darwin':
    hmmer_dirpath = os.path.join(busco_dirpath, 'hmmer-3.1b2-mac')
else:
    hmmer_dirpath = os.path.join(busco_dirpath, 'hmmer-3.1b2')

def all_required_binaries_exist(bin_dirpath, binary):
    if not os.path.isfile(os.path.join(bin_dirpath, binary)):
        return False
    return True

def do(contigs_fpaths, out_dirpath, clade):

    logger.print_timestamp()
    logger.info('Running BUSCO...')
    if qconfig.archaea:
        logger.info('BUSCO does not work with archaeal genomes.')
        return

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
            logger.info('Failed finding conservative genes.')
            return

    if not os.path.isdir(out_dirpath):
        os.makedirs(out_dirpath)

    log_path = os.path.join(out_dirpath, 'busco.log')
    err_path = os.path.join(out_dirpath, 'busco.err')
    open(log_path, 'w').close()
    open(err_path, 'w').close()
    results = []
    for index, contigs_fpath in enumerate(contigs_fpaths):
        if not os.path.isabs(contigs_fpath):
            current_dir = os.getcwd()
            contigs_fpath = os.path.join(current_dir, contigs_fpath)
        results.append(busco_v1.do(['-n', str(index), '-in', contigs_fpath,'-o', qutils.name_from_fpath(contigs_fpath), '-l', clade, '-m', 'genome', '-f', '-c', str(qconfig.max_threads)],
                                                         out_dirpath))
    logger.print_timestamp()
    # saving results
    for i, contigs_fpath in enumerate(contigs_fpaths):
        report = reporting.get(contigs_fpath)
        total, complete, part = results[i]
        report.add_field(reporting.Fields.CORE_COMPLETE, ('%.2f' % (float(complete)*100.0/total)))
        report.add_field(reporting.Fields.CORE_PART, ('%.2f' % (float(part)*100.0/total)))

    logger.info('Done.')

