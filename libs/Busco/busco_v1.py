#!/usr/bin/python

# BUSCO - Benchmarking sets of Universal Single-Copy Orthologs.

# Copyright (C) 2015 E. Zdobnov lab: F. Simao Neto
# <felipe.simao@unige.ch> based on code by R. Waterhouse.

# BUSCO is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# BUSCO is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# -------------------------------------------------------------------------------#
import Queue

import os
from os.path import join as join
import argparse
from collections import deque
import threading
import time
import platform
from libs import qconfig, qutils
from libs.log import get_logger

logger = get_logger(qconfig.LOGGER_DEFAULT_NAME)

import shutil
import shlex


busco_dirpath = join(qconfig.LIBS_LOCATION, 'Busco')

augustus_short_dirpath = join(busco_dirpath, 'augustus-3.0.3')

if platform.system() == 'Darwin':
    sed_cmd = "sed -i '' "
    hmmer_dirpath = join(busco_dirpath, 'hmmer-3.1b2-mac/src')
    augustus_dirpath = join(busco_dirpath, 'augustus-3.0.3/bin-mac')
else:
    sed_cmd = 'sed -i '
    hmmer_dirpath = join(busco_dirpath, 'hmmer-3.1b2/src')
    augustus_dirpath = join(busco_dirpath, 'augustus-3.0.3/bin')


def hmmer_fpath(fname):
    return join(hmmer_dirpath, fname)


def blast_fpath(fname):
    blast_dirpath = os.path.join(qconfig.LIBS_LOCATION, 'blast', qconfig.platform_name)
    return os.path.join(blast_dirpath, fname)


def august_fpath(fname):
    return join(augustus_dirpath, fname)


def do(f_args, output_dir):
    # ------------------------------------ Argument parser START ----------------------------------------#
    start_time = time.time()
    parser = argparse.ArgumentParser(
        description='Welcome to the Benchmarking set of Universal Single Copy Orthologs (BUSCO).\n\n For further usage information, please check the README file provided with this distrubution.',
        usage='busco_v1.py -in [SEQUENCE_FILE] -l [LINEAGE] -o [OUTPUT_NAME] [OTHER OPTIONS]')
    parser.add_argument('-g', '--genome', '-in', metavar='FASTA FILE', type=str,
                        help='Input file in fasta format.\nCan be a genome, proteome or transcriptome. Default analysis is run on the genome mode, for other files please specify the mode with (-m [MODE])\n')  # genome assembly file
    parser.add_argument('-c', '--cpu', metavar='N', type=str,
                        help='Number of threads/cores to use.')  # Number of available threads
    parser.add_argument('-a', '--abrev', '-o', metavar='output', type=str,
                        help='How to name output and temporary files.')  # Four letter abbreviation for use with genome assembly
    parser.add_argument('--ev', '-e', '-ev', metavar='N', type=float,
                        help='E-value cutoff for BLAST searches. (Default: 0.01)')  # evalue option
    parser.add_argument('-m', '--mode', metavar='mode', type=str,
                        help='which module to run the analysis to run, valid modes are \'all\'(genome assembly), \'OGS\' (gene set / proteome) and \'Trans\' (transcriptome).\n Defaults to \'all\'')
    parser.add_argument('-l', '--clade', '--lineage', metavar='lineage', type=str,
                        help='Which BUSCO lineage to be used.')  # lineage
    parser.add_argument('-f', action='store_true', default=False, dest='force',
                        help='Force rewrting of existing files. Must be used when output files with the provided name already exist.')
    parser.add_argument('-sp', '--species', default='generic', metavar='species', type=str,
                        help='Name of existing Augustus species gene finding metaparameters. (Default: generic)')
    parser.add_argument('-flank', '--flank', '-F', metavar='flanks', type=int,
                        help='Flanking sequence size for candidate regions. If not provided, flank size is calculated based on genome size with a range from 5 to 20 Kbp.')
    parser.add_argument('-Z', '--size', metavar='dbsize', type=int,
                        help='HMM library total size (Z). Important if using external datasets')
    parser.add_argument('--long', action='store_true', default=False, dest='long',
                        help='Optimization mode Augustus self-training (Default: Off) adds ~20h extra run time, but can improve results for some non-model organisms')
    parser.add_argument('-n', default='1', type=str, help='File index')

    args = vars(parser.parse_args(f_args))  #parse the arguments
    index = int(args['n'])
    assembly_name = args['abrev']
    mainout = join(output_dir, 'run_%s/' % assembly_name)  #final output directory
    summary_path = join(output_dir, 'short_summary_%s' % assembly_name)
    log_path = join(output_dir, 'busco.log')
    err_path = join(output_dir, 'busco.err')
    logger.info('  ' + qutils.index_to_str(index) + assembly_name)

    if os.path.isfile(summary_path):
        logger.info('  ' + qutils.index_to_str(index) + 'Using existing BUSCO files...')
        results = open(summary_path)
        total_buscos, part_buscos, complete_buscos = 0, 0, 0
        for line in results:
            if 'Complete' in line and 'Single' in line:
                complete_buscos = int(line.split()[0])
            elif 'Fragmented' in line:
                part_buscos = int(line.split()[0])
            elif 'Total' in line:
                total_buscos = int(line.split()[0])
        return total_buscos, complete_buscos, part_buscos
    if not os.path.exists(mainout) and assembly_name is not None:
        os.makedirs(mainout)
    else:
        if not args['force']:
            logger.info(
                'A run with that name already exists!\nIf are sure you wish to rewrite existing files please use the -f option')
            raise SystemExit

    target_species = args['species']

    ev_cut = 0.01  # default e-value cuttof
    try:
        if args['ev'] != ev_cut and args['ev'] is not None:
            logger.info('WARNING: You are using a custom e-value cutoff')
            ev_cut = args['ev']
    except:
        pass

    valid_clade_info = {'bacteria': 107114, 'eukaryota': 41317}
    # print(args['clade'])
    try:
        if args['clade'] != None:
            clade = args['clade']
            if clade.startswith(('b', 'B')):
                clade = 'bacteria'
                maxflank = 5000
            elif clade.startswith(('eu', 'Eu', 'E')):
                clade = 'eukaryota'
                maxflank = 20000
            clade_name = clade.strip('/').split('/')[-1].lower()
            if clade_name in valid_clade_info:
                Z = valid_clade_info[clade_name]
            else:
                try:
                    Z = args['dbsize'];
                except:
                    print('Please indicate the size of the custom HMM database using the (-Z integer)')
                    raise SystemExit
    except:
        print(
        'Please indicate the full path to a BUSCO clade: Eukaryota, Metazoa, Arthropoda, Vertebrata or Fungi\nExample: -l /path/to/Arthropoda')
        raise SystemExit

    cpus = 1  # 1 core default
    try:
        if args['cpu'] != cpus and args['cpu'] is not None:
            cpus = args['cpu']
    except:
        pass

    modes = ['all', 'blast', 'hmmer', 'augustus', 'parser', 'hmmer+', 'OGS', 'transcriptome', 'trans', 'ogs',
             'genome']  # valid modes
    mode = 'genome'  # unless otherwise specified, run on all (mode for genome assembly)
    try:
        if args['mode'] is not None and args['mode'] in modes:
            mode = args['mode']
            if mode == 'ogs':
                mode = 'OGS'
            elif mode == 'all' or mode == 'genome':
                mode = 'genome'
            elif mode == 'transcriptome':
                mode = 'trans'
    except:
        logger.info(
            'Error: Unknown mode specified * %s *, please check the documentation for valid modes.' % args['mode'])
        raise SystemExit
    flank = 5000
    if mode == 'genome' or mode == 'blast':
        size = os.path.getsize(args['genome']) / 1000
        flank = int(size / 50)
        if flank < 5000:
            flank = 5000
        elif flank > maxflank:
            flank = maxflank
    # ------------------------------------ Argument parser END ----------------------------------------#
    def startQueue(commands, cpus):
        exitFlag = 0
        queueLock = threading.Lock()
        threadList=[]
        for i in range(int(cpus)):
            threadList.append("Thread-%s" % str(i+1))
        workQueue = Queue.Queue(len(commands))

        class myThread (threading.Thread):
            def __init__(self, threadID, name, q):
                threading.Thread.__init__(self)
                self.threadID = threadID
                self.name = name
                self.q = q

            def run(self):
                process_data(self.name, self.q)#,self.alpha)

        def process_data(threadName, q):
            while not exitFlag:
                if not workQueue.empty():
                    data = q.get()
                    qutils.call_subprocess(shlex.split(data[0]), stdout=open(data[1], 'w'), stderr=open(err_path, 'a'))
                    q.task_done()
                time.sleep(1)

        threads = []
        threadID = 1
        needed = len(commands)
        mark = 0

        # Create new threads
        for tName in threadList:
            mark += 1
            thread = myThread(threadID, tName, workQueue)
            thread.start()
            threads.append(thread)
            threadID += 1
            if mark >= needed:
                break

        # Fill the queue
        queueLock.acquire()
        for word in commands:
            workQueue.put(word)
        queueLock.release()

        # Wait for queue to empty
        while not workQueue.empty():
            pass
        # Notify threads it's time to exit
        exitFlag = 1

        # Wait for all threads to complete
        for t in threads:
            t.join()
        exitFlag = 0

    def measuring(nested):
        if isinstance(nested, str):
            return ('0')
        scaffolds = list(nested.keys())
        if len(nested) == 1:
            total_len = [0]
            for hit in nested[scaffolds[0]]:
                total_len[0] += hit[1] - hit[0]
        elif len(nested) > 1:
            total_len = [0] * len(nested)
            for entry in range(0, len(scaffolds)):
                for hit in nested[scaffolds[entry]]:
                    total_len[entry] += hit[1] - hit[0]
        try:
            return total_len
        except:
            pass

    def extract(path, group, type):
        count = 0
        f = open('%saugustus/%s' % (path, group))
        if group.endswith(('.1', '.2', '.3')):
            extract_fpath = '%saugustus_proteins/%s.fas.%s' % (path, group[:-6], group[-1])
        else:
            extract_fpath = '%saugustus_proteins/%s.fas' % (path, group[:-4])
        check = 0
        out = open(extract_fpath, 'w')
        while True:
            line = f.readline()
            if not line:
                break
            if line.startswith('# start gene'):
                line = f.readline()
                line = line.split()
                places = [line[0], line[3], line[4]]
            elif line.startswith('# protein'):
                line = line.strip().split('[')
                count += 1
                out.write('>%s%s[%s:%s-%s]\n' % (type, count, places[0], places[1], places[2]))
                if line[1][-1] == ']':
                    line[1] = line[1][:-1]
                out.write(line[1])
                check = 1
            else:
                if line.startswith('# end'):
                    check = 0
                    out.write('\n')
                elif check == 1:
                    line = line.split()[1]
                    if line[-1] == ']':
                        line = line[:-1]
                    out.write(line)
        out.close()
        if os.path.getsize(extract_fpath) == 0:
            os.remove(extract_fpath)

    def disentangle(deck):
        structure = deque([deck.popleft()])
        try:
            while 1:
                temp = deck.popleft()
                start = temp[0];
                end = temp[1]
                for i in range(0, len(structure)):
                    ds = structure[i][0];
                    de = structure[i][1]
                    if start < ds and end < ds:  #fully before
                        if i == 0:  #first entry, just appendleft
                            structure.appendleft(temp)
                            break
                        else:
                            new = structure[0:i];
                            new.append(temp)
                            for z in range(i, len(structure)):
                                new.append(structure[z])
                            break
                    elif ds > start > end > ds:  #end overlaps inside, but the start is before
                        structure[i][0] = start
                        break
                    elif ds < start < de < end:  #start overlaps inside, but the end is after
                        structure[i][1] = end
                        break
                    elif start > de and end > de:  #fully after
                        if i == len(structure) - 1:  #only if its the last entry can it be safely added to structure
                            structure.append(temp)
                    elif start < ds and end > de:  #current structure is found fully inside the current entry
                        structure[i] = temp
        except:
            return structure


    def gargantua(deck):
        total = 0
        for entry in deck:
            total += entry[1] - entry[0]
        return total


    def shrink(number):
        number = number * 100
        if number >= 10:
            number = str(number)[:2]
        elif 10 > number > 0:
            number = str(number)[:3]
        return number

    #---------------------------BLAST steps START -------------------------------------------#

    #Make a blast database and run tblastn

    if mode == 'genome' or mode == 'blast' or mode == 'trans':
        logger.debug('  ' + qutils.index_to_str(index) + 'Running tBlastN')
        out_fpath = join(mainout, assembly_name)
        cmd = blast_fpath('makeblastdb') + (' -in %s -dbtype nucl -out %s' % (args['genome'], out_fpath))
        qutils.call_subprocess(shlex.split(cmd), stdout=open(log_path, 'a'), stderr=open(err_path, 'a'))
        logger.debug('')
        cmd = blast_fpath('tblastn') + (' -num_threads %s -query %s/ancestral -db %s -out %s_tblastn -outfmt 7' % (
            cpus, join(busco_dirpath, clade), out_fpath, out_fpath))
        qutils.call_subprocess(shlex.split(cmd), stdout=open(log_path, 'a'), stderr=open(err_path, 'a'))
        logger.debug('')

    coord_path = '%s_coordinates' % mainout

    #Get coordinates for a genome analysis
    if mode == 'genome' or mode == 'blast':
        logger.debug('  ' + qutils.index_to_str(index) + 'Getting coordinates for candidate regions!')
        f = open('%s_tblastn' % out_fpath)  #open input file
        out = open(coord_path, 'w')  #open Coordinates output file
        dic = {}
        coords = {}
        for i in f:
            if i.startswith('#'):
                pass
            else:
                line = i.strip().split()
                name = line[0]
                scaff = line[1]
                hitstart = int(line[6])
                hitend = int(line[7])
                postart = int(line[8])
                posend = int(line[9])
                e_val = float(line[10])
                sizer = int(line[3])
                if posend < postart:  #for minus-strand genes, invert coordinates for convenience
                    temp = posend
                    posend = postart
                    postart = temp
                if name not in dic.keys():  #create new entry in dictionary for current BUSCO
                    dic[name] = [scaff]
                    coords[name] = {}
                    coords[name][scaff] = [postart, posend, deque([[hitstart, hitend]]), sizer]
                elif scaff not in dic[name] and len(dic[name]) < 3:  #get just the top3 scoring regions
                    dic[name].append(scaff)
                    coords[name][scaff] = [postart, posend, deque([[hitstart, hitend]]), sizer]
                elif scaff in dic[name] and e_val < ev_cut:  #scaffold already checked, now update coordinates
                    if postart < coords[name][scaff][0] and coords[name][scaff][
                        0] - postart <= 50000:  #starts before, and withing 50kb of current position
                        coords[name][scaff][0] = postart
                        coords[name][scaff][2].append([hitstart, hitend])
                    if posend > coords[name][scaff][1] and posend - coords[name][scaff][
                        1] <= 50000:  #ends after and within 50 kbs
                        coords[name][scaff][1] = posend
                        coords[name][scaff][3] = hitend
                        coords[name][scaff][2].append([hitstart, hitend])
                    elif postart > coords[name][scaff][0] and postart < coords[name][scaff][
                        1]:  #starts inside current coordinates
                        if posend < coords[name][scaff][1]:
                            coords[name][scaff][2].append(
                                [hitstart, hitend])  #if ending inside, just add alignemnt positions to deque
                        elif posend > coords[name][scaff][1]:  #if ending after current coordinates, extend
                            coords[name][scaff][2][1] = posend
                            coords[name][scaff][2].append([hitstart, hitend])
        for i in coords:
            contest = {}
            maxi = 0
            for contig in coords[i]:
                sizer = disentangle(coords[i][contig][2])
                if contig not in contest:
                    contest[contig] = 0
                size = gargantua(sizer)
                contest[contig] = size
                if size > maxi:
                    maxi = size
            for contig in contest:
                #if contest[contig]>0.7*maxi:
                out.write(
                    '%s\t%s\t%s\t%s\n' % (
                        i, contig, max(0, coords[i][contig][0] - flank), coords[i][contig][1] + flank))
        out.close()

    #Get coordinates, candidate regions and translate sequences (transcriptome analysis)
    if mode == 'transcriptome' or mode == 'trans':
        logger.debug('  ' + qutils.index_to_str(index) + 'Getting coordinates for candidate transcripts!')
        f = open('%s_tblastn' % mainout)  #open input file
        dic = {}
        transdic = {}
        for i in f:  #get a dictionary of BUSCO matches vs candidate scaffolds
            if i.startswith('#'):
                pass
            else:
                line = i.strip().split()
                name = line[0]
                scaff = line[1]
                e_val = float(line[10])
                leng = int(line[3])
                if name not in dic.keys() and e_val <= ev_cut:
                    dic[name] = [scaff]
                    maxi = leng
                    transdic[scaff] = name
                elif e_val <= ev_cut and scaff not in dic[name] and len(dic[name]) < 3 and leng >= 0.7 * maxi:
                    dic[name].append(scaff)
                    transdic[scaff] = name

        scaff_list = []  #list of unique scaffolds with buscos matches
        for busco in dic:
            for scaff in dic[busco]:
                if scaff not in scaff_list:
                    scaff_list.append(scaff)
        logger.debug('  ' + qutils.index_to_str(index) + 'Extracting candidate transcripts!')
        f = open(args['genome'])
        check = 0
        for i in f:
            if i.startswith('>'):
                i = i.strip().split()
                i = i[0][1:]
                if i in scaff_list:
                    scaff_path = join(output_dir, assembly_name, i + assembly_name)
                    out = open('%s_.temp' % scaff_path, 'w')
                    out.write('>%s\n' % i)
                    check = 1
                else:
                    check = 0
            elif check == 1:
                out.write(i)
        out.close()
        if not os.path.exists('%stranslated_proteins' % mainout):
            os.makedirs('%stranslated_proteins' % mainout)
        files = os.listdir('.')
        lista = []
        for entry in files:
            if entry.endswith(assembly_name + '_.temp'):
                lista.append(entry)

        logger.debug('  ' + qutils.index_to_str(index) + 'Translating candidate transcripts !')
        for entry in lista:
            #logger.info(entry)a=input('press to continue')
            cmd = ('transeq -clean -frame 6 -trim -sequence %(scaffold)s -outseq %(translated_scaffold)s.fas' % {
                'scaffold': entry,
                'translated_scaffold': join(mainout, 'translated_proteins', entry.split(assembly_name)[0] + '_ts')})
            qutils.call_subprocess(shlex.split(cmd), stdout=open(log_path, 'a'), stderr=open(err_path, 'a'))
        f2 = open(join(busco_dirpath, clade, 'scores_cutoff'))  #open target scores file
        #Load dictionary of HMM expected scores and full list of groups
        score_dic = {}
        for i in f2:
            i = i.strip().split()
            try:
                score_dic[i[0]] = float(i[1])  #float [1] = mean value [2] = minimum value
            except:
                pass
        totalbuscos = len(list(score_dic.keys()))


    #---------------------------AUGUSTUS steps START -------------------------------------------#

    #Step-3
    #Extract candidate contigs/scaffolds from genome assembly (necessary because augustus doesn't handle multi-fasta files when running on a specific target region)
    if mode == 'genome' or mode == 'augustus':
        # target_species=species_list[0]
        logger.debug('  ' + qutils.index_to_str(index) + 'pre-Augustus scaffold extraction  ')
        coord = open(coord_path)
        dic = {}
        scaff_list = []
        for i in coord:
            i = i.strip().split()
            if len(i) != 2:
                dic[i[0]] = [i[1], i[2], i[3]]
                if i[1] not in scaff_list:
                    scaff_list.append(i[1])
        f = open(args['genome'])
        check = 0
        for i in f:
            if i.startswith('>'):
                i = i.split()
                i = i[0][1:]
                if i in scaff_list:
                    scaff_path = join(mainout, i + assembly_name)
                    out = open('%s_.temp' % scaff_path, 'w')
                    out.write('>%s\n' % (i))
                    check = 1
                else:
                    check = 0
            elif check == 1:
                out.write(i)
        out.close()

        #Step-4
        #Augustus search on candidate regions using the pre-built Block profiles (msa2prfl.pl)
        logger.info('  ' + qutils.index_to_str(index) + 'Running Augustus prediction')
        aug_out_path = join(mainout, 'augustus')
        if not os.path.exists(aug_out_path):
            os.makedirs(aug_out_path)
        dic = {}
        commands = []
        f = open(coord_path)
        for i in f:
            i = i.strip().split('\t')
            name = i[0]
            if name not in dic:
                dic[name] = [[i[1], i[2], i[3]]]  #scaffold,start and end
            elif name in dic:
                dic[name].append([i[1], i[2], i[3]])
        for i in dic:
            if len(dic[i]) > 1:
                for z in range(0, len(dic[i])):
                    command = august_fpath('augustus') + (' --proteinprofile=%(prot_profile)s --predictionStart=%(start_coord)s '
                    '--predictionEnd=%(end_coord)s --species=%(species)s \"%(scaffold)s\"' % {
                    'prot_profile': join(busco_dirpath, clade, 'prfl/' + i + '.prfl'), 'start_coord': dic[i][z][1], 'end_coord': dic[i][z][2],
                    'species': target_species, 'scaffold': join(mainout, dic[i][z][0] + args['abrev'] + '_.temp')})
                    out_aug = join(mainout, 'augustus/%s.out.%s' % (i, str(z+1)))
                    commands.append((command, out_aug))
            else:
                command = august_fpath('augustus') + (' --proteinprofile=%(prot_profile)s --predictionStart=%(start_coord)s '
                '--predictionEnd=%(end_coord)s --species=%(species)s \"%(scaffold)s\" ' % {
                'prot_profile': join(busco_dirpath, clade, 'prfl/' + i + '.prfl'), 'start_coord': dic[i][0][1], 'end_coord': dic[i][0][2],
                'species': target_species, 'scaffold': join(mainout, dic[i][0][0] + args['abrev'] + '_.temp')})
                out_aug = join(mainout, 'augustus/%s.out' % i)
                commands.append((command, out_aug))
        startQueue(commands, cpus)
        # ---------------------------AUGUSTUS steps END -------------------------------------------#
        # ---------------------------HMMER steps START -------------------------------------------#

    if mode == 'genome' or mode == 'hmmer':  # should be augustus\
        # STEP-1 EXTRACT AUGUSTUS PROTEINS
        logger.debug('  ' + qutils.index_to_str(index) + 'Extracting predicted proteins  ')
        files = os.listdir(join(mainout, 'augustus'))
        count = 0
        check = 0
        aug_prot_path = join(mainout, 'augustus_proteins')
        for i in files:
            cmd = sed_cmd + " \'1,3d\' %saugustus/%s" % (mainout, i)
            qutils.call_subprocess(shlex.split(cmd), stderr=open(err_path, 'a'))
        if not os.path.exists(aug_prot_path):
            os.makedirs(aug_prot_path)
        for i in files:
            extract(mainout, i, 'g')

    #Run HMMer (genome mode)
    if mode == 'genome' or mode == 'hmmer':
        logger.debug('  ' + qutils.index_to_str(index) + 'Running HMMER to confirm orthology of predicted proteins  ')
        hmmer_out_path = join(mainout, 'hmmer_output')
        files = os.listdir(aug_prot_path)
        if not os.path.exists(hmmer_out_path):
            os.makedirs(hmmer_out_path)
        for i in files:
            fpath = join(aug_prot_path, i)
            if i.endswith('.fas'):
                name = i[:-4]
                cmd = hmmer_fpath('hmmsearch') + (
                    ' --domtblout %(output_file)s.out -Z %(db_size)s --cpu %(cpu)s %(group_file)s.hmm %(input_file)s ' %
                    {'input_file': fpath, 'db_size': Z, 'cpu': cpus,
                     'group_file': join(busco_dirpath, clade, 'hmms', name),
                     'output_file': join(hmmer_out_path, name), 'erroutput': err_path})
                qutils.call_subprocess(shlex.split(cmd), stdout=open(os.devnull, 'w'), stderr=open(err_path, 'a'))
            elif i.endswith(('.1', '.2', '.3')):
                name = i[:-6]
                cmd = hmmer_fpath('hmmsearch') + (
                    ' --domtblout %(output_file)s -Z %(db_size)s --cpu %(cpu)s %(group_file)s.hmm %(input_file)s ' %
                    {'input_file': fpath, 'db_size': Z, 'cpu': cpus,
                     'group_file': join(busco_dirpath, clade, 'hmms', name),
                     'output_file': join(hmmer_out_path, '%s.out.%s' % (name, i[-1]))})
                qutils.call_subprocess(shlex.split(cmd), stdout=open(os.devnull, 'w'), stderr=open(err_path, 'a'))

    #Run HMMer (transcriptome mode)
    if mode == 'trans' or mode == 'transcriptome':
        logger.debug('  ' + qutils.index_to_str(index) + 'Running HMMER to confirm transcript orthology  ')
        trans_protein_path = join(mainout, 'translated_proteins')
        files = os.listdir(trans_protein_path)
        if not os.path.exists(hmmer_out_path):
            os.makedirs(hmmer_out_path)
        group = ''
        grouplist = []
        for i in files:
            if i.endswith('.fas'):
                fpath = join(trans_protein_path, i)
                f = open(trans_protein_path)
                name = i[:-7]
                group = transdic[name]
                if group not in grouplist:
                    grouplist.append(group)
                    cmd = hmmer_fpath('hmmsearch') + (
                        ' --domtblout %(output_file)s.out.1 -Z %(db_size)s --cpu %(cpu)s %(group_file)s.hmm %(input_file)s ' % \
                        {'input_file': fpath, 'db_size': Z, 'cpu': '1',
                         'group_file': join(busco_dirpath, clade, 'hmms', group),
                         'output_file': join(hmmer_out_path, group)})
                    qutils.call_subprocess(shlex.split(cmd), stdout=open(log_path, 'a'), stderr=open(err_path, 'a'))
                else:
                    grouplist.append(group)
                    cmd = hmmer_fpath('hmmsearch') + (
                        ' --domtblout %(output_file)s.out.%(count)s -Z %(db_size)s --cpu %(cpu)s %(group_file)s.hmm %(input_file)s ' % \
                        {'input_file': fpath, 'db_size': Z, 'cpu': '1',
                         'group_file': join(busco_dirpath, clade, 'hmms', group),
                         'output_file': join(hmmer_out_path, group),
                         'count': str(grouplist.count(group))})
                    qutils.call_subprocess(shlex.split(cmd), stdout=open(log_path, 'a'), stderr=open(err_path, 'a'))

    #OGS/Proteome module
    if mode == 'OGS':
        if not os.path.exists(hmmer_out_path):
            os.makedirs(hmmer_out_path)
        files = os.listdir(join(busco_dirpath, clade, 'hmms'))
        f2 = open(join(busco_dirpath, clade, 'scores_cutoff'))  #open target scores file
        #Load dictionary of HMM expected scores and full list of groups
        score_dic = {}
        for i in f2:
            i = i.strip().split()
            try:
                score_dic[i[0]] = float(i[1])  #[1] = mean value [2] = minimum value
            except:
                pass
        totalbuscos = len(list(score_dic.keys()))
        for i in files:
            name = i[:-4]
            if name in score_dic:
                cmd = hmmer_fpath('hmmsearch') + (
                    '--domtblout %(output_file)s.out -Z %(db_size)s --cpu %(cpu)s %(group_file)s.hmm %(input_file)s ' % \
                    {'input_file': args['genome'], 'db_size': Z, 'cpu': cpus,
                     'group_file': join(busco_dirpath, clade, 'hmms', name),
                     'output_file': join(hmmer_out_path, name)})
                qutils.call_subprocess(shlex.split(cmd), stdout=open(log_path, 'a'), stderr=open(err_path, 'a'))

    ###  *get list to be re-run
    if mode == 'genome' or mode == 'hmmer':
        logger.debug('  ' + qutils.index_to_str(index) + 'Parsing HMMER results  ')
        #Open the output file if no name was specified the default name will be used
        f2 = open(join(busco_dirpath, clade, 'scores_cutoff'))  #open target scores file
        #Load dictionary of HMM expected scores and full list of groups
        score_dic = {}
        for i in f2:
            i = i.strip().split()
            try:
                score_dic[i[0]] = float(i[1])
            except:
                pass
        totalbuscos = float(len(list(score_dic.keys())))
        f = open(coord_path)
        dic = {}
        for i in f:
            i = i.strip().split('\t')
            name = i[0]
            if name not in dic:
                dic[name] = [[i[1], i[2], i[3]]]  #scaffold,start and end
            elif name in dic:
                dic[name].append([i[1], i[2], i[3]])

    ####Make summary

    #Categorizing genes found in Complete multi-copy and partial hits
    leng_dic = {}
    sd_dic = {}
    complete = []
    frag = []
    done = []
    cc = []
    fcc = 0
    mcc = []
    unique = []
    if mode == 'genome' or mode != 'OGS' or mode == 'report' or mode == 'hmmer':
        temp = os.listdir(hmmer_out_path)
        files = []
        for i in temp:
            if i.endswith(('.out', '.1', '.2', '.3')):
                files.append(i)
        f = open(join(busco_dirpath, clade, 'lengths_cutoff'))
        for line in f:
            line = line.strip().split()
            leng_dic[line[0]] = float(line[3])
            sd_dic[line[0]] = float(line[2])
        for entry in files:
            f = open('%shmmer_output/%s' % (mainout, entry))
            hit_dic = {}
            for line in f:
                if line.startswith('#'):
                    pass
                else:
                    line = line.strip().split()
                    score = float(line[7])
                    group = line[3]
                    prot = line[0]
                    tlen = int(line[2])
                    qlen = int(line[5])
                    if tlen > 30 * qlen:
                        pass
                    else:
                        if prot not in hit_dic.keys() and score >= score_dic[group]:
                            hit_dic[prot] = [[int(line[15]), int(line[16])]]
                        elif score >= score_dic[group]:
                            hit_dic[prot].append([int(line[15]), int(line[16])])
            length = measuring(hit_dic)
            try:  #get maximum length of the putative gene in question
                if len(length) == 1:
                    length = length[0]
                else:
                    length = max(length) + 1
                sigma = abs(leng_dic[group] - length) / sd_dic[group]
                if sigma <= 2:
                    complete.append(entry)
                    cc.append(group)
                elif sigma > 2:
                    frag.append(entry)
            except:
                pass
        #check the multi hits
        for entry in complete:
            if entry.endswith('.out'):
                name = entry[:-4]
            else:
                name = entry[:-6]
            if name in done:
                if name not in mcc:
                    mcc.append(name)
            done.append(name)
        for i in cc:
            if i not in mcc:
                unique.append(i)
        for entry in frag:
            if entry.endswith('.out'):
                name = entry[:-4]
            else:
                name = entry[:-6]
            if name not in done and entry not in complete:
                done.append(name)
                fcc += 1

    if mode == 'OGS':
        complete = {}
        frag = {}
        done = []
        fcc = []
        temp = os.listdir(hmmer_out_path)
        files = []
        for i in temp:
            if i.endswith(('.out', '.1', '.2', '.3')):
                files.append(i)
        f = open(join(busco_dirpath, clade, 'lengths_cutoff'))
        for line in f:
            line = line.strip().split()
            leng_dic[line[0]] = float(line[3])
            sd_dic[line[0]] = float(line[2])
        for entry in files:
            f = open('%shmmer_output/%s' % (mainout, entry))
            hit_dic = {}
            for line in f:
                if line.startswith('#'):
                    pass
                else:
                    line = line.strip().split()
                    score = float(line[7])
                    group = line[3]
                    prot = line[0]
                    tlen = int(line[2])
                    qlen = int(line[5])
                    prediction = line[0]
                    if group not in complete:
                        complete[group] = []
                        frag[group] = []
                    if tlen > 30 * qlen:
                        pass
                    else:
                        if prot not in hit_dic.keys() and score >= score_dic[group]:
                            hit_dic[prot] = [[int(line[15]), int(line[16]), line[7]]]
                        elif score >= score_dic[group]:
                            hit_dic[prot].append([int(line[15]), int(line[16]), line[7]])
            lengths = measuring(hit_dic)
            try:  #get maximum length of the putative gene in question
                if len(lengths) == 1:
                    length = lengths[0]
                    sigma = abs(leng_dic[group] - length) / sd_dic[group]
                    if sigma <= 2:
                        complete[group].append([list(hit_dic.keys())[lengths.index(length)],
                                                hit_dic[list(hit_dic.keys())[lengths.index(length)]][0][2], length])
                    elif sigma > 2:
                        frag[group].append(list(hit_dic.keys())[lengths.index(length)])
                else:
                    for length in lengths:
                        #length=max(lengths)+1
                        sigma = abs(leng_dic[group] - length) / sd_dic[group]
                        if sigma <= 2:
                            complete[group].append([list(hit_dic.keys())[lengths.index(length)],
                                                    hit_dic[list(hit_dic.keys())[lengths.index(length)]][0][2],
                                                    length])
                        elif sigma > 2:
                            frag[group].append(list(hit_dic.keys())[lengths.index(length)])
            except:
                pass
        #check the multi hits
        for entry in complete:
            if len(complete[entry]) == 0:
                pass
            elif len(complete[entry]) == 1:  #complete
                cc.append(entry)
            elif len(complete[entry]) > 1:
                mcc.append(entry)
        for entry in frag:
            if len(complete[entry]) != 0:
                pass
            elif frag[entry] != []:
                fcc.append(entry)

    #summarize results, logger.info and write to output files
    summary = open(summary_path, 'w')

    summary.write(
        '#Summarized BUSCO benchmarking for file: %s\n#BUSCO was run in mode: %s\n\n' % (args['genome'], mode))
    if mode != 'OGS' and mode != 'trans':
        summary.write('Summarized benchmarks in BUSCO notation:\n\tC:%s%%[D:%s%%],F:%s%%,M:%s%%,n:%s\n\n' % (
            shrink((len(set(cc)) + len(set(mcc))) / totalbuscos), shrink(len(set(mcc)) / totalbuscos),
            shrink(fcc / totalbuscos), shrink((totalbuscos - (len(set(cc)) + fcc)) / totalbuscos), int(totalbuscos)))
    elif mode == 'OGS':
        summary.write('Summarized benchmarks in BUSCO notation:\n\tC:%s%%[D:%s%%],F:%s%%,M:%s%%,n:%s\n\n' % (
            shrink((len(set(cc)) + len(set(mcc))) / totalbuscos), shrink(len(set(mcc)) / totalbuscos),
            shrink(len(fcc) / totalbuscos),
            shrink((totalbuscos - (len(set(cc)) + len(set(mcc)) + len(fcc))) / totalbuscos),
            int(totalbuscos)))
    elif mode == 'trans':
        summary.write('Summarized benchmarks in BUSCO notation:\n\tC:%s%%[D:%s%%],F:%s%%,M:%s%%,n:%s\n\n' % (
            shrink(len(set(cc)) / totalbuscos), shrink(len(set(mcc)) / totalbuscos), shrink(fcc / totalbuscos),
            shrink((totalbuscos - (len(set(cc)) + fcc)) / totalbuscos), int(totalbuscos)))

    summary.write('Representing:\n')
    if mode != 'trans' and mode != 'OGS':
        summary.write('\t%s\tComplete Single-copy BUSCOs\n' % (len(set(cc))))
        summary.write('\t%s\tComplete Duplicated BUSCOs\n' % (len(set(mcc))))
    elif mode == 'OGS':
        summary.write('\t%s\tComplete Single-copy BUSCOs\n' % (len(set(cc)) + len(set(mcc))))
        summary.write('\t%s\tComplete Duplicated BUSCOs\n' % (len(set(mcc))))
    elif mode == 'trans':
        summary.write('\t%s\tComplete Single-copy BUSCOs\n' % (len(set(cc)) - len(set(mcc))))
        summary.write('\t%s\tComplete Duplicated BUSCOs\n' % (len(set(mcc))))
    if mode != 'OGS':
        summary.write('\t%s\tFragmented BUSCOs\n' % fcc)
        summary.write('\t%s\tMissing BUSCOs\n' % (int(totalbuscos) - (len(set(cc)) + fcc)))
    elif mode == 'OGS':
        summary.write('\t%s\tFragmented BUSCOs\n' % (len(fcc)))
        summary.write('\t%s\tMissing BUSCOs\n' % (int(totalbuscos) - (len(set(cc)) + len(set(mcc)) + len(fcc))))

    summary.write('\t%s\tTotal BUSCO groups searched\n' % int(totalbuscos))
    summary.close()
    summary = open(join(mainout, 'full_table_%s' % assembly_name), 'w')
    #write correct header
    if mode == 'genome' or mode == 'report':
        summary.write('#BUSCO_group\tStatus\tScaffold\tStart\tEnd\tBitscore\tLength\n')
    elif mode == 'OGS':
        summary.write('#BUSCO_group\tStatus\tGene\tBitscore\tLength\n')
    elif mode == 'trans' or mode == 'transcriptome':
        summary.write('#BUSCO_group\tStatus\tTranscript\tBitscore\tLength\n')

    temp = os.listdir('%shmmer_output' % mainout)
    done = []
    files = []
    for i in temp:
        if i.endswith(('.out', '.1', '.2', '.3')):
            files.append(i)

    for i in files:
        if i.endswith('.out'):
            name = i[:-4]
            marker = 0
        elif i.endswith(('.1', '.2', '.3')):
            name = i[:-6]
            marker = int(i[-1]) - 1
        f = open('%shmmer_output/%s' % (mainout, i))
        score = []
        hit_dic = {}
        for line in f:
            if line.startswith('#'):
                pass
            else:
                line = line.strip().split()
                score.append(float(line[7]))
                group = line[3]
                prot = line[0]
                tlen = int(line[2])
                qlen = int(line[5])
                prediction = line[0]
                if prot not in hit_dic.keys() and float(line[7]) >= score_dic[group]:
                    hit_dic[prot] = [[int(line[15]), int(line[16]), line[7]]]
                elif float(line[7]) >= score_dic[group]:
                    hit_dic[prot].append([int(line[15]), int(line[16]), line[7]])
        length = measuring(hit_dic)
        if mode == 'genome' or mode == 'report' or mode == 'hmmer':
            if hit_dic == {}:
                pass
            elif i in complete and name not in mcc:
                summary.write('%s\tComplete\t%s\t%s\t%s\t%s\t%s\n' % (
                    name, dic[group][marker][0], dic[group][marker][1], dic[group][marker][2], max(score),
                    max(length) + 1))
            elif i in complete and name in mcc:
                summary.write('%s\tDuplicated\t%s\t%s\t%s\t%s\t%s\n' % (
                    name, dic[group][marker][0], dic[group][marker][1], dic[group][marker][2], max(score),
                    max(length) + 1))
            elif i in frag and name not in cc and name not in done:
                summary.write('%s\tFragmented\t%s\t%s\t%s\t%s\t%s\n' % (
                    name, dic[group][marker][0], dic[group][marker][1], dic[group][marker][2], max(score),
                    max(length) + 1))
        elif mode == 'OGS':
            if hit_dic == {}:
                pass
            elif name in complete and name not in mcc and name in cc:
                summary.write('%s\tComplete\t%s\t%s\t%s\n' % (name, complete[name][0][0], max(score), max(length) + 1))
            elif name in mcc:
                for entry in complete[name]:
                    summary.write('%s\tDuplicated\t%s\t%s\t%s\n' % (name, entry[0], entry[1], entry[2] + 1))
            elif name in fcc and name not in cc:
                summary.write('%s\tFragmented\t%s\t%s\t%s\n' % (name, frag[name][0], max(score), max(length) + 1))
        elif mode == 'trans' or mode == 'Transcriptome':
            if hit_dic == {}:
                pass
            elif i in complete and name not in mcc:
                summary.write('%s\tComplete\t%s\t%s\t%s\n' % (name, dic[group][marker], max(score), max(length) + 1))
            elif i in complete and name in mcc:
                summary.write('%s\tDuplicated\t%s\t%s\t%s\n' % (name, dic[group][marker], max(score), max(length) + 1))
            elif i in frag and name not in cc and name not in done:
                summary.write('%s\tFragmented\t%s\t%s\t%s\n' % (name, dic[group][marker], max(score), max(length) + 1))
    summary.close()
    f.close()

    f = open(join(mainout, 'full_table_%s' % assembly_name), 'r')
    lista = []
    for i in f:
        i = i.strip().split()
        if i[0] not in lista:
            lista.append(i[0])
    f.close()
    out = open(join(mainout, 'missing_buscos_list_%s' % assembly_name), 'w')  #get final list of missing buscos
    f = open(join(mainout, 'full_table_%s' % assembly_name), 'a')
    for i in score_dic.keys():
        if i in lista:
            pass
        else:
            out.write(i + '\n')
            f.write('%s\tMissing\n' % (i))
    out.close()
    f.close()

    #######retraining

    if mode == 'genome' or mode == 'genome' or mode == 'hmmer':

        if not os.path.exists('%sselected' % mainout):
            os.makedirs('%sselected' % mainout)
        if not os.path.exists('%sgffs' % mainout):
            os.makedirs('%sgffs' % mainout)
        if not os.path.exists('%sgb' % mainout):
            os.makedirs('%sgb' % mainout)

        f = open(mainout + ('full_table_%s' % assembly_name))
        lista = []
        re_run = []
        for line in f:
            status = line.split()[1]
            if status == 'Complete':
                lista.append(line.split()[0])
            elif status == 'Missing' or status == 'Fragmented':
                re_run.append(line.split()[0])

        files = os.listdir('%shmmer_output' % mainout)
        chosen = []
        for i in files:
            if i.endswith(('.out', '.1', '.2', '.3')):
                if i.endswith('.out'):
                    name = i[:-4]
                else:
                    name = i[:-6]
                if name in lista:
                    cmd = 'cp %shmmer_output/%s %sselected/' % (mainout, i, mainout)
                    qutils.call_subprocess(shlex.split(cmd), stdout=open(log_path, 'a'), stderr=open(err_path, 'a'))
                    chosen.append(i)

        for entry in chosen:
            f = open('%sselected/%s' % (mainout, entry))
            out = open('%sgffs/%s' % (mainout, entry), 'w')
            choicy = ''
            for line in f:
                if line.startswith('#'):
                    pass
                elif choicy == '':
                    choicy = line.split()[0].split('[')[0]
            f.close()
            f = open('%saugustus/%s' % (mainout, entry))
            check = 0
            for line in f:
                if line.startswith('# start gene'):
                    name = line.strip().split()[-1]
                    if name == choicy:
                        check = 1
                elif line.startswith('#'):
                    check = 0
                elif check == 1:
                    out.write(line)
            out.close()
        f.close()
        if len(chosen) > 0:
            train_set_fpath = '%straining_set_%s' % (mainout, assembly_name)
            gff_files = ','.join('%sgffs/%s' % (mainout, entry) for entry in chosen)
            cmd = join(augustus_short_dirpath, 'scripts/gff2gbSmallDNA.pl %s %s 1000 %s' % (
                gff_files, args['genome'], train_set_fpath))
            qutils.call_subprocess(shlex.split(cmd), stdout=open(log_path, 'a'), stderr=open(err_path, 'a'))
            logger.debug('  ' + qutils.index_to_str(index) + 'Training augustus gene predictor')
            os.chdir(augustus_short_dirpath)
            prokaryotic = ''
            if clade == 'bacteria':
                prokaryotic = '--prokaryotic'
            cmd = join(augustus_short_dirpath,
                       'scripts/new_species.pl --species=%s --AUGUSTUS_CONFIG_PATH=%s %s' % (
                           assembly_name, join(augustus_short_dirpath, 'config/'),
                           prokaryotic))  # create new species config file from template
            qutils.call_subprocess(shlex.split(cmd), stdout=open(log_path, 'a'), stderr=open(err_path, 'a'))
            cmd = august_fpath('etraining') + (' --species=%s %straining_set_%s --stopCodonExcludedFromCDS=false  --AUGUSTUS_CONFIG_PATH=%s ' % (
                assembly_name, mainout, assembly_name, join(augustus_short_dirpath, 'config/')))  # train on new training set (complete single copy buscos
            return_code = qutils.call_subprocess(shlex.split(cmd), stdout=open(log_path, 'a'), stderr=open(err_path, 'a'))

        if args['long']:
            logger.info(
                '  ' + qutils.index_to_str(index) + 'Optimizing augustus metaparameters, this may take around 20 hours')
            cmd = join(augustus_short_dirpath,
                       'scripts/optimize_augustus.pl --species=%s %straining_set_%s --AUGUSTUS_CONFIG_PATH=%s' % (
                           assembly_name, mainout, assembly_name,
                           join(augustus_short_dirpath, 'config/')))
            qutils.call_subprocess(shlex.split(cmd), stdout=open(log_path, 'a'),
                                   stderr=open(err_path, 'a'))  # train on new training set (complete single copy buscos)
            cmd = august_fpath('etraining') + (' --species=%s %straining_set_%s' % (
                assembly_name, mainout, assembly_name))  # train on new training set (complete single copy buscos)
            qutils.call_subprocess(shlex.split(cmd), stdout=open(log_path, 'a'), stderr=open(err_path, 'a'))
        if len(re_run) > 0 and len(chosen) > 0 and return_code == 0:
            logger.info('  ' + qutils.index_to_str(
                index) + 'Re-running failed predictions with different constraints, total number %s  ' % len(re_run))
            target_species = assembly_name
            commands = []
            hammers = []
            seds = []
            ripped = []
            for item in re_run:
                if item not in dic:  # no coordinates found
                    pass
                elif len(dic[item]) > 1:  # more than one target coordinate
                    count = 0
                    for entry in dic[item]:
                        count += 1
                        command = august_fpath('augustus') + (
                            ' --proteinprofile=%(clade)s/%(prot_profile)s.prfl --predictionStart=%(start_coord)s --predictionEnd=%(end_coord)s --species=%(species)s \"%(scaffold)s\"' %
                            {'prot_profile': 'prfl/' + item, 'start_coord': entry[1], 'end_coord': entry[2],
                             'clade': join(busco_dirpath, clade),
                             'species': target_species, 'scaffold': join(mainout, entry[0] + assembly_name + '_.temp')})
                        out_aug = join(mainout, 'augustus/%s.out.%s' % (item, count))
                        commands.append((command, out_aug))
                        command = sed_cmd + ' \'1,3d\' %s' % (mainout + 'augustus/' + item + '.out.' + str(count))
                        log_out = os.devnull
                        seds.append((command, log_out))
                        ripped.append(item + '.out.' + str(count))
                        command = hmmer_fpath('hmmsearch') + (
                            ' --domtblout %(output_file)s -Z %(db_size)s --cpu %(cpu)s %(group_file)s.hmm %(input_file)s' %
                            {'input_file': join(aug_prot_path, '%s.fas.%s' % (item, count)), 'db_size': Z,
                             'cpu': '1',
                             'group_file': join(busco_dirpath, clade, 'hmms', item),
                             'output_file': join(hmmer_out_path, '%s.out.%s' % (item, count))})
                        hammers.append((command, log_out))
                elif len(dic[item]) == 1:
                    entry = dic[item][0]
                    command = august_fpath('augustus') + (
                        ' --proteinprofile=%(clade)s/%(prot_profile)s.prfl --predictionStart=%(start_coord)s --predictionEnd=%(end_coord)s -'
                        '-species=%(species)s \"%(scaffold)s\" ' %
                        {'prot_profile': 'prfl/' + item, 'start_coord': entry[1], 'end_coord': entry[2],
                         'clade': join(busco_dirpath, clade),
                         'species': target_species, 'scaffold': join(mainout, entry[0] + assembly_name + '_.temp')})
                    out_aug = join(mainout, 'augustus/%s.out' % item)
                    log_out = os.devnull
                    commands.append((command, out_aug))
                    command = sed_cmd + ' \'1,3d\' %s' % (mainout + 'augustus/' + item + '.out')
                    seds.append((command, log_out))
                    ripped.append(item + '.out')
                    command = hmmer_fpath('hmmsearch') + (
                        ' --domtblout %(output_file)s.out -Z %(db_size)s  --cpu %(cpu)s %(group_file)s.hmm %(input_file)s.fas' %
                        {'input_file': mainout + 'augustus_proteins/' + item, 'db_size': Z, 'cpu': '1',
                         'group_file': join(busco_dirpath, clade, 'hmms', item),
                         'output_file': join(hmmer_out_path, item)})
                    hammers.append((command, log_out))

            ###retraining and running over
            startQueue(commands, cpus)
            startQueue(seds, cpus)
            for entry in ripped:
                extract(mainout, entry, 'p')
            startQueue(hammers, cpus)

    shutil.rmtree(join(augustus_short_dirpath, 'config/species', assembly_name), ignore_errors=True)

    #parse results and write final summary
    #Categorizing genes found in Complete multi-copy and partial hits
    leng_dic = {}
    sd_dic = {}
    complete = []
    frag = []
    done = []
    cc = []
    fcc = 0
    mcc = []
    unique = []

    if mode == 'genome' or mode == 'genome' or mode == 'hmmer':
        temp = os.listdir(hmmer_out_path)
        files = []
        for i in temp:
            if i.endswith(('.out', '.1', '.2', '.3')):
                files.append(i)
        f = open(join(busco_dirpath, clade, 'lengths_cutoff'))
        for line in f:
            line = line.strip().split()
            leng_dic[line[0]] = float(line[3])
            sd_dic[line[0]] = float(line[2])
        for entry in files:
            f = open(join(hmmer_out_path, entry))
            hit_dic = {}
            for line in f:
                if line.startswith('#'):
                    pass
                else:
                    line = line.strip().split()
                    score = float(line[7])
                    group = line[3]
                    prot = line[0]
                    tlen = int(line[2])
                    qlen = int(line[5])
                    if tlen > 30 * qlen:
                        pass
                    else:
                        if prot not in hit_dic.keys() and score >= score_dic[group]:
                            hit_dic[prot] = [[int(line[15]), int(line[16])]]
                        elif score >= score_dic[group]:
                            hit_dic[prot].append([int(line[15]), int(line[16])])
            length = measuring(hit_dic)
            try:  #get maximum length of the putative gene in question
                if len(length) == 1:
                    length = length[0]
                else:
                    length = max(length) + 1
                sigma = abs(leng_dic[group] - length) / sd_dic[group]
                if sigma <= 2:
                    complete.append(entry)
                    cc.append(group)
                elif sigma > 2:
                    frag.append(entry)
            except:
                pass
        #check the multi hits
        for entry in complete:
            if entry.endswith('.out'):
                name = entry[:-4]
            else:
                name = entry[:-6]
            if name in done:
                if name not in mcc:
                    mcc.append(name)
            done.append(name)
        for i in cc:
            if i not in mcc:
                unique.append(i)
        for entry in frag:
            if entry.endswith('.out'):
                name = entry[:-4]
            else:
                name = entry[:-6]
            if name not in done and entry not in complete:
                done.append(name)
                fcc += 1

                #summarize results, logger.info and write to output files
        summary = open(summary_path, 'w')

        logger.info('  %sComplete BUSCOs found: %s (%s duplicated), partially recovered: %s,\
                    total groups: %s' % (qutils.index_to_str(index), len(set(cc)), len(mcc), fcc, int(totalbuscos)))

        summary.write(
            '#Summarized BUSCO benchmarking for file: %s\n#BUSCO was run in mode: %s\n\n' % (args['genome'], mode))
        summary.write('Summarized benchmarks in BUSCO notation:\n\tC:%s%%[D:%s%%],F:%s%%,M:%s%%,n:%s\n\n' % (
            shrink(len(set(cc)) / totalbuscos), shrink(len(set(mcc)) / totalbuscos), shrink(fcc / totalbuscos),
            shrink((totalbuscos - (len(set(cc)) + fcc)) / totalbuscos), int(totalbuscos)))

        summary.write('Representing:\n')
        summary.write('\t%s\tComplete Single-Copy BUSCOs\n' % (len(set(cc))))
        summary.write('\t%s\tComplete Duplicated BUSCOs\n' % (len(set(mcc))))
        summary.write('\t%s\tFragmented BUSCOs\n' % fcc)
        summary.write('\t%s\tMissing BUSCOs\n' % (int(totalbuscos) - (len(set(cc)) + fcc)))

        summary.write('\t%s\tTotal BUSCO groups searched\n' % int(totalbuscos))
        summary.close()
        summary = open(join(mainout, 'full_table_%s' % assembly_name), 'w')
        summary.write('#BUSCO_group\tStatus\tScaffold\tStart\tEnd\tBitscore\tLength\n')

        temp = os.listdir(hmmer_out_path)
        done = []
        files = []
        for i in temp:
            if i.endswith(('.out', '.1', '.2', '.3')):
                files.append(i)

        for i in files:
            if i.endswith('.out'):
                name = i[:-4]
                marker = 0
            elif i.endswith(('.1', '.2', '.3')):
                name = i[:-6]
                marker = int(i[-1]) - 1
            f = open('%shmmer_output/%s' % (mainout, i))
            score = []
            hit_dic = {}
            for line in f:
                if line.startswith('#'):
                    pass
                else:
                    line = line.strip().split()
                    score.append(float(line[7]))
                    group = line[3]
                    prot = line[0]
                    tlen = int(line[2])
                    qlen = int(line[5])
                    prediction = line[0]
                    if prot not in hit_dic.keys() and float(line[7]) >= score_dic[group]:
                        hit_dic[prot] = [[int(line[15]), int(line[16]), line[7]]]
                    elif float(line[7]) >= score_dic[group]:
                        hit_dic[prot].append([int(line[15]), int(line[16]), line[7]])
            length = measuring(hit_dic)
            if hit_dic == {}:
                pass
            elif i in complete and name not in mcc:
                summary.write('%s\tComplete\t%s\t%s\t%s\t%s\t%s\n' % (
                    name, dic[group][marker][0], dic[group][marker][1], dic[group][marker][2], max(score),
                    max(length) + 1))
            elif i in complete and name in mcc:
                summary.write('%s\tDuplicated\t%s\t%s\t%s\t%s\t%s\n' % (
                    name, dic[group][marker][0], dic[group][marker][1], dic[group][marker][2], max(score),
                    max(length) + 1))
            elif i in frag and name not in cc and name not in done:
                summary.write('%s\tFragmented\t%s\t%s\t%s\t%s\t%s\n' % (
                    name, dic[group][marker][0], dic[group][marker][1], dic[group][marker][2], max(score),
                    max(length) + 1))
        summary.close()

        f = open(join(mainout, 'full_table_%s' % assembly_name), 'r')
        lista = []
        for i in f:
            i = i.strip().split()
            if i[0] not in lista:
                lista.append(i[0])
        out = open(join(mainout, 'missing_buscos_list_%s' % assembly_name),
                   'w')  #get final list of missing buscos
        f = open(join(mainout, 'full_table_%s' % assembly_name), 'a')
        for i in score_dic.keys():
            if i in lista:
                pass
            else:
                out.write('%s\n' % i)
                f.write('%s\tMissing\n' % (i))
        out.close()
        f.close()
        if not qconfig.debug:
            shutil.rmtree(mainout, ignore_errors=True)

    if mode == 'OGS':
        return totalbuscos, len(set(cc)), len(fcc)
    else:
        return totalbuscos, len(set(cc)), fcc