#!/usr/bin/python

import sys
import os
import shutil
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

if len(sys.argv) != 3:
    print("Usage: " + sys.argv[0] + "  <.ini file> <output dir>")
    sys.exit()
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
iniFileName = sys.argv[1]
iniFile = open(iniFileName, 'r')
line = iniFile.readline()
#print(line.split()[0])
assert(line.split()[0] == 'collections')
collection_num = int(line.split()[1])
line = iniFile.readline()
line = iniFile.readline()
assert (line.strip() == 'collection_names')
collection_names = []
for i in range(collection_num):
    line = iniFile.readline().strip()
    collection_names.append(line)




line = iniFile.readline()
line = iniFile.readline()
assert(line.split()[0] == 'datasets')
datasets_num = int(line.split()[1])
colors = ["red", "blue", "green", "yellow", "black", "grey"]
assert(collection_num <= colors.__len__())

assemblies = []

line = iniFile.readline()
line = iniFile.readline()
assert (line.strip() == 'dataset_names')
dataset_names = []
for i in range(datasets_num):
    line = iniFile.readline().strip()
    dataset_names.append(line)



for i in range(collection_num):
    line = iniFile.readline()
    line = iniFile.readline()
#    print(line)
    assert(line.split()[0] == 'collection')
    assert(int(line.split()[1]) == i + 1)
    assemblies.append([])
    for j in range (datasets_num):
        line = iniFile.readline().strip();
        assemblies[i].append(line);

#print(collection_num, datasets_num)
#for i in range(collection_num):
#    for j in range (datasets_num):
#        print(assemblies[i][j])
#    print('\n')

line = iniFile.readline()
line = iniFile.readline()
assert(line.split()[0] == 'creatures')
creatures_num = int(line.split()[1])
assert(creatures_num <= datasets_num)

creatures_index = []
for i in range(datasets_num):
    creatures_index.append([-1])
for i in range (creatures_num):
    line = iniFile.readline()
    datasets = line.split()
    for j in datasets:
        creatures_index[int(j) - 1] = i

for i in range(datasets_num):
    assert(creatures_index[i] in range ( creatures_num ))


line = iniFile.readline()
line = iniFile.readline()
assert(line.split()[0] == 'min_contig_length')
min_len = int(line.split()[1])

line = iniFile.readline()
line = iniFile.readline()
assert(line.split()[0] == 'reference')
reference_num = int(line.split()[1])
assert(reference_num == creatures_num)
references = []
for i in range (reference_num):
    line = iniFile.readline()
    references.append(line.strip());


line = iniFile.readline()
line = iniFile.readline()
assert(line.split()[0] == 'genes')
genes_num = int(line.split()[1])
assert(genes_num == creatures_num or genes_num == 0)
genes = []
for i in range (genes_num):
    line = iniFile.readline()
    genes.append(line.strip());


line = iniFile.readline()
line = iniFile.readline()
assert(line.split()[0] == 'operons')
operons_num = int(line.split()[1])
assert(operons_num == creatures_num or operons_num == 0)
operons = []
for i in range (operons_num):
    line = iniFile.readline()
    operons.append(line.strip());

line = iniFile.readline()
line = iniFile.readline()
assert(line.split()[0] == 'metrics')
metrics_num = int(line.split()[1])

metrics = []
for i in range (metrics_num):
    line = iniFile.readline()
    metrics.append(line.strip());



folders = []
for i in range (datasets_num):
    folders.append('tmp_res'+str(i))

    quast_line = './quast.py -R ' + references[creatures_index[i]] + ' -M ' + str(min_len)
    if (operons_num != 0):
        quast_line += ' -O ' + operons[creatures_index[i]]
    if (genes_num != 0):
        quast_line += ' -G ' + genes[creatures_index[i]]

    for j in range (collection_num):
        shutil.rmtree(folders[i] + '_' + str(j), True)
        run_line = quast_line + ' -o ' + folders[i] + '_' + str(j) + (' ' + assemblies[j][i])
        print(quast_line)
        os.system(run_line)

output_dir = sys.argv[2];
os.system('mkdir -p ' + output_dir)
for metric in metrics:
    results = []
    print 'drawing... '  + metric
    ymax = 0;
    for i in range (datasets_num):
        #    print (i)
        #    print(columns)
        results.append([])
        for j in range(collection_num):
            resultsFileName = folders[i] + '_' + str(j) + "/transposed_report.tsv"
            resultsFile = open(resultsFileName, 'r')
            columns = map(lambda s: s.strip(), resultsFile.readline().split('\t'))
            values = map(lambda s: s.strip(), resultsFile.readline().split('\t'))
            #        print(values)
#            print (values)
            if values[columns.index(metric)].split()[0] == 'None' :
                metr_res = 0
            else:
                metr_res = float(values[columns.index(metric)].split()[0])
            ymax = max(ymax, metr_res)
            results[i].append(metr_res);
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    #ax = fig.add_subplot(1,1,1, aspect='equal')
    plt.xticks(range(1, datasets_num + 1) , dataset_names,  size='small')
    title = metric
    #for j in range(collection_num):
    #    title += colors[j] + "  for " + str(j) + " \n"
    plt.title(title)

    #ax.set_xticks(range(1,datasets_num + 1))
    #ax.set_xticklabels(assemblies[0])
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height*1.0])

    for j in range(collection_num):
        to_plot = []
        arr = range(1, datasets_num + 1)
        for i in range(datasets_num):
            to_plot.append(results[i][j])
            arr[i] += 0.07 *  (j - (collection_num-1) * 0.5)
        ax.plot( arr, to_plot, 'ro', color=colors[j])
    plt.xlim([0,datasets_num + 1])
    plt.ylim([0, math.ceil(ymax *  1.05)])
    #    ax.plot(range(1, datasets_num + 1), to_plot, 'ro', color=colors[j])
    legend = []
    for j in range(collection_num):
        legend.append(collection_names[j])



    ax.legend(legend, loc = 'center left', bbox_to_anchor = (1.0, 0.5))
    #plt.legend(legend, font='small', loc=(1.1,0.5))


    F = plt.gcf()


    DPI = F.get_dpi()
#    print "DPI:", DPI
    DefaultSize = F.get_size_inches()
    F.set_size_inches(2*DefaultSize)
    plt.savefig(output_dir +'/' + metric+'.jpg')

