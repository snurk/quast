#!/usr/bin/python

import sys
import os
import shutil
import re
import getopt
import subprocess
import quast

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
iniFileName = 'comparison.ini'
iniFile = open(iniFileName, 'r')
line = iniFile.readline()
#print(line.split()[0])
assert(line.split()[0] == 'collections')
collection_num = int(line.split()[1])
line = iniFile.readline()
line = iniFile.readline()
assert(line.split()[0] == 'datasets')
datasets_num = int(line.split()[1])

assemblies = []

for i in range(collection_num):
    line = iniFile.readline()
    line = iniFile.readline()
    print(line)
    assert(line.split()[0] == 'collection')
    assert(int(line.split()[1]) == i + 1)
    assemblies.append([])
    for j in range (datasets_num):
        line = iniFile.readline().strip();
        assemblies[i].append(line);

print(collection_num, datasets_num)
for i in range(collection_num):
    for j in range (datasets_num):
        print(assemblies[i][j])
    print('\n')

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
assert(line.split()[0] == 'reference')
reference_num = int(line.split()[1])
assert(reference_num == creatures_num)
references = []
for i in range (reference_num):
    line = iniFile.readline()
    references.append(line.strip());

for i in range (datasets_num):
    quast_line = './quast.py -R ' + references[creatures_index[i]] + ' -o tmp_res' + str(i)
    for j in range (collection_num):
        quast_line += (' ' + assemblies[j][i])
    print(quast_line)
    os.system(quast_line)