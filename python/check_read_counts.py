import os
import re

os.chdir("Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/WGS_Northwestern/")

# Run count_fastq_lines.sh first, then this code to get any FATSQ with discordant read counts

def stringToList(string):
    listRes = list(string.split(" "))
    return listRes

a = []
used = set()
# with open('read_counts.txt') as file:
# with open('read_counts_Ashley.txt') as file: #
# with open('read_counts_HLHS.txt') as file:
# with open('read_counts_LMNA.txt') as file:
with open('read_counts_LVNC.txt') as file:
  for line in file:
    # print(line.rstrip())
    line=line.rstrip()
    sample=line.split('/',2)[1]
    string=line.split('/',2)[2]
    file=string[:string.index(" has")]
    numbers=re.findall(r'has \d+ lines', line)
    a_tmp = stringToList(sample)
    a.append(a_tmp)
    A=re.findall(r'\d+', str(numbers))[0]
    B=re.findall(r'\d+', str(numbers))[1]
    if A != B:
      print(file + ' of ' + sample + ' = ' + 'Not OK')
    # else:
    #   print(file + ' of ' + sample + ' = ' + 'OK')
    # all.samples = set(list(k[0] for k in a))
    # list(set([*i for i in a]))
  unique_samples = list(map(list, set(map(tuple, a))))
  print(unique_samples)   
  len(unique_samples)

