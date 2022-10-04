import os
import re

os.chdir("Z:/ResearchHome/Groups/sapkogrp/projects/CAB/common/WGS_Northwestern/")


with open('read_counts.txt') as file:
  for line in file:
    # print(line.rstrip())
    line=line.rstrip()
    sample=line.split('/',2)[1]
    numbers=re.findall(r'has \d+ lines', line)
    A=re.findall(r'\d+', str(numbers))[0]
    B=re.findall(r'\d+', str(numbers))[1]
    if A != B:
      print(sample + ' = ' + 'Not OK')



