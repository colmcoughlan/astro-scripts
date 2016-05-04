import os
import sys

if(len(sys.argv)!=2):
        print("\tError: Takes 1 argument.")
        print("\tUseage: make_concat_list <directory>")
        sys.exit()
else:
        inputdir = str(sys.argv[1])
        print '\tReading from:', inputdir

concatlist=[]
for file in os.listdir(inputdir):
    if file.endswith((".MS",".ms",'.ndppp_prep_cal')):
        concatlist.append(file)

with open('concat_list.txt', 'w') as f:
        f.write(str(concatlist)+'\n')

print('\t'+str(len(concatlist))+' files identified')
