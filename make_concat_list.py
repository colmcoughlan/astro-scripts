import os
import sys
import fnmatch

if(len(sys.argv)<4):
        print("\tError: Takes at least 3 arguments.")
        print("\tUseage: make_concat_list <stem> <ending> <directory(s)>")
        sys.exit()
else:
	stem = str(sys.argv[1])
	ending = str(sys.argv[2])
	ndirs = len(sys.argv)-3
	inputdirs=[[] for _ in range(ndirs)]

	for i in range(ndirs):
        	inputdirs[i] = str(sys.argv[i+3])
        print '\tReading from:', inputdirs

concatlist=[]
for inputdir in inputdirs:
	for file in os.listdir(inputdir):
	    if fnmatch.fnmatch(file, stem+'*'+ending):
		concatlist.append(inputdir+'/'+file)

with open('concat_list.txt', 'w') as f:
        f.write(str(concatlist)+'\n')

print('\t'+str(len(concatlist))+' files identified')
