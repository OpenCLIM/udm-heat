import subprocess
from os.path import isfile, join, isdir
from os import listdir, remove, mkdir

data_path = '/data'

if isdir(join(data_path,'temp')) is False:
    mkdir(join(data_path,'temp'))
else:
    files = [f for f in listdir(join(data_path,'temp')) if isfile(join(data_path,'temp',f))]
    for file in files:
        remove(join(data_path,'temp',file))


if isdir(join(data_path,'outputs')) is False:
    mkdir(join(data_path,'outputs'))

# empty output dir
files = [f for f in listdir(join(data_path,'outputs')) if isfile(join(data_path,'outputs',f))]
for file in files:
    remove(join(data_path,'outputs',file))

#get list of input files to loop through

files = [f for f in listdir(join(data_path,'inputs','layers')) if isfile(join(data_path,'inputs','layers',f))]
print('files to loop through')

for file in files:
    print(file)
    file_name = file.split('.')[0]
    subprocess.run(["gdalwarp", "-te", "0", "12000", "660000", "1212000", "-tr", "12000", "12000", "-r", "sum", join(data_path,'inputs', 'layers', file), join(data_path, 'temp','%s-%s.vrt'%(file_name,'temp') )])

    subprocess.run(["gdal_translate", "-of", "AAIGrid",join(data_path,"temp","%s-%s.vrt"%(file_name,'temp')), join(data_path,"outputs","gb-2017-detached-dph-12km-sum.asc")]) 
