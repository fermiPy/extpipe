import os
import sys
import glob

prefix = sys.argv[1]

files = sorted(glob.glob(os.path.join(prefix,'*')))
files = [f for f in files if os.path.isdir(f)]


files_per_job = 100
njob = len(files)/files_per_job+1

for i in range(njob):
    job_files = ' '.join(files[i*files_per_job:(i+1)*files_per_job])

    outfile = os.path.join(prefix,'table_%s_%02i.fits'%(prefix,i))
    outlog = os.path.splitext(outfile)[0] + '.log'
    
    cmd = 'bsub -W 500 -R rhel60 '
    cmd += '-oo %s haloanalysis-aggregate %s --output=%s'%(outlog, job_files,
                                                           outfile)

    #print(cmd)
    os.system(cmd)    
    #print(i,len(job_files))

sys.exit(0)
                  
for i in range(24):
    
    cmd = 'bsub -W 500 -R rhel60 '
    cmd += '-oo table_%s_%02i.log haloanalysis-aggregate %s/3fgl_j%02i* --output=table_%s_%02i.fits'%(prefix,
                                                                                                      i,prefix,
                                                                                                      i,prefix,i)

    print(cmd)
    os.system(cmd)
    
#bsub -W 500 -R rhel60 -oo table_$1_00_04.log haloanalysis-aggregate $1/3fgl_j0[0-4]* --output=table_$1_00_04.fits 
#bsub -W 500 -R rhel60 -oo table_$1_05_09.log haloanalysis-aggregate $1/3fgl_j0[5-9]* --output=table_$1_05_09.fits 
#bsub -W 500 -R rhel60 -oo table_$1_10_14.log haloanalysis-aggregate $1/3fgl_j1[0-4]* --output=table_$1_10_14.fits 
#bsub -W 500 -R rhel60 -oo table_$1_15_19.log haloanalysis-aggregate $1/3fgl_j1[5-9]* --output=table_$1_15_19.fits 
#bsub -W 500 -R rhel60 -oo table_$1_20_24.log haloanalysis-aggregate $1/3fgl_j2[0-4]* --output=table_$1_20_24.fits
