import os
import subprocess as sp

username = 'mcanyelles'

def read_IDs(username=username):
    command = "qstat -u %s" % username

    lines = sp.getoutput(command).split('\n')
    
    if lines != []:
        IDs = []
        for i in lines[5:]:
            if i[86:87] == 'R' and 'scan' in i:
                IDs.append(i[:6])

        return IDs
    
    else :
        return None

def read_paths(ID):
    command = 'qstat -f %s.kirk.uab.es' % ID
    lines = sp.getoutput(command).split('\n')
    
    for line in range(len(lines)):
        if 'Job_Name' in lines[line]:
            job_name = lines[line][15:]
        if 'init_work_dir' in lines[line]:
            work_dir = lines[line][20:] + ((lines[line + 1].replace(' ', '')).replace('\t', '')).replace('\n', '')

    return job_name, work_dir

def print_PES(ID, job_name, work_dir):
    
    print(ID, ': ', job_name)
    PES = work_dir + '/PES.plt'
    with open(PES) as f:
        print(f.read())
    

def main(username=username):
    IDs = read_IDs()
    if IDs != None:
        for ID in IDs:
            job_name, work_dir = read_paths(ID)
            print_PES(ID, job_name, work_dir)

if __name__ == '__main__':
    main()
