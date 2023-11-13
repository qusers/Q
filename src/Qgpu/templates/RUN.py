import multiprocessing
import shutil
import os

def worker(i):
    os.chdir('298.0/{:02d}'.format(i))
$runfiles
    shutil.move('out/qfep.out', '../../results/qfep-{:02d}.out'.format(i))
    os.chdir('../../')
    os.system('tar -cvf raw_data-{:02d}.tar.gz 298.0/{:02d}'.format(i,i))
    shutil.rmtree('298.0/{:02d}'.format(i))
    return

if __name__ == '__main__':
    jobs = []
    for i in range($processors):
        p = multiprocessing.Process(target=worker,args=(i,))
        jobs.append(p)
        p.start()