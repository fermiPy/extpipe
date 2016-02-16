import os

def collect_dirs(dirs):

    o = []
    for d in dirs:
        if not os.path.isdir(d):
            continue
        
        for subdir in os.listdir(d):

            subdir = os.path.join(d,subdir)
            if not os.path.isdir(subdir):
                continue

            o += [subdir]

        o += [d]

    return o
