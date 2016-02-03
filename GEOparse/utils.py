import os
from errno import EEXIST
from sys import stderr, stdout

def mkdir_p(path_to_dir):
    try:
        os.makedirs(path_to_dir)
    except OSError as e: # Python >2.5
        if e.errno == EEXIST and os.path.isdir(path_to_dir):
            stderr.write("Directory %s already exists. Skipping.\n" % path_to_dir)
        else:
            raise e
