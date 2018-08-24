#!/usr/bin/env python2.7

import sys
import glob
import subprocess 
ps = []
for file in glob.glob("tests_protoDC2_lc/*.py"):
    print "\n\n=================\n"
    print "./"+file+" "+sys.argv[1]+" --nonX"
    p1 = subprocess.Popen("./"+file+" "+sys.argv[1]+" --nonX", shell=True)
    p1.wait()
