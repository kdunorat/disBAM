import os
cmd = 'grep ^FREADS | cut -f 2- stats.txt > teste.txt'
os.cmd(cmd)