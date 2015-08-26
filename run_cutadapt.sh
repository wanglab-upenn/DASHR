#!/bin/bash
# export PYTHONPATH=$PYTHONPATH:$HOME/lib/python/:$HOME/lib64/python/:$HOME/lib/python2.6/site-packages/
export PYTHONPATH=$PYTHONPATH:$HOME/lib/python/:$HOME/lib64/python/:/home/pry/lib/python/
#/genomics/share/python2.6.5/bin/python $HOME/bin/cutadapt $@
#python2.6 $HOME/bin/cutadapt $@
$HOME/bin/cutadapt-1.8.1/bin/cutadapt $@


