#!/usr/bin/python3.7

import os
import subprocess
import logging
from multiprocessing import Pool
import sys
import configparser

def sub_sh(sh):
    logging.info("Start runing step6 synteny : %s.."%sh)
    stdout=subprocess.run(['sh',sh],stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    logging.debug('{} stdout:\n'.format(sh) +stdout.stdout.decode('utf-8'))
    logging.debug('{} stderr:\n'.format(sh) +stdout.stderr.decode('utf-8'))
    logging.info("step6 synteny : %s has done."%sh)


def synteny(sp_list,config,out,diamond,seqkit,pep_postfix,args):
    logging.info("Writing synteny script files...")
    sp_num=len(sp_list)
    op=''
    sh_list=[]
    step6_sh = open(os.sep.join([out,'step6.sh']),"w")
    for k,v in options:
        op=' '.join((op,' '.join((k,v))))
    for i in sp_list:
        SH = open(os.sep.join([out,'sample_sh',i+".sh"]),"w")
        SH.write('python -m jcvi.compara.catalog ortholog --no_strip_names --cscore=.99 '+os.sep.join([args.i,i])+" "+os.sep.join([args.i,i])+"\n")
        SH.close()
        sh_list.append(os.sep.join([out,'sample_sh',i+".sh"]))

    if not args.only_script:
        try:
            subprocess.run(jcvi,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        except FileNotFoundError:
            logging.error("jcvi can not run.")
            sys.exit(0)

        sh_pool=Pool(int(p))
        sh_pool.map(sub_sh,sh_list)
        sh_pool.close()
        sh_pool.join()
        logging.info("ALL synteny  has done.")
        logging.info("Start sorting results..")
        stdout=subprocess.run(['sh',os.sep.join([out,'summary.sh'])],stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    step1_sh.writelines("\n".join(sh_list))
    step1_sh.write(os.sep.join([out,'summary.sh']))
    step1_sh.close()
