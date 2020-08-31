#!/usr/bin/python3.7

import os
import subprocess
import logging
import sys

def run_dolloparsimony(step1out,step2out,step3out,dolloparsimony,args):
    logging.info("Writing step3 script files...")
    step3_sh=open(os.sep.join([step3out,'step3.sh']),"w")
    step3_sh.write(' '.join([dolloparsimony,args.tree,os.sep.join([step1out,'all_sample2fa.list']),os.sep.join([step2out,'allmcl.out.OGs.group']),os.sep.join([step3out,'dollop.out'])+"\n"]))
    step3_sh.write("less "+os.sep.join([step3out,'dollop.out','summary.tree'])+"|cut -d '[' -f2|cut -d ']' -f1  >"+os.sep.join([step3out,'dollop.out.scv']))
    step3_sh.close()

    if not args.only_script:
        logging.info("Start runing step3 dollop: %s.."%(os.sep.join([step3out,'step3.sh'])))
        stdout=subprocess.run(['sh',os.sep.join([step3out,'step3.sh'])],stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logging.debug('{} stdout:\n'.format('step2.sh') +stdout.stdout.decode('utf-8'))
        logging.debug('{} stderr:\n'.format('step2.sh') +stdout.stderr.decode('utf-8'))
        logging.info("step3 dollop has done.")
