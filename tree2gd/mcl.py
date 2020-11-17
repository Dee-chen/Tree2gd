#!/usr/bin/python3.7
import os
import subprocess
import logging
import sys

def run_phymcl(step1out,phymcl,args,options,step2out):
    logging.info("Writing step2 script files...")
    op=''
    for k,v in options:
        op=' '.join((op,' '.join((k,v))))
    step2_sh=open(os.sep.join([step2out,'step2.sh']),"w")
    step2_sh.write(phymcl+ " -in " +os.sep.join([step1out,'all_blastp.out'])+" -threads "+ str(args.t) + " -length " + os.sep.join([step1out,'all_sample.faa.length']) +
                        " -species " + os.sep.join([step1out,'all_sample2fa.list']) + " -tree " + os.sep.join([step2out,"phymcl.input.tree"]) + " -out " + os.sep.join([step2out,'allmcl.out'])+"\n")
    step2_sh.close()
    if not args.only_script:
        try:
            subprocess.run(phymcl,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        except FileNotFoundError:
            logging.error("phymcl can not run.")
            sys.exit(0)
        logging.info("Start runing step2 MCL: %s.."%(os.sep.join([step2out,'step2.sh'])))
        stdout=subprocess.run(['sh',os.sep.join([step2out,'step2.sh'])],stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logging.debug('{} stdout:\n'.format('step2.sh') +stdout.stdout.decode('utf-8'))
        logging.debug('{} stderr:\n'.format('step2.sh') +stdout.stderr.decode('utf-8'))
        logging.info("step2 mcl has done.")
