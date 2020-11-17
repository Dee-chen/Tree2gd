#!/usr/bin/python3.7

import os
import subprocess
import logging
from multiprocessing import Pool
import sys
import configparser

def sub_sh(sh):
    logging.info("Start runing step1 blastp: %s.."%sh)
    stdout=subprocess.run(['sh',sh],stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    logging.debug('{} stdout:\n'.format(sh) +stdout.stdout.decode('utf-8'))
    logging.debug('{} stderr:\n'.format(sh) +stdout.stderr.decode('utf-8'))
    logging.info("step1 blastp: %s has done."%sh)


def all2all_diamond(sp_list,config,out,diamond,seqkit,pep_postfix,args):
    logging.info("Writing step1 script files...")
    sp_num=len(sp_list)
    op=''
    sh_list=[]
    p=config.get("diamond","-p") if config.has_option("diamond","-p") and config.get("diamond","-p")!="" else args.t
    if int(p)>args.t:
        logging.warning("The number of diamond threads %s is greater than the total program thread %s, it will be set to %s."%(str(p),str(args.t),str(args.t)))
        p=args.t
    if not config.has_section('diamond'):
        config.add_section('diamond')
    config.set("diamond","-p",str(p))
    options=config.items("diamond")
    summary_sh = open(os.sep.join([out,'summary.sh']),"w")
    step1_sh = open(os.sep.join([out,'step1.sh']),"w")
    for k,v in options:
        op=' '.join((op,' '.join((k,v))))
    for i in sp_list:
        SH = open(os.sep.join([out,'sample_sh',i+".sh"]),"w")
        SH.write(diamond + ' makedb -d '+os.sep.join([args.i,i+pep_postfix])+' --in '+os.sep.join([args.i,i+pep_postfix])+"\n")
        for j in sp_list:
            SH.write(diamond+' blastp -f 6 --more-sensitive '+op+' -d '+os.sep.join([args.i,i+pep_postfix])+' -q '+os.sep.join([args.i,j+pep_postfix])+' -o '+os.sep.join([out,'allsample_blast',i+"-"+j+".blast\n"]))
            summary_sh.write("cat "+os.sep.join([out,'allsample_blast',i+"-"+j+".blast"])+" >> "+os.sep.join([out,'all_blastp.out'])+"\n")
        summary_sh.writelines(["cat "+os.sep.join([args.i,i+pep_postfix])+ " >> " +os.sep.join([out,'all_sample.pep.faa'])+"\n",
                        seqkit+" fx2tab -l " + os.sep.join([args.i,i+pep_postfix])+"|awk \'{print $1\"\\t\"$NF}\' >> "+os.sep.join([out,'all_sample.faa.length'])+"\n",
                        seqkit+" fx2tab -l " + os.sep.join([args.i,i+pep_postfix])+"|awk \'{print $1\"\\t"+i+"\"}\' >> "+os.sep.join([out,'all_sample2fa.list'])+"\n"])
        SH.close()
        sh_list.append(os.sep.join([out,'sample_sh',i+".sh"]))
    summary_sh.close()

    if not args.only_script:
        try:
            subprocess.run(diamond,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        except FileNotFoundError:
            logging.error("diamond can not run.")
            sys.exit(0)
        try:
            subprocess.run(seqkit,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        except FileNotFoundError:
            logging.error("seqkit can not run.")
            sys.exit(0)
        sh_pool=Pool(args.t//int(p))
        sh_pool.map(sub_sh,sh_list)
        sh_pool.close()
        sh_pool.join()
        logging.info("ALL step1 blastp has done.")
        logging.info("Start sorting results..")
        stdout=subprocess.run(['sh',os.sep.join([out,'summary.sh'])],stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    step1_sh.writelines("\n".join(sh_list))
    step1_sh.write(os.sep.join([out,'summary.sh']))
    step1_sh.close()
