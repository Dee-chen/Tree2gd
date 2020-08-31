#!/usr/bin/python3.7

import argparse
#import threading
import time
import os
import configparser
import logging
import coloredlogs
import sys
import subprocess
import multiprocessing

pwd= os.getcwd()
parser = argparse.ArgumentParser(description='Tree2GD:A pipeline for WGD')
parser.add_argument('-l', type=str, required=True, metavar='list',help='The species list.')
parser.add_argument('-i', type=str, required=True, metavar='input_dir',help='The CDS and pep dir.')
parser.add_argument('-tree', type=str, required=True, metavar='phytree.nwk',help='The phytree file.')
parser.add_argument('-t', type=int, default='1', help='Thread num.default:1')
parser.add_argument('-step', type=str, default='12345', help='which step you needs')
parser.add_argument('-log', type=str, help='log file name.if not will print on stdout')
parser.add_argument('-only_script', action='store_true', default=False,help='Only generate scripts, not run automatically.')
parser.add_argument('-o', type=str, default=os.sep.join([pwd,'output']),help='The output dir.default:'+os.sep.join([pwd,'output']))
args = parser.parse_args()

if not args.log:
        coloredlogs.install(
                fmt='%(asctime)s: %(levelname)s\t%(message)s',
                level='debug'.upper(), stream=sys.stdout
        )
else:
        print('Logs will be written to {}'.format(args.log))
        logging.basicConfig(
                    filename=args.log,
                    filemode='a',
                    format='%(asctime)s: %(levelname)s\t%(message)s',
                    datefmt='%H:%M:%S',
                    level='debug'.upper()
                    )
home_dir=os.path.dirname(os.path.abspath(__file__))
cf=configparser.ConfigParser()
#cf["software"] = {'diamond': 'diamond','makeblastdb': 'makeblastdb'}
cf.read("config.ini")
sections = cf.sections()

try:
    sp_file=open(args.l)
    sp_list=sp_file.read().splitlines()
except FileNotFoundError:
    logging.error('The species list :%s  do not exist' %(args.l))

if not os.path.exists(args.o):
    logging.info('The output dir :%s does not exist,will create it.'%(args.o))
    os.mkdir(args.o)

if '1' in args.step:
    step1out=os.sep.join([args.o,'step1.blastp'])
    diamond=cf.get("software", "diamond") if cf.has_option("software", "diamond") else "diamond"
    seqkit=cf.get("software", "seqkit") if cf.has_option("software", "seqkit") else "seqkit"
    #try:
    #    makeblastdb = cf.get("software", "makeblastdb")
    #except configparser.NoOptionError:
    #    logging.error("makeblastdb didn't set in config file")
        #sys.exit(0)

    if not os.path.exists(step1out):
        logging.info('The Step1 output dir :%s does not exist,will create it.'%(step1out))
        os.mkdir(step1out)
    if not os.path.exists(os.sep.join([step1out,'sample_sh'])):
        logging.info('The Step1 scripts dir :%s does not exist,will create it.'%(os.sep.join([step1out,'sample_sh'])))
        os.mkdir(os.sep.join([step1out,'sample_sh']))
    if not os.path.exists(os.sep.join([step1out,'allsample_blast'])):
        logging.info('The Step1 blastp result dir :%s does not exist,will create it.'%(os.sep.join([step1out,'allsample_blast'])))
        os.mkdir(os.sep.join([step1out,'allsample_blast']))

    from tree2gd.blastp import *
    if not os.path.exists(args.i):
        logging.warning('The input dir:%s does not exist,program can not run automatically,will change to only generate scripts mode.'%(args.i))
        args.only_script=True
    pep_postfix=cf.get("postfix", "pep") if cf.has_option("postfix", "pep") else ".pep"
    for sp in sp_list:
        if not os.path.exists(os.sep.join([args.i,sp+pep_postfix])):
            logging.warning('The input pep file :%s does not exist,program can not run automatically,will change to only generate scripts mode.'%(os.sep.join([args.i,sp+pep_postfix])))
            args.only_script=True
    all2all_diamond(sp_list,cf,step1out,diamond,seqkit,pep_postfix,args)
    logging.info('pep all2all diamond done.')

if '2' in args.step:
    step1out=os.sep.join([args.o,'step1.blastp'])
    step2out=os.sep.join([args.o,'step2.MCL'])
    phymcl=cf.get("software", "phymcl") if cf.has_option("software", "phymcl") else "phymcl"
    if not os.path.exists(step2out):
        logging.info('The Step2 output dir :%s does not exist,will create it.'%(step2out))
        os.mkdir(step2out)
    from tree2gd.mcl import *
    run_phymcl(step1out,phymcl,args,cf.items('phymcl'),step2out)


if '3' in args.step:
    step1out=os.sep.join([args.o,'step1.blastp'])
    step2out=os.sep.join([args.o,'step2.MCL'])
    step3out=os.sep.join([args.o,'step3.dollop'])
    dolloparsimony=cf.get("software", "dolloparsimony") if cf.has_option("software", "dolloparsimony") else home_dir+"/software/dolloparsimony"
    if not os.path.exists(step3out):
        logging.info('The Step3 output dir :%s does not exist,will create it.'%(step3out))
        os.mkdir(step3out)
    from tree2gd.dollop import *
    run_dolloparsimony(step1out,step2out,step3out,dolloparsimony,args)

if '4' in args.step:
    step1out=os.sep.join([args.o,'step1.blastp'])
    step2out=os.sep.join([args.o,'step2.MCL'])
    step4out=os.sep.join([args.o,'step4.WGD'])
    if not os.path.exists(step4out):
        logging.info('The Step4 output dir :%s does not exist,will create it.'%(step4out))
        os.mkdir(step4out)
    from tree2gd.wgd import *
    run_tree2gd(step1out,args,cf,step2out,step4out)


if '5' in args.step:
    step1out=os.sep.join([args.o,'step1.blastp'])
    step4out=os.sep.join([args.o,'step4.WGD'])
    step5out=os.sep.join([args.o,'step5.KaKs'])
    gene_pairs_idmap=[]
    if not os.path.exists(step5out):
        logging.info('The Step5 output dir :%s does not exist,will create it.'%(step5out))
        os.mkdir(step5out)
    cds_postfix=cf.get("postfix", "cds") if cf.has_option("postfix", "cds") else ".cds"
    for i in sp_list:
        src_path = os.sep.join([args.i,i+cds_postfix])
        with open(src_path, 'r') as sf, open(os.sep.join([step5out,"all_cds.fa"]), 'a+') as df:
            df.write(sf.read()+"\n")
    with open(os.sep.join([step4out,'Tree2GD_out','gd.gene_pairs.txt']),"r") as gene_pairs,open(os.sep.join([step5out,"all_gene_pairs.list"]), 'a+') as out:
        for l in gene_pairs:
            if not "#" in l:
                out.write(l.strip().split("\t")[4]+"\t"+l.strip().split("\t")[5]+"\n")
                gene_pairs_idmap.append([l.strip().split("\t")[2],l.strip().split("\t")[3],l.strip().split("\t")[4],l.strip().split("\t")[5]])
    from tree2gd.kaks import *
    run_kaks(sp_list,step1out,args,cf,step4out,step5out,gene_pairs_idmap)
