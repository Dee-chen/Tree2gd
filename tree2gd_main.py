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
import shutil
from Bio import Phylo
from tree2gd.seq import read_fasta_file

def err_exit():
    sys.exit('\033[1;31;47m!!The Tree2gd program exited abnormally, please check the log file !!\033[0m')

def cmd_check(cmd):
    try:
        p=subprocess.run(cmd,check=True,stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
        logging.debug(p.stdout.decode("utf-8").strip())
        pass
    except FileNotFoundError:
        logging.error("The program %s path exists, please check the settings"%(cmd))
        err_exit()
    except subprocess.CalledProcessError:
        logging.error("The program %s not working properly, please check the settings"%(cmd))
        err_exit()
    logging.info("The program %s working properly [OK]"%(cmd))

def main():
    pwd= os.getcwd()
    parser = argparse.ArgumentParser(description='Tree2GD:A pipeline for WGD')
    parser.add_argument('-i', type=str, required=True, metavar='input_dir',help='The CDS and pep dir.')
    parser.add_argument('-tree', type=str, required=True, metavar='phytree.nwk',help='The phytree file.')
    parser.add_argument('-t', type=int,metavar='t',default='1', help='Thread num.default:1')
    parser.add_argument('-o', type=str,metavar='outputdir',default=os.sep.join([pwd,'output']),help='The output dir.default:'+os.sep.join([pwd,'output']))
    parser.add_argument('--step', type=str,metavar='123456',default='123456', help='Which steps you need.default:123456(Choose from numbers: such as \'234\')')
    parser.add_argument('--log', type=str, metavar='logfile',help='log file name,or log will print on stdout')
    parser.add_argument('--config', type=str,metavar='config.ini', help='config.ini configuration file, leave it blank to run with default parameters and the program\'s own software version.')
    parser.add_argument('--debug', action='store_true', default=False, help='The log file will contain the output of each software itself, which is convenient for finding errors (-log is required)')
    parser.add_argument('--only_script', action='store_true', default=False,help='Only generate scripts, not run automatically.')
    parser.add_argument('--cds2tree', action='store_true', default=False,help='Use cds sequence to construct gene tree.')
    args = parser.parse_args()

    if not args.log:
        coloredlogs.install(
                fmt='%(asctime)s: %(levelname)s\t%(message)s',
                evel='info'.upper(), stream=sys.stdout
            )
    else:
        if args.debug:
            l='debug'
        else:
            l='info'
        print('Logs will be written to {}'.format(args.log))
        if os.path.exists(args.log):
            os.remove(args.log)
        logging.basicConfig(
                    filename=args.log,
                    filemode='a',
                    format='%(asctime)s: %(levelname)s\t%(message)s',
                    datefmt='%H:%M:%S',
                    level=l.upper()
                    )
    home_dir=os.path.dirname(os.path.abspath(__file__))
    cf=configparser.ConfigParser()

    if not args.config:
        logging.info('The user does not provide a configuration file, it will run with the default parameters and the software version of the program.')
        cf.add_section('software')
        #cf.read(os.sep.join(home_dir,"config.ini")
    else:
        cf.read(args.config)
    if not cf.has_section('software'):
        cf.add_section('software')
    abs_out=os.path.abspath(args.o)
    sp_list=[]
    logging.info("Start checking the input tree file..")
    try:
        tree_file=Phylo.read(args.tree,'newick')
        for i in tree_file.get_terminals():
            sp_list.append(i.name)
    except FileNotFoundError:
        logging.error('Checking Newicktree file :[ERR] -- The tree file :%s  does not exist' %(args.tree))
        err_exit()
    except Phylo.NewickIO.NewickError as err:
        logging.error('Checking Newicktree file :[ERR] -- %s'%(err))
        err_exit()
    logging.info("Checking Newicktree file :[OK]")
    logging.info("Start checking the input fasta file..")
    pep_postfix=cf.get("postfix", "pep") if cf.has_option("postfix", "pep") else ".pep"
    cds_postfix=cf.get("postfix", "cds") if cf.has_option("postfix", "cds") else ".cds"
    for sp in sp_list:
        seqid_list=[]
        if not os.path.exists(os.sep.join([args.i,sp+pep_postfix])):
            logging.error('The input pep file :%s does not exist.'%(os.sep.join([args.i,sp+pep_postfix])))
            err_exit()
        for s in read_fasta_file(os.sep.join([args.i,sp+pep_postfix])):
            seqid,seq = s.name,s.seq
            if '.' in s.seq:
                logging.error("The %s sequence of %s contains the character \".\", diamond will not run normally, please replace it with the \"*\" character."%(seqid,sp+pep_postfix))
                err_exit()
            seqid_list.append(s.name)
        logging.info("Checking pep file %s:[OK]"%(sp+pep_postfix))
        if  '5' in args.step or args.cds2tree:
            logging.info("Need to perform kaks calculation or use cds for tree structure, start to check the correspondence between cds file %s and pep :"%(sp+cds_postfix))
            for s in read_fasta_file(os.sep.join([args.i,sp+cds_postfix])):
                seqid,seq = s.name,s.seq
                #if '.' in s.seq:
                #    logging.error("The %s sequence of %s contains the character \".\", please replace it with the \"*\" character."%(seqid)%(sp+pep_postfix))
                #    sys.exit('\033[1;31;47m!!The Tree2gd program exited abnormally, please check the log file !!\033[0m')
                if s.name not in seqid_list:
                    logging.error("The %s sequence id in cdsfile:%s  does not have corresponding sequence in pep file."%(s.name,sp+cds_postfix))
                    err_exit()
            logging.info("Checking cds file %s:[OK]"%(sp+cds_postfix))
    logging.info("Checking all input file :[OK]")

    logging.info("Start checking required software..")
    if '1' in args.step:
        diamond=cf.get("software", "diamond") if cf.has_option("software", "diamond") and cf.get("software", "diamond")!="" else os.sep.join([home_dir,"software", "diamond"])
        seqkit=cf.get("software", "seqkit") if cf.has_option("software", "seqkit") and cf.get("software", "seqkit")!="" else os.sep.join([home_dir,"software", "seqkit"])
        cmd_check([diamond,"help"])
        cmd_check(seqkit)
    if '2' in args.step:
        phymcl=cf.get("software", "phymcl") if cf.has_option("software", "phymcl") and cf.get("software", "phymcl")!="" else os.sep.join([home_dir,"software", "PhyloMCL"])
        cmd_check(phymcl)
    if '3' in args.step:
        dolloparsimony=cf.get("software", "dolloparsimony") if cf.has_option("software", "dolloparsimony") and cf.get("software", "dolloparsimony")!=""  else os.sep.join([home_dir,"software", "dolloparsimony"])
        cmd_check(dolloparsimony)
    if '4' in args.step:
        tree2gd=cf.get("software", "tree2gd") if cf.has_option("software", "tree2gd") and cf.get("software", "tree2gd")!="" else os.sep.join([home_dir,"software", "Tree2GD"])
        cmd_check(tree2gd)
        cf.set("software", "tree2gd",tree2gd)
        muscle=cf.get("software", "muscle") if cf.has_option("software", "muscle") and cf.get("software", "muscle")!="" else os.sep.join([home_dir,"software", "muscle"])
        cmd_check([muscle,"-version"])
        cf.set("software", "muscle",muscle)
        iqtree=cf.get("software", "iqtree") if cf.has_option("software", "iqtree") and cf.get("software", "iqtree")!="" else os.sep.join([home_dir,"software", "iqtree"])
        cmd_check([iqtree,"-h"])
        cf.set("software","iqtree",iqtree)
        if args.cds2tree:
            Epal2nal=cf.get("software", "Epal2nal") if cf.has_option("software", "Epal2nal") and cf.get("software", "Epal2nal")!="" else os.sep.join([home_dir,"software", "Epal2nal.pl"])
            cmd_check(Epal2nal)
            cf.set("software","Epal2nal",Epal2nal)
    if '5' in args.step:
        KaKs_Calculator=cf.get("software", "KaKs_Calculator") if cf.has_option("software", "KaKs_Calculator") and cf.get("software", "KaKs_Calculator")!="" else os.sep.join([home_dir,"software", "KaKs_Calculator"])
        cmd_check(KaKs_Calculator)
        cf.set("software", "KaKs_Calculator",KaKs_Calculator)
        ParaAT=cf.get("software", "ParaAT") if cf.has_option("software", "ParaAT") and cf.get("software", "ParaAT")!="" else os.sep.join([home_dir,"software", "ParaAT.pl"])
        cmd_check(ParaAT)
        cf.set("software", "ParaAT",ParaAT)
        ParaAT_aligncmd=cf.get("software", "ParaAT_aligncmd") if cf.has_option("software", "ParaAT_aligncmd") and cf.get("software", "ParaAT_aligncmd")!="" else os.sep.join([home_dir,"software", "muscle"])
        cf.set("software", "ParaAT_aligncmd",ParaAT_aligncmd)
        Epal2nal=cf.get("software", "Epal2nal") if cf.has_option("software", "Epal2nal") and cf.get("software", "Epal2nal")!="" else os.sep.join([home_dir,"software", "Epal2nal.pl"])
        cmd_check(Epal2nal)
        cf.set("software", "Epal2nal",Epal2nal)


    if not os.path.exists(abs_out):
        logging.info('The output dir :%s does not exist,will create it.'%(abs_out))
        os.mkdir(abs_out)

    if '1' in args.step:
        step1out=os.sep.join([abs_out,'step1.blastp'])
        if not os.path.exists(step1out):
            logging.info('The Step1 output dir :%s does not exist,will create it.'%(step1out))
        else:
            logging.warning('The Step1 output dir :%s has exist,will backup it as %s_backup and create new.'%(step1out,step1out))
            os.rename(step1out,step1out+"_backup")
        os.mkdir(step1out)
        if not os.path.exists(os.sep.join([step1out,'sample_sh'])):
            logging.info('The Step1 scripts dir :%s does not exist,will create it.'%(os.sep.join([step1out,'sample_sh'])))
            os.mkdir(os.sep.join([step1out,'sample_sh']))
        if not os.path.exists(os.sep.join([step1out,'allsample_blast'])):
            logging.info('The Step1 blastp result dir :%s does not exist,will create it.'%(os.sep.join([step1out,'allsample_blast'])))
            os.mkdir(os.sep.join([step1out,'allsample_blast']))

        from tree2gd.blastp import all2all_diamond
        all2all_diamond(sp_list,cf,step1out,diamond,seqkit,pep_postfix,args)
        logging.info('pep all2all diamond done.')

    if '2' in args.step:
        step1out=os.sep.join([abs_out,'step1.blastp'])
        step2out=os.sep.join([abs_out,'step2.MCL'])
        if not os.path.exists(step2out):
            logging.info('The Step2 output dir :%s does not exist,will create it.'%(step2out))
        else:
            logging.warning('The Step2 output dir :%s has exist,will backup it as %s_backup and create new.'%(step2out,step2out))
            os.rename(step2out,step2out+"_backup")
        os.mkdir(step2out)
        if not cf.has_section('phymcl'):
            cf.add_section('phymcl')
        from tree2gd.mcl import run_phymcl
        run_phymcl(step1out,phymcl,args,cf.items('phymcl'),step2out)


    if '3' in args.step:
        step1out=os.sep.join([abs_out,'step1.blastp'])
        step2out=os.sep.join([abs_out,'step2.MCL'])
        step3out=os.sep.join([abs_out,'step3.dollop'])
        if not os.path.exists(step3out):
            logging.info('The Step3 output dir :%s does not exist,will create it.'%(step3out))
        else:
            logging.warning('The Step3 output dir :%s has exist,will backup it as %s_backup and create new.'%(step3out,step3out))
            os.rename(step3out,step3out+"_backup")
        os.mkdir(step3out)
        from tree2gd.dollop import run_dolloparsimony
        run_dolloparsimony(step1out,step2out,step3out,dolloparsimony,args)

    if '4' in args.step:
        step1out=os.sep.join([abs_out,'step1.blastp'])
        step2out=os.sep.join([abs_out,'step2.MCL'])
        step4out=os.sep.join([abs_out,'step4.WGD'])

        if not os.path.exists(step4out):
            logging.info('The Step4 output dir :%s does not exist,will create it.'%(step4out))
        else:
            logging.warning('The Step4 output dir :%s has exist,will backup it as %s_backup and create new.'%(step4out,step4out))
            os.rename(step4out,step4out+"_backup")
        os.mkdir(step4out)
        if args.cds2tree:
            for i in sp_list:
                src_path = os.sep.join([args.i,i+cds_postfix])
                with open(src_path, 'r') as sf, open(os.sep.join([step1out,"all_cds.fa"]), 'a+') as df:
                    lines=sf.readlines()
                    for line in lines:
                        df.write(line.strip().split(" ")[0]+"\n")
        from tree2gd.wgd import run_tree2gd
        run_tree2gd(step1out,args,cf,step2out,step4out)
        shutil.rmtree(os.sep.join([step4out,'all_fa']))
        shutil.rmtree(os.sep.join([step4out,'muscle_out']))

    if '5' in args.step:
        step1out=os.sep.join([abs_out,'step1.blastp'])
        step4out=os.sep.join([abs_out,'step4.WGD'])
        step5out=os.sep.join([abs_out,'step5.KaKs'])
        gene_pairs_idmap=[]
        if not os.path.exists(step5out):
            logging.info('The Step5 output dir :%s does not exist,will create it.'%(step5out))
        else:
            logging.warning('The Step5 output dir :%s has exist,will backup it as %s_backup and create new.'%(step5out,step5out))
            os.rename(step5out,step5out+"_backup")
        os.mkdir(step5out)
        cds_postfix=cf.get("postfix", "cds") if cf.has_option("postfix", "cds") else ".cds"
        with open(os.sep.join([step1out,'all_sample.pep.faa']), 'r') as pep,open(os.sep.join([step5out,"all_pep.fa"]), 'a+') as pep_out:
            lines=pep.readlines()
            for line in lines:
                pep_out.write(line.strip().split(" ")[0]+"\n")
        for i in sp_list:
            src_path = os.sep.join([args.i,i+cds_postfix])
            with open(src_path, 'r') as sf, open(os.sep.join([step5out,"all_cds.fa"]), 'a+') as df:
                lines=sf.readlines()
                for line in lines:
                    df.write(line.strip().split(" ")[0]+"\n")
        with open(os.sep.join([step4out,'Tree2GD_out','gd.gene_pairs.txt']),"r") as gene_pairs,open(os.sep.join([step5out,"all_gene_pairs.list"]), 'a+') as out:
            for l in gene_pairs:
                if not "#" in l:
                    out.write(l.strip().split("\t")[4]+"\t"+l.strip().split("\t")[5]+"\n")
                    gene_pairs_idmap.append([l.strip().split("\t")[2],l.strip().split("\t")[3],l.strip().split("\t")[4],l.strip().split("\t")[5]])
        from tree2gd.kaks import run_kaks
        run_kaks(sp_list,step1out,args,cf,step4out,step5out,gene_pairs_idmap)
        shutil.rmtree(os.sep.join([step5out,"all_gene_pairs_kaks"]))

    if '6' in args.step:
        step2out=os.sep.join([abs_out,'step2.MCL'])
        step3out=os.sep.join([abs_out,'step3.dollop'])
        step4out=os.sep.join([abs_out,'step4.WGD'])
        step5out=os.sep.join([abs_out,'step5.KaKs'])
        step6out=os.sep.join([abs_out,'step6.plot_summary'])
        if not os.path.exists(step6out):
            logging.info('The Step6 output dir :%s does not exist,will create it.'%(step6out))
        else:
            logging.warning('The Step6 output dir :%s has exist,will backup it as %s_backup and create new.'%(step6out,step6out))
            os.rename(step6out,step6out+"_backup")
        os.mkdir(step6out)
        shutil.copyfile(os.sep.join([step2out,"allmcl.out.OGs.tabular"]),os.sep.join([step6out,"allmcl.out.OGs.tabular"]))
        shutil.copyfile(os.sep.join([step3out,"dollop.out.scv"]),os.sep.join([step6out,"dollop.out.scv"]))
        os.mkdir(os.sep.join([step6out,"Tree2GD_out"]))
        shutil.copyfile(os.sep.join([step4out,"Tree2GD_out","gd.gene_pairs.txt"]),os.sep.join([step6out,"Tree2GD_out","gd.gene_pairs.txt"]))
        shutil.copyfile(os.sep.join([step4out,"Tree2GD_out","GDtype_stat.txt"]),os.sep.join([step6out,"Tree2GD_out","GDtype_stat.txt"]))
        shutil.copyfile(os.sep.join([step4out,"Tree2GD_out","Phtree.nwk"]),os.sep.join([step6out,"Tree2GD_out","Phtree.nwk"]))
        shutil.copyfile(os.sep.join([step4out,"Tree2GD_out","summarytable.txt"]),os.sep.join([step6out,"Tree2GD_out","summarytable.txt"]))
        shutil.copytree(os.sep.join([step5out,"sp_kaks_out"]),os.sep.join([step6out,"sp_kaks_out"]))
        shutil.copyfile(os.sep.join([home_dir,"software","Tree2GD_draw.R"]),os.sep.join([step6out,"Tree2GD_draw.R"]))
        from tree2gd.plot import run_plot
        run_plot(step6out,args,cf)
        os.remove(os.sep.join([step6out,"Rplots.pdf"]))


def test():
    home_dir=os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(description='test cmd for Tree2GD')
    parser.add_argument('-t', type=int, default='1', help='Thread num.default:1')
    parser.add_argument('--config',default=os.sep.join([home_dir,"example_data","config.ini"]),type=str, help='config.ini configuration file, leave it blank to run with default parameters and the program\'s own software version.')
    args = parser.parse_args()
    subprocess.run("Tree2gd -i "+os.sep.join([home_dir,"example_data"])+" -tree "+os.sep.join([home_dir,"example_data","example.tree"])+" --config "+args.config+" -o ./Tree2gd_test_out -t "+str(args.t),shell=True)

if __name__ == '__main__':
    main()
