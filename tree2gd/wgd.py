#!/usr/bin/python3.7
import os
import subprocess
import logging
import sys
from tree2gd.seq import read_fasta_file
#from seq import read_fasta_file
import configparser
import subprocess
from multiprocessing import Pool


def mcl2fasta(minimal_taxa,outdir,step2out,step1out,args):
        logging.info("Reading mcl output file")
        clusterDICT = {} #key is seqID, value is clusterID
        gene_idmap={}
        fa_list=[]
        count = 0
        if outdir[-1] != "/": outdir += "/"
        with open(os.sep.join([step1out,'all_sample2fa.list']),"rU") as idmap_file:
            for l in idmap_file:
                sp = l.strip().split('\t')
                gene_idmap[sp[0]]=sp[1]
        logging.info("gene_idmap read done.")
        with open(os.sep.join([step2out,'allmcl.out.OGs.group']),"rU") as infile:
                for line in infile:
                        spls=[]
                        if len(line) < 3: continue #ignore empty lines
                        genes = line.strip().split('\t')
                        for gene in genes:
                            spls.append(gene_idmap[gene])
                        if len(set(spls)) >= minimal_taxa:
                                count += 1
                                clusterID = str(count)
                                for seqID in genes:
                                        clusterDICT[seqID] = clusterID
        logging.info("clusters with at least "+str(minimal_taxa)+" taxa read")

        logging.info("Reading the fasta file "+os.sep.join([step1out,'all_sample.pep.faa']))

        for s in read_fasta_file(os.sep.join([step1out,'all_sample.pep.faa'])):
                seqid,seq = s.name,s.seq
                try:
                        clusterID = clusterDICT[seqid]
                        with open(outdir+"cluster"+clusterID+".fa","a") as outfile:
                            outfile.write(">"+seqid+"\n"+seq+"\n")
                            fa_list.append("cluster"+str(clusterID))
                except:
                        pass # Those seqs that did not go in a cluster with enough taxa
                        # will not be in clusterDICT
        fa_list=list(set(fa_list))
        if args.cds2tree:
            logging.info("Choose to use cds sequence for gene tree construction and sort out the corresponding cds sequence files.")
            for s in read_fasta_file(os.sep.join([step1out,'all_cds.fa'])):
                    seqid,seq = s.name,s.seq
                    try:
                            clusterID = clusterDICT[seqid]
                            with open(outdir+"cluster"+clusterID+".cds","a") as outfile:
                                outfile.write(">"+seqid+"\n"+seq+"\n")
                    except:
                            pass
        return fa_list

def fa2tree(muscle,iqtree,fa,step4out,iqtreeop,args,Epal2nal):
    muscle_stdout=subprocess.run([muscle,'-in',os.sep.join([step4out,'all_fa',fa+".fa"]),'-out',os.sep.join([step4out,'muscle_out',fa+".pep.aln"]),],stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    logging.debug('{} stdout:\n'.format(muscle) +muscle_stdout.stdout.decode('utf-8'))
    logging.debug('{} stderr:\n'.format(muscle) +muscle_stdout.stderr.decode('utf-8'))
    if args.cds2tree:
        Epal2nal_stdout=subprocess.run(Epal2nal+" "+os.sep.join([step4out,'muscle_out',fa+".pep.aln"])+" "+os.sep.join([step4out,'all_fa',fa+".cds"])+" -output fasta "+" > "+os.sep.join([step4out,'muscle_out',fa+".cds.aln"]),shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logging.debug('{} stdout:\n'.format(Epal2nal) +Epal2nal_stdout.stdout.decode('utf-8'))
        logging.debug('{} stderr:\n'.format(Epal2nal) +Epal2nal_stdout.stderr.decode('utf-8'))
        iqtree_stdout=subprocess.run(iqtree+" -s "+os.sep.join([step4out,'muscle_out',fa+".cds.aln"])+" "+iqtreeop+" -pre "+os.sep.join([step4out,'iqtree_out',fa+".iqtree"]),shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logging.debug('{} stdout:\n'.format(iqtree) +iqtree_stdout.stdout.decode('utf-8'))
        logging.debug('{} stderr:\n'.format(iqtree) +iqtree_stdout.stderr.decode('utf-8'))
    else:
        iqtree_stdout=subprocess.run(iqtree+" -s "+os.sep.join([step4out,'muscle_out',fa+".pep.aln"])+" "+iqtreeop+" -pre "+os.sep.join([step4out,'iqtree_out',fa+".iqtree"]),shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logging.debug('{} stdout:\n'.format(iqtree) +iqtree_stdout.stdout.decode('utf-8'))
        logging.debug('{} stderr:\n'.format(iqtree) +iqtree_stdout.stderr.decode('utf-8'))
def run_fa2tree(args):
    fa2tree(args[0],args[1],args[2],args[3],args[4],args[5],args[6])

def run_tree2gd(step1out,args,config,step2out,step4out):
    logging.info("Writing step4 script files...")
    step4_sh=open(os.sep.join([step4out,'step4.sh']),"w")
    iqtreeop=''
    if not config.has_section('iqtree'):
        config.add_section('iqtree')
    if not config.has_option("iqtree", "-bb"):
        config.set("iqtree", "-bb","1000")
    if config.get("iqtree", "-bb")=="0":
        config.remove_option("iqtree", "-bb")
    if not config.has_option("iqtree", "-m"):
        if args.cds2tree:
            config.set("iqtree", "-m","HKY")
        else:
            config.set("iqtree", "-m","JTT+G4")
    for k,v in config.items("iqtree"):
        iqtreeop=' '.join((iqtreeop,' '.join((k,v))))
    tree2gd=config.get("software", "tree2gd")
    muscle=config.get("software", "muscle")
    iqtree=config.get("software", "iqtree")
    if args.cds2tree:
        Epal2nal=config.get("software", "Epal2nal")
    else:
        Epal2nal="Epal2nal.pl"
    step4_sh.writelines(["for F in "+os.sep.join([step4out,'all_fa/'])+"*.fa;do "+muscle+" -in $F -out "+os.sep.join([step4out,'muscle_out/'])+"`basename $F`.aln;done\n",
                            "for F in "+os.sep.join([step4out,'muscle_out/'])+"*.aln;do "+iqtree+iqtreeop+" -s $F -pre "+os.sep.join([step4out,'iqtree_out/'])+"`basename $F`.iqtree;done \n",
                            "ls "+os.sep.join([step4out,'iqtree_out/'])+"*.contree > "+os.sep.join([step4out,'tree.list.tmp'])+"\n",
                            "cat -n "+os.sep.join([step4out,'tree.list.tmp'])+" > "+ os.sep.join([step4out,'tree.list'])+"\n",
                            tree2gd+" "+args.tree+" "+os.sep.join([step1out,'all_sample2fa.list'])+" "+os.sep.join([step4out,'tree.list'])+" "+os.sep.join([step4out,'Tree2GD_out/'])
                            ])
    step4_sh.close()
    if not args.only_script:
        if not os.path.exists(os.sep.join([step4out,'all_fa/'])):
            os.mkdir(os.sep.join([step4out,'all_fa/']))
        minimal_taxa=config.get("mcl2fasta", "min_taxa") if config.has_option("mcl2fasta", "min_taxa") and config.get("mcl2fasta", "min_taxa")!="" else '4'
        if int(minimal_taxa)<4:
            logging.warning("The mcl2fasta min_taxa cannot be set less than 4, otherwise a valid tree cannot be obtained ,it will be changed to the default 4!")
            minimal_taxa="4"
        fa_list=mcl2fasta(int(minimal_taxa),os.sep.join([step4out,'all_fa/']),step2out,step1out,args)
        logging.info("FASTA file selection has been completed ,according to MCL results.")
        try:
            subprocess.run(muscle,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        except FileNotFoundError:
            logging.error("muscle can not run.")
            sys.exit(0)
        try:
            subprocess.run(iqtree,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        except FileNotFoundError:
            logging.error("iqtree can not run.")
            sys.exit(0)
        try:
            subprocess.run(tree2gd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        except FileNotFoundError:
            logging.error("tree2gd can not run.")
            sys.exit(0)
        sh_list=[]
        if not os.path.exists(os.sep.join([step4out,'muscle_out/'])):
            os.mkdir(os.sep.join([step4out,'muscle_out/']))
        if not os.path.exists(os.sep.join([step4out,'iqtree_out/'])):
            os.mkdir(os.sep.join([step4out,'iqtree_out/']))
        for fa in fa_list:
            tmp=(muscle,iqtree,fa,step4out,iqtreeop,args,Epal2nal)
            sh_list.append(tmp)
        sh_pool=Pool(args.t)
        logging.info("Start runing step4 fa2tree.")
        sh_pool.map(run_fa2tree,sh_list)
        sh_pool.close()
        sh_pool.join()
        logging.info("step4 fa2tree has done.")
        logging.info("All sequences have completed the construction of the tree.")
        logging.info("Start WGD calculation.")
        tree2gdop=""
        if not config.has_section('tree2gd'):
            config.add_section('tree2gd')
        if not config.has_option("iqtree", "-bb"):
            config.set("tree2gd", "--bp","0")
        for k,v in config.items("tree2gd"):
            tree2gdop=' '.join((tree2gdop,'='.join((k,v))))
        subprocess.run("ls "+os.sep.join([step4out,"iqtree_out","*.treefile"])+" > "+os.sep.join([step4out,'tree.list.tmp']),shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        subprocess.run("cat -n "+os.sep.join([step4out,'tree.list.tmp'])+">"+os.sep.join([step4out,'tree.list']),shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        tree2gd_out=subprocess.run(tree2gd+" "+args.tree+" "+os.sep.join([step1out,'all_sample2fa.list'])+" "+os.sep.join([step4out,'tree.list'])+" "+os.sep.join([step4out,'Tree2GD_out/'])+" "+tree2gdop,shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logging.debug('{} stdout:\n'.format(tree2gd) +tree2gd_out.stdout.decode('utf-8'))
        logging.debug('{} stderr:\n'.format(tree2gd) +tree2gd_out.stderr.decode('utf-8'))
        logging.info("step4 Tree2GD has done.")
