#!/usr/bin/python3.7
import os
import subprocess
import logging
import sys
from tree2gd.seq import read_fasta_file
import configparser
import subprocess
from multiprocessing import Pool
import re

def run_kaks(sp_list,step1out,args,config,step4out,step5out,gene_pairs_idmap):
    logging.info("Writing step5 script files...")
    op=''
    if not config.has_section('ParaAT'):
        config.add_section('ParaAT')
    if not config.has_option("ParaAT", "-mod"):
        config.set("ParaAT", "-mod","YN")
    if not config.has_option("ParaAT", "-m"):
        config.set("ParaAT", "-m","muscle")
    for k,v in config.items("ParaAT"):
        op=' '.join((op,' '.join((k,v))))
    step5_sh=open(os.sep.join([step5out,'step5.sh']),"w")
    KaKs_Calculator=config.get("software", "KaKs_Calculator")
    ParaAT=config.get("software", "ParaAT")
    ParaAT_aligncmd=config.get("software", "ParaAT_aligncmd")
    Epal2nal=config.get("software", "Epal2nal")
    step5_sh.writelines(['perl '+os.sep.join([os.path.dirname(os.path.abspath(__file__)),"software","ParaAT.pl"])+" -kaks "+KaKs_Calculator+op+" -agcmd "+ParaAT_aligncmd+" -k -h "+os.sep.join([step5out,"all_gene_pairs.list"])
                            +" -n "+os.sep.join([step5out,"all_cds.fa"])+" -a "+os.sep.join([step5out,'all_sample.pep.faa'])+" -p "+str(args.t)+" -k -f axt -o "+os.sep.join([step5out,"all_gene_pairs_kaks"])+"\n",
                            "cd "+os.sep.join([step5out,"all_gene_pairs_kaks"])+"\n",
                            "for i in `ls *.axt`;do axt2one-line.py $i ${i}.one-line;done\n",
                            "ls *.one-line|while read id;do calculate_4DTV_correction.pl $id >${id%%one-line}4dtv;done\n"
                            "for i in `ls *.4dtv`;do awk 'NR>1{print $1\"\t\"$2}' $i >>all-4dtv.txt;done\n",
                            "for i in `ls *.kaks`;do awk 'NR>1{print $1\"\t\"$3\"\t\"$4\"\t\"$5}' $i >>all-kaks.txt;done\n",
                            "sort all-4dtv.txt|uniq >all-4dtv.results\n",
                            "sort all-kaks.txt|uniq >all-kaks.results\n",
                            "join -1 1 -2 1 all-kaks.results all-4dtv.results >all-results.txt\n",
                            "sed -i '1i\Seq\tKa\tKs\tKa/Ks\t4dtv_corrected' all-results.txt\n"
                            ])
    step5_sh.close()
    if not args.only_script:
        logging.info("Start KaKs calculation.")
        stdout=subprocess.run('perl '+ParaAT+' -kaks '+KaKs_Calculator+" "+op+' -h '+os.sep.join([step5out,"all_gene_pairs.list"])+' -n '+os.sep.join([step5out,"all_cds.fa"])+' -a '+os.sep.join([step5out,'all_pep.fa'])+' -agcmd '+ParaAT_aligncmd+' -p2n '+Epal2nal+' -p '+str(args.t)+" -f axt -k -o "+os.sep.join([step5out,"all_gene_pairs_kaks"]),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        logging.info("ParaAT done.")
        logging.debug('{} stdout:\n'.format("ParaAT.pl") +stdout.stdout.decode('utf-8'))
        logging.debug('{} stderr:\n'.format("ParaAT.pl") +stdout.stderr.decode('utf-8'))
        if not os.path.exists(os.sep.join([step5out,"sp_kaks_out"])):
            os.mkdir(os.sep.join([step5out,"sp_kaks_out"]))
        os.chdir(os.sep.join([step5out,"all_gene_pairs_kaks"]))
        sh_pool=Pool(args.t)
        for pair in gene_pairs_idmap:
            sh_pool.apply_async(sub_sh,(pair,),callback=write_result)
        sh_pool.close()
        sh_pool.join()
        logging.info("ALL step5 kaks has done.")


def sub_sh(pair):
    out=[]
    axt2oneline((pair[2]+"-"+pair[3]+".cds_aln.axt"),(pair[2]+"-"+pair[3]+".one-line"))
    stdout=subprocess.run("perl "+os.sep.join([os.path.dirname(os.path.abspath(__file__)),"..","software","calculate_4DTV_correction.pl"])+" "+pair[2]+"-"+pair[3]+".one-line"+" > "+pair[2]+"-"+pair[3]+".4dtv",shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    with open(pair[2]+"-"+pair[3]+".cds_aln.axt.kaks") as kaks_file,open(pair[2]+"-"+pair[3]+".4dtv") as DTV_file:
        kaks=kaks_file.readlines()[1].split("\t")
        dtv=DTV_file.readlines()[1].split("\t")
        out=[pair[1],pair[0],pair[2],pair[3],kaks[2],kaks[3],kaks[4],dtv[1]]
    return out

def write_result(input):
    with open(os.sep.join(["..","sp_kaks_out",input[0]+".kaks_4dtv.result"]),"a+") as out:
        out.write("\t".join([input[0],input[1],input[2],input[3],input[4],input[5],input[6],input[7]])+"\n")

def axt2oneline(reader,out):#copyed
    mydict = {}
    writer = open(out,"w")
    with open(reader) as fh:
            for line in fh:
                    if re.search(r"\d+",line):
                        header = line.strip()
                        mydict[header] = []
                    else:
                        mydict[header].append(line.strip())
    for key,value in mydict.items():
        seqs = "".join(value)
        length = len(seqs)
        mid = int(length/2)
        seq1 = seqs[0:mid]
        seq2 = seqs[mid:]
        writer.write(key + "\n" + "\n".join([seq1,seq2]))
    writer.close()
