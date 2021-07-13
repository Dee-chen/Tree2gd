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
    if not config.has_section('KaKs_Calculator'):
        config.add_section('KaKs_Calculator')
    if not config.has_option("KaKs_Calculator", "-m"):
        config.set("KaKs_Calculator", "-m","YN")
    for k,v in config.items("KaKs_Calculator"):
        op=' '.join((op,' '.join((k,v))))
    step5_sh=open(os.sep.join([step5out,'step5.sh']),"w")
    KaKs_Calculator=config.get("software", "KaKs_Calculator")
    Kaks_aligncmd=config.get("software", "muscle")
    Epal2nal=config.get("software", "Epal2nal")
    step5_sh.writelines([   "perl ParaAT.pl -kaks "+KaKs_Calculator+op+" -h "+os.sep.join([step5out,"all_gene_pairs.list"])
                            +" -n "+os.sep.join([step5out,"all_cds.fa"])+" -a "+os.sep.join([step5out,'all_pep.fa'])+" -p "+str(args.t)+" -k -f axt -o "+os.sep.join([step5out,"all_gene_pairs_kaks"])+"\n",
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
        pepmap={}
        cdsmap={}
        logging.info("Readfile")
        for s in read_fasta_file(os.sep.join([step5out,'all_pep.fa'])):
                pepmap[s.name]=s.seq
        for s in read_fasta_file(os.sep.join([step5out,'all_cds.fa'])):
                cdsmap[s.name]=s.seq
        if not os.path.exists(os.sep.join([step5out,"sp_kaks_out"])):
            os.mkdir(os.sep.join([step5out,"sp_kaks_out"]))
        if not os.path.exists(os.sep.join([step5out,"all_gene_pairs_kaks"])):
            os.mkdir(os.sep.join([step5out,"all_gene_pairs_kaks"]))
        logging.info("Make dir")
        os.chdir(os.sep.join([step5out,"all_gene_pairs_kaks"]))
        logging.info("Start gene pairs kaks calculation.")
        a=0
        pairs_num=len(gene_pairs_idmap)
        b=pairs_num//20
        arg_list=[]

        for pair in gene_pairs_idmap:
            tmp=(pair,pepmap,cdsmap,KaKs_Calculator,Kaks_aligncmd,Epal2nal,op)
            arg_list.append(tmp)
        sh_pool=Pool(args.t)
        sh_pool.map(run_ks,arg_list)
            #a=a+1
            #sh_pool.apply_async(sub_sh,args=(pair,pepmap,cdsmap,KaKs_Calculator,Kaks_aligncmd,Epal2nal,op),callback=write_result)
            #if(a%b==0):
            #    logging.info("%s \%(%s pairs) completed."%((a//b)*5)%(a))
        sh_pool.close()
        sh_pool.join()
        logging.info("Gene pairs kaks calculation finish.")
        logging.info("Start summarize the Ks results of each species.")
        for pair in gene_pairs_idmap:
            out=sum_Ks(pair)
            write_result(out)
            a=a+1
            if(a%b==0):
                logging.info("%d %% (%d pairs) completed."%((a//b)*5,a))

        logging.info("ALL step5 kaks has done.")

def run_ks(args):
    sub_sh(args[0],args[1],args[2],args[3],args[4],args[5],args[6])

def sum_Ks(pair):
    out=[]
    with open(pair[2]+"-"+pair[3]+".cds_aln.axt.kaks") as kaks_file,open(pair[2]+"-"+pair[3]+".4dtv") as DTV_file:
        kaks=kaks_file.readlines()[1].split("\t")
        dtv=DTV_file.readlines()[1].split("\t")
        out=[pair[1],pair[0],pair[2],pair[3],kaks[2],kaks[3],kaks[4],dtv[1]]
    return out

def sub_sh(pair,pepmap,cdsmap,KaKs_Calculator,Kaks_aligncmd,Epal2nal,op):
    out=[]
    home_dir=os.path.dirname(os.path.abspath(__file__))
    try:
        pepseq1 = pepmap[pair[2]]
        pepseq2 = pepmap[pair[3]]
        cdsseq1 = cdsmap[pair[2]]
        cdsseq2 = cdsmap[pair[3]]
        with open(pair[2]+"-"+pair[3]+".pep","a") as outfile:
            outfile.write(">"+pair[2]+"\n"+pepseq1+"\n")
            outfile.write(">"+pair[3]+"\n"+pepseq2+"\n")
        with open(pair[2]+"-"+pair[3]+".cds","a") as outfile:
            outfile.write(">"+pair[2]+"\n"+cdsseq1+"\n")
            outfile.write(">"+pair[3]+"\n"+cdsseq2+"\n")

    except:
        logging.info("err")
        pass
    stdout=subprocess.run(Kaks_aligncmd+" -quiet  -in "+pair[2]+"-"+pair[3]+".pep  -out "+pair[2]+"-"+pair[3]+".pep_aln >> msg.msa",shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    Epal2nal_stdout=subprocess.run(Epal2nal+" "+pair[2]+"-"+pair[3]+".pep_aln "+pair[2]+"-"+pair[3]+".cds -output fasta -nogap -nomismatch > "+pair[2]+"-"+pair[3]+".cds_aln",shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    Fasta2AXT(pair[2]+"-"+pair[3]+".cds_aln",pair[2]+"-"+pair[3]+".cds_aln.axt")
    kaks_stdout=subprocess.run(KaKs_Calculator+op+" -i "+pair[2]+"-"+pair[3]+".cds_aln.axt -o  "+pair[2]+"-"+pair[3]+".cds_aln.axt.kaks >> msg.kaks",shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    axt2oneline((pair[2]+"-"+pair[3]+".cds_aln.axt"),(pair[2]+"-"+pair[3]+".one-line"))
    stdout=subprocess.run("perl "+os.sep.join([home_dir,"software","calculate_4DTV_correction.pl"])+" "+pair[2]+"-"+pair[3]+".one-line"+" > "+pair[2]+"-"+pair[3]+".4dtv",shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)


def write_result(input):
    with open(os.sep.join(["..","sp_kaks_out",input[0]+".kaks_4dtv.result"]),"a+") as out:
        out.write("\t".join([input[0],input[1],input[2],input[3],input[4],input[5],input[6],input[7]])+"\n")


def Fasta2AXT(input,output):
    id=""
    seq=""
    for s in read_fasta_file(input):
        if id=="":
            id=s.name
        else:
            id=id+"-"+s.name
        if seq=="":
            seq=s.seq
        else:
            seq=seq+'\n'+s.seq
    with open(output,"w") as outfile:
        outfile.write(id+"\n"+seq+"\n")

def axt2oneline(reader,out):#copyed
    mydict = {}
    writer = open(out,"w")
    with open(reader) as fh:
            for line in fh:
                    if re.search("-",line):
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
