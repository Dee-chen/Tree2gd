#!/usr/bin/python3.7
import os
import subprocess
import logging
import sys
import configparser
import subprocess
from multiprocessing import Pool
import re

from pyecharts import options as opts
from pyecharts.charts import Bar
from pyecharts.faker import Faker


def run_plot(step6out,args,cf):
    logging.info("Start step6 summary plot...")
    os.chdir(step6out)
    if not args.only_script:
        stdout=subprocess.run('Rscript '+'Tree2GD_draw.R',shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        logging.info("R plot done.")
        logging.debug('{} stdout:\n'.format("Tree2GD_draw.R") +stdout.stdout.decode('utf-8'))
        logging.debug('{} stderr:\n'.format("Tree2GD_draw.R") +stdout.stderr.decode('utf-8'))
        logging.info("Start html plot ...")
        if not os.path.exists(os.sep.join([step6out,"ks_html_plot"])):
            os.mkdir(os.sep.join([step6out,"ks_html_plot"]))
        file_names=[]
        file_names = os.listdir("html_plot_in")
        os.chdir(os.sep.join([step6out,"html_plot_in"]))
        sh_pool=Pool(args.t)
        sh_pool.map(html_plot,file_names)
        sh_pool.close()
        sh_pool.join()
        logging.info("ALL step6 plot has done.")

def html_plot(file):
    sp=file.split(".")
    with open(file) as f:
        l=f.readlines()
        i=0
        ks={}
        for line in l:
            if i==0:
                head=line.strip().split("\t")
                i=1
            else:
                t=line.strip().split("\t")
                name=t[0]
                del t[0]
                ks[name]=t
        bar=Bar(init_opts=opts.InitOpts(page_title=sp[0]+"-ks.htmlpolt",width="1500px",height="700px"))
        bar.add_xaxis(head)
        for key,value in ks.items():
            bar.add_yaxis(key,value,category_gap=0,stack="stack1",is_large=True)
        bar.set_series_opts(label_opts=opts.LabelOpts(is_show=False))
        bar.set_global_opts(
        title_opts=opts.TitleOpts(title=sp[0]+".ks.polt",subtitle="Plot by Tree2gd v1.0"),
        toolbox_opts=opts.ToolboxOpts(is_show=True,pos_top='10%',orient='vertical',pos_left='right',feature=opts.ToolBoxFeatureOpts(
        save_as_image=opts.ToolBoxFeatureSaveAsImageOpts(title="Save as png"),
        restore=opts.ToolBoxFeatureRestoreOpts(title="Recovery"),
        data_view=opts.ToolBoxFeatureDataViewOpts(title="Data",lang=['Data View','Close','Refresh']),
        data_zoom=opts.ToolBoxFeatureDataZoomOpts(is_show=False),
        magic_type=opts.ToolBoxFeatureMagicTypeOpts(line_title="Switch to line plot",bar_title="Switch to histogram plot",stack_title="Switch to stack plot",tiled_title="Switch to tiled plot"),
        brush=opts.ToolBoxFeatureBrushOpts(rect_title="Rectangle selection",polygon_title="Circle selection",line_x_title="X line selection",line_y_title="Y line selection",keep_title="Keep selection",clear_title="Clear selection"),
        )),
        legend_opts=opts.LegendOpts(pos_left = 'right',),
        datazoom_opts=[opts.DataZoomOpts(range_start=0,range_end=60),opts.DataZoomOpts(type_="inside")],
        xaxis_opts=opts.AxisOpts(name="Ks"),)
        bar.render(os.sep.join(["..","ks_html_plot",sp[0]+".ks.plot.html"]))
