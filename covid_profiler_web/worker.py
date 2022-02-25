from celery import Celery
import subprocess
import json
import sqlite3
import pathogenprofiler as pp
import os
import sys
from flask import Flask
from pymongo import MongoClient
import covid_profiler as cp


celery = Celery('tasks', broker='redis://localhost:6379/0')

@celery.task
def run_phylogeny(file,uniq_id,working_dir="/tmp/"):
    cp.log("This is the worker. Running %s" % uniq_id)
    count=0
    for l in open(file):
        if ">" in l:
            count+=1
    cp.log(f"THE COUNT IS {count}")
    with open(f"{working_dir}/{uniq_id}.numseqs.txt","w") as O:
        O.write(str(count))
    if count<=500:
        pp.run_cmd("covid_profiler_align_fasta.py --fasta %s --working-dir %s --out %s" % (file,working_dir,uniq_id))
    return True

@celery.task
def run_profile(uniq_id,storage_dir,fasta=None,R1 = None, R2 = None):
    cp.log("This is the worker. Running %s" % uniq_id)
    if fasta:
        pp.run_cmd("covid-profiler.py profile --fasta %s --prefix %s --dir %s" % (fasta,uniq_id,storage_dir))
    elif R1 and not R2:
        pp.run_cmd("covid-profiler.py profile -1 %s --prefix %s --dir %s" % (R1,uniq_id,storage_dir))
    elif R1 and R2:
        pp.run_cmd("covid-profiler.py profile -1 %s -2 %s --prefix %s --dir %s" % (R1,R2,uniq_id,storage_dir))
    else:
        sys.stderr.write("ERROR!!! Check file inputs to profile worker!")
    pp.run_cmd("zip -j %s/%s.zip %s/%s*" % (storage_dir, uniq_id, storage_dir, uniq_id))
    results = json.load(open("%s/%s.results.json" % (storage_dir,uniq_id)))

    if R1:
        pp.run_cmd("bcftools view %s/%s.vcf.gz > %s/%s.vcf" % (storage_dir, uniq_id, storage_dir, uniq_id))
        for l in pp.cmd_out("bedtools genomecov -ibam %s/%s.bam -d | datamash mean 3" % (storage_dir,uniq_id)):
            cp.log(l)
            results["mean_depth"] = round(float(l.strip()),2)
    results["num_variants"] = len(results["variants"])

    client = MongoClient()
    db = client.test_database
    db.profiler_results.find_one_and_update(
        {"_id":uniq_id},
        {
            "$set": {"results":results,"status":"done"}
        }
    )

    return True

@celery.task
def run_primer_conservation(primerF,primerR,probe,uniq_id,save_dir):
    pp.run_cmd("primer_analysis.py --fp %s --rp %s --probe %s --dir %s --out %s --write-json" % (primerF,primerR,probe,save_dir,uniq_id))

    return True
