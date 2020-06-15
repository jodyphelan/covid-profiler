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
def profile(uniq_id,storage_dir,fasta=None,R1 = None, R2 = None):
    cp.log("This is the worker. Running %s" % uniq_id)
    if fasta:
        pp.run_cmd("covid-profiler.py profile --fasta %s --prefix %s --dir %s" % (fasta,uniq_id,storage_dir))
    elif R1 and not R2:
        pp.run_cmd("covid-profiler.py profile -1 %s --prefix %s --dir %s" % (R1,uniq_id,storage_dir))
    elif R1 and R2:
        pp.run_cmd("covid-profiler.py profile -1 %s -2 %s --prefix %s --dir %s" % (R1,R2,uniq_id,storage_dir))
    else:
        sys.stderr.write("ERROR!!! Check file inputs to profile worker!")

    results = json.load(open("%s/%s.results.json" % (storage_dir,uniq_id)))

    print("Updating database")
    print(uniq_id)
    print(results)
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
def profile_primer(primerF,primerR,probe,uniq_id,save_dir):
    pp.run_cmd("covid-profiler.py primer --primerF %s --primerR %s --probe %s --out %s/%s.csv" % (primerF,primerR,probe,save_dir,uniq_id))
    pp.run_cmd("covid_plot_primers.R %s/%s.csv %s %s %s" % (save_dir,uniq_id,primerF,primerR,probe))
    pp.run_cmd("rm %s/%s.csv" % (save_dir,uniq_id))
    client = MongoClient()
    db = client.test_database
    db.primer_results.find_one_and_update(
        {"_id":uniq_id},
        {
            "$set": {"status":"done"}
        }
    )

    return True
