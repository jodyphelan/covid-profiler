from celery import Celery
import subprocess
import json
import sqlite3
import pathogenprofiler as pp
import os
import sys
from flask import Flask
from pymongo import MongoClient

celery = Celery('tasks', broker='pyamqp://guest@localhost//')


@celery.task
def profile(fasta,uniq_id,storage_dir):
    pp.run_cmd("covid-profiler.py profile --fasta %s --prefix %s --dir %s" % (fasta,uniq_id,storage_dir))

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
    pp.run_cmd("covid-profiler.py primer --msa /tmp/gisaid.filtered.fasta --primerF %s --primerR %s --probe %s --out /tmp/%s.csv" % (primerF,primerR,probe,uniq_id))
    pp.run_cmd("Rscript /Users/jodyphelan/github/covid-profiler/scripts/plot_primers.R /tmp/%s.csv %s %s %s" % (uniq_id,primerF,primerR,probe))
    pp.run_cmd("mv /tmp/%s.svg %s/" % (uniq_id,save_dir))

    client = MongoClient()
    db = client.test_database
    db.primer_results.find_one_and_update(
        {"_id":uniq_id},
        {
            "$set": {"status":"done"}
        }
    )

    return True
