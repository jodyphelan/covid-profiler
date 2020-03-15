from celery import Celery
import subprocess
import json
import sqlite3
import pathogenprofiler as pp
import os
import sys


celery = Celery('tasks', broker='pyamqp://guest@localhost//')



@celery.task
def profile(fasta,uniq_id,db,storage_dir):
    pp.run_cmd("covid-profiler.py profile --fasta %s > %s/results.txt" % (fasta,storage_dir))
    pp.run_cmd("covid-profiler.py position_isolate --fasta %s > %s/node.txt" % (fasta,storage_dir))
    results = {}
    for l in open("%s/results.txt" % (storage_dir)):
        row = l.strip().split()
        results["type"] = row[1]
    results["branch"] = open("%s/node.txt" % storage_dir).readline().strip()

    conn = sqlite3.connect(db)
    c = conn.cursor()
    c.execute("UPDATE results SET result = ? where id = ?", (json.dumps(results),uniq_id,))
    conn.commit()

    return True
