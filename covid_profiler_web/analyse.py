import sys
from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for, Response
)
from werkzeug.exceptions import abort
import subprocess as sp
#from aseantb.auth import login_required
from covid_profiler_web.db import get_db, get_mongo_db
from covid_profiler_web.worker import run_profile, run_phylogeny, run_primer_conservation
import uuid
from werkzeug.utils import secure_filename
import os
from flask import current_app as app
import datetime
import re
import covid_profiler as cp
import csv
import json
import gzip

bp = Blueprint('analyse', __name__)

def check_file_type(filename):
    for l in sp.Popen(f"file {filename}",shell=True,stdout=sp.PIPE).stdout:
        row = l.decode().strip().split()
        if row[1]!="ASCII" and row[1]!="gzip":
            return "unknown"
        if row[1]=="ASCII":
            first_char = open(filename).readline()[0]
        else:
            first_char =  gzip.open(filename).readline().decode()[0]
        if first_char==">":
            return "fasta"
        elif first_char=="@":
            return "fastq"
        else:
            return "unknown"


def run_sample(mongo,user_id,uniq_id,sample_name,type,f1,f2=None):
    cp.log("Running sample with id:%s, file1:%s and file2:%s\n" % (uniq_id,f1,f2 ))
    filename1 = secure_filename(f1.filename)
    server_fname1 = os.path.join(app.config["UPLOAD_FOLDER"], filename1)
    f1.save(server_fname1)
    if f2:
        filename2 = secure_filename(f2.filename)
        server_fname2 = os.path.join(app.config["UPLOAD_FOLDER"], filename2)
        f2.save(server_fname2)
    else:
        server_fname2 = None
    filetype = check_file_type(server_fname1)
    print(filetype)
    if filetype!=type:
        return "Error! File type doesn't look like %s, please check." % type
    mongo.db.profiler_results.insert_one({
        "_id":uniq_id,"user_id":user_id,"sample_name":sample_name,"results":{},
        "created":datetime.datetime.now().strftime("%d-%M-%Y %H:%M:%S"),
        "data_type": type, "status":"processing"
    })
    if type=="fasta":
        run_profile.delay(fasta=server_fname1,uniq_id=uniq_id,storage_dir=app.config["APP_ROOT"]+url_for('static', filename='results'))
    elif type=="fastq":
        run_profile.delay(R1=server_fname1, R2=server_fname2, uniq_id=uniq_id,storage_dir=app.config["APP_ROOT"]+url_for('static', filename='results'))
    else:
        sys.stderr.write("ERROR!!! Not fastq or fasta?")


@bp.route('/analyse/profile',methods=('GET', 'POST'))
def profile():
    mongo = get_mongo_db()
    if request.method == 'POST':
        cp.log("Recieved POST request with files %s\n" % (", ".join([request.files[x].filename for x in request.files])) )
        error=None
        # username = g.user['username'] if g.user else 'private'
        username = 'private'
        uniq_id = str(uuid.uuid4())
        sample_name = uniq_id
        if request.files['file1'].filename=="":
            error = "No files found, please try again!"
        if "fasta_submit" in request.form and error==None:
            error = run_sample(mongo,username,uniq_id,sample_name,"fasta",request.files['file1'])
            if error==None:
                return redirect(url_for('results.run_result', sample_id=uniq_id))
        elif "fastq_submit" in request.form and error==None:
            error = run_sample(mongo,username,uniq_id,sample_name,"fastq",request.files['file1'],request.files['file2'])
            if error==None:
                return redirect(url_for('results.run_result', sample_id=uniq_id))

        flash(error)
    return render_template('analyse/profile.html')


@bp.route('/analyse')
def analyse():
    mongo = get_mongo_db()
    return render_template('analyse/analyse.html')


@bp.route('/analyse/phylogeny',methods=('GET', 'POST'))
def phylogeny():
    mongo = get_mongo_db()
    if request.method == 'POST':
        uniq_id = str(uuid.uuid4())
        if request.files['file1'].filename=="":
            error = "No files found, please try again!"
        filename1 = secure_filename(request.files['file1'].filename)
        server_fname1 = os.path.join(app.config["UPLOAD_FOLDER"], filename1)
        request.files['file1'].save(server_fname1)
        run_phylogeny.delay(server_fname1,uniq_id,working_dir=app.config["APP_ROOT"]+url_for('static', filename='results'))

        return redirect(url_for('analyse.phylogeny_result', run_id=uniq_id))


    return render_template('analyse/phylogeny.html')


@bp.route('/analyse/phylogeny/results/<uuid:run_id>')
def phylogeny_result(run_id):
    mongo = get_mongo_db()
    file_prefix = "%s/%s" % (app.config["APP_ROOT"]+url_for('static', filename='results'), run_id)
    status = {
        "msa":"OK" if os.path.isfile(file_prefix+".aln") else "Processing",
        "vcf": "OK" if os.path.isfile(file_prefix+".vcf.gz") else "Processing",
        "variant_info": "OK" if os.path.isfile(file_prefix+".variant_info.csv") else "Processing",
        "phylogeny": "OK" if os.path.isfile(file_prefix+".aln.iqtree") else "Processing",
        "num_seqs": int(open(file_prefix+".numseqs.txt").readline().strip()) if os.path.isfile(file_prefix+".numseqs.txt") else "Processing",
    }
    print(status)
    csv_data = []
    if status["variant_info"]=="OK":
        for row in csv.DictReader(open(file_prefix+".variant_info.csv")):
            csv_data.append(row)

    return render_template('analyse/phylogeny_result.html',run_id=run_id, status=status, csv_data=csv_data)


@bp.route('/analyse/primers',methods=('GET', 'POST'))
def primers():
    mongo = get_mongo_db()
    if request.method == 'POST':
        uniq_id = str(uuid.uuid4())
        primerF = request.form["primerF"]
        primerR = request.form["primerR"]
        probe = request.form["probe"]
        save_dir = app.root_path+"/"+app.static_url_path+"/results/"
        run_primer_conservation.delay(primerF=primerF,primerR=primerR,probe=probe,uniq_id=uniq_id, save_dir=save_dir)
        return redirect(url_for('analyse.primer_result', run_id=str(uniq_id)))


    return render_template('analyse/diagnostics.html')


@bp.route('/analyse/primers/<uuid:run_id>',methods=('GET', 'POST'))
def primer_result(run_id):
    mongo = get_mongo_db()
    json_file = app.root_path+"/"+app.static_url_path+"/results/%s.plots.json" % str(run_id)
    data = json.load(open(json_file)) if os.path.isfile(json_file) else None

    return render_template('analyse/primer_result.html',run_id=str(run_id),data = data)
