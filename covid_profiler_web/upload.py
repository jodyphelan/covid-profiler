from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for, Response
)
from werkzeug.exceptions import abort
import subprocess
#from aseantb.auth import login_required
from covid_profiler_web.db import get_db, get_mongo_db
from covid_profiler_web.worker import profile
import uuid
from werkzeug.utils import secure_filename
import os
from flask import current_app as app
import datetime
import re
bp = Blueprint('upload', __name__)

def run_sample(mongo,user_id,uniq_id,sample_name,type,f1,f2=None):
    filename1 = secure_filename(f1.filename)
    server_fname1 = os.path.join(app.config["UPLOAD_FOLDER"], filename1)
    f1.save(server_fname1)
    if f2:
        filename2 = secure_filename(f2.filename)
        server_fname2 = os.path.join(app.config["UPLOAD_FOLDER"], filename2)
        f2.save(server_fname2)
    else:
        server_fname2 = None

    mongo.db.profiler_results.insert_one({
        "_id":uniq_id,"user_id":user_id,"sample_name":sample_name,"results":{},
        "created":datetime.datetime.now().strftime("%d-%M-%Y %H:%M:%S"),
        "status":"processing"
    })
    if type=="fasta":
        profile.delay(fasta=server_fname1,uniq_id=uniq_id,storage_dir=app.config["UPLOAD_FOLDER"])
    elif type=="fastq":
        profile.delay(R1=server_fname1, R2=server_fname2, uniq_id=uniq_id,storage_dir=app.config["UPLOAD_FOLDER"])
    else:
        sys.stderr.write("ERROR!!! Not fastq or fasta?")

@bp.route('/upload',methods=('GET', 'POST'))
def upload():
    mongo = get_mongo_db()
    if request.method == 'POST':
        error=None
        # username = g.user['username'] if g.user else 'private'
        username = 'private'
        uniq_id = str(uuid.uuid4())
        sample_name = uniq_id
        if request.files['file1'].filename=="":
            error = "No files found, please try again!"
        if "fasta_submit" in request.form and error==None:
            run_sample(mongo,username,uniq_id,sample_name,"fasta",request.files['file1'])
            return redirect(url_for('results.run_result', sample_id=uniq_id))
        elif "fastq_submit" in request.form and error==None:
            run_sample(mongo,username,uniq_id,sample_name,"fastq",request.files['file1'],request.files['file2'])
            return redirect(url_for('results.run_result', sample_id=uniq_id))

        flash(request.form)
        flash(error)
    return render_template('upload/upload.html')
