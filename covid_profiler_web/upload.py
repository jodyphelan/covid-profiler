from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for, Response
)
from werkzeug.exceptions import abort
import subprocess
#from aseantb.auth import login_required
from covid_profiler_web.db import get_db
from covid_profiler_web.worker import profile
import uuid
from werkzeug.utils import secure_filename
import os
from flask import current_app as app

import re
bp = Blueprint('upload', __name__)

def run_sample(db,user_id,uniq_id,sample_name,f1):
    filename1 = secure_filename(f1.filename)
    f1.save(os.path.join(app.config["UPLOAD_FOLDER"], filename1))
    server_fname1 = "/tmp/%s" % filename1
    db.execute("INSERT INTO results (id,user_id,sample_name,result) VALUES (?,?,?,?)",(uniq_id,user_id,sample_name,"{}"))
    db.commit()
    profile.delay(fasta=server_fname1,uniq_id=uniq_id,db=app.config["DATABASE"],storage_dir=app.config["UPLOAD_FOLDER"])

@bp.route('/upload',methods=('GET', 'POST'))
def upload():
    db = get_db()
    if request.method == 'POST':
        print(request.form)
        error=None
        # username = g.user['username'] if g.user else 'private'
        username = 'private'
        if "single_sample_submit" in request.form:
            uniq_id = str(uuid.uuid4())
            if "sample_name" in request.form:
                sample_name = request.form["sample_name"] if request.form["sample_name"]!="" else uniq_id
            else:
                sample_name = uniq_id
            if request.files['file1'].filename=="":
                error = "No file found for read 1, please try again!"
            if error==None:
                run_sample(db,username,uniq_id,sample_name,request.files['file1'])
                return redirect(url_for('results.run_result', sample_id=uniq_id))
        flash(error)
    return render_template('upload/upload.html')
