from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for, Response
)
from werkzeug.exceptions import abort
import json
# from covid_profiler_web.auth import login_required
from covid_profiler_web.db import get_db, get_mongo_db
from covid_profiler_web.worker import profile_primer
import tbprofiler as tbp
bp = Blueprint('diagnostics', __name__)
import os
from flask import current_app as app
import uuid

@bp.route('/diagnostics',methods=('GET', 'POST'))
def diagnostics():
    mongo = get_mongo_db()
    if request.method == 'POST':
        uniq_id = str(uuid.uuid4())
        primerF = request.form["primerF"]
        primerR = request.form["primerR"]
        probe = request.form["probe"]
        save_dir = app.root_path+"/"+app.static_url_path+"/results/"
        mongo.db.primer_results.insert_one({"_id": str(uniq_id),"status": "processing","primerF":primerF, "primerR":primerR, "probe":probe})
        profile_primer.delay(primerF=primerF,primerR=primerR,probe=probe,uniq_id=uniq_id, save_dir=save_dir)
        return redirect(url_for('diagnostics.primer_result', run_id=str(uniq_id)))


    return render_template('diagnostics/diagnostics.html')


@bp.route('/diagnostics/<uuid:run_id>',methods=('GET', 'POST'))
def primer_result(run_id):
    mongo = get_mongo_db()
    parameters = mongo.db.primer_results.find_one({'_id':str(run_id)})
    print({'_id':str(run_id)})
    print(parameters)
    return render_template('diagnostics/primer_result.html',run_id=str(run_id),parameters=parameters)
