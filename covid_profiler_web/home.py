from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for, current_app
)
from werkzeug.exceptions import abort
import json
#from aseantb.auth import login_required
from covid_profiler_web.db import get_db

bp = Blueprint('home', __name__)


@bp.route('/',methods=('GET', 'POST'))
def index():
    db = get_db()
    newick = open(current_app.config["UPLOAD_FOLDER"]+"/covid_public.tree").readline().strip()
    if request.method == 'POST':
        return redirect(url_for('results.run_result', sample_id=request.form["sample_id"]))
    return render_template('home/index.html',tree = newick)
