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
    tmp = db.execute("SELECT * FROM tree ORDER BY id DESC LIMIT 1" ).fetchone()
    meta = {}
    for row in db.execute("SELECT * FROM tree_data" ).fetchall():
        meta[row["id"]] = dict(row)
    tree = {"newick":tmp["newick"], "created":tmp["created"], "meta": json.dumps(meta)}

    if request.method == 'POST':
        return redirect(url_for('results.run_result', sample_id=request.form["sample_id"]))
    return render_template('home/index.html',tree = tree)
