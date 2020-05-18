from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for, current_app
)
from werkzeug.exceptions import abort
import json
#from aseantb.auth import login_required
from covid_profiler_web.db import get_db, get_mongo_db

bp = Blueprint('home', __name__)


@bp.route('/',methods=('GET', 'POST'))
def index():
    db = get_db()
    mongo = get_mongo_db()
    tree_text = mongo.db.tree.find_one()["tree"]
    meta = mongo.db.meta.find_one()
    del meta["_id"]
    meta = json.dumps(meta)

    tree = {"newick":tree_text, "created":"NA", "meta": meta}


    if request.method == 'POST':
        return redirect(url_for('results.run_result', sample_id=request.form["sample_id"]))
    return render_template('home/index.html',tree = tree)




@bp.route('/phylocanvas')
def phylocanvas():
    db = get_db()
    mongo = get_mongo_db()

    meta = mongo.db.meta.find_one()
    del meta["_id"]
    meta = json.dumps(meta)

    return render_template('results/phylocanvas.html',meta = meta)
