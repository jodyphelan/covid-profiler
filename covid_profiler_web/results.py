from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for, Response
)
from werkzeug.exceptions import abort
import json
# from covid_profiler_web.auth import login_required
from covid_profiler_web.db import get_db
import tbprofiler as tbp
bp = Blueprint('results', __name__)
import os
from flask import current_app as app


@bp.route('/results/<uuid:sample_id>',methods=('GET', 'POST'))
def run_result(sample_id):
    db = get_db()

    tmp = db.execute("SELECT * FROM results WHERE id = ?", (str(sample_id),) ).fetchone()
    if tmp == None:
        error = "Run does not exist"
        abort(404)
    run = dict(tmp)
    run["result"] = json.loads(tmp["result"])
    tmp = db.execute("SELECT * FROM tree ORDER BY id DESC LIMIT 1" ).fetchone()
    meta = {}
    for row in db.execute("SELECT * FROM tree_data" ).fetchall():
        meta[row["id"]] = dict(row)
    tree = {"newick":tmp["newick"], "created":tmp["created"], "meta": json.dumps(meta)}
    return render_template('results/run_result.html',run=run, tree = tree)


@bp.route('/mutations/position/<int:position>')
def mutation(position):
    db = get_db()

    tmp = db.execute("SELECT * FROM mutations WHERE position = ?", (position,)).fetchone()
    if tmp == None:
        error = "Mutation does not exist"
        abort(404)
    mutation = dict(tmp)
    print(mutation)
    tmp = db.execute("SELECT * FROM tree ORDER BY id DESC LIMIT 1" ).fetchone()
    meta = {}
    for row in db.execute("SELECT * FROM tree_data" ).fetchall():
        meta[row["id"]] = dict(row)
    tree = {"newick":tmp["newick"], "created":tmp["created"], "meta": json.dumps(meta)}
    tmp = json.dumps(mutation)
    mutation["json_string"] = tmp
    return render_template('results/mutations.html',mutation=mutation, tree = tree)

@bp.route('/mutations/')
def mutation_table():
    db = get_db()
    mutations = db.execute("SELECT * FROM mutations").fetchall()
    return render_template('results/mutation_table.html',mutations=mutations)


@bp.route('/sra',methods=('GET', 'POST'))
def sra():
    return result_table(request,"public")


def result_table(request,user):
    db = get_db()
    if request.method=='POST':
        if "search_strains_button" in request.form:
            sql_query = "select id,sample_name,created,status,lineage,drtype from results WHERE user_id = '%s'" % user
            filters = []
            key_values = list(request.form.lists())
            for key,values in list(request.form.lists()):
                if values==[""]:
                    continue
                elif key=="sample_name":
                    filters.append("( %s )" % (" OR ".join(["sample_name = '%s'" % (run_id.strip()) for run_id in values[0].split(",")])))
                elif key=="project_id":
                    pass
                    # filters.append("( %s )" % (" OR ".join(["project_id = '%s'" % (run_id.strip()) for run_id in values[0].split(",")])))
                elif key=="drtype":
                    filters.append("( %s )" % (" OR ".join(["drtype = '%s'" % (drtype) for drtype in values])))
                elif key=="lineage":
                    filters.append("( %s )" % (" OR ".join(["lineage LIKE 'lineage%s%%'" % (lineage.strip().replace("lineage","")) for lineage in values[0].split(",")])))
                else:
                    pass

            if len(filters)>0:
                sql_query = sql_query + " AND %s" % (" AND ".join(filters))
            tmp = db.execute(sql_query).fetchall()
            return render_template('results/result_table.html',results = tmp, user=user)
        else:
            # print(request.form)
            if request.form["button"]=="download":
                ids = list(json.loads(request.form["ids"]).keys())
                cmd = "select * from full_results where id in ( %s )" % ", ".join(["'%s'" % x for x in ids])
                data = db.execute(cmd).fetchall()
                fieldnames = [x["name"] for x in  db.execute("PRAGMA table_info(full_results)").fetchall()]
                csv_text = ",".join(fieldnames) + "\n"
                for row in data:
                    csv_text = csv_text + ",".join(['"%s"' % row[c] if (row[c]!=None and row[c]!="") else '"-"' for c in fieldnames]) + "\n"
                return Response(csv_text,mimetype="text/csv",headers={"Content-disposition": "attachment; filename=result.csv"})
            elif request.form["button"]=="delete":
                ids = list(json.loads(request.form["ids"]).keys())
                cmd = "DELETE FROM results WHERE id in ( %s )" % ", ".join(["'%s'" % x for x in ids])
                db.execute(cmd)
                db.commit()
    return render_template('results/result_table.html', user=user)
