import os

from flask import Flask


def create_app(test_config=None):
    # create and configure the app
    app = Flask(__name__, instance_relative_config=True)
    app.config.from_mapping(
        SECRET_KEY='dev',
        DATABASE=os.path.join(app.instance_path, 'flaskr.sqlite'),
        UPLOAD_FOLDER="/tmp",
        APP_ROOT=os.path.dirname(os.path.abspath(__file__)),
        MONGO_URI = "mongodb://localhost:27017/test_database"
    )

    if test_config is None:
        # load the instance config, if it exists, when not testing
        app.config.from_pyfile('config.py', silent=True)
    else:
        # load the test config if passed in
        app.config.from_mapping(test_config)

    # ensure the instance folder exists
    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass

    from . import db
    db.init_app(app)

    from . import home
    app.register_blueprint(home.bp)

    from . import results
    app.register_blueprint(results.bp)

    from . import analyse
    app.register_blueprint(analyse.bp)

    from . import diagnostics
    app.register_blueprint(diagnostics.bp)

    # from . import auth
    # app.register_blueprint(auth.bp)
    #
    # from . import users
    # app.register_blueprint(users.bp)
    #
    # from . import mutations
    # app.register_blueprint(mutations.bp)

    return app
