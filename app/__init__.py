from flask import Flask
from flask_sqlalchemy import SQLAlchemy


# Initialize the app
app = Flask(__name__)
app.config.from_pyfile("config.py")

db = SQLAlchemy(app)

# Register blueprints
from .api import api
from .frontend import frontend

app.register_blueprint(api, url_prefix="/api")
app.register_blueprint(frontend, url_prefix="/")
