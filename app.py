from flask import Flask
from config import Config
from flask_sqlalchemy import SQLAlchemy
from flask_login import LoginManager
import os
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

logger.info("This is an info message")
logger.error("This is an error message")

print("Hello world", flush=True)
print("Hello world", flush=True)
app = Flask(__name__)
app.config.from_object(Config)
app.config['UPLOAD_FOLDER'] = 'uploads'  # Temporary upload folder
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)  # Ensure folder exists

db = SQLAlchemy(app)
login_manager = LoginManager(app)
login_manager.login_view = 'login'

from routes import *

if __name__ == '__main__':
    app.run(debug=True)
