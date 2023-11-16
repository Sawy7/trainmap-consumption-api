import logging
import sys
import os
logging.basicConfig(stream=sys.stderr)
sys.path.insert(0, os.path.dirname(__file__))
from main import app as application
