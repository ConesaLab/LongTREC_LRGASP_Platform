from app import app
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

logger.info("This is an info message")
logger.error("This is an error message")

print("Hello world", flush=True)
print("Hello world", flush=True)

application = app  # Some servers expect 'application' variable
