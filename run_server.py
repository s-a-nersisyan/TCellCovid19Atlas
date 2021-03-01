import sys

from app import app

if __name__ == "__main__":
    port = sys.argv[1]
    app.run(host="0.0.0.0", port=port, debug=True)
