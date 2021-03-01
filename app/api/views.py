from flask import jsonify

from . import api


@api.route("/", methods=["GET"])
def hello_world_json():
    return jsonify({"Hello": "world!"})
