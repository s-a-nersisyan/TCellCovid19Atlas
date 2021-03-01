from . import frontend


@frontend.route("/", methods=["GET"])
def hello_world():
    return "Hello world!"
