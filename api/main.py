import os
import json
from flask import Flask, request, jsonify, make_response
from flask_restful import Api, Resource
from flasgger import Swagger, swag_from
from werkzeug.exceptions import abort
from tconsumption import Consumption

PROJECT_ROOT = os.path.dirname(os.path.realpath(__file__))

app = Flask(__name__)
app.config['JSON_SORT_KEYS'] = False
app.url_map.strict_slashes = False

swagger_config = {
    "headers": [
    ],
    "specs": [
        {
            "endpoint": 'apispec',
            "route": '/apispec.json',
            "rule_filter": lambda rule: True,  # all in
            "model_filter": lambda tag: True,  # all in
        }
    ],
    "static_url_path": "/flasgger_static",
    "swagger_ui": True,
    "specs_route": "/apidocs"
}

api = Api(app)
swagger = Swagger(app, template_file="swag_definitions/template.yml", config=swagger_config)

def validation_error(err, data, schema):
    """
    Custom validation error handler which produces 404 Not Found
    response in case validation fails instead of 400 Bad Request
    """
    response = make_response(jsonify({
        "status": 400,
        "mistake_path": err.json_path,
        "validator": f"{err.validator} == {err.validator_value}",
        "message": err.message
    }), 400)
    abort(response)

class ConsumptionCall(Resource):
    @swag_from("swag_definitions/consumption.yml", validation=True, validation_error_handler=validation_error)
    def post(self):
        post_json = request.get_json()

        c = Consumption()
        c.load(
            post_json["coordinates"],
            post_json["station_orders"],
            post_json["velocity_ways"]
        )
        c.run()

        return {
            "force_values": c.series["force_values"],
            "exerted_force_values": c.series["exerted_force_values"],
            "dist_values": c.series["dist_values"],
            "acceleration_values": c.series["acceleration_values"],
            "velocity_values": c.series["velocity_values"]
        }

api.add_resource(ConsumptionCall, '/consumption')

if __name__ == "__main__":
    app.run(debug=True)
