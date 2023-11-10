import os
from flask import Flask, request, jsonify, make_response
from flask_restful import Api, Resource
from flasgger import Swagger, swag_from
from werkzeug.exceptions import abort
from tconsumption import Consumption
import numpy as np

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

        # Params
        for k in post_json["params"]:
            c.params[k] = post_json["params"][k]
        c.params["power_limit"] *= 1000 # kW to W

        # Variable params
        for k in post_json["variable_params"]:
            c.variable_params[k] = post_json["variable_params"][k]
        if c.variable_params["Compensation polynomial"] is not None:
            c.variable_params["Compensation polynomial"] = np.poly1d(c.variable_params["Compensation polynomial"])

        c.load(
            post_json["rail_definition"]["coordinates"],
            post_json["rail_definition"]["station_orders"],
            post_json["rail_definition"]["velocity_ways"]
        )
        c.run()
        
        if post_json["output_options"]["energy_in_kwh"]:
            exerted_energy = [x/3600000 for x in c.series["energy_from_exerted_force"]]
        else:
            exerted_energy = c.series["energy_from_exerted_force"]

        return {
            "force_values": c.series["force_values"].tolist(),
            "exerted_force_values": c.series["exerted_force_values"].tolist(),
            "dist_values": c.series["dist_values"].tolist(),
            "acceleration_values": c.series["acceleration_values"].tolist(),
            "velocity_values": c.series["velocity_values"].tolist(),
            "exerted_energy": exerted_energy,
            "elevation_values": c.series["elevation_values"]
        }

api.add_resource(ConsumptionCall, '/consumption')

if __name__ == "__main__":
    app.run(debug=True)
