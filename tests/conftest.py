import datetime
import json

import pytest

from tilepy.include.CampaignDefinition import ObservationParameters


def load_test_cases():
    """Reads test cases from a JSON file."""
    with open("tests/test_cases.json", "r") as f:
        return json.load(f)  # Returns a list of dictionaries


@pytest.fixture(scope="module")
def all_cases():
    """Fixture that provides all test cases as a list."""
    return load_test_cases()


@pytest.fixture(scope="module", params=load_test_cases())
def parsed_obs_parameters(request):
    """
    Parses ObservationParameters from the ini file and adds additional arguments.
    """
    test_case = request.param

    # Initialize ObservationParameters
    obspar = ObservationParameters()

    # Load the config (.ini) file
    obspar.from_configfile(test_case["cfgFile"])

    # Add other parsed arguments from JSON
    obspar.add_parsed_args(
        test_case["skymap"],
        datetime.datetime.fromisoformat(test_case["obsTime"]),
        test_case["datasetDir"],
        test_case["galcatName"],
        test_case["outDir"],
        test_case["pointingsFile"],
    )

    return obspar
