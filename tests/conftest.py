import datetime
import json
from pathlib import Path

import pytest

from tilepy.include.CampaignDefinition import ObservationParameters


def load_test_cases():
    """Reads test cases from a JSON file."""
    data_path = Path(__file__).resolve().parent / "test_cases.json"
    with data_path.open(encoding="utf-8") as f:
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
        test_case.get("igrfcoeffs", None),
    )

    return obspar
