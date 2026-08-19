import datetime
import json
from pathlib import Path

import pytest
import requests

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


@pytest.fixture
def skip_if_skymap_unreachable(parsed_obs_parameters):
    """Skip the test if the skymap URL can't be fetched, instead of failing
    later with a confusing error once the download degrades downstream results."""
    url = parsed_obs_parameters.skymap
    try:
        response = requests.get(url, timeout=15, stream=True)
        response.close()
    except requests.exceptions.RequestException as exc:
        pytest.skip(f"Skymap URL unreachable, skipping test: {url} ({exc})")
    if response.status_code >= 400:
        pytest.skip(
            f"Skymap URL returned status {response.status_code}, skipping test: {url}"
        )
