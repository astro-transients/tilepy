# tests/test_scheduler.py
import os
import shutil

import pytest

from tilepy.include.ObservationScheduler import GetSchedule


@pytest.fixture
def setup_and_teardown(parsed_obs_parameters):
    """Fixture to setup and teardown output directory"""
    # Ensure the output directory exists and is clean
    if os.path.exists(parsed_obs_parameters.outDir):
        # Remove the entire output directory if it exists
        shutil.rmtree(parsed_obs_parameters.outDir)

    # Ensure the output directory exists
    os.makedirs(parsed_obs_parameters.outDir, exist_ok=True)

    # Yield parsed_obs_parameters to the test function
    yield parsed_obs_parameters

    # Teardown: Remove the output directory or files after the test is done (optional)
    if os.path.exists(parsed_obs_parameters.outDir):
        shutil.rmtree(
            parsed_obs_parameters.outDir
        )  # Ensure directory is removed after test


def test_get_schedule(setup_and_teardown):
    """Tests GetSchedule using config from parsed_obs_parameters (defined in conftest.py)"""
    parsed_obs_parameters = (
        setup_and_teardown  # This gets the value yielded from the fixture
    )

    # Run GetSchedule with the parsed parameters
    GetSchedule(parsed_obs_parameters)

    # Check if the file has been created
    # assert os.path.exists(output_file), f"Expected output file {output_file} not found."

    # Optionally, check the content of the file
    # with open(output_file, 'r') as f:
    #     content = f.read()

    # Replace the expected_content with the actual expected content for comparison
    # expected_content = "Expected content here"  # Adjust accordingly
    # assert content == expected_content, f"Content of {output_file} does not match the expected content."
