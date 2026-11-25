"""Tests of the ObservationParameters constructor."""

from tilepy.include.CampaignDefinition import ObservationParameters


def test_constructor_keeps_the_values_given():
    obspar = ObservationParameters(
        shape="polygon", numberSides=6, alphaR=0.3, betaR=0.7
    )

    assert obspar.shape == "polygon"
    assert obspar.numberSides == 6
    assert obspar.alphaR == 0.3
    assert obspar.betaR == 0.7
