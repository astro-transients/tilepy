# Reproduce form
# We then need add "pybtex" in the ../pyproject.toml, to be installed by sphinx

"""Pybtex style to display all authors in a reference.
Reproduced from
https://github.com/m4opt/m4opt/blob/main/m4opt/utils/pybtex/styles/short_alpha.py

"""

from pybtex.plugin import register_plugin
from pybtex.style.formatting.alpha import Style as AlphaStyle
from pybtex.style.template import FieldIsMissing, join, node, sentence


@node
def full_names(children, context, role):
    """Return formatted names."""

    assert not children

    try:
        persons = context["entry"].persons[role]
    except KeyError:
        raise FieldIsMissing(role, context["entry"])

    style = context["style"]
    formatted_names = [
        style.format_name(person, style.abbreviate_names) for person in persons
    ]
    return join(sep=", ", sep2=" and ", last_sep=", and ")[formatted_names].format_data(
        context
    )


class ShortAlphaStyle(AlphaStyle):
    def format_names(self, role, as_sentence=True):
        formatted_names = full_names(role)
        if as_sentence:
            return sentence[formatted_names]
        else:
            return formatted_names


register_plugin("pybtex.style.formatting", "short_alpha", ShortAlphaStyle)
