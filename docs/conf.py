# Unless otherwise indicated, all code in this project is licensed under the two-clause BSD license.
# Copyright (c) 2007-2025 by the Sphinx team (see AUTHORS file).
# All rights reserved.
# See full license details at: https://github.com/sphinx-doc/sphinx/blob/master/LICENSE.rst


import os
import re
import sys
from importlib import import_module

import astropy_sphinx_theme  # required for theme detection by Sphinx, even if unused  # noqa

DEFAULT_PROJECT = "tilepy"
DEFAULT_AUTHOR = "Halim, Monica, Fabian"
DEFAULT_RELEASE = "1.0"

sys.path.insert(0, os.path.abspath("../tilepy/include"))

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

try:
    try:
        import tomllib  # Python 3.11+
    except ImportError:
        import tomli as tomllib  # pip install tomli for <3.11

    pyproject_path = os.path.join(os.path.dirname("__file__"), "pyproject.toml")
    with open(pyproject_path, "rb") as f:
        meta = tomllib.load(f).get("project", {})

    project = meta.get("name", DEFAULT_PROJECT)
    authors = meta.get("authors", [{"name": DEFAULT_AUTHOR}])
    author = authors[0].get("name", DEFAULT_AUTHOR) if authors else DEFAULT_AUTHOR

except Exception:
    project = DEFAULT_PROJECT
    author = DEFAULT_AUTHOR

# Try to get the version from the installed package (setuptools_scm)

try:
    import_module(project)
    package = sys.modules[project]
    release = package.__version__

    # Use a regex to extract just the short X.Y version. X.Y.Z part (e.g., 2.7.6)
    match = re.match(r"^(\d+\.\d+\.\d+)", release)
    version = match.group(1) if match else release

except Exception:
    release = DEFAULT_RELEASE
    version = DEFAULT_RELEASE

# Only include dev docs in dev version.
dev = "dev" in release

copyright = f"2023, {author}"


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
html_title = "TilePy"

rst_epilog = r"""
.. |TilePy| replace:: *Tilepy*
.. |TilepyGitHub| replace:: `tilepy <https://github.com/astro-transients/tilepy>`__
.. |TilepyDocs| replace:: `tilepy documentation <https://readthedocs.org/projects/tilepy/badge/?version=implement-readthedocs>`__
"""

# Extensions for Sphinx
extensions = [
    "sphinx.ext.autodoc",  # Automatic documentation for Python objects (functions, classes, etc.)
    "sphinx.ext.napoleon",  # Support for Google-style and NumPy-style docstrings
    "sphinx.ext.intersphinx",  # Cross-references to the docs of other projects
    "sphinx.ext.todo",  # Use .. todo:: directives in your docs
    "sphinx.ext.viewcode",  # Add links to highlighted source code in the API docs
    "sphinxcontrib.bibtex",  # Scientific references (for citations)
    "sphinx.ext.autosummary",  # Automatically generate summaries for modules, functions, classes, etc.
    "sphinx.ext.extlinks",  # Create custom external links (e.g., to GitHub issues or web pages)
    "sphinx.ext.coverage",  # Display code coverage information in the documentation
    "sphinx.ext.inheritance_diagram",  # Generate inheritance diagrams for classes
    "sphinx.ext.graphviz",  # Integrate Graphviz diagrams and graphs into the documentation
    # other extensions...
]

# # Intersphinx mappings for cross-referencing other Python library documentation
# intersphinx_mapping = {
#     "numpy": ("https://numpy.org/doc/stable/", None),
#     "astropy": ("https://docs.astropy.org/en/stable/", None),
# }


todo_include_todos = True  # Show TODO items in the built documentation

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
templates_path = ["_templates"]

# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# Theme and Customizations
# html_theme = "alabaster"
html_theme = "bootstrap-astropy"
html_static_path = ["_static"]

html_theme_options = {
    "github_url": "https://github.com/weizmannk/tilepy",
    "use_edit_page_button": True,
}


# Documentation site title (optional; defaults to "<project> v<release> documentation")
html_title = "TilePy"

htmlhelp_basename = project + "doc"

# Prefixes that are ignored for sorting the Python module index
modindex_common_prefix = ["tilepy."]


html_context = {
    "default_mode": "light",
    "to_be_indexed": ["stable", "latest"],
    "is_development": dev,
    "github_user": "weizmannk",
    "github_repo": "tilepy",
    "github_version": "implement-readthedocs",
    "doc_path": "docs",
}

# -- Options for plot_directive -----------------------------------------------
plot_include_source = True
plot_formats = [("svg", 300), ("pdf", 300)]
plot_html_show_formats = False


# -- Options for the edit_on_github extension ---------------------------------
edit_on_github_project = "weizmannk/tilepy"
edit_on_github_branch = "implement-readthedocs"

edit_on_github_source_root = ""
edit_on_github_doc_root = "docs"

# Generate the URL for editing on GitHub
edit_on_github_url = f"https://github.com/{edit_on_github_project}/edit/{edit_on_github_branch}/"  # Link to GitHub editor


# -- API documentation options -----------------------------------------------
autodoc_typehints = "description"  # Show type hints in the API documentation (optional)

# -- BibTeX for scientific references ----------------------------------------
bibtex_bibfiles = ["refs.bib"]  # For bibliography
bibtex_default_style = "plain"
