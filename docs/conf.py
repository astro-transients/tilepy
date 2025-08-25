# -*- coding: utf-8 -*-
# Copyright (C) 2024, tilepy developers
# Licensed under the GNU license - see ../LICENSE.rst
#
# Inspired by Sphinx_Astropy, Sphinx, m4opt.

# FIXME: ImportError: numpy.core.multiarray failed to import
# Une a numpy<2 for 'healpy==1.16.2',
# SO we need to fix it in the pyprojet.toml

import os
import re
import sys
from importlib import import_module

import astropy_sphinx_theme  # required for theme detection by Sphinx, even if unused  # noqa

# docs/_pybtex/short_alpha.py:
sys.path.insert(0, os.path.abspath("_pybtex"))
import short_alpha  # noqa: F401, E402

try:
    from sphinx_astropy.conf.v2 import *  # noqa
except ImportError:
    print(
        "ERROR: the documentation requires the sphinx-astropy package to be installed"
    )
    sys.exit(1)

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
html_title = "tilepy"

rst_epilog = r"""
.. |tilepy| replace:: *tilepy*

.. |tilepyDocs| replace:: `tilepy documentation <https://readthedocs.org/projects/tilepy/badge/?version=latest>`__

.. |tilepyGitHub| image:: https://img.shields.io/badge/GitHub-tilepy-9400D3?logo=github
   :target: https://github.com/astro-transients/tilepy
   :alt: tilepy on GitHub

.. |tilepyEmail| image:: https://img.shields.io/badge/Email-astro.tilepy@gmail.com-0078D4?logo=gmail
   :target: mailto:astro.tilepy@gmail.com
   :alt: Contact tilepy Team

.. |Forum| image:: https://img.shields.io/badge/Forum-Colibri-4B286D?logo=discourse
   :target: https://forum.astro-colibri.science/c/instrumentation-and-tools/tilepy
   :alt: tilepy Forum

.. |email| image:: /_static/email.svg
   :width: 16px
   :align: middle
"""

highlight_language = "python3"


# --- Additional/local Sphinx extensions ---
# Add only those extensions that are NOT already included in sphinx_astropy.conf.v2

try:
    extensions  # noqa: F405
except NameError:
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
        "sphinx_automodapi.automodapi",
        "sphinxcontrib.jquery",
        "sphinx_copybutton",
        "sphinx_astropy.ext.intersphinx_toggle",
        "sphinx_automodapi.smart_resolver",
        "pytest_doctestplus.sphinx.doctestplus",
        "matplotlib.sphinxext.plot_directive",
        "sphinx_astropy.ext.generate_config",
        "sphinx_astropy.ext.missing_static",
        "sphinx_astropy.ext.changelog_links",
        # other extensions...
    ]

# Remove 'numpydoc' if present (we want to use 'napoleon' instead)
if "numpydoc" in extensions:
    extensions.remove("numpydoc")

additional_extensions = [
    "sphinx.ext.napoleon",
    "sphinxcontrib.typer",
    "sphinxcontrib.bibtex",
    "sphinx_astropy.ext.edit_on_github",
    "sphinx.ext.mathjax",
    "nbsphinx",  # To display jupyter notebook
    "sphinx_gallery.gen_gallery",
    "sphinx_design",
]
for ext in additional_extensions:
    if ext not in extensions:
        extensions.append(ext)

# Activate TODO
todo_include_todos = True


# -- Options for nbsphinx -----------------------------------------------------

examples_dirs = ["./tutorials"]
gallery_dirs = ["auto_tutorials"]

sphinx_gallery_conf = {
    "examples_dirs": examples_dirs,
    "gallery_dirs": gallery_dirs,
    "capture_repr": ("_repr_html_", "__repr__"),
}

# Intersphinx mappings for cross-referencing other Python library documentation
intersphinx_mapping = {
    # "numpy": ("https://numpy.org/doc/stable/", None),
    "astropy": ("https://docs.astropy.org/en/stable/", None),
    "ligo.skymap": ("https://lscsoft.docs.ligo.org/ligo.skymap/", None),
    "userguide": ("https://emfollow.docs.ligo.org/userguide/", None),
}

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "_templates"]


# Documentation site title (optional; defaults to "<project> v<release> documentation")
html_title = "tilepy"

htmlhelp_basename = project + "doc"

# Prefixes that are ignored for sorting the Python module index
modindex_common_prefix = ["tilepy."]

html_context = {
    "default_mode": "light",
    "to_be_indexed": ["stable", "latest"],
    "is_development": dev,
    "github_user": "astro-transients",
    "github_repo": "tilepy",
    "github_version": "master",
    "doc_path": "docs",
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# Theme and Customizations
try:
    html_theme  # noqa: F405
except NameError:
    html_theme = "pydata_sphinx_theme"


html_static_path = ["_static"]

html_theme_options = {
    "github_url": "https://github.com/astro-transients/tilepy",
    "use_edit_page_button": True,
}

# Ignore issue/PR links; set retries/timeouts for external link checking
linkcheck_ignore = [
    r"https://github\.com/astro-transients/tilepy/(?:issues|pull)/\d+",
]
linkcheck_retry = 5
linkcheck_timeout = 180
linkcheck_anchors = False

# -- Options for the edit_on_github extension ---------------------------------
edit_on_github_project = "astro-transients/tilepy"
edit_on_github_branch = "master"

edit_on_github_source_root = ""
edit_on_github_doc_root = "docs"

# Generate the URL for editing on GitHub
edit_on_github_url = f"https://github.com/{edit_on_github_project}/edit/{edit_on_github_branch}/"  # Link to GitHub editor


# If we want to use latex format
latex_documents = [
    ("index", project + ".tex", project + " Documentation", author, "manual")
]

# -- Options for plot_directive -----------------------------------------------
plot_include_source = True
plot_formats = [("svg", 300), ("pdf", 300)]
plot_html_show_formats = False

# -- API documentation options -----------------------------------------------
autodoc_typehints = "description"  # Show type hints in the API documentation (optional)

# -- BibTeX for scientific references ----------------------------------------
bibtex_bibfiles = ["refs.bib"]  # For bibliography
bibtex_default_style = "short_alpha"
# bibtex_default_style = "plain"
# bibtex_author_limit = 3
