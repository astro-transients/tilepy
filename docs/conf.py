import os
import sys

sys.path.insert(0, os.path.abspath("../tilepy/include"))

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "tilepy"
copyright = "2023, Halim, Monica, Fabian"
author = "Halim, Monica, Fabian"
release = "1.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
# html_title = "TilePy"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",  # Support for Google-style and NumPy-style docstrings
    "sphinx.ext.intersphinx",  # Cross-references to the docs of other projects
    "sphinx.ext.todo",  # Use .. todo:: directives in your docs
    "sphinx.ext.viewcode",  # Add links to highlighted source code in the API docs
    "sphinxcontrib.bibtex",  # Scientific references (for citations)
    # other extensions...
]

# # Intersphinx mappings for cross-referencing other Python library documentation
# intersphinx_mapping = {
#     "numpy": ("https://numpy.org/doc/stable/", None),
#     "astropy": ("https://docs.astropy.org/en/stable/", None),
# }

todo_include_todos = True  # Show TODO items in the built documentation

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "alabaster"
html_static_path = ["_static"]


# Documentation site title (optional; defaults to "<project> v<release> documentation")
html_title = "tilepy Documentation"

# -- API documentation options -----------------------------------------------
autodoc_typehints = "description"  # Show type hints in the API documentation (optional)

# -- BibTeX for scientific references ----------------------------------------
bibtex_bibfiles = ["refs.bib"]  # To bibliography
