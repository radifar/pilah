# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------

project = u"pilah"
copyright = u"2024, Muhammad Radifar & Enade Perdana Istyastono"
author = u"Muhammad Radifar"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ["myst_parser", "sphinxcontrib.video", "sphinx_design"]
myst_enable_extensions = ["colon_fence"]
# autoapi_dirs = ["../src"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "venv", ".venv"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "furo"

html_static_path = ["_static"]

html_theme_options = {
    "sidebar_hide_name": True,
    "light_logo": "PiLAH_web.png",
    "dark_logo": "PiLAH_web.png",
}

