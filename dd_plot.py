import matplotlib
import matplotlib.pyplot as pyplot
import numpy

TEXT_WIDTH_PT = 506.295 # Get this from LaTeX using \the\textwidth

def figure_size(figure_size_scale):
  inches_per_pt = 1.0 / 72.27 # Convert pt to inch
  golden_mean = (numpy.sqrt(5.0) - 1.0) / 2.0 # Aesthetic ratio (you could change this)
  figure_width = TEXT_WIDTH_PT * inches_per_pt * figure_size_scale
  figure_height = figure_width * golden_mean
  return [figure_width, figure_height]

def initialize_figure(plot_style, figure_size_scale):
  publication_with_latex = {
    "pgf.texsystem": "pdflatex", # change this if using xetex or lautex
    # "text.usetex": True, # use LaTeX to write all text
    "font.family": "serif",
    "font.serif": [], # blank entries should cause plots to inherit fonts from the document
    "font.sans-serif": [],
    "font.monospace": [],
    "axes.labelsize": 8, # LaTeX default is 10pt font.
    "font.size": 8,
    "legend.fontsize": 8, # Make the legend/label fonts a little smaller
  }
  pyplot.style.use(plot_style)
  matplotlib.rcParams.update(publication_with_latex)
  matplotlib.rcParams['savefig.dpi'] = 125
  matplotlib.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath,amssymb,amsfonts}"]
  matplotlib.rcParams['figure.figsize'] = figure_size(figure_size_scale)
