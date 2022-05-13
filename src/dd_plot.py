import argparse

import matplotlib
import matplotlib.pyplot as pyplot
import numpy

# TEXT_WIDTH_PT = 506.295 # Get this from LaTeX using \the\textwidth
TEXT_WIDTH_PT = 656.295 # Get this from LaTeX using \the\textwidth
TEMP_FOLDER = 'temp'

def figure_size(figure_size_scale):
    inches_per_pt = 1.0 / 72.27 # Convert pt to inch
    # golden_mean = (numpy.sqrt(5.0) - 1.0) / 2.0 # Aesthetic ratio (you could change this)
    golden_mean = 9/21
    figure_width = TEXT_WIDTH_PT * inches_per_pt * figure_size_scale
    figure_height = figure_width * golden_mean
    return [figure_width, figure_height]

def initialize_figure(plot_style, figure_size_scale):
    publication_with_latex = {
        "pgf.texsystem": "pdflatex", # change this if using xetex or lautex
        "text.usetex": True, # use LaTeX to write all text
        "font.family": "serif",
        "axes.labelsize": 8, # LaTeX default is 10pt font.
        "font.size": 8,
        "legend.fontsize": 8, # Make the legend/label fonts a little smaller
        "savefig.dpi": 125,
        "text.latex.preamble": r"\usepackage{amsmath,amssymb,amsfonts}",
        "figure.figsize": figure_size(figure_size_scale)
    }
    pyplot.style.use(plot_style)
    matplotlib.rcParams.update(publication_with_latex)

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', metavar = 'FILE', help = 'data to plot')
    parser.add_argument('--export', choices = ['pdf', 'png', 'eps'], help = 'export plot to file')
    return parser

def plot(name, export):
    if export == 'pdf':
        pyplot.savefig(r'{0}/{1}.pdf'.format(TEMP_FOLDER, name),
                       bbox_inches = 'tight', pad_inches = 0.05)
    elif export == 'eps':
        pyplot.savefig(r'{0}/{1}.eps'.format(TEMP_FOLDER, name),
                       bbox_inches = 'tight', pad_inches = 0.05)
    elif export == 'png':
        pyplot.savefig(r'{0}/{1}.png'.format(TEMP_FOLDER, name),
                       bbox_inches = 'tight', pad_inches = 0.05)
    else:
        pyplot.show()
