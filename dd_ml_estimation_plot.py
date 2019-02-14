# Tool for plotting the results from dd_ml_estimation.
# Run: python -B ./dd_ml_estimation_plot.py FILE [--export {pdf}]

import argparse
import math
import matplotlib.pyplot as pyplot
import numpy
import os

import dd_plot

PLOT_STYLE = 'default'
FIGURE_SIZE_SCALE = 0.6
X_LABELS, Y_LABELS = 11, 6

PASTOR_SATORRAS_PATTERN = [False, True, False]
CHUNG_LU_PATTERN = [False, False, True]

def read_data(data):
  values = [[float(value) if value != '' else None for value in parameter.split(',')] for parameter in data.strip().split(' ')]
  if (numpy.equal(values[0][:len(PASTOR_SATORRAS_PATTERN)], None) == PASTOR_SATORRAS_PATTERN).all():
    pyplot.xlabel(r'$p$')
    pyplot.ylabel(r'$r$', rotation = 0)
    x, _, y, z = zip(*values)
  elif (numpy.equal(values[0][:len(CHUNG_LU_PATTERN)], None) == CHUNG_LU_PATTERN).all():
    pyplot.xlabel(r'$p$')
    pyplot.ylabel(r'$q$', rotation = 0)
    x, y, _, z = zip(*values)
  else:
    raise Exception('Unidentified type of data: {0}'.format(numpy.equal(values[0], None)))

  x_range, y_range = numpy.unique(x), numpy.unique(y)
  matrix = numpy.zeros((len(y_range), len(x_range)))
  for px, py, pz in zip(x, y, z):
    matrix[len(y_range) - 1 - numpy.searchsorted(y_range, py), numpy.searchsorted(x_range, px)] = pz
  return x_range, y_range, matrix

def plot_data(file, filename, export):
  name = os.path.splitext(filename.strip())[0]
  dd_plot.initialize_figure(PLOT_STYLE, FIGURE_SIZE_SCALE)
  x, y, data = read_data(file.readlines()[0])
  image = pyplot.imshow(data, cmap = 'BuPu', interpolation = 'nearest', aspect = 'auto')
  # alternatives: pcolormesh + cmap {'YlGnBu', 'Pastel1'} + interpolation {'nearest', 'bilinear'} + aspect {'auto'}

  labels_x = numpy.linspace(round(min(x)), round(max(x)), num = X_LABELS, endpoint = True)
  pyplot.gca().set_xticks(numpy.linspace(0, len(x) - 1, num = X_LABELS, endpoint = True))
  pyplot.gca().set_xticklabels([(round(float(label), 2)) for label in labels_x])
  labels_y = numpy.flip(numpy.linspace(round(min(y)), round(max(y)), num = Y_LABELS, endpoint = True))
  pyplot.gca().set_yticks(numpy.linspace(0, len(y) - 1, num = Y_LABELS, endpoint = True))
  pyplot.gca().set_yticklabels([(round(float(label), 1)) for label in labels_y])
  pyplot.gcf().colorbar(image)

  if export == 'pdf':
    pyplot.savefig('temp/{0}.pdf'.format(name), bbox_inches = 'tight', pad_inches = 0.05)
  else:
    pyplot.show()

parser = argparse.ArgumentParser()
parser.add_argument('filename', metavar = 'FILE', help = 'data to plot')
parser.add_argument('--export', choices = ['pdf'], help = 'export plot to file')
args = parser.parse_args()
with open('temp/' + args.filename) as file:
  plot_data(file, args.filename, args.export)
