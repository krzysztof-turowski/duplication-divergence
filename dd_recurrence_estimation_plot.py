# Tool for plotting the results from dd_recurrence_estimation.
# Run: python -B ./dd_recurrence_estimation_plot.py -B FILE [--export {pdf}]

import argparse
import matplotlib.pyplot as pyplot
import numpy
import os
import sys

import dd_plot

PLOT_STYLE = 'ggplot'
FIGURE_SIZE_SCALE = 0.5

PASTOR_SATORRAS_PATTERN = [False, True, False]
CHUNG_LU_PATTERN = [False, False, True]

def read_data(data):
  values = numpy.array([[[float(value) if value != '' else None for value in interval.split(',')]
      for interval in parameter.split(';')] for parameter in data.strip().split(' ')])
  if (numpy.equal(values[0][0], None) == PASTOR_SATORRAS_PATTERN).all():
    pyplot.xlabel(r'$p$')
    pyplot.ylabel(r'$r$', rotation = 0)
    if values.shape[1] == 3:
      return values[:, 0, 0].astype(float), values[:, 0, 2].astype(float), values[:, 1, 2].astype(float), values[:, 2, 2].astype(float)
    elif values.shape[1] == 1:
      return values[:, 0, 0].astype(float), None, values[:, 0, 2].astype(float), None
    else:
      raise Exception('Unidentified shape of data: {0}'.format(values.shape))
  elif (numpy.equal(values[0][0], None) == CHUNG_LU_PATTERN).all():
    pyplot.xlabel(r'$p$')
    pyplot.ylabel(r'$q$', rotation = 0)
    if values.shape[1] == 3:
      return values[:, 0, 0].astype(float), values[:, 0, 1].astype(float), values[:, 1, 1].astype(float), values[:, 2, 1].astype(float)
    elif values.shape[1] == 1:
      return values[:, 0, 0].astype(float), None, values[:, 0, 1].astype(float), None
    else:
      raise Exception('Unidentified shape of data: {0}'.format(values.shape))
  else:
    raise Exception('Unidentified type of data: {0}'.format(numpy.equal(values[0], None)))

def plot_parameter(X, Y_min, Y, Y_max, label, color):
  pyplot.plot(X, Y, color = color, marker = None, alpha = 0.75, label = label)
  if Y_min is not None and Y_max is not None:
    pyplot.fill_between(X, Y_min, Y_max, color = color, alpha = 0.25)

def plot_data(file, filename, export):
  name = os.path.splitext(filename.strip())[0]
  dd_plot.initialize_figure(PLOT_STYLE, FIGURE_SIZE_SCALE)
  data = file.readlines()
  plot_parameter(*read_data(data[0]), label = 'Degree', color = 'red')
  plot_parameter(*read_data(data[1]), label = 'Degree squared', color = 'magenta')
  plot_parameter(*read_data(data[2]), label = 'Open triangles', color = 'blue')
  plot_parameter(*read_data(data[3]), label = 'Triangles', color = 'green')
  pyplot.legend(loc = 'upper left')
  pyplot.gca().get_xaxis().set_major_locator(pyplot.MaxNLocator(8))
  pyplot.gca().get_yaxis().set_major_locator(pyplot.MaxNLocator(8))
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
