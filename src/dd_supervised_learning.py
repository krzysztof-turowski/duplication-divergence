import matplotlib.pyplot as pyplot

import dd_plot

PLOT_STYLE = 'ggplot'
FIGURE_SIZE_SCALE = 0.5

alpha = [0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1]
arxiv_sort_by_p_uv_sum_precision = [
    0.553, 0.603, 0.704, 0.808, 0.889, 0.953, 0.976, 0.989, 0.995, 0.998]
arxiv_p_uv_threshold_precision = [
    0.551, 0.597, 0.710, 0.808, 0.892, 0.954, 0.977, 0.988, 0.996, 0.998]
simplewiki_sort_by_p_uv_sum_precision = [
    0.572, 0.631, 0.748, 0.845, 0.914, 0.964, 0.982, 0.991, 0.997, 0.998]
simplewiki_p_uv_threshold_precision = [
    0.571, 0.628, 0.747, 0.844, 0.916, 0.964, 0.982, 0.991, 0.997, 0.998]
collegemsg_sort_by_p_uv_sum_precision = [
    0.513, 0.519, 0.561, 0.625, 0.713, 0.842, 0.911, 0.955, 0.982, 0.992]
collegemsg_p_uv_threshold_precision = [
    0.517, 0.541, 0.574, 0.625, 0.712, 0.840, 0.916, 0.953, 0.983, 0.992]

def plot_data(A, B, filename):
    dd_plot.initialize_figure(PLOT_STYLE, FIGURE_SIZE_SCALE)
    pyplot.ylabel(r'$\theta$')
    pyplot.xlabel(r'$\alpha$')
    pyplot.ylim(0.5, 1.0)
    pyplot.semilogx(
        alpha, A, color = 'red', marker = None, alpha = 0.8, linestyle = '-',
        label = r'\texttt{sort-by-$p_{uv}$-sum}, $|C| = 1$')
    pyplot.semilogx(
        alpha, B, color = 'blue', marker = None, alpha = 0.8, linestyle = ':',
        label = r'\texttt{$p_{uv}$-threshold}, $\tau = 0.5$')
    pyplot.legend(loc = 'lower right')
    dd_plot.plot(filename, 'pdf')
    pyplot.close()

plot_data(arxiv_p_uv_threshold_precision, arxiv_sort_by_p_uv_sum_precision, 'arxiv')
plot_data(simplewiki_p_uv_threshold_precision, simplewiki_sort_by_p_uv_sum_precision, 'simplewiki')
plot_data(collegemsg_p_uv_threshold_precision, collegemsg_sort_by_p_uv_sum_precision, 'collegemsg')
