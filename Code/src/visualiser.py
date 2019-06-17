import numpy as np
import matplotlib
import networkx as nx
import matplotlib.pyplot as plt
import os.path

"""
This class defines methods to visualise project components in an intuitive
way. This includes:
    - networks
    - traffic assignment matrices
    - flows

For some very nice examples of matrix visualisation on how to create annonated
heatmaps (helper functions in this code are taken from it):
https://matplotlib.org/gallery/images_contours_and_fields/image_annotated_heatmap.html#sphx-glr-gallery-images-contours-and-fields-image-annotated-heatmap-py
"""

# Some helper functions for matrix visualisation
def heatmap(data, row_labels, col_labels, 
            ax=None, title=None, xlabel=None, ylabel=None, 
            show_cbar=True, cbar_kw={}, cbarlabel="", **kwargs
    ):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Arguments:
        data       : A 2D numpy array of shape (N,M)
        row_labels : A list or array of length N with the labels for the rows
        col_labels : A list or array of length M with the labels for the columns
    Optional arguments:
        ax         : A matplotlib.axes.Axes instance to which the heatmap
                    is plotted. If not provided, use current axes or
                    create a new one.
        title      : Figure's title.
        xlabel     : X-axis title.
        ylabel     : Y-axis title.
        show_cbar  : Whether to show colorbar.
        cbar_kw    : A dictionary with arguments to
                    :meth:`matplotlib.Figure.colorbar`.
        cbarlabel  : The label for the colorbar
    All other arguments are directly passed on to the imshow call.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar if required
    if show_cbar:
        cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
        cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    if xlabel:
        ax.set_xlabel(xlabel, labelpad=8)
    ax.set_yticklabels(row_labels)
    if ylabel:
        ax.set_ylabel(ylabel, labelpad=8)

    # Figure's title
    if title:
        ax.set_title(title)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=False, bottom=True, labeltop=False, labelbottom=True)

    # Rotate the tick labels and set their alignment.
    # plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
    #         rotation_mode="anchor")
    plt.setp(ax.get_xticklabels(), rotation=45, ha='center')#, position=(0,-0.1))

    # Turn spines off and create white grid.
    for _, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                    textcolors=["black", "white"],
                    threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Arguments:
        im         : The AxesImage to be labeled.
    Optional arguments:
        data       : Data used to annotate. If None, the image's data is used.
        valfmt     : The format of the annotations inside the heatmap.
                    This should either use the string format method, e.g.
                    "$ {x:.2f}", or be a :class:`matplotlib.ticker.Formatter`.
        textcolors : A list or array of two color specifications. The first is
                    used for values below a threshold, the second for those
                    above.
        threshold  : Value in data units according to which the colors from
                    textcolors are applied. If None (the default) uses the
                    middle of the colormap as separation.

    Further arguments are passed on to the created text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center", verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[im.norm(data[i, j]) > threshold])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts

class Vis():
    """
    Visualiser class
    """

    net = None
    last_fig = None # Last figure plotted

    def __init__(self, net):
        self.net = net

    def visualise_network(self, title=None, figsize=None):
        # Calculate grid dimensions
        # TODO: Include title and other options
        # TODO: Better use of network properties to determine figures' size
        if not figsize:
            side = min(len(self.net.nodes),10)
            figsize = (side,side)
        self.last_fig, _ = plt.subplots(figsize=figsize)
        nx.draw_kamada_kawai(self.net.G, with_labels=True)

    def visualise_ss_am(self, mat='A', figsize=None, show=True):
        """Prints specified single-step assignment matrix as a heatmap.
        
        Keyword Arguments:
            mat {str} -- Assignment matrix to be visualised, 'A', 'P', or 
                'A_path' (default: {'A'})
            figsize {tuple} -- Size of figure to be saved (default: {None})
            show {bool} -- Whether to plot figure or not (defult: {True})
        """

        if mat == 'A':
            cols = self.net.od_pairs
            xlabel = 'OD pairs'
            cmap = 'Reds'
        elif mat == 'P':
            cols = self.net.origins
            xlabel = 'Origins'
            cmap = 'Reds'
        elif mat == 'A_path':
            cols = list(map(
                lambda i: 'p_{i}'.format(i=i),
                range(len(self.net.paths))
            ))
            xlabel = 'Paths'
            valfmt = '{x: d}'
            show_cbar = False
            cmap = 'binary'
            print('List of paths:')
            for (i,p) in enumerate(self.net.paths):
                print('\tp_{i} = {p}'.format(i=i, p=p))
        else:
            raise ValueError(
                '''Assignment matrix {m} is not recognised, possible options are:
                    - A (for OD flows)
                    - P (for O flows)
                    - A_ms (for path flows)'''.format(m=mat))

        vals = getattr(self.net, mat)
        title = 'Assignment matrix ' + mat
        rows = self.net.links
        valfmt = '{x: .2f}'
        show_cbar = True
        textcolors = ['black', 'white']
        if not figsize:
            w = min(len(rows),10)
            h = min(len(cols),10)
            figsize = (h,w)
        self.last_fig, ax = plt.subplots(figsize=figsize)
        im = heatmap(
            vals, rows, cols, ax=ax, title=title, xlabel=xlabel, ylabel='Links',
            show_cbar=show_cbar, cmap=cmap, cbarlabel="traffic fraction per link"
        )
        annotate_heatmap(im, valfmt=valfmt, textcolors=textcolors)

        # fig.tight_layout()
        if show:
            plt.show()

    def visualise_ms_am(self, mat='A_ms', figsize=None, show=True):
        """Prints specified multi-step assignment matrix as a heatmap.
        
        Keyword Arguments:
            mat {str} -- Assignment matrix to be visualised (default: {'A_ms'})
            figsize {tuple} -- Size of figure to be saved (default: {None})
            show {bool} -- Whether to plot figure or not (default: {True})
        """
        rows = self.net.links
        textcolors = ['black', 'white']
        valfmt = '{x: .2f}'
        if mat == 'A_ms':
            if not figsize:
                w = len(self.net.A_ms) * min(10, len(self.net.links))
                h = min(10, len(self.net.od_pairs)+2)
                figsize = (h,w)
            self.last_fig, axs = plt.subplots(len(self.net.A_ms), 1, figsize=figsize)

            cols = self.net.od_pairs
            show_cbar = False
            cmap = 'Reds'
            for t in range(len(self.net.A_ms)):
                # if t > 0:
                #     show_cbar = False
                vals = self.net.A_ms[t]
                title = 'Assignment matrix A_ms[{t}]'.format(t=t)
                im = heatmap(
                    vals, rows, cols, ax=axs[t], title=title,
                    show_cbar=show_cbar, cmap=cmap, cbarlabel="traffic fraction per link"
                )
                annotate_heatmap(im, valfmt=valfmt, textcolors=textcolors)
            # TODO: take care of xlabel and ylabel
        elif mat == 'P_ms':
            if not figsize:
                w = len(self.net.P_ms) * min(10, len(self.net.links))
                h = min(10, len(self.net.origins)+2)
                figsize = (h,w)
            self.last_fig, axs = plt.subplots(len(self.net.P_ms), 1, figsize=figsize)

            cols = self.net.origins
            show_cbar = False
            cmap = 'Reds'
            for t in range(len(self.net.P_ms)):
                # if t > 0:
                #     show_cbar = False
                vals = self.net.P_ms[t]
                title = 'Assignment matrix A_ms[{t}]'.format(t=t)
                im = heatmap(
                    vals, rows, cols, ax=axs[t], title=title,
                    show_cbar=show_cbar, cmap=cmap, cbarlabel="traffic fraction per link"
                )
                annotate_heatmap(im, valfmt=valfmt, textcolors=textcolors)
            # TODO: take care of xlabel and ylabel
        if show:
            self.last_fig.show()



    def visualise_flows(self, flows='link'):
        pass

    def save_figure(self, fname):
        """Save last figure plotted with this instance.
        
        Arguments:
            fname {str} -- String containing path to filename
        """
        # Avoid overwriting
        file_exists = os.path.isfile(fname)
        i=1
        next_fname = fname + '({i})'.format(i=i)
        while file_exists:
            fname = next_fname
            file_exists = os.path.isfile(fname)
            lbrack = next_fname.rfind('(')
            rbrack = next_fname.rfind(')')
            next_fname = next_fname[:lbrack+1] + str(i) + next_fname[rbrack:]
            i += 1
        self.last_fig.save_figure(fname, format='pdf')