import itertools

import matplotlib as mpl
from matplotlib import gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from seaborn.matrix import (Grid, _matrix_mask, _convert_colors, dendrogram,
                            heatmap)
from seaborn.utils import despine


def merge_dicts(*dict_args):
    """Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result


class ClusterGrid(Grid):
    def __init__(self,
                 data,
                 pivot_kws=None,
                 z_score=None,
                 standard_scale=None,
                 figsize=None,
                 row_colors=None,
                 col_colors=None,
                 row_ratios=None,
                 col_ratios=None,
                 row_cluster=True,
                 col_cluster=True,
                 colorbar=True,
                 mask=None):
        """Grid object for organizing clustered heatmap input on to axes"""

        if isinstance(data, pd.DataFrame):
            self.data = data
        else:
            self.data = pd.DataFrame(data)

        self.data2d = self.format_data(self.data, pivot_kws, z_score,
                                       standard_scale)

        self.mask = _matrix_mask(self.data2d, mask)

        if figsize is None:
            width, height = 10, 10
            figsize = (width, height)

        self.row_colors, self.row_color_labels = \
            self._preprocess_colors(data, row_colors, axis=0)
        self.col_colors, self.col_color_labels = \
            self._preprocess_colors(data, col_colors, axis=1)

        self.col_cluster = col_cluster
        self.row_cluster = row_cluster

        self.dendrogram_row = None
        self.dendrogram_col = None

        self.colorbar = colorbar

        self._setup_axes(figsize, row_ratios, col_ratios)

    def _setup_axes(self, figsize, row_ratios, col_ratios):
        """Setup the different axes that make up the ClusterGrid.
           Axes that are not required under the passed arguments,
           (missing side_colors or disabled dendrograms for example)
           are ommitted and returned empty (as None).
        """

        # Setup figure.
        self.fig = plt.figure(figsize=figsize)

        # Setup axes grid.
        height_ratios = self._dim_ratios(
            figsize=figsize, axis=0, ratios=col_ratios)

        width_ratios = self._dim_ratios(
            figsize=figsize, axis=1, ratios=row_ratios)

        self.gs = gridspec.GridSpec(
            len(height_ratios),
            len(width_ratios),
            wspace=0.01,
            hspace=0.01,
            width_ratios=width_ratios,
            height_ratios=height_ratios)

        # Set dendrogram axes.
        ax_index = slice(0, 2) if self.colorbar else 0

        if self.row_cluster:
            self.ax_row_dendrogram = self.fig.add_subplot(
                self.gs[-1, ax_index], axisbg="white")
        else:
            self.ax_row_dendrogram = None

        if self.col_cluster:
            self.ax_col_dendrogram = self.fig.add_subplot(
                self.gs[ax_index, -1], axisbg="white")
        else:
            self.ax_col_dendrogram = None

        # Set heatmap axis.
        self.ax_heatmap = self.fig.add_subplot(self.gs[-1, -1])

        # Set colorbar axis.
        if self.colorbar:
            self.cax = self.fig.add_subplot(self.gs[0, 0])
        else:
            self.cax = None

        # Set color axes.
        if self.row_colors is not None:
            self.ax_row_colors = self.fig.add_subplot(self.gs[-1, -2])
        else:
            self.ax_row_colors = None

        if self.col_colors is not None:
            self.ax_col_colors = self.fig.add_subplot(self.gs[-2, -1])
        else:
            self.ax_col_colors = None

    def _dim_ratios(self, figsize, axis, ratios=None):
        """Get the proportions of the figure taken up by each axes
        """

        # Merge default ratios with any overridden ratios
        default_ratios = {'dendrogram': min(2. / figsize[axis], .2),
                          'side_colors': 0.05,
                          'heatmap': 0.8}

        ratios = merge_dicts(default_ratios, ratios or {})

        # Determine which clustering/colors to use
        if axis == 0:
            side_colors = self.col_colors
            side_cluster = self.col_cluster
        else:
            side_colors = self.row_colors
            side_cluster = self.row_cluster

        # Add the colorbar
        if self.colorbar:
            colorbar_width = .8 * ratios['dendrogram']
            colorbar_height = .2 * ratios['dendrogram']

            if axis == 0:
                dim_ratios = [colorbar_width, colorbar_height]
            else:
                dim_ratios = [colorbar_height, colorbar_width]
        else:
            dim_ratios = [ratios['dendrogram']]

        if side_colors is not None:
            # Add room for the colors
            dim_ratios += [ratios['side_colors']]

        # Add the ratio for the heatmap itself
        dim_ratios += [ratios['heatmap']]

        return dim_ratios

    def _preprocess_colors(self, data, colors, axis):
        """Preprocess {row/col}_colors to extract labels and convert colors."""
        labels = None

        if colors is not None:
            if isinstance(colors, (pd.DataFrame, pd.Series)):
                # Ensure colors match data indices
                if axis == 0:
                    colors = colors.ix[data.index]
                else:
                    colors = colors.ix[data.columns]

                # Replace na's with background color
                colors = colors.fillna('white')

                # Extract color values and labels from frame/series
                if isinstance(colors, pd.DataFrame):
                    labels = list(colors.columns)
                    colors = colors.T.values
                else:
                    labels = [colors.name]
                    colors = colors.values

            colors = _convert_colors(colors)

        return colors, labels

    def format_data(self, data, pivot_kws, z_score=None, standard_scale=None):
        """Extract variables from data or use directly."""

        # Either the data is already in 2d matrix format, or need to do a pivot
        if pivot_kws is not None:
            data2d = data.pivot(**pivot_kws)
        else:
            data2d = data

        if z_score is not None and standard_scale is not None:
            raise ValueError(
                'Cannot perform both z-scoring and standard-scaling on data')

        if z_score is not None:
            data2d = self.z_score(data2d, z_score)
        if standard_scale is not None:
            data2d = self.standard_scale(data2d, standard_scale)
        return data2d

    @staticmethod
    def z_score(data2d, axis=1):
        """Standarize the mean and variance of the data axis

        Parameters
        ----------
        data2d : pandas.DataFrame
            Data to normalize
        axis : int
            Which axis to normalize across. If 0, normalize across rows, if 1,
            normalize across columns.

        Returns
        -------
        normalized : pandas.DataFrame
            Noramlized data with a mean of 0 and variance of 1 across the
            specified axis.
        """
        if axis == 1:
            z_scored = data2d
        else:
            z_scored = data2d.T

        z_scored = (z_scored - z_scored.mean()) / z_scored.std()

        if axis == 1:
            return z_scored
        else:
            return z_scored.T

    @staticmethod
    def standard_scale(data2d, axis=1):
        """Divide the data by the difference between the max and min

        Parameters
        ----------
        data2d : pandas.DataFrame
            Data to normalize
        axis : int
            Which axis to normalize across. If 0, normalize across rows, if 1,
            normalize across columns.
        vmin : int
            If 0, then subtract the minimum of the data before dividing by
            the range.

        Returns
        -------
        standardized : pandas.DataFrame
            Noramlized data with a mean of 0 and variance of 1 across the
            specified axis.

        >>> import numpy as np
        >>> d = np.arange(5, 8, 0.5)
        >>> ClusterGrid.standard_scale(d)
        array([ 0. ,  0.2,  0.4,  0.6,  0.8,  1. ])
        """
        # Normalize these values to range from 0 to 1
        if axis == 1:
            standardized = data2d
        else:
            standardized = data2d.T

        subtract = standardized.min()
        standardized = (standardized - subtract) / (
            standardized.max() - standardized.min())

        if axis == 1:
            return standardized
        else:
            return standardized.T

    @staticmethod
    def color_list_to_matrix_and_cmap(colors, ind, axis=0):
        """Turns a list of colors into a numpy matrix and matplotlib colormap

        These arguments can now be plotted using heatmap(matrix, cmap)
        and the provided colors will be plotted.

        Parameters
        ----------
        colors : list of matplotlib colors
            Colors to label the rows or columns of a dataframe.
        ind : list of ints
            Ordering of the rows or columns, to reorder the original colors
            by the clustered dendrogram order
        axis : int
            Which axis this is labeling

        Returns
        -------
        matrix : numpy.array
            A numpy array of integer values, where each corresponds to a color
            from the originally provided list of colors
        cmap : matplotlib.colors.ListedColormap

        """
        # check for nested lists/color palettes.
        # Will fail if matplotlib color is list not tuple
        if any(issubclass(type(x), list) for x in colors):
            all_colors = set(itertools.chain(*colors))
            n = len(colors)
            m = len(colors[0])
        else:
            all_colors = set(colors)
            n = 1
            m = len(colors)
            colors = [colors]
        color_to_value = dict((col, i) for i, col in enumerate(all_colors))

        matrix = np.array([color_to_value[c] for color in colors
                           for c in color])

        shape = (n, m)
        matrix = matrix.reshape(shape)
        matrix = matrix[:, ind]
        if axis == 0:
            # row-side:
            matrix = matrix.T

        cmap = mpl.colors.ListedColormap(all_colors)
        return matrix, cmap

    def savefig(self, *args, **kwargs):
        if 'bbox_inches' not in kwargs:
            kwargs['bbox_inches'] = 'tight'
        self.fig.savefig(*args, **kwargs)

    def plot_dendrograms(self, row_cluster, col_cluster, metric, method,
                         row_linkage, col_linkage):
        # Plot the row dendrogram
        if row_cluster:
            self.dendrogram_row = dendrogram(
                self.data2d,
                metric=metric,
                method=method,
                label=False,
                axis=0,
                ax=self.ax_row_dendrogram,
                rotate=True,
                linkage=row_linkage)

            despine(ax=self.ax_row_dendrogram, bottom=True, left=True)

        # Plot the column dendrogram
        if col_cluster:
            self.dendrogram_col = dendrogram(
                self.data2d,
                metric=metric,
                method=method,
                label=False,
                axis=1,
                ax=self.ax_col_dendrogram,
                linkage=col_linkage)

            despine(ax=self.ax_col_dendrogram, bottom=True, left=True)

    def plot_colors(self, xind, yind, **kws):
        """Plots color labels between the dendrogram and the heatmap

        Parameters
        ----------
        heatmap_kws : dict
            Keyword arguments heatmap
        """
        # Remove any custom colormap and centering
        kws = kws.copy()
        kws.pop('cmap', None)
        kws.pop('center', None)
        kws.pop('vmin', None)
        kws.pop('vmax', None)
        kws.pop('xticklabels', None)
        kws.pop('yticklabels', None)

        if self.row_colors is not None:
            matrix, cmap = self.color_list_to_matrix_and_cmap(
                self.row_colors, yind, axis=0)

            # Get row_color labels
            if self.row_color_labels is not None:
                row_color_labels = self.row_color_labels
            else:
                row_color_labels = False

            heatmap(
                matrix,
                cmap=cmap,
                cbar=False,
                ax=self.ax_row_colors,
                xticklabels=row_color_labels,
                yticklabels=False,
                **kws)

            # Adjust rotation of labels
            if row_color_labels is not False:
                plt.setp(self.ax_row_colors.get_xticklabels(), rotation=90)

        if self.col_colors is not None:
            matrix, cmap = self.color_list_to_matrix_and_cmap(
                self.col_colors, xind, axis=1)

            # Get col_color labels
            if self.col_color_labels is not None:
                col_color_labels = self.col_color_labels
            else:
                col_color_labels = False

            heatmap(
                matrix,
                cmap=cmap,
                cbar=False,
                ax=self.ax_col_colors,
                xticklabels=False,
                yticklabels=col_color_labels,
                **kws)

            # Adjust rotation of labels, place on right side
            if col_color_labels is not False:
                self.ax_col_colors.yaxis.tick_right()
                plt.setp(self.ax_col_colors.get_yticklabels(), rotation=0)

    def plot_matrix(self, colorbar_kws, xind, yind, **kws):
        self.data2d = self.data2d.iloc[yind, xind]
        self.mask = self.mask.iloc[yind, xind]

        # Try to reorganize specified tick labels, if provided
        xtl = kws.pop("xticklabels", True)
        try:
            xtl = np.asarray(xtl)[xind]
        except (TypeError, IndexError):
            pass
        ytl = kws.pop("yticklabels", True)
        try:
            ytl = np.asarray(ytl)[yind]
        except (TypeError, IndexError):
            pass

        heatmap(
            self.data2d,
            ax=self.ax_heatmap,
            cbar=self.colorbar,
            cbar_ax=self.cax,
            cbar_kws=colorbar_kws,
            mask=self.mask,
            xticklabels=xtl,
            yticklabels=ytl,
            **kws)
        self.ax_heatmap.yaxis.set_ticks_position('right')
        self.ax_heatmap.yaxis.set_label_position('right')

    def plot(self, metric, method, colorbar_kws, row_linkage, col_linkage,
             **kws):
        colorbar_kws = {} if colorbar_kws is None else colorbar_kws
        self.plot_dendrograms(
            self.row_cluster,
            self.col_cluster,
            metric,
            method,
            row_linkage=row_linkage,
            col_linkage=col_linkage)
        try:
            xind = self.dendrogram_col.reordered_ind
        except AttributeError:
            xind = np.arange(self.data2d.shape[1])
        try:
            yind = self.dendrogram_row.reordered_ind
        except AttributeError:
            yind = np.arange(self.data2d.shape[0])

        self.plot_colors(xind, yind, **kws)
        self.plot_matrix(colorbar_kws, xind, yind, **kws)
        return self


def clustermap(data,
               pivot_kws=None,
               method='average',
               metric='euclidean',
               z_score=None,
               standard_scale=None,
               figsize=None,
               cbar_kws=None,
               row_cluster=True,
               col_cluster=True,
               row_linkage=None,
               col_linkage=None,
               row_colors=None,
               col_colors=None,
               row_ratios=None,
               col_ratios=None,
               colorbar=True,
               mask=None,
               **kwargs):
    """Plot a hierarchically clustered heatmap of a pandas DataFrame

    Parameters
    ----------
    data: pandas.DataFrame
        Rectangular data for clustering. Cannot contain NAs.
    pivot_kws : dict, optional
        If `data` is a tidy dataframe, can provide keyword arguments for
        pivot to create a rectangular dataframe.
    method : str, optional
        Linkage method to use for calculating clusters.
        See scipy.cluster.hierarchy.linkage documentation for more information:
        http://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
    metric : str, optional
        Distance metric to use for the data. See
        scipy.spatial.distance.pdist documentation for more options
        http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html
    z_score : int or None, optional
        Either 0 (rows) or 1 (columns). Whether or not to calculate z-scores
        for the rows or the columns. Z scores are: z = (x - mean)/std, so
        values in each row (column) will get the mean of the row (column)
        subtracted, then divided by the standard deviation of the row (column).
        This ensures that each row (column) has mean of 0 and variance of 1.
    standard_scale : int or None, optional
        Either 0 (rows) or 1 (columns). Whether or not to standardize that
        dimension, meaning for each row or column, subtract the minimum and
        divide each by its maximum.
    figsize: tuple of two ints, optional
        Size of the figure to create.
    cbar_kws : dict, optional
        Keyword arguments to pass to ``cbar_kws`` in ``heatmap``, e.g. to
        add a label to the colorbar.
    {row,col}_cluster : bool, optional
        If True, cluster the {rows, columns}.
    {row,col}_linkage : numpy.array, optional
        Precomputed linkage matrix for the rows or columns. See
        scipy.cluster.hierarchy.linkage for specific formats.
    {row,col}_colors : list-like or pandas DataFrame/Series, optional
        List of colors to label for either the rows or columns. Useful to
        evaluate whether samples within a group are clustered together. Can
        use nested lists or DataFrame for multiple color levels of labeling.
        If given as a DataFrame or Series, labels for the colors are extracted
        from the DataFrames column names or from the name of the Series.
        DataFrame/Series colors are also matched to the data by their
        index, ensuring colors are drawn in the correct order.
    mask : boolean array or DataFrame, optional
        If passed, data will not be shown in cells where ``mask`` is True.
        Cells with missing values are automatically masked. Only used for
        visualizing, not for calculating.
    kwargs : other keyword arguments
        All other keyword arguments are passed to ``sns.heatmap``

    Returns
    -------
    clustergrid : ClusterGrid
        A ClusterGrid instance.

    Notes
    -----
    The returned object has a ``savefig`` method that should be used if you
    want to save the figure object without clipping the dendrograms.

    To access the reordered row indices, use:
    ``clustergrid.dendrogram_row.reordered_ind``

    Column indices, use:
    ``clustergrid.dendrogram_col.reordered_ind``

    Examples
    --------

    Plot a clustered heatmap:

    .. plot::
        :context: close-figs

        >>> import seaborn as sns; sns.set()
        >>> flights = sns.load_dataset("flights")
        >>> flights = flights.pivot("month", "year", "passengers")
        >>> g = sns.clustermap(flights)

    Don't cluster one of the axes:

    .. plot::
        :context: close-figs

        >>> g = sns.clustermap(flights, col_cluster=False)

    Use a different colormap and add lines to separate the cells:

    .. plot::
        :context: close-figs

        >>> cmap = sns.cubehelix_palette(as_cmap=True, rot=-.3, light=1)
        >>> g = sns.clustermap(flights, cmap=cmap, linewidths=.5)

    Use a different figure size:

    .. plot::
        :context: close-figs

        >>> g = sns.clustermap(flights, cmap=cmap, figsize=(7, 5))

    Standardize the data across the columns:

    .. plot::
        :context: close-figs

        >>> g = sns.clustermap(flights, standard_scale=1)

    Normalize the data across the rows:

    .. plot::
        :context: close-figs

        >>> g = sns.clustermap(flights, z_score=0)

    Use a different clustering method:

    .. plot::
        :context: close-figs

        >>> g = sns.clustermap(flights, method="single", metric="cosine")

    Add colored labels on one of the axes:

    .. plot::
        :context: close-figs

        >>> season_colors = (sns.color_palette("BuPu", 3) +
        ...                  sns.color_palette("RdPu", 3) +
        ...                  sns.color_palette("YlGn", 3) +
        ...                  sns.color_palette("OrRd", 3))
        >>> g = sns.clustermap(flights, row_colors=season_colors)

    """
    plotter = ClusterGrid(
        data,
        pivot_kws=pivot_kws,
        figsize=figsize,
        row_colors=row_colors,
        col_colors=col_colors,
        z_score=z_score,
        standard_scale=standard_scale,
        row_ratios=row_ratios,
        col_ratios=col_ratios,
        row_cluster=row_cluster,
        col_cluster=col_cluster,
        colorbar=colorbar,
        mask=mask)

    return plotter.plot(
        metric=metric,
        method=method,
        colorbar_kws=cbar_kws,
        row_linkage=row_linkage,
        col_linkage=col_linkage,
        **kwargs)
