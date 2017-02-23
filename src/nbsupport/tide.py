import matplotlib.pyplot as plt
import matplotlib.patches as mpl_patches
import seaborn as sns


def plot(data,
         ax=None,
         overall_efficiency=None,
         legend=True,
         legend_kws=None,
         label_kws=None):

    # Check data columns.
    # assert all([c in data.columns for c in ['shift', 'percentage', 'pvalue']])

    data = data.sort_values(by='shift')

    if ax is None:
        _, ax = plt.subplots()

    color_map = {'p < 0.001': (0.98, 0.1, 0.11), 'p > 0.001': 'black'}

    palette = {row.shift: color_map['p < 0.001']
               if row.pvalue < 0.001 else color_map['p > 0.001']
               for row in data.itertuples()}
    palette[0] = (0.89, 0.59, 0.59)

    sns.barplot(data=data, x='shift', y='percentage', palette=palette, ax=ax)
    sns.despine()

    ax.set_xlabel('\u2190 Deletions / Insertions \u2192')
    ax.set_ylabel('Gene editing (percentage)')

    ax.set_ylim(0, 100)

    if overall_efficiency is not None:
        ax.annotate(
            xy=(0.05, 0.94),
            xycoords='axes fraction',
            s='Gene editing = {}%'.format(overall_efficiency),
            **(label_kws or {}))

    if legend:
        _add_legend(color_map, ax=ax, **(legend_kws or {}))

    _add_percentages(data, ax, **(label_kws or {}))

    return ax


def _add_legend(color_dict, ax, **kwargs):
    patches = [mpl_patches.Patch(
        color=color, label=label) for label, color in color_dict.items()]
    return ax.legend(handles=patches, **kwargs)


def _add_percentages(data, ax, max_pval=0.001, **kwargs):
    for row, patch in zip(data.itertuples(), ax.patches):
        if row.pvalue < max_pval:
            ax.text(
                x=patch.get_x() + (patch.get_width() / 2),
                y=patch.get_height() + 1.5,
                ha='center',
                s='{}%'.format(row.percentage),
                **kwargs)
