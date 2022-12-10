import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def make_sensitivity_heatmap(data, centre_point, xlabel, ylabel, ax=None, 
                             title=None, minmax=None, transform=None, cmap='coolwarm', 
                             reverse_y=True, fmt='{:.2f}'):
    """Make a custom heatmap plot to show the sensitivity results
    data.
    """

    # See documentation example:
    # https://matplotlib.org/stable/gallery/images_contours_and_fields/image_annotated_heatmap.html

    # Transform data
    if callable(transform):
        plot_data = transform(data)
    elif transform is None:
        plot_data = data
    elif transform == 'diff':
        centre_value = data.loc[centre_point]
        plot_data = data - centre_value
    elif transform == 'fraction':
        centre_value = data.loc[centre_point]
        plot_data = (data - centre_value) / centre_value
    elif transform == 'pct':
        centre_value = data.loc[centre_point]
        plot_data = 100 * (data - centre_value) / centre_value
    else:
        raise ValueError("invalid plot type")
    
    if minmax is None:
        minmax = (np.min(plot_data.values), np.max(plot_data.values))

    if reverse_y:
        plot_data = plot_data.iloc[::-1, :]

    pct_diff_idx = np.round(100 * (plot_data.index - centre_point[0]) 
                    / centre_point[0]).astype('int')
    pct_diff_cols = np.round(100 * (plot_data.columns - centre_point[1]) 
                     / centre_point[1]).astype('int')

    if ax is None:
        ax = plt.gca()

    im = ax.imshow(plot_data, aspect='auto', vmin=minmax[0], vmax=minmax[1], 
                   cmap=cmap)

    # Show all ticks and label them with the respective list entries
    ax.set_xticks(np.arange(len(plot_data.columns)), labels=pct_diff_cols)
    ax.set_yticks(np.arange(len(plot_data.index)), labels=pct_diff_idx)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title is not None:
        ax.set_title(title)

    # Loop over data dimensions and create text annotations.
    for i, idx in enumerate(plot_data.index):
        for j, col in enumerate(plot_data.columns):
            text = ax.text(j, i, fmt.format(plot_data.loc[idx, col]), 
                           ha="center", va="center")
    #fig.colorbar()