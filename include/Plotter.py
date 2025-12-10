"""
mplhep-based plotting module for Galapago-FWK.

Provides CMS-style plots using matplotlib and mplhep.
"""

import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import os


# Set CMS style by default
plt.style.use(hep.style.CMS)


def root_to_numpy(hist):
    """
    Convert ROOT TH1 histogram to numpy arrays.

    Args:
        hist: ROOT TH1F/TH1D histogram

    Returns:
        tuple: (bin_contents, bin_edges, bin_errors)
    """
    nbins = hist.GetNbinsX()
    contents = np.array([hist.GetBinContent(i) for i in range(1, nbins + 1)])
    errors = np.array([hist.GetBinError(i) for i in range(1, nbins + 1)])
    edges = np.array([hist.GetBinLowEdge(i) for i in range(1, nbins + 2)])
    return contents, edges, errors


def plot_histogram(hist, ax=None, label=None, color=None, histtype='step',
                   yerr=False, **kwargs):
    """
    Plot a ROOT histogram using mplhep.

    Args:
        hist: ROOT TH1 histogram
        ax: matplotlib axes (creates new if None)
        label: Legend label
        color: Line/fill color
        histtype: 'step', 'fill', or 'errorbar'
        yerr: Show error bars
        **kwargs: Additional arguments passed to hep.histplot

    Returns:
        matplotlib axes
    """
    if ax is None:
        fig, ax = plt.subplots()

    contents, edges, errors = root_to_numpy(hist)

    if histtype == 'errorbar':
        # Data points with error bars
        centers = (edges[:-1] + edges[1:]) / 2
        ax.errorbar(centers, contents, yerr=errors, fmt='o',
                    color=color, label=label, markersize=4, **kwargs)
    else:
        hep.histplot(contents, edges, ax=ax, label=label, color=color,
                     histtype=histtype, yerr=errors if yerr else None, **kwargs)

    return ax


def plot_data_mc(histoDATA, histoMC, output_path, xlabel='', ylabel='Events',
                 data_label='Data', mc_label='MC', mc_color='#1f77b4',
                 lumi=None, year=None, xlim=None, ylim=None,
                 ratio=True, logy=False, formats=['png', 'pdf']):
    """
    Create a Data/MC comparison plot with optional ratio panel.

    Args:
        histoDATA: ROOT histogram for data
        histoMC: ROOT histogram for MC
        output_path: Output file path (without extension)
        xlabel: X-axis label
        ylabel: Y-axis label
        data_label: Legend label for data
        mc_label: Legend label for MC
        mc_color: Color for MC histogram
        lumi: Luminosity value for label (fb^-1)
        year: Year for CMS label
        xlim: X-axis limits (tuple)
        ylim: Y-axis limits (tuple)
        ratio: Show ratio panel
        logy: Logarithmic y-axis
        formats: Output formats

    Returns:
        matplotlib figure
    """
    data_contents, edges, data_errors = root_to_numpy(histoDATA)
    mc_contents, _, mc_errors = root_to_numpy(histoMC)
    centers = (edges[:-1] + edges[1:]) / 2

    if ratio:
        fig, (ax_main, ax_ratio) = plt.subplots(2, 1, figsize=(10, 10),
                                                 gridspec_kw={'height_ratios': [3, 1]},
                                                 sharex=True)
        fig.subplots_adjust(hspace=0.05)
    else:
        fig, ax_main = plt.subplots(figsize=(10, 8))
        ax_ratio = None

    # Main panel - MC as filled histogram
    hep.histplot(mc_contents, edges, ax=ax_main, label=mc_label,
                 color=mc_color, histtype='fill', alpha=0.7)
    hep.histplot(mc_contents, edges, ax=ax_main,
                 color=mc_color, histtype='step', linewidth=1.5)

    # Main panel - Data as points
    ax_main.errorbar(centers, data_contents, yerr=data_errors, fmt='o',
                     color='black', label=data_label, markersize=5, capsize=2)

    # Styling
    ax_main.set_ylabel(ylabel)
    ax_main.legend(loc='upper right', frameon=False)

    if logy:
        ax_main.set_yscale('log')
        ax_main.set_ylim(bottom=0.1)
    else:
        ax_main.set_ylim(bottom=0)

    if ylim:
        ax_main.set_ylim(ylim)

    # CMS label
    if lumi and year:
        hep.cms.label(ax=ax_main, data=True, lumi=lumi, year=year, loc=0)
    elif lumi:
        hep.cms.label(ax=ax_main, data=True, lumi=lumi, loc=0)
    else:
        hep.cms.label(ax=ax_main, data=True, loc=0)

    # Ratio panel
    if ratio and ax_ratio is not None:
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio_vals = np.where(mc_contents > 0, data_contents / mc_contents, 0)
            ratio_errs = np.where(mc_contents > 0, data_errors / mc_contents, 0)

        ax_ratio.errorbar(centers, ratio_vals, yerr=ratio_errs, fmt='o',
                         color='black', markersize=5, capsize=2)
        ax_ratio.axhline(1, color='gray', linestyle='--', linewidth=1)
        ax_ratio.set_ylabel('Data / MC')
        ax_ratio.set_xlabel(xlabel)
        ax_ratio.set_ylim(0.5, 1.5)

        if xlim:
            ax_ratio.set_xlim(xlim)
    else:
        ax_main.set_xlabel(xlabel)
        if xlim:
            ax_main.set_xlim(xlim)

    # Save
    os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else '.', exist_ok=True)
    for fmt in formats:
        fig.savefig(f"{output_path}.{fmt}", dpi=150, bbox_inches='tight')
        print(f"Saved: {output_path}.{fmt}")

    return fig


def plot_stacked(histograms, labels, colors, output_path, histoDATA=None,
                 xlabel='', ylabel='Events', data_label='Data',
                 lumi=None, year=None, ratio=True, logy=False,
                 formats=['png', 'pdf']):
    """
    Create a stacked histogram plot with optional data overlay.

    Args:
        histograms: List of ROOT histograms to stack
        labels: List of labels for each histogram
        colors: List of colors for each histogram
        output_path: Output file path (without extension)
        histoDATA: Optional data histogram to overlay
        xlabel: X-axis label
        ylabel: Y-axis label
        data_label: Legend label for data
        lumi: Luminosity for CMS label
        year: Year for CMS label
        ratio: Show ratio panel
        logy: Logarithmic y-axis
        formats: Output formats

    Returns:
        matplotlib figure
    """
    # Convert all histograms
    all_contents = []
    edges = None
    for hist in histograms:
        contents, edges, _ = root_to_numpy(hist)
        all_contents.append(contents)

    if ratio and histoDATA:
        fig, (ax_main, ax_ratio) = plt.subplots(2, 1, figsize=(10, 10),
                                                 gridspec_kw={'height_ratios': [3, 1]},
                                                 sharex=True)
        fig.subplots_adjust(hspace=0.05)
    else:
        fig, ax_main = plt.subplots(figsize=(10, 8))
        ax_ratio = None

    # Stacked histograms
    hep.histplot(all_contents, edges, ax=ax_main, label=labels,
                 color=colors, histtype='fill', stack=True, alpha=0.8)

    # Data overlay
    if histoDATA:
        data_contents, _, data_errors = root_to_numpy(histoDATA)
        centers = (edges[:-1] + edges[1:]) / 2
        ax_main.errorbar(centers, data_contents, yerr=data_errors, fmt='o',
                        color='black', label=data_label, markersize=5, capsize=2)

    # Styling
    ax_main.set_ylabel(ylabel)
    ax_main.legend(loc='upper right', frameon=False)

    if logy:
        ax_main.set_yscale('log')
        ax_main.set_ylim(bottom=0.1)
    else:
        ax_main.set_ylim(bottom=0)

    # CMS label
    if lumi:
        hep.cms.label(ax=ax_main, data=histoDATA is not None, lumi=lumi, year=year, loc=0)
    else:
        hep.cms.label(ax=ax_main, data=histoDATA is not None, loc=0)

    # Ratio panel
    if ratio and ax_ratio is not None and histoDATA:
        data_contents, _, data_errors = root_to_numpy(histoDATA)
        mc_total = np.sum(all_contents, axis=0)
        centers = (edges[:-1] + edges[1:]) / 2

        with np.errstate(divide='ignore', invalid='ignore'):
            ratio_vals = np.where(mc_total > 0, data_contents / mc_total, 0)
            ratio_errs = np.where(mc_total > 0, data_errors / mc_total, 0)

        ax_ratio.errorbar(centers, ratio_vals, yerr=ratio_errs, fmt='o',
                         color='black', markersize=5, capsize=2)
        ax_ratio.axhline(1, color='gray', linestyle='--', linewidth=1)
        ax_ratio.set_ylabel('Data / MC')
        ax_ratio.set_xlabel(xlabel)
        ax_ratio.set_ylim(0.5, 1.5)
    elif ax_ratio is None:
        ax_main.set_xlabel(xlabel)

    # Save
    os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else '.', exist_ok=True)
    for fmt in formats:
        fig.savefig(f"{output_path}.{fmt}", dpi=150, bbox_inches='tight')
        print(f"Saved: {output_path}.{fmt}")

    return fig


def plot_simple(hist, output_path, xlabel='', ylabel='Events', label=None,
                color='#1f77b4', histtype='step', lumi=None, year=None,
                logy=False, formats=['png', 'pdf']):
    """
    Create a simple single histogram plot.

    Args:
        hist: ROOT histogram
        output_path: Output file path (without extension)
        xlabel: X-axis label
        ylabel: Y-axis label
        label: Legend label
        color: Histogram color
        histtype: 'step' or 'fill'
        lumi: Luminosity for CMS label
        year: Year for CMS label
        logy: Logarithmic y-axis
        formats: Output formats

    Returns:
        matplotlib figure
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    contents, edges, errors = root_to_numpy(hist)
    hep.histplot(contents, edges, ax=ax, label=label, color=color,
                 histtype=histtype, yerr=errors)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if label:
        ax.legend(loc='upper right', frameon=False)

    if logy:
        ax.set_yscale('log')
        ax.set_ylim(bottom=0.1)
    else:
        ax.set_ylim(bottom=0)

    if lumi:
        hep.cms.label(ax=ax, data=False, lumi=lumi, year=year, loc=0)
    else:
        hep.cms.label(ax=ax, data=False, loc=0)

    # Save
    os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else '.', exist_ok=True)
    for fmt in formats:
        fig.savefig(f"{output_path}.{fmt}", dpi=150, bbox_inches='tight')
        print(f"Saved: {output_path}.{fmt}")

    return fig
