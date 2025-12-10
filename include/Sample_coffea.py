"""
Simplified Coffea-based Sample module for CMS NanoAOD analysis.
Uses coffea, awkward arrays, and hist for analysis.
"""

import awkward as ak
import numpy as np
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
import hist
from hist import Hist
import subprocess
import hashlib
import os
from typing import List, Optional, Union

# XRootD redirector for CMS data access
XROOTD_REDIRECTOR = "root://cms-xrd-global.cern.ch/"


def queryDAS(dataset: str, limit: int = 0) -> List[str]:
    """
    Query CMS DAS for files in a dataset.

    Args:
        dataset: DAS dataset name (e.g., '/Muon0/Run2023C-PromptNanoAODv12-v1/NANOAOD')
        limit: Maximum number of files to return (0 = all)

    Returns:
        List of xrootd file paths
    """
    query = f"file dataset={dataset}"
    cmd = ["dasgoclient", "--query", query]
    if limit > 0:
        cmd.extend(["--limit", str(limit)])

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        files = [f for f in result.stdout.strip().split('\n') if f]
        xrootd_files = [XROOTD_REDIRECTOR + f for f in files]
        print(f"DAS query returned {len(xrootd_files)} files for {dataset}")
        return xrootd_files
    except subprocess.CalledProcessError as e:
        print(f"DAS query failed: {e.stderr}")
        return []
    except FileNotFoundError:
        print("ERROR: dasgoclient not found. Please setup CMS environment (cmsenv)")
        return []


def getCachedDASFiles(dataset: str, cache_dir: str = ".das_cache") -> List[str]:
    """Get DAS files with local caching to avoid repeated queries."""
    os.makedirs(cache_dir, exist_ok=True)
    cache_key = hashlib.md5(dataset.encode()).hexdigest()
    cache_file = os.path.join(cache_dir, f"{cache_key}.txt")

    if os.path.exists(cache_file):
        with open(cache_file, 'r') as f:
            files = [line.strip() for line in f if line.strip()]
            print(f"Using cached DAS query for {dataset} ({len(files)} files)")
            return files

    files = queryDAS(dataset)
    if files:
        with open(cache_file, 'w') as f:
            f.write('\n'.join(files))
    return files


########################################################################################
################################### Sample Class #######################################
########################################################################################
class Sample:
    """
    Simplified Sample class for NanoAOD analysis with coffea.

    Usage:
        # Create sample from DAS dataset
        data = Sample(
            name='Muon0_Run2024I',
            dataset='/Muon0/Run2024I-PromptNanoAODv12-v1/NANOAOD',
            isdata=True,
            xsection=1.0,
            file_limit=2
        )

        # Apply selections
        data.addSelection('nMuon > 1')
        data.addSelection('Muon_pt[:, 0] > 25')

        # Create histogram
        h = data.getHist(
            var='DiMuon_mass',
            bins=100,
            range=(40, 150),
            lumi=35.9
        )
    """

    def __init__(self, name: str, dataset: str, isdata: bool = False,
                 xsection: float = 1.0, file_limit: int = 0, label: str = None,
                 color: str = '#000000'):
        """
        Initialize Sample.

        Args:
            name: Sample name
            dataset: DAS dataset path
            isdata: True if data, False if MC
            xsection: Cross section in pb (only for MC)
            file_limit: Limit number of files to load (0 = all)
            label: Label for plotting (defaults to name)
            color: Color for plotting (hex string)
        """
        self.name = name
        self.dataset = dataset
        self.isData = isdata
        self.xSection = xsection
        self.label = label or name
        self.color = color
        self.file_limit = file_limit

        # Storage
        self.events = None
        self.file_paths = []
        self.selections = []
        self.sum_genweight = 0.0
        self.nevents = 0

        # Load events
        self._loadEvents()
        self._setupWeights()

        print(f"Loaded {self.name}: {self.nevents} events from {len(self.file_paths)} files")

    def _loadEvents(self):
        """Load events from DAS dataset using coffea."""
        # Get file list from DAS
        self.file_paths = getCachedDASFiles(self.dataset)

        if self.file_limit > 0:
            self.file_paths = self.file_paths[:self.file_limit]

        if not self.file_paths:
            raise RuntimeError(f"No files found for dataset {self.dataset}")

        # Load events
        all_events = []
        for fpath in self.file_paths:
            try:
                events = NanoEventsFactory.from_root(
                    fpath,
                    schemaclass=NanoAODSchema,
                ).events()
                all_events.append(events)
            except Exception as e:
                print(f"Warning: Failed to load {fpath}: {e}")
                continue

        # Concatenate
        if all_events:
            self.events = ak.concatenate(all_events, axis=0)
            self.nevents = len(self.events)
        else:
            raise RuntimeError(f"Could not load any events for {self.name}")

    def _setupWeights(self):
        """Setup event weights."""
        if not self.isData:
            # MC: sum of genWeights for normalization
            self.sum_genweight = float(ak.sum(self.events.genWeight))
            lum_weight = self.xSection / self.sum_genweight

            # Add weight fields
            self.events = ak.with_field(
                self.events,
                self.events.genWeight * lum_weight,
                'eventWeight'
            )
        else:
            # Data: unit weight
            self.sum_genweight = float(self.nevents)
            self.events = ak.with_field(
                self.events,
                ak.ones_like(self.events.event),
                'eventWeight'
            )

    def addSelection(self, cut: str):
        """
        Add a selection cut.

        Args:
            cut: Selection string in awkward/numpy syntax
                 Examples: 'events.nMuon > 1', 'events.Muon.pt[:, 0] > 25'
        """
        self.selections.append(cut)
        return self

    def clearSelections(self):
        """Clear all selections."""
        self.selections = []
        return self

    def _applySelections(self):
        """Apply all selection cuts and return filtered events."""
        events = self.events

        for cut in self.selections:
            # Evaluate cut
            mask = eval(cut, {'events': events, 'ak': ak, 'np': np})
            events = events[mask]

        return events

    def getHist(self, var: str, bins: int, range: tuple,
                lumi: float = 1.0, weight: str = None) -> Hist:
        """
        Create 1D histogram.

        Args:
            var: Variable to plot (use 'events.Muon.pt[:, 0]' syntax)
            bins: Number of bins
            range: (min, max) tuple
            lumi: Integrated luminosity in fb^-1
            weight: Optional weight expression (default: eventWeight * lumi)

        Returns:
            hist.Hist object
        """
        # Apply selections
        events = self._applySelections()

        # Get variable values
        values = eval(var, {'events': events, 'ak': ak, 'np': np})

        # Get weights
        if weight is None:
            if self.isData:
                weights = ak.ones_like(values)
            else:
                weights = events.eventWeight * lumi
        else:
            weights = eval(weight, {'events': events, 'ak': ak, 'np': np})

        # Create histogram
        h = Hist(
            hist.axis.Regular(bins, range[0], range[1], name=var),
            storage=hist.storage.Weight()
        )

        # Fill
        h.fill(ak.flatten(values), weight=ak.flatten(weights))

        return h

    def getHist2D(self, varx: str, vary: str,
                  binsx: int, rangex: tuple,
                  binsy: int, rangey: tuple,
                  lumi: float = 1.0) -> Hist:
        """
        Create 2D histogram.

        Args:
            varx: X variable
            vary: Y variable
            binsx: Number of X bins
            rangex: (xmin, xmax)
            binsy: Number of Y bins
            rangey: (ymin, ymax)
            lumi: Integrated luminosity in fb^-1

        Returns:
            hist.Hist object (2D)
        """
        # Apply selections
        events = self._applySelections()

        # Get values
        xvalues = eval(varx, {'events': events, 'ak': ak, 'np': np})
        yvalues = eval(vary, {'events': events, 'ak': ak, 'np': np})

        # Get weights
        if self.isData:
            weights = ak.ones_like(xvalues)
        else:
            weights = events.eventWeight * lumi

        # Create histogram
        h = Hist(
            hist.axis.Regular(binsx, rangex[0], rangex[1], name=varx),
            hist.axis.Regular(binsy, rangey[0], rangey[1], name=vary),
            storage=hist.storage.Weight()
        )

        # Fill
        h.fill(
            ak.flatten(xvalues),
            ak.flatten(yvalues),
            weight=ak.flatten(weights)
        )

        return h

    def getYield(self, lumi: float = 1.0) -> tuple:
        """
        Get event yield after selections.

        Args:
            lumi: Integrated luminosity in fb^-1

        Returns:
            (yield, uncertainty) tuple
        """
        events = self._applySelections()

        if self.isData:
            weights = ak.ones_like(events.event)
        else:
            weights = events.eventWeight * lumi

        total = float(ak.sum(weights))
        # Simple uncertainty: sqrt(sum of weights squared)
        unc = float(np.sqrt(ak.sum(weights**2)))

        return (total, unc)

    def printInfo(self):
        """Print sample information."""
        print("=" * 50)
        print(f"Sample: {self.name}")
        print(f"Dataset: {self.dataset}")
        print(f"Type: {'Data' if self.isData else 'MC'}")
        if not self.isData:
            print(f"Cross section: {self.xSection} pb")
            print(f"Sum genWeight: {self.sum_genweight:.2e}")
        print(f"Events loaded: {self.nevents}")
        print(f"Files: {len(self.file_paths)}")
        print(f"Selections: {len(self.selections)}")
        for i, sel in enumerate(self.selections, 1):
            print(f"  {i}. {sel}")
        print("=" * 50)


########################################################################################
################################ Helper Functions ######################################
########################################################################################

def plotComparison(h_data: Hist, h_mc: Hist,
                   xlabel: str, ylabel: str = 'Events',
                   data_label: str = 'Data', mc_label: str = 'MC',
                   output: str = 'plot.png', ratio: bool = True):
    """
    Plot data/MC comparison using mplhep.

    Args:
        h_data: Data histogram
        h_mc: MC histogram
        xlabel: X-axis label
        ylabel: Y-axis label
        data_label: Label for data
        mc_label: Label for MC
        output: Output filename
        ratio: Show ratio panel
    """
    import matplotlib.pyplot as plt
    import mplhep as hep

    plt.style.use(hep.style.CMS)

    if ratio:
        fig, (ax, rax) = plt.subplots(
            2, 1,
            figsize=(10, 10),
            gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05}
        )
    else:
        fig, ax = plt.subplots(figsize=(10, 8))

    # Main plot
    hep.histplot(h_mc, ax=ax, label=mc_label, histtype='fill', alpha=0.7)
    hep.histplot(h_data, ax=ax, label=data_label, histtype='errorbar',
                 color='black', markersize=8)

    ax.set_ylabel(ylabel, fontsize=14)
    ax.legend(fontsize=12)
    ax.set_yscale('log')

    if ratio:
        ax.set_xlabel('')
        ax.tick_params(labelbottom=False)

        # Ratio plot
        ratio_vals = h_data.values() / h_mc.values()
        ratio_err = np.sqrt(h_data.variances()) / h_mc.values()

        bin_centers = (h_data.axes[0].edges[:-1] + h_data.axes[0].edges[1:]) / 2

        rax.errorbar(bin_centers, ratio_vals, yerr=ratio_err,
                     fmt='o', color='black', markersize=6)
        rax.axhline(1, color='gray', linestyle='--')
        rax.set_xlabel(xlabel, fontsize=14)
        rax.set_ylabel('Data / MC', fontsize=12)
        rax.set_ylim(0.5, 1.5)
        rax.grid(True, alpha=0.3)
    else:
        ax.set_xlabel(xlabel, fontsize=14)

    # CMS label
    hep.cms.label(ax=ax, data=True, lumi='35.9', year='2024', loc=0)

    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {output}")
    plt.close()


def combineHistograms(histograms: List[Hist]) -> Hist:
    """
    Combine multiple histograms by summing them.

    Args:
        histograms: List of hist.Hist objects

    Returns:
        Combined histogram
    """
    if not histograms:
        raise ValueError("Empty histogram list")

    result = histograms[0].copy()
    for h in histograms[1:]:
        result = result + h

    return result
