#!/usr/bin/env python3
"""
Harvesting script for CMS Run 3 Muon analysis.
Creates dimuon invariant mass histograms comparing data and MC.

Usage:
    python harvesting_Muon0_Summer23.py [options]

Examples:
    # Run with default DAS samples
    python harvesting_Muon0_Summer23.py

    # Run with custom dat file (e.g., local EOS paths)
    python harvesting_Muon0_Summer23.py --dat dat/Samples_cern_Muon_provisionalRun3.dat

    # Test mode with limited files
    python harvesting_Muon0_Summer23.py --test --nfiles 2
"""

import ROOT as r
from ROOT import gROOT
import os
import argparse

import include.Sample as Sample
import include.helper as helper
import include.Plotter as Plotter

# Determine working directory from script location
WORKPATH = os.path.dirname(os.path.abspath(__file__)) + '/'


def setup_root():
    """Setup ROOT for batch processing."""
    gROOT.SetBatch(1)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Create dimuon mass histograms from CMS Run 3 data',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--dat', type=str,
                        default='dat/Samples_DAS_Muon_Run3.dat',
                        help='Sample definition file')
    parser.add_argument('--lumi', type=float, default=35.9,
                        help='Integrated luminosity (fb^-1)')
    parser.add_argument('--output', '-o', type=str, default='plots',
                        help='Output directory for plots')
    parser.add_argument('--test', action='store_true',
                        help='Test mode: use limited files')
    parser.add_argument('--nfiles', type=int, default=1,
                        help='Number of files to use in test mode')
    parser.add_argument('--year', type=str, default='2024',
                        help='Year for CMS label')
    parser.add_argument('--ratio', action='store_true', default=True,
                        help='Show ratio panel in plots')
    parser.add_argument('--no-ratio', action='store_false', dest='ratio',
                        help='Hide ratio panel in plots')
    return parser.parse_args()


def load_samples(datfile, data_names, mc_names, file_limit=0):
    """
    Load data and MC samples from dat file.

    Args:
        datfile: Path to sample definition file
        data_names: List of data sample names to load
        mc_names: List of MC sample names to load
        file_limit: Limit files per sample (0 = all)

    Returns:
        tuple: (data_tree, mc_tree)
    """
    data_file = helper.selectSamples(datfile, data_names, 'DATA')
    mc_file = helper.selectSamples(datfile, mc_names, 'MC')

    treeDATA = Sample.Tree(fileName=data_file, name='DATA', isdata=1, file_limit=file_limit)
    treeMC = Sample.Tree(fileName=mc_file, name='MC', isdata=0, file_limit=file_limit)

    return treeDATA, treeMC


def create_dimuon_mass_histogram(tree, lumi, name_suffix=''):
    """
    Create dimuon invariant mass histogram.

    Args:
        tree: Sample.Tree object
        lumi: Integrated luminosity
        name_suffix: Suffix for histogram name

    Returns:
        ROOT TH1F histogram
    """
    cut = '(nMuon > 1 && Muon_charge[0]*Muon_charge[1] < 0.0)'

    return tree.getTH1F(
        lumi=lumi,
        name=f'DiMuon_mass{name_suffix}',
        var='DiMuon_mass',
        nbin=100,
        xmin=40.0,
        xmax=150.0,
        cut=cut,
        xlabel='Dimuon mass (GeV)'
    )


def get_sample_color(tree):
    """Extract color from tree's first block's first sample."""
    if tree.blocks and tree.blocks[0].samples:
        return tree.blocks[0].samples[0].color
    return '#000000'


def main():
    """Main execution."""
    args = parse_args()

    # Setup
    setup_root()
    helper.ensureDirectory(args.output)

    # Sample configuration
    datfile = WORKPATH + args.dat
    data_samples = ['Muon0_Run2024I_v1']
    mc_samples = ['DYJetsToLL_M-50']

    # Variable definitions for dimuon mass calculation
    definitions_file = 'definitions/defs_MuonVariables.json'

    print(f"{'='*60}")
    print(f"Galapago-FWK: Dimuon Mass Analysis (mplhep)")
    print(f"{'='*60}")
    print(f"Dat file:    {args.dat}")
    print(f"Luminosity:  {args.lumi} fb^-1")
    print(f"Output:      {args.output}")
    print(f"Test mode:   {args.test}")
    print(f"Year:        {args.year}")
    print(f"{'='*60}")

    # Load samples
    print("\nLoading samples...")
    treeDATA, treeMC = load_samples(
        datfile,
        data_samples, mc_samples,
        file_limit=args.nfiles if args.test else 0
    )

    # Get MC color from dat file
    mc_color = get_sample_color(treeMC)
    print(f"MC color: {mc_color}")

    # Apply variable definitions
    treeDATA.setDefinitions(config=definitions_file)
    treeMC.setDefinitions(config=definitions_file)

    # Create histograms
    print("\nCreating histograms...")
    histoDATA = create_dimuon_mass_histogram(treeDATA, args.lumi, '_data')
    histoMC = create_dimuon_mass_histogram(treeMC, args.lumi, '_mc')

    # Print yields
    print(f"\nYields in mass window [40, 150] GeV:")
    print(f"  Data: {histoDATA.Integral():.0f} events")
    print(f"  MC:   {histoMC.Integral():.1f} events (scaled to {args.lumi} fb^-1)")

    # Create plot using mplhep
    print("\nCreating comparison plot with mplhep...")
    output_path = os.path.join(args.output, 'DiMuon_mass')

    Plotter.plot_data_mc(
        histoDATA=histoDATA,
        histoMC=histoMC,
        output_path=output_path,
        xlabel='Dimuon mass (GeV)',
        ylabel='Events / GeV',
        data_label='Data',
        mc_label='Drell-Yan',
        mc_color=mc_color,
        lumi=args.lumi,
        year=args.year,
        ratio=args.ratio,
        formats=['png', 'pdf']
    )

    print("\nDone!")


if __name__ == "__main__":
    main()
