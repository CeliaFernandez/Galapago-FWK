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
import numpy as np

import include.Sample as Sample
import include.helper as helper
import include.Canvas as Canvas

# Determine working directory from script location
WORKPATH = os.path.dirname(os.path.abspath(__file__)) + '/'


def setup_style():
    """Load CMS TDR plot style."""
    gROOT.LoadMacro(WORKPATH + 'include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()


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
    parser.add_argument('--output', '-o', type=str, default='.',
                        help='Output directory for plots')
    parser.add_argument('--test', action='store_true',
                        help='Test mode: use limited files')
    parser.add_argument('--nfiles', type=int, default=1,
                        help='Number of files to use in test mode')
    return parser.parse_args()


def load_samples(datfile, lumi, data_names, mc_names, file_limit=0):
    """
    Load data and MC samples from dat file.

    Args:
        datfile: Path to sample definition file
        lumi: Integrated luminosity
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


def plot_comparison(histoMC, histoDATA, output_dir):
    """
    Create and save Data/MC comparison plot.

    Args:
        histoMC: MC histogram
        histoDATA: Data histogram
        output_dir: Directory to save plots
    """
    canvas = Canvas.Canvas(
        'DiMuon_mass', 'png,pdf',
        0.16, 0.72, 0.56, 0.82, 1
    )
    canvas.addHisto(histoMC, 'HIST', 'MC (DY)', 'l', r.kBlue, True, 0)
    canvas.addHisto(histoDATA, 'P, SAME', 'Data', 'p', r.kBlack, True, 1)
    canvas.save(1, 0, 0, '', '', outputDir=output_dir, isPrivate=False, xlog=False)

    print(f"Plot saved to {output_dir}/DiMuon_mass.png")


def main():
    """Main execution."""
    args = parse_args()

    # Setup
    setup_style()
    helper.ensureDirectory(args.output)

    # Sample configuration
    datfile = WORKPATH + args.dat
    data_samples = ['Muon0_Run2024I_v1']
    mc_samples = ['DYJetsToLL_M-50']

    # Variable definitions for dimuon mass calculation
    definitions_file = 'definitions/defs_MuonVariables.json'

    print(f"{'='*60}")
    print(f"Galapago-FWK: Dimuon Mass Analysis")
    print(f"{'='*60}")
    print(f"Dat file:    {args.dat}")
    print(f"Luminosity:  {args.lumi} fb^-1")
    print(f"Output:      {args.output}")
    print(f"Test mode:   {args.test}")
    print(f"{'='*60}")

    # Load samples
    print("\nLoading samples...")
    treeDATA, treeMC = load_samples(
        datfile, args.lumi,
        data_samples, mc_samples,
        file_limit=args.nfiles if args.test else 0
    )

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

    # Create plot
    print("\nCreating comparison plot...")
    plot_comparison(histoMC, histoDATA, args.output)

    print("\nDone!")


if __name__ == "__main__":
    main()
