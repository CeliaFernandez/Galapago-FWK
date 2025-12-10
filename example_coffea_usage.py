#!/usr/bin/env python3
"""
Example usage of the simplified coffea-based Sample class.
Demonstrates how to create dimuon mass plots with data and MC.
"""

from include.Sample_coffea import Sample, plotComparison

def main():
    print("=" * 60)
    print("Coffea Sample Example: Dimuon Mass Analysis")
    print("=" * 60)

    # Configuration
    lumi = 35.9  # fb^-1

    # Create Data sample
    print("\nLoading Data sample...")
    data = Sample(
        name='Muon0_Run2024I',
        dataset='/Muon0/Run2024I-PromptNanoAODv12-v1/NANOAOD',
        isdata=True,
        file_limit=2,  # Use 2 files for testing
        label='Data',
        color='#000000'
    )

    # Create MC sample
    print("\nLoading MC sample...")
    mc = Sample(
        name='DYJetsToLL',
        dataset='/DYJetsToLL_M-50/Run3Summer23-PromptNanoAODv12-v1/NANOAOD',
        isdata=False,
        xsection=6077.22,  # pb
        file_limit=2,
        label='Drell-Yan',
        color='#f18f01'
    )

    # Apply selections
    print("\nApplying selections...")

    # Dimuon selection
    data.addSelection('events.nMuon >= 2')
    data.addSelection('events.Muon.charge[:, 0] * events.Muon.charge[:, 1] < 0')  # Opposite charge

    mc.addSelection('events.nMuon >= 2')
    mc.addSelection('events.Muon.charge[:, 0] * events.Muon.charge[:, 1] < 0')

    # Print info
    data.printInfo()
    mc.printInfo()

    # Get yields
    print("\nYields after selections:")
    data_yield, data_unc = data.getYield(lumi)
    mc_yield, mc_unc = mc.getYield(lumi)
    print(f"Data: {data_yield:.0f} ± {data_unc:.0f}")
    print(f"MC:   {mc_yield:.1f} ± {mc_unc:.1f}")

    # Create histograms
    print("\nCreating dimuon mass histograms...")

    # Note: You'll need to define DiMuon_mass if it's not already in the events
    # For now, using Muon_pt as an example
    h_data = data.getHist(
        var='events.Muon.pt[:, 0]',  # Leading muon pT
        bins=100,
        range=(0, 200),
        lumi=lumi
    )

    h_mc = mc.getHist(
        var='events.Muon.pt[:, 0]',
        bins=100,
        range=(0, 200),
        lumi=lumi
    )

    # Plot comparison
    print("\nCreating comparison plot...")
    plotComparison(
        h_data=h_data,
        h_mc=h_mc,
        xlabel='Leading Muon $p_T$ [GeV]',
        ylabel='Events / 2 GeV',
        data_label='Data 2024',
        mc_label='DY $\\rightarrow \\mu\\mu$',
        output='plots/coffea_muon_pt.png',
        ratio=True
    )

    print("\n" + "=" * 60)
    print("Done!")
    print("=" * 60)


if __name__ == '__main__':
    main()
