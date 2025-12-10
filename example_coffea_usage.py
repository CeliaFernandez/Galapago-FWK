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
        name='Muon_Run2022G_NanoAODv15',
        dataset='/Muon/Run2022G-NanoAODv15-v1/NANOAOD',
        isdata=True,
        file_limit=1,  # Use 2 files for testing
        label='Data',
        color='#000000'
    )

    # Create MC sample
    print("\nLoading MC sample...")
    mc = Sample(
        name='DYJetsToLL',
        dataset='/DYto2Mu-2Jets_Bin-MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v6/NANOAODSIM',
        isdata=False,
        xsection=6077.22,  # pb
        file_limit=1,
        label='Drell-Yan',
        color='#f18f01'
    )

    # Apply selections
    print("\nApplying selections...")

    # Dimuon selection
    data.addSelection('ak.num(events.Muon) >= 2')
    data.addSelection('events.Muon.charge[:, 0] * events.Muon.charge[:, 1] < 0')  # Opposite charge

    mc.addSelection('ak.num(events.Muon) >= 2')
    mc.addSelection('events.Muon.charge[:, 0] * events.Muon.charge[:, 1] < 0')

    # Print info
    data.printInfo()
    mc.printInfo()

    # Compute dimuon mass using auxiliar function
    print("\nComputing dimuon invariant mass...")
    data.computeDimuonMass()
    mc.computeDimuonMass()

    # Get yields
    print("\nYields after selections:")
    data_yield, data_unc = data.getYield(lumi)
    mc_yield, mc_unc = mc.getYield(lumi)
    print(f"Data: {data_yield:.0f} ± {data_unc:.0f}")
    print(f"MC:   {mc_yield:.1f} ± {mc_unc:.1f}")

    # Create histograms
    print("\nCreating dimuon mass histograms...")

    # Dimuon invariant mass (Z peak region)
    h_data = data.getHist(
        var='events.DiMuon_mass',
        bins=100,
        range=(60, 120),
        lumi=lumi
    )

    h_mc = mc.getHist(
        var='events.DiMuon_mass',
        bins=100,
        range=(60, 120),
        lumi=lumi
    )

    # Plot comparison (absolute yields)
    print("\nCreating comparison plot...")
    plotComparison(
        h_data=h_data,
        h_mc=h_mc,
        xlabel='Dimuon Invariant Mass $m_{\\mu\\mu}$ [GeV]',
        ylabel='Events / 0.6 GeV',
        data_label='Data 2024',
        mc_label='DY $\\rightarrow \\mu\\mu$',
        output='plots/coffea_dimuon_mass.png',
        ratio=True
    )

    # Create normalized histograms for shape comparison
    print("\nCreating normalized histograms...")
    h_data_norm = data.getHist(
        var='events.DiMuon_mass',
        bins=100,
        range=(60, 120),
        lumi=lumi,
        normalize=True
    )

    h_mc_norm = mc.getHist(
        var='events.DiMuon_mass',
        bins=100,
        range=(60, 120),
        lumi=lumi,
        normalize=True
    )

    # Plot normalized comparison
    print("\nCreating normalized comparison plot...")
    plotComparison(
        h_data=h_data_norm,
        h_mc=h_mc_norm,
        xlabel='Dimuon Invariant Mass $m_{\\mu\\mu}$ [GeV]',
        ylabel='Normalized to Unity',
        data_label='Data 2024',
        mc_label='DY $\\rightarrow \\mu\\mu$',
        output='plots/coffea_dimuon_mass_normalized.png',
        ratio=True,
        normalized=True
    )

    print("\n" + "=" * 60)
    print("Done!")
    print("=" * 60)


if __name__ == '__main__':
    main()
