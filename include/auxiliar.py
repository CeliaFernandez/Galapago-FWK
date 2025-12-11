"""
Auxiliary functions for Sample_coffea module.
Contains helper functions for computing derived variables.
"""

import awkward as ak


def computeDimuonMass(sample):
    """
    Compute dimuon invariant mass from two leading muons.
    Adds 'DiMuon_mass' field to events.
    Uses 4-momentum: mass = sqrt((p1 + p2)^2)

    Args:
        sample: Sample instance with events containing Muon collection

    Returns:
        The sample instance (for method chaining)

    Note:
        Events with fewer than 2 muons will have DiMuon_mass = -999
    """
    import numpy as np

    # Create mask for events with at least 2 muons
    has_dimuon = ak.num(sample.events.Muon) >= 2

    # Initialize mass array with -999 (invalid value)
    mass = ak.where(
        has_dimuon,
        -999.0,  # placeholder, will be overwritten
        -999.0   # events without 2 muons
    )

    # For events with >= 2 muons, compute the mass
    # Get muon 4-vectors (using coffea vector operations)
    mu1 = sample.events.Muon[has_dimuon][:, 0]
    mu2 = sample.events.Muon[has_dimuon][:, 1]

    # Calculate invariant mass using 4-vector addition
    dimuon = mu1 + mu2
    dimuon_mass = dimuon.mass

    # Fill in the computed masses
    mass = ak.where(has_dimuon, dimuon_mass, -999.0)

    # Add to events
    sample.events = ak.with_field(sample.events, mass, 'DiMuon_mass')

    return sample
