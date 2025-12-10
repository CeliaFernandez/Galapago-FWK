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
    """
    # Get muon 4-vectors (using coffea vector operations)
    mu1 = sample.events.Muon[:, 0]
    mu2 = sample.events.Muon[:, 1]

    # Calculate invariant mass using 4-vector addition
    dimuon = mu1 + mu2
    mass = dimuon.mass

    # Add to events
    sample.events = ak.with_field(sample.events, mass, 'DiMuon_mass')

    return sample
