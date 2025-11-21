//
// LIBRARY FOR SNT SCOUTING DIMUON ANALYSIS
// @fernance (25 Jul 2023)
//


// Function to make dimuon vertices out of the associated muons
// OS muon pairs
auto MakeDimuonVertices(ROOT::RVec<int> SVOv, ROOT::RVec<int> SV, ROOT::RVec<int> ch) {
  int ovidx;
  int vidx;
  ROOT::RVec<int> dmuidxs;
  for (int m = 0; m < ch.size(); m++) {
    ovidx = -1;
    vidx = -1;
    if (SVOv.at(m) > -1) {
      ovidx = SVOv.at(m);
    } else if (SV.at(m) > -1) {
      vidx = SV.at(m);
    }
    for (int mm = 0; mm < ch.size(); mm++) {
      if (m==mm) continue;
      if (ch.at(m)*ch.at(mm) > 0) continue;
        if ( !(ovidx > -1 && ovidx == SVOv.at(mm)) && vidx > -1 && vidx == SV.at(mm)) {
        if (!(std::find(dmuidxs.begin(), dmuidxs.end(), m) != dmuidxs.end() || std::find(dmuidxs.begin(), dmuidxs.end(), mm) != dmuidxs.end())) {
          dmuidxs.push_back(m);
          dmuidxs.push_back(mm);
        }
      }
    }
  }
  return dmuidxs;
}

// Function to compute invariant mass from vector of TLorentzVectors
auto GetVectorMass(ROOT::RVec<TLorentzVector> x) {
  ROOT::RVec<double> mass;
  for (auto m:x)
    mass.push_back(m.M());
  return mass;
}

// Function to create a new RVec from selected RVec indexes
template <typename TYPE>
inline auto GetElements(ROOT::RVec<TYPE> v, ROOT::RVec<int> idx) {
  ROOT::RVec<TYPE> out;
  for (auto x:idx)
    out.push_back(v.at(x));
  return out;
}


