"""
    PrecursorChromatogram

Type that represents an chromatogram for a `Precursor`

### Fields

- rts::Vector{Float32} -- Retention times for scans targeting this precursor
- last_scan_idx::Vector{Int64} --  Keeps track of the last_scan from which a transition has been added (should change to a Ref in the future)
- transitions::UnorderedDictionary{String, Vector{Float32}} -- Dictionary with keys for each transition, and values for their intensities in 
each matched scan. Those vectors should all have the same length and correspond to the `rts` field. 
- best_psm::NamedTuple{(:rt, :scan_idx, :name, :mz, :intensity), Tuple{Float32, Int64, Vector{String}, Vector{Float32}, Vector{Float32} -- List of 
transitions, mz's, and their intensities for the best scoring PSM for the precursor.  

### Examples

- PrecursorChromatogram() = PrecursorChromatogram(Float32(0), 0, Float32[], UnorderedDictionary{String, Vector{Float32}}(), (rt = Float32(0), scan_idx = 0, name = String[], mz = Float32[], intensity = Float32[]))

### GetterMethods

- getRTs(pc::PrecursorChromatogram) = pc.rts
- getLastScanIdx(pc::PrecursorChromatogram) = pc.last_scan_idx
- getTransitions(pc::PrecursorChromatogram) = pc.transitions
- getBestPSM(pc::PrecursorChromatogram) = pc.best_psm

### Methods

- initMatchedPrecursors(bestPSMs::DataFrame) -- Initializes a dictionary container of `PrecursorChromatogram` from a DataFrame of best PSMs. 
"""
struct PrecursorChromatogram
    rts::Vector{Float32}
    scan_idxs::Vector{Int64}
    transitions::UnorderedDictionary{String, Vector{Float32}}
    best_psm::NamedTuple{(:rt, :scan_idx, :name, :mz, :intensity), Tuple{Float32, Int64, Vector{String}, Vector{Float32}, Vector{Float32}}}
end

PrecursorChromatogram() = PrecursorChromatogram(Float32[], Int64[], UnorderedDictionary{String, Vector{Float32}}(), (rt = Float32(0), scan_idx = 0, name = String[], mz = Float32[], intensity = Float32[]))
getRTs(pc::PrecursorChromatogram) = pc.rts
getLastScanIdx(pc::PrecursorChromatogram) = pc.scan_idxs[end]
getTransitions(pc::PrecursorChromatogram) = pc.transitions
getBestPSM(pc::PrecursorChromatogram) = pc.best_psm