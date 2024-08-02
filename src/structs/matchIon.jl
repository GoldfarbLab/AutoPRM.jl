"""
    FragmentMatch

Type that represents a match between a fragment ion and a mass spectrum peak

### Fields

- transition::Transition -- Represents a fragment ion
- intensity::Float32 -- Intensity of the matching peak
- match_mz::Float32 -- M/Z of the matching empirical peak. 
    NOT the transition mass and may differ from getMZ(transition) by some error
- count::UInt8 -- Number of matches (may want to count the number of matches if the spectrum is researched at 
    different mass offsets, such as in a cross correlation score)
- peak_int::Int64 -- Index of the matching peak in the mass spectrum. 

### Examples

- `FragmentMatch(transition::Transition, intensity::Float32, mass::Float32, count::UInt8, peak_ind::Int64) --
    default internal constructor
- `FragmentMatch()` -- constructor for null/empty precursor

### GetterMethods

- getMZ(f::FragmentMatch) = getMZ(f.transition)
- getLow(f::FragmentMatch) = getLow(f.transition)
- getHigh(f::FragmentMatch) = getHigh(f.transition)
- getPrecID(f::FragmentMatch) = getPrecID(f.transition)
- getCharge(f::FragmentMatch) = getCharge(f.transition)
- getIsotope(f::FragmentMatch) = getIsotope(f.transition)
- getIonType(f::FragmentMatch) = getIonType(f.transition)
- getInd(f::FragmentMatch) = getInd(f.transition)
"""
mutable struct FragmentMatch{T<:AbstractFloat}
    transition::Transition
    intensity::T
    match_mz::T
    count::UInt8
    peak_ind::Int64
    scan_idx::UInt32
    ms_file_idx::UInt32
end
export FragmentMatch

FragmentMatch() = FragmentMatch(Transition(), Float64(0), Float64(0), UInt8(0), 0, UInt32(0), UInt32(0))
getMZ(f::FragmentMatch) = getMZ(f.transition)
getLow(f::FragmentMatch) = getLow(f.transition)
getHigh(f::FragmentMatch) = getHigh(f.transition)
getPrecID(f::FragmentMatch) = getPrecID(f.transition)
getCharge(f::FragmentMatch) = getCharge(f.transition)
getIsotope(f::FragmentMatch) = getIsotope(f.transition)
getIonType(f::FragmentMatch) = getIonType(f.transition)
getInd(f::FragmentMatch) = getInd(f.transition)
getPeakInd(f::FragmentMatch) = f.peak_ind
getIntensity(f::FragmentMatch) = f.intensity
getCount(f::FragmentMatch) = f.count
getMSFileID(f::FragmentMatch) = f.ms_file_idx

