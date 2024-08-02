function importScripts(package_root::String)

    [include(joinpath(package_root, "src", "structs", jl_file)) for jl_file in ["precursor.jl",
                                                            "precursorChromatogram.jl",
                                                            "precursorContainers.jl",                                                   
                                                            "PSM.jl",
                                                            "surveyPrecursor.jl",
                                                            "matchIon.jl",
                                                            "FastXTandem.jl",
                                                            "buildPrecursorTable.jl"]]
    [include(joinpath(package_root, "src", "utils", jl_file)) for jl_file in ["applyMods.jl",
                                                            "getMS1PeakHeights.jl",
                                                            "binaryRangeQuery.jl",
                                                            "LFQ.jl",                                                   
                                                            "matchpeaks.jl",
                                                            "precursorChromatogram.jl",
                                                            "searchRAW.jl",
                                                            "plotPRM.jl",
                                                            "selectTransitions.jl",
                                                            "getBestPSMs.jl",
                                                            "getBestTransitions.jl",
                                                            "initTransitions.jl",
                                                            "getScanPairs.jl",
                                                            "parEstimation.jl",
                                                            "writeTables.jl",
                                                            "buildPrecursorTable.jl"]]

end