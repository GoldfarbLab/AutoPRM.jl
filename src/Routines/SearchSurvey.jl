function SearchSurvey(params_path::String)

    [include(joinpath(package_root, 
    "src", "Routines","PRM","IS-PRM-SURVEY", jl_file)
    ) for jl_file in ["plotPRM.jl"]]

[include(joinpath(package_root, 
    "src", "utils", jl_file)
    ) for jl_file in ["writeTables.jl"]]

	[include(joinpath(package_root, 
    "src", "structs", jl_file)
    ) for jl_file in ["buildPrecursorTable.jl"]]


    params = JSON.parse(read(params_path, String));
    MS_DATA_DIR = params["ms_data_dir"];
	PRECURSOR_LIST_PATH = params["peptide_list_path"];
    #SPEC_LIB_DIR = params["library_folder"];
    MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) for file in filter(file -> isfile(joinpath(MS_DATA_DIR, file)) && match(r"\.arrow$", file) != nothing, readdir(MS_DATA_DIR))]
    println("Processing: "*string(length(MS_TABLE_PATHS))*" files")


    function parse_mods(fixed_mods)
        fixed_mods_parsed = Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}()
        for mod in fixed_mods
            push!(fixed_mods_parsed, (p=Regex(mod[1]), r = mod[2]))
        end
        return fixed_mods_parsed
    end

    #Parse argments
    params = (
    right_precursor_tolerance = Float64(params["right_precursor_tolerance"]),
    left_precursor_tolerance = Float64(params["left_precursor_tolerance"]),
    precursor_rt_tolerance = Float64(params["precursor_rt_tolerance"]),
    b_ladder_start = Int64(params["b_ladder_start"]),
    y_ladder_start = Int64(params["y_ladder_start"]),
    precursor_charges = [UInt8(charge) for charge in params["precursor_charges"]],
    precursor_isotopes = [UInt8(isotope) for isotope in params["precursor_isotopes"]],
    transition_charges = [UInt8(charge) for charge in params["transition_charges"]],
    transition_isotopes = [UInt8(isotope) for isotope in params["transition_isotopes"]],
    fragment_match_ppm = Float64(params["fragment_match_ppm"]),
    minimum_fragment_count = UInt8(params["minimum_fragment_count"]),
    fragments_to_select = UInt8(params["fragments_to_select"]),
    precursort_rt_window = Float64(params["precursor_rt_window"]),
    max_variable_mods = Int(params["max_variable_mods"]),
    fixed_mods = parse_mods(params["fixed_mods"]),
    variable_mods = parse_mods(params["variable_mods"]),
    modification_masses = Dict{String, Float64}(k => Float64(v) for (k, v) in params["modification_masses"]),
    ms_file_conditions = params["ms_file_conditions"]
    )

    ##########
    #Read Precursor Table
    ##########
    @time begin

        ptable = PrecursorTable()
        buildPrecursorTable!(ptable, 
                            params[:fixed_mods], 
                            params[:variable_mods], 
                            params[:max_variable_mods], 
                            PRECURSOR_LIST_PATH)
        addPrecursors!(
                            ptable, 
                            params[:precursor_charges], 
                            params[:precursor_isotopes], 
                            params[:modification_masses]
                            )
        addTransitions!(
                            ptable, 
                            params[:transition_charges],
                            params[:transition_isotopes],
                            params[:b_ladder_start],
                            params[:y_ladder_start],
                            params[:fragment_match_ppm]
                        )
        ##########
        #Search Survey Runs
        ##########
        MS_TABLES = Dict{UInt32, Arrow.Table}()
        combined_scored_psms = makePSMsDict(FastXTandem())
        combined_fragment_matches = Dict{UInt32, Vector{FragmentMatch}}()
        MS_RT = Dict{UInt32, Vector{Float32}}()
            for (ms_file_idx, MS_TABLE_PATH) in enumerate(MS_TABLE_PATHS)
        
                MS_TABLE = Arrow.Table(MS_TABLE_PATH)
                MS_RT[UInt32(ms_file_idx)] = MS_TABLE[:retentionTime]
                scored_psms, fragment_matches = SearchRAW(
                                                        MS_TABLE, 
                                                        ptable, 
                                                        selectTransitions, 
                                                        params[:right_precursor_tolerance],
                                                        params[:left_precursor_tolerance],
                                                        params[:transition_charges],
                                                        params[:transition_isotopes],
                                                        params[:b_ladder_start],
                                                        params[:y_ladder_start],
                                                        params[:fragment_match_ppm],
                                                        UInt32(ms_file_idx)
                                                        )
                for key in keys(combined_scored_psms)
                    append!(combined_scored_psms[key], scored_psms[key])
                end
                combined_fragment_matches[UInt32(ms_file_idx)] = fragment_matches
            end
        
        ##########
        #Get Best PSMs for Each Peptide
        ##########
            best_psms = surveyGetBestPSMs(combined_scored_psms, ptable, MS_RT, params[:minimum_fragment_count])
        ##########
        #Get MS1 Peak Heights
        ##########
            #First key is ms_file_idx (identifier of the ms file), second key is pep_idx (peptide id)
            ms1_peak_heights = UnorderedDictionary{UInt32, UnorderedDictionary{UInt32, Float64}}()
            #Peak heights are zero to begin with
            precursor_idxs = unique(best_psms[!,:precursor_idx])
            for (ms_file_idx, MS_TABLE_PATH) in enumerate(MS_TABLE_PATHS)
                MS_TABLE = Arrow.Table(MS_TABLE_PATH)
                insert!(ms1_peak_heights, 
                        UInt32(ms_file_idx), 
                        UnorderedDictionary(precursor_idxs, zeros(Float64, length(precursor_idxs)))
                        )
        
                getMS1PeakHeights!( ms1_peak_heights[ms_file_idx],
                                    ptable,
                                    MS_TABLE[:retentionTime], 
                                    MS_TABLE[:mz_array], 
                                    MS_TABLE[:intensity_array], 
                                    MS_TABLE[:msOrder],
                                    best_psms[!,:retention_time], 
                                    best_psms[!,:precursor_idx], 
                                    best_psms[!,:ms_file_idx],
                                    Float32(0.25), 
                                    params[:right_precursor_tolerance], 
                                    params[:left_precursor_tolerance],
                                    UInt32(ms_file_idx))
            end
        
            #Add MS1 Heights to the best_psms DataFrame 
            transform!(best_psms, AsTable(:) => ByRow(psm -> ms1_peak_heights[psm[:ms_file_idx]][psm[:precursor_idx]]) => :ms1_peak_height)
        ##########
        #Get Chromatograms for the best precursors in each file. 
        ##########
            precursor_chromatograms = UnorderedDictionary{UInt32, UnorderedDictionary{UInt32, PrecursorChromatogram}}()
            for (ms_file_idx, MS_TABLE_PATH) in enumerate(MS_TABLE_PATHS)
                MS_TABLE = Arrow.Table(MS_TABLE_PATH)
                insert!(precursor_chromatograms, UInt32(ms_file_idx), initPrecursorChromatograms(best_psms, UInt32(ms_file_idx)) |> (best_psms -> fillPrecursorChromatograms!(best_psms, 
                                                                                                                            combined_fragment_matches[UInt32(ms_file_idx)], 
                                                                                                                            MS_TABLE, 
                                                                                                                            params[:precursor_rt_tolerance],
                                                                                                                            UInt32(ms_file_idx))
                                                                                                                )
                        ) 
            end 
        
            #Names and charges for the "n" most intense fragment ions for each precursor
            #display(getBestPSM(precursor_chromatograms[psm[:ms_file_idx]][psm[:precursor_idx]]))
            #println(typeof(getBestPSM(precursor_chromatograms[psm[:ms_file_idx]][psm[:precursor_idx]])))
        
        
            transform!(best_psms, AsTable(:) => ByRow(psm -> surveyGetBestTransitions(getBestPSM(precursor_chromatograms[psm[:ms_file_idx]][psm[:precursor_idx]]),
                                                                                maximum_fragment_count = params[:fragments_to_select])) => :best_transitions)
            transform!(best_psms, AsTable(:) => ByRow(psm -> getBestPSM(precursor_chromatograms[psm[:ms_file_idx]][psm[:precursor_idx]])[:name][psm[:best_transitions]]) => :transition_names)
            transform!(best_psms, AsTable(:) => ByRow(psm -> getBestPSM(precursor_chromatograms[psm[:ms_file_idx]][psm[:precursor_idx]])[:mz][psm[:best_transitions]]) => :transition_mzs)
        
            function sumTopN(intensities::Vector{T}) where {T<:AbstractFloat}
                topN = min(3, length(intensities))
                return sum(sort(intensities, rev = true)[1:topN])
            end
        
            function NthIntensity(intensities::Vector{T}) where {T<:AbstractFloat}
                topN = min(3, length(intensities))
                return sort(intensities, rev = true)[topN]
            end
            transform!(best_psms, AsTable(:) => ByRow(psm -> sumTopN(getBestPSM(precursor_chromatograms[psm[:ms_file_idx]][psm[:precursor_idx]])[:intensity])) => :sumTopN)
            transform!(best_psms, AsTable(:) => ByRow(psm -> NthIntensity(getBestPSM(precursor_chromatograms[psm[:ms_file_idx]][psm[:precursor_idx]])[:intensity])) => :NthIntensity)
        
        ##########
        #Apply conditions to MS Files. 
        ##########
        MS_FILE_ID_TO_NAME = Dict(
                                    zip(
                                    [UInt32(i) for i in 1:length(MS_TABLE_PATHS)], 
                                    [splitpath(filepath)[end] for filepath in MS_TABLE_PATHS]
                                    )
                                )
        transform!(best_psms, AsTable(:) => ByRow(psm -> MS_FILE_ID_TO_NAME[psm[:ms_file_idx]]) => :file_name)
        
            MS_FILE_ID_TO_CONDITION = Dict(
                                            zip(
                                            [key for key in keys(MS_FILE_ID_TO_NAME)], 
                                            ["NONE" for key in keys(MS_FILE_ID_TO_NAME) ]
                                            )
                                        )
        
        for (file_id, file_name) in MS_FILE_ID_TO_NAME 
            for (condition, value) in params[:ms_file_conditions]
                if occursin(condition, file_name)
                    MS_FILE_ID_TO_CONDITION[file_id] = condition
                end
            end
        end
        
        transform!(best_psms, AsTable(:) => ByRow(psm -> MS_FILE_ID_TO_CONDITION[psm[:ms_file_idx]]) => :condition)
        ##########
        #Write Method Files
        ##########
        #Get best_psm for each peptide across all ms_file
        surveyWriteIAPIMethod(best_psms, joinpath(MS_DATA_DIR, "precursors_summary.csv"))
        best_psms = combine(sdf -> sdf[argmax(sdf.hyperscore), :], groupby(best_psms, :pep_idx)) 
    
        filter!(row -> (row.total_ions > params[:minimum_fragment_count]), best_psms)
    
        surveyWriteTransitionList(best_psms, joinpath(MS_DATA_DIR, "transition_list.csv"))
        surveyWriteIAPIMethod(best_psms, joinpath(MS_DATA_DIR, "iapi_method.csv"))
        
        println(" Scored "*string(size(best_psms)[1])*" precursors")
        ##########
        #Make Plots
        ##########
        if true==true
        #if ARGS["make_plots"] != "false"
            @time for (ms_file_idx, MS_TABLE_PATH) in enumerate(MS_TABLE_PATHS)
                MS_TABLE = Arrow.Table(MS_TABLE_PATH)
                #sample_name = split(MS_FILE_ID_TO_NAME[ms_file_idx], ".")[1]
                sample_name = split(basename(MS_FILE_ID_TO_NAME[ms_file_idx]), ".")[1]
                println("sample name ", sample_name)
                best_spectra_path = ""
                precursor_chromatograms_path = ""
                if Sys.iswindows()
                    best_spectra_path = joinpath(":c","figures","best_spectra", sample_name)
                    precursor_chromatograms_path = joinpath(":c","figures","precursor_chromatogram", sample_name)
                else
                    best_spectra_path = joinpath("./","figures","best_spectra", sample_name)
                    precursor_chromatograms_path = joinpath("./","figures","precursor_chromatogram", sample_name)
                end
        
                plotAllBestSpectra(precursor_chromatograms[ms_file_idx], 
                                    ptable, 
                                    MS_TABLE,
                                    best_spectra_path,
                                    join(sample_name*"_best_spectra"*".pdf"))
                
                plotAllFragmentIonChromatograms(precursor_chromatograms[ms_file_idx], 
                                                ptable,
                                                precursor_chromatograms_path,
                                                join(sample_name*"_precursor_chromatograms"*".pdf"))
            end
        end
    end
end