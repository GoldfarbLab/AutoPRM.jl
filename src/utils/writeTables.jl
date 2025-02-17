function surveyWriteTransitionListSureQuant(best_psms::DataFrame, f_out::String)

    for precursor_z in unique(best_psms[!, :precursor_charge])
        for heavy_AA in ["Hlys", "Harg"]
            sub_psms = filter(row -> row.precursor_charge == precursor_z && occursin(heavy_AA, row.sequence), best_psms)
            open(f_out * "_" * heavy_AA * "_" * string(precursor_z) * ".csv", "w") do io
                # loop over data and write each line
                write(io, join(["Compound","m/z","Group ID"],",")*"\n")
                for row in eachrow(sub_psms)
                    for i in range(1, length(row[:transition_names]))
                        data = push!([row[:transition_names][i], 
                                    row[:transition_mzs][i],
                                    row[:precursor_mz]])
                        write(io, join(data,",")*"\n")
                    end
                end
            end
        end
    end
end

function surveyWriteSILsSureQuant(best_psms::DataFrame, f_out::String)
    for precursor_z in unique(best_psms[!, :precursor_charge])
        for heavy_AA in ["Hlys", "Harg"]
            sub_psms = filter(row -> row.precursor_charge == precursor_z && occursin(heavy_AA, row.sequence), best_psms)
            open(f_out * "_" * heavy_AA * "_" * string(precursor_z) * ".csv", "w") do io
                # loop over data and write each line
                write(io, join(["Compound","m/z","Intensity Threshold"],",")*"\n")
                for row in eachrow(sub_psms)
                    data = push!([row[:sequence], 
                                row[:precursor_mz],
                                "10000"])
                    write(io, join(data,",")*"\n")
                end
            end
        end
    end
end

function surveyWriteTransitionList(best_psms::DataFrame, f_out::String)
    open(f_out, "w") do io
        # loop over data and write each line
        write(io, join(["protein_name","sequence","precursor_charge","precursor_isotope","transition_names"],",")*"\n")
        for row in eachrow(best_psms)

            #data = join(append!([row[:proteinNames]*","*row[:sequence]], row[:names]),",")
            data = push!([row[:protein_names], 
                            row[:sequence],
                            row[:precursor_charge],
                            row[:precursor_isotope]], #replace(row[:sequence], r"\[(.*?)\]" => "")
                            join(row[:transition_names],";"))
            write(io, join(data,",")*"\n")
        end
    end
end

function writeTransitionList(best_psms::DataFrame, f_out::String)
    open(f_out, "w") do io
        # loop over data and write each line
        for row in eachrow(best_psms)

            #data = join(append!([row[:proteinNames]*","*row[:sequence]], row[:names]),",")
            data = append!([row[:protein_name], 
                            row[:sequence]], #replace(row[:sequence], r"\[(.*?)\]" => "")
                            row[:transition_names])
            write(io, join(data,",")*"\n")
        end
    end
end

function surveyWriteIAPIMethod(best_psms::DataFrame, f_out::String)
    open(f_out, "w") do io
        # loop over data and write each line
        write(io, join(["protein_name","sequence","precursor_mz","precursor_charge", "retention_time","precursor_intensity","hyperscore","NthIntensity","sumTopN","file_name","condition","transition_mz"],",")*"\n")
        for row in eachrow(best_psms)

            #data = join(append!([row[:proteinNames]*","*row[:sequence]], row[:names]),",")
            data = [row[:protein_names], row[:sequence], row[:precursor_mz], row[:precursor_charge], row[:retention_time], row[:ms1_peak_height], row[:hyperscore], row[:NthIntensity], row[:sumTopN], row[:file_name], row[:condition]]#, row[:transition_mzs])
            write(io, join(data,",")*","*join(row[:transition_mzs], ";")*"\n")
        end
    end
end


function writeIAPIMethod(best_psms::DataFrame, f_out::String)
    open(f_out, "w") do io
        # loop over data and write each line
        write(io, join(["protein_name","sequence","precursor_mz","precursor_intensity","condition","transition_mz"],",")*"\n")
        for row in eachrow(best_psms)

            #data = join(append!([row[:proteinNames]*","*row[:sequence]], row[:names]),",")
            data = append!([row[:protein_name], row[:sequence], row[:precursor_mz], row[:retention_time], row[:ms1_peak_height], row[:condition]], row[:transition_mzs])
            write(io, join(data,",")*"\n")
        end
    end
end