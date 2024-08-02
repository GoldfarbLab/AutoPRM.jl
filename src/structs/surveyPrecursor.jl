struct SurveyPrecursor{T<:AbstractFloat,U<:Unsigned}
    Compound::String
    mz::T
    z::U
    intensity_threshold::T
end