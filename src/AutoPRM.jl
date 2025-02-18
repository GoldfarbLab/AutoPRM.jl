module AutoPRM 
__precompile__(false)
##########
#Load Dependencies 
##########
using ArgParse, Arrow
using Base.Iterators
using Combinatorics, CSV
using DataFrames, Dictionaries
using JSON
using Missings
using PDFmerger, Plots, PrettyPrinting, Printf, ProgressBars
using RobustModels
using SpecialFunctions, Statistics, Suppressor
using Tables, Test
using DataFramesMeta

const package_root = dirname(@__DIR__);
include("utils/importScripts.jl")
importScripts(package_root)
include(joinpath(package_root, "src","Routines","BuildSurvey.jl"))
include(joinpath(package_root, "src","Routines","SearchSurvey.jl"))
include(joinpath(package_root, "src","Routines","SearchPRM.jl"))
export SearchPRM, SearchSurvey, BuildSurvey

end
