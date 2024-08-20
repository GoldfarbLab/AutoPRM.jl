# AutoPRM: Method Preparation and Data Analysis for Internal Standard Triggered Parallel Reaction Monitoring (IS-PRM). 

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://nwamsley1.github.io/Titus.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://nwamsley1.github.io/Titus.jl/dev/)
[![Build Status](https://github.com/nwamsley1/Titus.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/nwamsley1/Titus.jl/actions/workflows/CI.yml?query=branch%3Amain)

[![Coverage](https://codecov.io/gh/nwamsley1/AutoPRM.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/nwamsley1/AutoPRM.jl)
## Aims
  Titus aims to simplify preparation and analysis of IS-PRM experiments [1]. Supports analysis of survey runs to characterize internal standard peptides, preparation of methods, peak area ratio estimation and protein-level quantitation.

## Installation
1) AutoPRM requires Julia 1.10. Download [julia](https://pages.github.com/) and add it to the PATH. 
2) Open an instance of the julia REPL
3) Type ";" to activate the shell from within the REPL. Then, navigate to the desired directory and clone the AutoPRM.jl repository.
```
shell> git clone https://github.com/nwamsley1/AutoPRM.jl.git
```
and not move into the package directory
```
shell> cd AutoPRM.jl
```
4) Return to julia by hitting the backspace key. Activate the julia package manager by typing "]" into the REPL and enter the following:
```
(@v1.10) pkg> activate
(@v1.10) pkg> develop ./
(@v1.10) pkg> add ./
```

## Usage
AutoPRM exports three methods "SearchPRM", "SearchSurvey", and "BuildSurvey". Each takes a single argument, that is a file path to a .json parameters files (examples included). To access these methods
import the AutoPRM package from within the Julia REPL. 
```
julia> using AutoPRM
```
### File Conversion
`AutoPRM` requires Thermo .raw files be converted to an Apache .arrow format with a specific column specification. Use [PioneerConverter](https://github.com/nwamsley1/PioneerConverter) to convert .raw files. 

### BuildSurvey
BuildSurvey constructs a table of precursors, m/z ratios, charge states, and intensity thresholds needed to build a SurveyMethod to analyze a new plate or aliquot of SIL peptides. This table can be directly loaded
into the Thernmo method editor. 

##### Peptide List File
BuildSurvey requires takes a tab-delimited petide list file with at least two columns, "PROTEIN" and "PRECURSOR". The PROTEIN column is the Gene name, Uniprot ID, etc. that you want to associate with the PRECURSOR.
The PRECURSOR is the peptide sequence using one of the 20 standard single-letter amino-acid abbreviasions. On Mac/Linux, you can inspect an example file "HUMAN_IMMUNOLOGY_PEPTIDE_LIST_SEP_2023.txt" with the `head` command:
```
╭─n.t.wamsley@3225-AD-00020.local ~/TEST_DATA/SUREQUANT_JUL24/HUMAN_IMMUNOLOGY  
╰─➤  head HUMAN_IMMUNOLOGY_PEPTIDE_LIST_SEP_2023.txt 
PROTEIN	PEPTIDE	VENDOR	ID	WELL
ARG1	GGVEEGPTVLR	Thermo	NR102342.1	A1
ARG1	VMEETLSYLLGR	Thermo	NR102342.2	A2
AXL	APLQGTLLGYR	Thermo	NR102342.3	A3
AXL	TATITVLPQQPR	Thermo	NR102342.4	A4
B2M	VNHVTLSQPK	Thermo	NR102342.5	A5
B2M	VEHSDLSFSK	Thermo	NR102342.6	A6
CCL20	QLANEGCDINAIIFHTK	Thermo	NR102342.7	A7
CCL20	LSVCANPK	Thermo	NR102342.8	A8
CD163	EAEFGQGTGPIWLNEVK	Thermo	NR102342.9	A9
```

##### Parameter File
The parameter file it a text file in the `.json` format. It contains instructions for building the survey run. It often makes sense to split a single peptide panel into multiple survey runs, for example, one per charge state. However, for a short list of peptides, we can prepare a single run. 
In this cas all the SIL pepetides in the mix are heavy, so we specify K[Hlys] and R[Harg] as fixed modifications. The `peptide_list_path` argument is an absolute path to the peptide list file (above). Example: 
```
╭─n.t.wamsley@3225-AD-00020.local ~/TEST_DATA/SUREQUANT_JUL24/HUMAN_IMMUNOLOGY  
╰─➤  cat params/build_survey_params.json
{
        "fixed_mods":[
                     ["C","C[Carb]"],
                     ["K$","K[Hlys]"],
                     ["R$","R[Harg]"]
        ],
        "modification_masses":
        {
        "Carb":57.021464,
        "Harg":10.008269,
        "Hlys":8.014199
        },
    "charges": [2, 3, 4],
"peptide_list_path":"/Users/n.t.wamsley/TEST_DATA/SUREQUANT_JUL24/HUMAN_IMMUNOLOGY/HUMAN_IMMUNOLOGY_PEPTIDE_LIST_SEP_2023.txt"
}
```
##### Run BuildSurvey
Using the list of SIL peptides, we will start Julia, import the `AuroPRM` library, and run the "BuildSurvey" method providing a path to the parameter file. This will build a folder `SurveyMethod` containing `SurveyMethod.csv` that we can import into the Survey method file in the Thermo method editor. 

```
╭─n.t.wamsley@3225-AD-00020.local ~/TEST_DATA/SUREQUANT_JUL24/HUMAN_IMMUNOLOGY  
╰─➤  julia
julia> using AutoPRM
[ Info: Precompiling AutoPRM [161ade51-a71e-4c5f-b262-f1262fd3217d]
[ Info: Skipping precompilation since __precompile__(false). Importing AutoPRM [161ade51-a71e-4c5f-b262-f1262fd3217d].
┌ Warning: Replacing docs for `AutoPRM.PSM :: Union{}` in module `AutoPRM`
└ @ Base.Docs docs/Docs.jl:243
julia> BuildSurvey("params/build_survey_params.json")
"/Users/n.t.wamsley/TEST_DATA/SUREQUANT_JUL24/HUMAN_IMMUNOLOGY/SurveyMethod/SurveyMethod.csv"
julia> exit()
╭─n.t.wamsley@3225-AD-00020.local ~/TEST_DATA/SUREQUANT_JUL24/HUMAN_IMMUNOLOGY  
╰─➤  head SurveyMethod/SurveyMethod.csv 
Compound,Formula,Adduct,m/z,z,intensity_threshold
ARG1_GGVEEGPTVLR[Harg],,,562.3026884000001,2,10000.0
ARG1_GGVEEGPTVLR[Harg],,,375.20421773333334,3,10000.0
ARG1_GGVEEGPTVLR[Harg],,,281.6549824,4,10000.0
ARG1_VMEETLSYLLGR[Harg],,,710.8726283999999,2,10000.0
ARG1_VMEETLSYLLGR[Harg],,,474.25084439999995,3,10000.0
ARG1_VMEETLSYLLGR[Harg],,,355.9399523999999,4,10000.0
AXL_APLQGTLLGYR[Harg],,,599.8445284000001,2,10000.0
AXL_APLQGTLLGYR[Harg],,,400.2321110666667,3,10000.0
AXL_APLQGTLLGYR[Harg],,,300.4259024,4,10000.0
```


### SearchSurvey
##### Convert Survey Run to Arrow
Move all survey runs into the 'survey_runs' folder. From the `PioneerConverter` directory, convert the `.raw` survey runs to `.arrow` format
```
╭─n.t.wamsley@3225-AD-00020.local ~/TEST_DATA/SUREQUANT_JUL24/NRF2_HP/PioneerConverter  ‹master*› 
╰─➤  dotnet run ../../HUMAN_IMMUNOLOGY/survey_runs/IMMUNOSIL_80nMSIL_1ugPEP_METHTEST_08022024_01.raw 
Converting: IMMUNOSIL_80nMSIL_1ugPEP_METHTEST_08022024_01
batchSize: 10000
n_threads: 2
Starting Conversion For: IMMUNOSIL_80nMSIL_1ugPEP_METHTEST_08022024_01
Execution Time: 6844 ms for IMMUNOSIL_80nMSIL_1ugPEP_METHTEST_08022024_01
```
##### Parameter File
The parameter file cat params/search_survey_params.json contains the instructions for how to search the survey run. Make sure to change the ms_data_dir to the directory containing the .arrow formatted survey run. Also set peptide_list_path to HUMAN_IMMUNOLOGY_PEPTIDE_LIST_SEP_2023.txt, the same precursor list used to build the survey run
```
╭─n.t.wamsley@3225-AD-00020.local ~/TEST_DATA/SUREQUANT_JUL24/HUMAN_IMMUNOLOGY  
╰─➤  cat params/search_survey_params.json
{
    "right_precursor_tolerance": 0.001,
    "left_precursor_tolerance": 0.001,
    "precursor_rt_tolerance": 0.3,
    "b_ladder_start": 3,
    "y_ladder_start": 4,
    "precursor_charges": [2, 3, 4],
    "precursor_isotopes": [0],
    "transition_charges": [1, 2],
    "transition_isotopes": [0],
    "fragment_match_ppm": 40,
    "minimum_fragment_count": 5,
    "fragments_to_select": 5,
    "precursor_rt_window": 0.3,
    "max_variable_mods": 2,
    "fixed_mods":[
                     ["C","C[Carb]"],
                     ["K$","K[Hlys]"],
                     ["R$","R[Harg]"]
    ],
    "variable_mods":
    [],
    "modification_masses":
        {
        "Carb":57.021464,
        "Harg":10.008269,
        "Hlys":8.014199
        },
    "ms_file_conditions":
        {
            "_35NCE_":"35NCE",
            "_40NCE_":"40NCE",
            "GAPDH":"GAPDH"
        },
    "ms_data_dir": "/Users/n.t.wamsley/TEST_DATA/SUREQUANT_JUL24/HUMAN_IMMUNOLOGY/survey_runs/arrow_out",
    "peptide_list_path": "/Users/n.t.wamsley/TEST_DATA/SUREQUANT_JUL24/HUMAN_IMMUNOLOGY/HUMAN_IMMUNOLOGY_PEPTIDE_LIST_SEP_2023.txt"
}
```
##### Run SearchSurvey
Searching the survey generates output in the same folder as the `.arrow` raw files. Outputs are `iapi_method.csv`, `transition_list.csv`, `precursors_summary.csv`, and a folder `figures`. The `iapi_method.csv` is an input for the Thermo IAPI SureQuant method. It specifies
the precursor m/z ratios, expected fragment ions, etc. To search an IS-PRM run, `transition_list.csv` is required as an input. `precursors_summary.csv` reports the best psm for each precursor accross all the survey runs searched. This can be useful for
combining the results of many surveys at different NCE or FAIMS CV values. The `figures` folder contains chromatograms and annotated spectra for the SIL peptides from the survey run. 
```
julia> SearchSurvey("params/search_survey_params.json")

shell> head survey_runs/arrow_out/iapi_method.csv
protein_name,sequence,precursor_mz,precursor_charge,retention_time,precursor_intensity,hyperscore,NthIntensity,sumTopN,file_name,condition,transition_mz
ARG1,GGVEEGPTVLR[Harg],562.3026884000001,2,39.35995,7.8305648e7,51.076622009277344,9.583797e6,3.9467732e7,IMMUNOSIL_80nMSIL_1ugPEP_METHTEST_08022024_01.arrow,NONE,652.3972;781.4386;910.4792;595.376;214.11655
ARG1,VMEETLSYLLGR[Harg],710.8726283999999,2,61.726166,1.5517284e7,54.347801208496094,1.9420418e6,8.438351e6,IMMUNOSIL_80nMSIL_1ugPEP_METHTEST_08022024_01.arrow,NONE,718.4056;932.5331;1061.574;831.4879;631.37427
AXL,APLQGTLLGYR[Harg],599.8445284000001,2,51.279736,4.1923796e7,54.990875244140625,4.6274795e6,1.9375808e7,IMMUNOSIL_80nMSIL_1ugPEP_METHTEST_08022024_01.arrow,NONE,789.456;917.5154;1030.5983;518.30206;732.4348
AXL,TATITVLPQQPR[Harg],667.8869284000001,2,45.642937,2.6542056e7,56.771541595458984,5.0537445e6,1.9140718e7,IMMUNOSIL_80nMSIL_1ugPEP_METHTEST_08022024_01.arrow,NONE,635.3557;948.5567;748.43976;847.5089;274.14142
...

shell> head survey_runs/arrow_out/transition_list.csv
protein_name,sequence,precursor_charge,precursor_isotope,transition_names
ARG1,GGVEEGPTVLR[Harg],2,0,y6+1;y7+1;y8+1;y5+1;b3+1
ARG1,VMEETLSYLLGR[Harg],2,0,y6+1;y8+1;y9+1;y7+1;y5+1
AXL,APLQGTLLGYR[Harg],2,0,y7+1;y8+1;y9+1;y4+1;y6+1
AXL,TATITVLPQQPR[Harg],2,0,y5+1;y8+1;y6+1;y7+1;b3+1
...

shell> head survey_runs/arrow_out/precursors_summary.csv
protein_name,sequence,precursor_mz,precursor_charge,retention_time,precursor_intensity,hyperscore,NthIntensity,sumTopN,file_name,condition,transition_mz
FOXP3,HNLSLHK[Hlys],428.7475883999999,2,24.704,1.0979155e6,25.879701614379883,40209.586,186015.69,IMMUNOSIL_80nMSIL_1ugPEP_METHTEST_08022024_01.arrow,NONE,719.4211;565.30206;702.3848;365.18796;605.3863
CD8A,AAEGLDTQR[Harg],485.74501339999995,2,27.910803,2.7493786e7,44.79069137573242,3.184588e6,1.2822804e7,IMMUNOSIL_80nMSIL_1ugPEP_METHTEST_08022024_01.arrow,NONE,699.37195;828.4157;529.2649;642.3484;272.124
...

```


|Name                |Default| Short        |Description                    |
 |--------------------|-------|-------------|--------------------|
 |params_json||mandatory|Path to a .json file with the parameters (see Configuration)
 |data_dir||mandatory|"Path to a folder with .arrow MS data tables"
 |precursor_list||mandatory|"Path to a tab delimited table of precursors"
 |--make_plots|true|-p|"Whether to make plots. Defaults to `true`"
 |--print_params|false|-s|"Whether to print the parameters from the json. Defaults to `false`"
 
#### Precursor List Example
```
ABCB6	YYNAESYEVER
ABCB6	IDGQDISQVTQASLR
ABCB6	ALNVLVPIFYR
ABHD4	YVSLPNQNK
ADD2	VNVADEVQR
.
.
.
```
#### Output
Writes transition_list.csv written into the user supplied `data_dir` (folder that contains the survey run). See below for example. 
#### Transition List Example
```
protein_name,sequence,precursor_charge,precursor_isotope,transition_names
ABCB6,YYNAESYEVER[Harg],2,0,y10+1;y7+1;y8+1;y9+1;y6+1
ABCB6,ALNVLVPIFYR[Harg],2,0,y6+1;y8+1;y7+1;b4+1;y9+1
ABHD4,YVSLPNQNK[Hlys],2,0,b4+1;y7+1;y5+2;y3+1;y5+1
.
.
.
```
### IS-PRM Analysis
###### POSIX
```
julia --threads 24 ./src/Routines/PRM/IS-PRM/routine.jl ./data/example_config/IS-PRM-TEST.json ./data/parquet ./data/parquet/transition_list.csv
```
###### Windows
```
julia --threads 24 .\\src\\Routines\\PRM\\IS-PRM\\routine.jl .\\data\\example_config\\IS-PRM-TEST.json .\\data\\parquet .\\data\\parquet\\transition_list.csv
```
|Name                |Default| Short        |Description                    |
 |--------------------|-------|-------------|--------------------|
 |params_json||mandatory|Path to a .json file with the parameters (see Configuration)
 |data_dir||mandatory|"Path to a folder with .arrow MS data tables"
 |transition_list||mandatory|"Path to a tab delimited table of transitions"
 |--make_plots|true|-p|"Whether to make plots. Defaults to `true`"
 |--print_params|false|-s|"Whether to print the parameters from the json. Defaults to `false`"

## Example Outputs

## Configuration files

### IS-PRM-Survey
```
{
    "right_precursor_tolerance": 0.001,
    "left_precursor_tolerance": 0.001,
    "precursor_rt_tolerance": 0.3,
    "b_ladder_start": 3,
    "y_ladder_start": 3,
    "precursor_charges": [2, 3, 4],
    "precursor_isotopes": [0],
    "transition_charges": [1, 2],
    "transition_isotopes": [0],
    "fragment_match_ppm": 40,
    "minimum_fragment_count": 5,
    "fragments_to_select": 5,
    "precursor_rt_window": 0.3,
    "max_variable_mods": 2,
    "fixed_mods":[
                     ["C","C[Carb]"],
                     ["K$","K[Hlys]"],
                     ["R$","R[Harg]"]
    ],
    "variable_mods":
    [],
    "modification_masses":
        {
        "Carb":57.021464,
        "Harg":10.008269,
        "Hlys":8.014199
        },
    "ms_file_conditions":
        {
            "_35NCE_":"35NCE",
            "_40NCE_":"40NCE",
            "GAPDH":"GAPDH"
        }
}
```

### IS-PRM
```
{
    "right_precursor_tolerance": 0.001,
    "left_precursor_tolerance": 0.001,
    "precursor_rt_tolerance": 0.3,
    "b_ladder_start": 3,
    "y_ladder_start": 3,
    "precursor_charges": [2, 3, 4],
    "precursor_isotopes": [0],
    "transition_charges": [1, 2],
    "transition_isotopes": [0],
    "fragment_match_ppm": 40,
    "minimum_fragment_count": 5,
    "fragments_to_select": 5,
    "precursor_rt_window": 0.3,
    "max_variable_mods": 2,
    "fixed_mods":[
                     ["C","C[Carb]"],
                     ["K$","K[Hlys]"],
                     ["R$","R[Harg]"]
    ],
    "variable_mods":
    [],
    "modification_masses":
        {
        "Carb":57.021464,
        "Harg":10.008269,
        "Hlys":8.014199
        },
    "ms_file_conditions":
        {
            "_35NCE_":"35NCE",
            "_40NCE_":"40NCE",
            "GAPDH":"GAPDH"
        }
}
```
## References
<a id="1">[1]</a> 
Gallien S, Kim SY, Domon B. Large-Scale Targeted Proteomics Using Internal Standard Triggered-Parallel Reaction Monitoring (IS-PRM). Mol Cell Proteomics. 2015 Jun;14(6):1630-44. doi: 10.1074/mcp.O114.043968. Epub 2015 Mar 9. PMID: 25755295; PMCID: PMC4458725.
<br>
<a id="1">[2]</a> 
Cox J, Hein MY, Luber CA, Paron I, Nagaraj N, Mann M. Accurate proteome-wide label-free quantification by delayed normalization and maximal peptide ratio extraction, termed MaxLFQ. Mol Cell Proteomics. 2014 Sep;13(9):2513-26. doi: 10.1074/mcp.M113.031591. Epub 2014 Jun 17. PMID: 24942700; PMCID: PMC4159666
<br>
<a id="1">[3]</a> 
Stopfer LE, Flower CT, Gajadhar AS, Patel B, Gallien S, Lopez-Ferrer D, White FM. High-Density, Targeted Monitoring of Tyrosine Phosphorylation Reveals Activated Signaling Networks in Human Tumors. Cancer Res. 2021 May 1;81(9):2495-2509. doi: 10.1158/0008-5472.CAN-20-3804. Epub 2021 Jan 28. PMID: 33509940; PMCID: PMC8137532.
<br>
<a id="1">[4]</a> 
Wamsley et al. Targeted proteomic quantitation of NRF2 signaling and predictive biomarkers in HNSCC. https://doi.org/10.1101/2023.03.13.532474 
