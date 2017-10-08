using MendelGWAS, MendelBase, SnpArrays

using Compat
import Compat: view
using DataFrames                        # From package DataFrames.
using Distributions                     # From package Distributions.
using GLM                               # From package GLM.
using Plots                             # From package Plots.
using StatsBase                         # From package StatsBase.

#
# Used to create different dictionaries so that testing different
# functions don't accidentally use a pre-modified dictionary.
function make_keyword_dic()
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    keyword["distribution"] = "" # Binomial(), Gamma(), Normal(), Poisson(), etc.
    keyword["link"] = ""         # LogitLink(), IdentityLink(), LogLink(), etc.
    keyword["lrt_threshold"] = 5e-8
    keyword["maf_threshold"] = 0.01
    keyword["manhattan_plot_file"] = ""
    keyword["regression"] = ""   # Linear, Logistic, or Poisson
    keyword["regression_formula"] = ""
    keyword["pcs"] = 0
    return keyword
end

@testset "add_pcs!" begin
    # correctness of pca() has already been tested in SnpArrays package
    # will only test basic behavior of pca in MendelGWAS
    keyword = make_keyword_dic()
    process_keywords!(keyword, "gwas 1 Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
       locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
       read_external_data_files(keyword)
    @test typeof(pedigree_frame) == DataFrame
    @test size(pedigree_frame) == (2200, 7)
    @test pedigree_frame[1, 6] == 0.201823285
    keyword["pcs"] = 0 #pcs should be ≥ 1
    @test_throws(ArgumentError, 
        MendelGWAS.add_pcs!(pedigree_frame, snpdata, keyword))
    keyword["pcs"] = 0.5
    @test_throws(MethodError, 
        MendelGWAS.add_pcs!(pedigree_frame, snpdata, keyword))
    keyword["pcs"] = Inf
    @test_throws(MethodError, 
        MendelGWAS.add_pcs!(pedigree_frame, snpdata, keyword))
    keyword["pcs"] = NaN
    @test_throws(MethodError, 
        MendelGWAS.add_pcs!(pedigree_frame, snpdata, keyword))
    keyword["pcs"] = 1
    MendelGWAS.add_pcs!(pedigree_frame, snpdata, keyword)
    @test size(pedigree_frame) == (2200, 7 + 1)
    keyword["pcs"] = 6
    MendelGWAS.add_pcs!(pedigree_frame, snpdata, keyword)
    @test size(pedigree_frame) == (2200, 7 + 6)
    keyword["pcs"] = 100
    MendelGWAS.add_pcs!(pedigree_frame, snpdata, keyword)
    @test size(pedigree_frame) == (2200, 7 + 100)
end

@testset "assign_method" begin
    #gwas 1 example
    keyword = make_keyword_dic()
    process_keywords!(keyword, "gwas 1 Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
       locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
       read_external_data_files(keyword)
    (fast_method, lhs, rhs) = assign_method!(keyword, pedigree_frame)
    @test fast_method == true
    @test string(rhs) == "Sex"
    @test string(lhs) == "Trait"
    @test keyword["distribution"] == Normal()
    @test keyword["link"] == IdentityLink()
    @test keyword["regression"] == "linear"

    #gwas 2 example
    keyword = make_keyword_dic()
    process_keywords!(keyword, "gwas 2 Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
       locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
       read_external_data_files(keyword)
    (fast_method, lhs, rhs) = assign_method!(keyword, pedigree_frame)
    @test fast_method == true
    @test string(rhs) == "Sex + BMI"
    @test string(lhs) == "Case_Control"
    @test keyword["distribution"] == Binomial()
    @test keyword["link"] == LogitLink()
    @test keyword["regression"] == "logistic"

    # non-existant example 3
    keyword = make_keyword_dic()
    pedigree_frame = readtable("test1.txt")
    keyword["regression"] = "Poisson"
    keyword["regression_formula"] = "BMI ~ Sex & Case_Control"
    (fast_method, lhs, rhs) = assign_method!(keyword, pedigree_frame)
    @test fast_method == false # because rhs contains interaction term
    @test string(rhs) == "Sex & Case_Control"
    @test string(lhs) == "BMI"
    @test keyword["distribution"] == Poisson()
    @test keyword["link"] == LogLink()
    @test keyword["regression"] == "poisson"

    # are errors are being thrown properly?
    keyword = make_keyword_dic()
    pedigree_frame = readtable("test1.txt")
    # check if regression type present
    @test_throws(ArgumentError, assign_method!(keyword, pedigree_frame))
    keyword["regression"] = "Linear"
    # check if "~" is present in regression formula 
    @test_throws(ArgumentError, assign_method!(keyword, pedigree_frame))
    keyword["regression_formula"] = "hello + look + at + me"
    @test_throws(ArgumentError, assign_method!(keyword, pedigree_frame))
    # Left regression formula lacking
    keyword["regression_formula"] = " ~ look + at + me"
    @test_throws(ArgumentError, assign_method!(keyword, pedigree_frame))
    # Left regression formula contains more than 1 
    keyword["regression_formula"] = "look + at ~ me"
    @test_throws(ArgumentError, assign_method!(keyword, pedigree_frame))
    # Left regression formula contains bad characters
    keyword["regression_formula"] = "look at+me&you^* ~ hi + there"
    @test_throws(ArgumentError, assign_method!(keyword, pedigree_frame))
    # Left regression formula not recognized
    keyword["regression_formula"] = "look ~ hi + there"
    @test_throws(ArgumentError, assign_method!(keyword, pedigree_frame))
end

@testset "change_sex_desig!" begin
    keyword = make_keyword_dic()
    process_keywords!(keyword, "gwas 1 Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
       locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
       read_external_data_files(keyword)

    @test size(pedigree_frame) == (2200, 7)
    @test typeof(pedigree_frame[:Sex]) == 
        DataArrays.DataArray{Int64,1}
    @test pedigree_frame[1, 5] == 2
    @test pedigree_frame[2, 5] == 1
    @test pedigree_frame[3, 5] == 2
    @test pedigree_frame[2200, 5] == 2

    sex_is_1_or_2 = true
    # assign males = true, females = false
    sex_vector_prechange = trues(size(pedigree_frame, 1)) 
    for i in 1:size(pedigree_frame, 1)
        if pedigree_frame[i, 5] == 2; sex_vector_prechange[i] = false; end
        if pedigree_frame[i, 5] != 1 && pedigree_frame[i, 5] != 2
            sex_is_1_or_2 = false
        end
    end
    @test sex_is_1_or_2 == true

    pedigree_frame = MendelGWAS.change_sex_desig!(person, 
        keyword, pedigree_frame)
    
    @test size(pedigree_frame) == (2200, 7)
    @test typeof(pedigree_frame[:Sex]) == 
        DataArrays.DataArray{Float64,1}
    @test pedigree_frame[1, 7] == 1.0
    @test pedigree_frame[2, 7] == -1.0
    @test pedigree_frame[3, 7] == 1.0
    @test pedigree_frame[2200, 7] == 1.0

    sex_is_neg1_or_1 = true 
    sex_vector_postchange = trues(size(pedigree_frame, 1))
    for i in 1:size(pedigree_frame, 1)
        # 1 is female and -1 is male, after changing
        if pedigree_frame[i, 7] == 1; sex_vector_postchange[i] = false; end
        if pedigree_frame[i, 7] != -1 && pedigree_frame[i, 7] != 1
            sex_is_neg1_or_1 = false
        end
    end
    @test sex_is_neg1_or_1 == true
    @test sex_vector_prechange == sex_vector_postchange
end

@testset "change_case_desig!" begin
    keyword = make_keyword_dic()
    process_keywords!(keyword, "gwas 2 Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
       locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
       read_external_data_files(keyword)
    lhs, rhs = parse("Case_Control"), parse("Sex + BMI") #formula for gwas 2
    case_label = keyword["affected_designator"]

    @test size(pedigree_frame) == (2200, 5)
    @test typeof(pedigree_frame[:Case_Control]) == 
        DataArrays.DataArray{String,1}
    @test isna(pedigree_frame[1, 3])
    @test isna(pedigree_frame[2, 3])
    @test pedigree_frame[3, 3] == "NaN"
    @test isna(pedigree_frame[4, 3])
    @test isna(pedigree_frame[5, 3])
    @test pedigree_frame[6, 3] == "  " #there will be an extra space
    @test pedigree_frame[7, 3] == "1"
    @test pedigree_frame[8, 3] == "0"
    @test pedigree_frame[2200, 3] == "1"

    pedigree_frame = MendelGWAS.change_case_desig!(lhs, case_label, 
        person, pedigree_frame)
    
    @test size(pedigree_frame) == (2200, 5)
    @test typeof(pedigree_frame[:Case_Control]) == 
        DataArrays.DataArray{Float64,1}
    @test isna(pedigree_frame[1, 5])
    @test isna(pedigree_frame[2, 5])
    @test isna(pedigree_frame[3, 5])
    @test isna(pedigree_frame[4, 5])
    @test isna(pedigree_frame[5, 5])
    @test isna(pedigree_frame[6, 5])
    @test pedigree_frame[7, 5] == 1.0
    @test pedigree_frame[2197, 5] == 1.0
    @test pedigree_frame[2198, 5] == 0.0
    @test pedigree_frame[2199, 5] == 0.0
    @test pedigree_frame[2200, 5] == 1.0
end

@testset "Missing data assignment in model construction" begin
    keyword = make_keyword_dic()
    process_keywords!(keyword, "gwas 2 Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
       locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
       read_external_data_files(keyword)
    case_label = keyword["affected_designator"]
    lhs, rhs = parse("Case_Control"), parse("Sex + BMI")
    fm = Formula(lhs, rhs) 
    model = ModelFrame(fm, pedigree_frame)

    @test typeof(model) <: ModelFrame
    @test size(model.df) == (2196, 3)
    @test model.df[1, 1] == "NaN"
    @test model.df[2, 1] == "  " #there will be an extra space for some reason
    @test model.df[3, 1] == "1"

    pedigree_frame = MendelGWAS.change_case_desig!(lhs, case_label, 
        person, pedigree_frame)
    model = ModelFrame(fm, pedigree_frame)

    @test typeof(model) <: ModelFrame
    @test size(model.df) == (2194, 3)
    @test model.df[1, 1] == 1.0
    @test model.df[2, 1] == 0.0
    @test model.df[3, 1] == 1.0
    @test model.df[2194, 1] == 1.0
end

@testset "basics" begin
    
end

@testset "correctness" begin
    
end

