using MendelGWAS, MendelBase, SnpArrays

using Compat
import Compat: view
using DataFrames                        # From package DataFrames.
using Distributions                     # From package Distributions.
using GLM                               # From package GLM.
using Plots                             # From package Plots.
using StatsBase                         # From package StatsBase.
using RDatasets

# Used to create different dictionaries so that testing different
# functions don't accidentally use an existing dictionary.
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
    (fast_method, lhs, rhs) = 
        MendelGWAS.assign_method!(keyword, pedigree_frame)
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
    (fast_method, lhs, rhs) = 
        MendelGWAS.assign_method!(keyword, pedigree_frame)
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
    (fast_method, lhs, rhs) = 
        MendelGWAS.assign_method!(keyword, pedigree_frame)
    @test fast_method == false # because rhs contains interaction term
    @test string(rhs) == "Sex & Case_Control"
    @test string(lhs) == "BMI"
    @test keyword["distribution"] == Poisson()
    @test keyword["link"] == LogLink()
    @test keyword["regression"] == "poisson"

    # The rest tests of errors are being thrown properly
    keyword = make_keyword_dic()
    pedigree_frame = readtable("test1.txt")
    # check if regression type present
    @test_throws(ArgumentError, 
        MendelGWAS.assign_method!(keyword, pedigree_frame))
    keyword["regression"] = "Linear"
    # check if "~" is present in regression formula 
    @test_throws(ArgumentError, 
        MendelGWAS.assign_method!(keyword, pedigree_frame))
    keyword["regression_formula"] = "hello + look + at + me"
    @test_throws(ArgumentError, 
        MendelGWAS.assign_method!(keyword, pedigree_frame))
    # Left regression formula lacking
    keyword["regression_formula"] = " ~ look + at + me"
    @test_throws(ArgumentError, 
        MendelGWAS.assign_method!(keyword, pedigree_frame))
    # Left regression formula contains more than 1 
    keyword["regression_formula"] = "look + at ~ me"
    @test_throws(ArgumentError, 
        MendelGWAS.assign_method!(keyword, pedigree_frame))
    # Left regression formula contains bad characters
    keyword["regression_formula"] = "look at+me&you^* ~ hi + there"
    @test_throws(ArgumentError, 
        MendelGWAS.assign_method!(keyword, pedigree_frame))
    # Left regression formula not recognized
    keyword["regression_formula"] = "look ~ hi + there"
    @test_throws(ArgumentError, 
        MendelGWAS.assign_method!(keyword, pedigree_frame))
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

@testset "run_fast_regression" begin
    # This function basically calls the regress function in MendelBase, 
    # which already contains its own tests. So we will just test some 
    # basic behavior.

    #
    # gwas example 1, as is
    #
    keyword = make_keyword_dic()
    process_keywords!(keyword, "gwas 1 Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
       locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
       read_external_data_files(keyword)
    # remember to change males to -1 and females to 1
    # otherwise the result will be wrong.
    pedigree_frame = MendelGWAS.change_sex_desig!(person, keyword, pedigree_frame)
    people = person.people
    snps = snpdata.snps
    dosage = zeros(people)
    copy!(dosage, view(snpdata.snpmatrix, :, 1); impute = true)

    lhs, rhs = parse("Trait"), parse("Sex")
    fm = Formula(lhs, rhs) 
    model = ModelFrame(fm, pedigree_frame)
    complete = completecases(model.df)
    regression_type = lowercase(keyword["regression"])
    io = keyword["output_unit"]

    (X, y, cases, complete, residual_base, base_estimate, 
        base_loglikelihood) = MendelGWAS.run_fast_regression(model, fm, 
        regression_type, lhs, io)

    @test size(X) == (2200, 2)
    @test sum(X[:, 1]) == 2200 #first column are 1's
    @test X[1, 2] == 1.0 #female
    @test X[2, 2] == -1.0 #male
    @test X[3, 2] == 1.0 
    @test X[2200, 2] == 1.0 
    @test size(y) == (2200,)
    @test y == model.df[:Trait]
    @test cases == 2200
    @test size(residual_base) == (2200,)
    @test residual_base == y - (X * base_estimate)
    @test signif(base_estimate[1], 6) == 0.142482
    #negative sign (-0.0129378) below if males = -1 and females = 1
    @test signif(base_estimate[2], 6) == -0.0129378 
    @test signif(base_loglikelihood, 6) == -1258.18

    #
    # gwas example 2, with additional missing data. 
    #
    keyword = make_keyword_dic()
    process_keywords!(keyword, "gwas 2 Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
       locus_frame, phenotype_frame, pedigree_frame, 
       snp_definition_frame) = read_external_data_files(keyword)
    lhs, rhs = parse("Case_Control"), parse("Sex + BMI")
    case_label = keyword["affected_designator"]
    # remember to change males to -1 and females to 1
    # and to exclude the first few missing data
    pedigree_frame = MendelGWAS.change_sex_desig!(person, keyword, 
        pedigree_frame)
    pedigree_frame = MendelGWAS.change_case_desig!(lhs, case_label,
        person, pedigree_frame)
    people = person.people
    snps = snpdata.snps
    dosage = zeros(people)
    copy!(dosage, view(snpdata.snpmatrix, :, 1); impute = true)

    fm = Formula(lhs, rhs) 
    model = ModelFrame(fm, pedigree_frame)
    complete = completecases(model.df)
    regression_type = lowercase(keyword["regression"])
    io = keyword["output_unit"]

    (X, y, cases, complete, residual_base, base_estimate, 
        base_loglikelihood) = MendelGWAS.run_fast_regression(model, fm, 
        regression_type, lhs, io)
    
    @test size(X) == (2194, 3)
    @test sum(X[:, 1]) ≈ 2194 #first column are 1's
    @test X[1, 2] ≈ -1.0 #male
    @test X[2, 2] ≈ 1.0 #female
    @test X[3, 2] ≈ 1.0 
    @test X[2194, 2] ≈ 1.0 
    @test size(y) == (2194,)
    @test y == model.df[:Case_Control]
    @test cases == 2194
    # logistic model returns no residual base
    @test typeof(residual_base) == String 
    # numbers below are slightly different because 
    # some data originally present were removed to test for missing data
    @test signif(base_estimate[1], 6) == -0.687241 
    @test signif(base_estimate[2], 6) == 0.012089
    @test signif(base_estimate[3], 6) == -0.00919787  
    @test signif(base_loglikelihood, 6) == -1312.1
end

@testset "use_glm_package" begin
    # This function basically calls glm() on some previously constructed
    # models and formula, so we will test by plugging in 2 basic examples
    # from the glm package.
    keyword = make_keyword_dic()
    process_keywords!(keyword, "gwas 2 Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
       locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
       read_external_data_files(keyword)

    keyword["link"] = IdentityLink()
    keyword["distribution"] = Normal()
    data = DataFrame(X=[1,2,3], Y=[2,4,7])
    fm = @formula(Y ~ X)
    io = keyword["output_unit"]
    model = ModelFrame(fm, data)
    test1 = MendelGWAS.use_glm_package(keyword, model, fm, io)

    @test signif(coef(test1)[1], 6) == -0.666667
    @test signif(coef(test1)[2], 6) == 2.5
    @test signif(stderr(test1)[1], 6) == 0.62361
    @test signif(stderr(test1)[2], 6) == 0.288675
    @test signif(deviance(test1), 6) == 0.166667
    @test dof_residual(test1) == 1.0

    keyword = make_keyword_dic()
    process_keywords!(keyword, "gwas 2 Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
       locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
       read_external_data_files(keyword)
    keyword["link"] = ProbitLink()
    keyword["distribution"] = Binomial()
    # keyword["output_unit"] = "test3.txt"
    data = DataFrame(X=[1,2,3], Y=[1,0,1])
    fm = @formula(Y ~ X)
    io = keyword["output_unit"]
    model = ModelFrame(fm, data)
    test2 = MendelGWAS.use_glm_package(keyword, model, fm, io)

    @test signif(coef(test2)[1], 6) == 0.430727
    @test signif(coef(test2)[2], 6) == -3.64399e-19
    @test signif(stderr(test2)[1], 6) == 1.98019
    @test signif(stderr(test2)[2], 6) == 0.91665
    @test signif(deviance(test2), 6) == 3.81909
    @test dof_residual(test2) == 1.0
end

@testset "run_score_test" begin
    
end

@testset "Does all the functions come together: gwas_option" begin
    keyword = make_keyword_dic()
    process_keywords!(keyword, "gwas 1 Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
       locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
       read_external_data_files(keyword)

    test1 = MendelGWAS.gwas_option(person, snpdata, pedigree_frame, keyword)
    @test test1 == false

    keyword = make_keyword_dic()
    process_keywords!(keyword, "gwas 1 Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
       locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
       read_external_data_files(keyword)
    test2 = MendelGWAS.gwas_option(person, snpdata, pedigree_frame, keyword)
    @test test2 == false
end

