using MendelGWAS, MendelBase, SnpArrays

using Compat
import Compat: view
using DataFrames                        # From package DataFrames.
using Distributions                     # From package Distributions.
using GLM                               # From package GLM.
using Plots                             # From package Plots.
using StatsBase                         # From package StatsBase.

info("Unit tests for MendelGWAS")

@testset "add_pcs! test" begin
    # correctness of pca() has already been tested in snparrays_test.jl
    # will only test basic behavior of pca in MendelGWAS

    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    keyword["distribution"] = "" # Binomial(), Gamma(), Normal(), Poisson(), etc.
    keyword["link"] = ""         # LogitLink(), IdentityLink(), LogLink(), etc.
    keyword["lrt_threshold"] = 5e-8
    keyword["maf_threshold"] = 0.01
    keyword["manhattan_plot_file"] = ""
    keyword["regression"] = ""   # Linear, Logistic, or Poisson
    keyword["regression_formula"] = ""
    keyword["pcs"] = 0
    process_keywords!(keyword, "gwas 1 Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
           locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
           read_external_data_files(keyword)

    @test size(pedigree_frame) == (2200, 7)
    keyword["pcs"] = 0 #pcs should be ≥ 1
    @test_throws(ArgumentError, MendelGWAS.add_pcs!(pedigree_frame, snpdata, keyword))
    keyword["pcs"] = 0.5
    @test_throws(MethodError, MendelGWAS.add_pcs!(pedigree_frame, snpdata, keyword))
    keyword["pcs"] = Inf
    @test_throws(MethodError, MendelGWAS.add_pcs!(pedigree_frame, snpdata, keyword))
    keyword["pcs"] = NaN
    @test_throws(MethodError, MendelGWAS.add_pcs!(pedigree_frame, snpdata, keyword))
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

@testset "gwas basics" begin
    #side assignment
    
end

@testset "gwas_option correctness" begin
    
end

@testset "gwas error messages" begin
    
end

