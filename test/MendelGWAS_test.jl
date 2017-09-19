using MendelGWAS, MendelBase, SnpArrays

using Compat
import Compat: view
using DataFrames                        # From package DataFrames.
using Distributions                     # From package Distributions.
using GLM                               # From package GLM.
using Plots                             # From package Plots.
using StatsBase                         # From package StatsBase.

info("Unit tests for MendelGWAS")



@testset "GWAS constructors" begin
    
end


@testset "add_pcs! test" begin
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    keyword["distribution"] = "" # Binomial(), Gamma(), Normal(), Poisson(), etc.
    keyword["link"] = ""         # LogitLink(), IdentityLink(), LogLink(), etc.
    keyword["lrt_threshold"] = 5e-8
    keyword["maf_threshold"] = 0.01
    keyword["manhattan_plot_file"] = ""
    keyword["regression"] = ""   # Linear, Logistic, or Poisson
    keyword["regression_formula"] = ""
    keyword["pcs"] = 6
    process_keywords!(keyword, "gwas 1 Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
           locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
           read_external_data_files(keyword)

    @test size(pedigree_frame) == (2200, 7)

    MendelGWAS.add_pcs!(pedigree_frame, snpdata, keyword)

    @test size(pedigree_frame) == (2200, 7 + keyword["pcs"])
    "hello"
end
