# A script to make the metadata.csv file located in inst/extdata of the package.
# See ?AnnotationHubData::makeAnnotationHubMetadata for a description of the
# metadata.csv file, expected fields and data types. This
# AnnotationHubData::makeAnnotationHubMetadata() function can be used to
# validate the metadata.csv file before submitting the package.
# See https://bioconductor.org/packages/release/bioc/vignettes/HubPub/inst/doc/CreateAHubPackage.html
# for more details about the content of the metadata file.

suppressPackageStartupMessages({
    library(dplyr)
})

df <- dplyr::bind_rows(
    ## Mouse ARCHS4
    data.frame(
        Title = "",
        Description = "",
        BiocVersion = "3.17",
        Genome = "",
        SourceType = "",  ## check AnnotationHubData::getValidSourceTypes()
        SourceUrl = "",
        SourceVersion = "",  ## ARCHS4 version
        Species = "Mus musculus",
        TaxonomyId = "10090",
        Coordinate_1_based = NA,
        DataProvider = "",
        Maintainer = "Panagiotis Papasaikas <panagiotis.papasaikas@fmi.ch>",
        RDataClass = "SummarizedExperiment",
        DispatchClass = "", ## AnnotationHub::DispatchClassList()
        Location_Prefix = "",
        RDataPath = ""
    ),
    ## Human ARCHS4
    data.frame(
        Title = "",
        Description = "",
        BiocVersion = "3.17",
        Genome = "",
        SourceType = "",  ## check AnnotationHubData::getValidSourceTypes()
        SourceUrl = "",
        SourceVersion = "",  ## ARCHS4 version
        Species = "Homo sapiens",
        TaxonomyId = "9606",
        Coordinate_1_based = NA,
        DataProvider = "",
        Maintainer = "Panagiotis Papasaikas <panagiotis.papasaikas@fmi.ch>",
        RDataClass = "SummarizedExperiment",
        DispatchClass = "", ## AnnotationHub::DispatchClassList()
        Location_Prefix = "",
        RDataPath = ""
    )

)

write.csv(
    filelist,
    file = "../extdata/metadata.csv",
    row.names = FALSE
)
