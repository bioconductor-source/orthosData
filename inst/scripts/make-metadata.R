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
    ## Mouse Context Encoder
    data.frame(
        Title = "ContextEncoder_Mouse",
        Description = "Context (gene lcpms) encoder of mouse contrast experiments compiled from ARCHS4",
        BiocVersion = "3.17",
        Genome = "GRCm39",
        SourceType = "HDF5",  ## check AnnotationHubData::getValidSourceTypes()
        SourceUrl = "https://maayanlab.cloud/archs4/",
        SourceVersion = "ARCHS4 v2.1.2",  ## ARCHS4 version
        Species = "Mus musculus",
        TaxonomyId = "10090",
        Coordinate_1_based = NA,
        DataProvider = "Ma'ayan Laboratory (https://labs.icahn.mssm.edu/maayanlab/)",
        Maintainer = "Panagiotis Papasaikas <panagiotis.papasaikas@fmi.ch>",
        RDataClass = "keras.engine.functional.Functional",
        DispatchClass = "keras::load_model_hdf5; require keras", ## AnnotationHub::DispatchClassList()
        Location_Prefix = "",
        RDataPath = "TRAINED_MODELS/ContextEncoder_ARCHS4_v212_Mouse"
    ),
    ## Human Context Encoder
    data.frame(
        Title = "ContextEncoder_Human",
        Description = "Context (gene lcpms) encoder of human contrast experiments compiled from ARCHS4",
        BiocVersion = "3.17",
        Genome = "GRCh38",
        SourceType = "HDF5",  ## check AnnotationHubData::getValidSourceTypes()
        SourceUrl = "https://maayanlab.cloud/archs4/",
        SourceVersion = "ARCHS4 v2.1.2",  ## ARCHS4 version
        Species = "Homo sapiens",
        TaxonomyId = "9606",
        Coordinate_1_based = NA,
        DataProvider = "Ma'ayan Laboratory (https://labs.icahn.mssm.edu/maayanlab/)",
        Maintainer = "Panagiotis Papasaikas <panagiotis.papasaikas@fmi.ch>",
        RDataClass = "keras.engine.functional.Functional",
        DispatchClass = "keras::load_model_hdf5; require keras", ## AnnotationHub::DispatchClassList()
        Location_Prefix = "",
        RDataPath = "TRAINED_MODELS/ContextEncoder_ARCHS4_v212_Mouse"
    ),
    ## Mouse Contrast Encoder
    data.frame(
        Title = "DeltaEncoder_Mouse",
        Description = "Contrast (gene lcpm deltas i.e LFCs) encoder of mouse contrast experiments compiled from ARCHS4",
        BiocVersion = "3.17",
        Genome = "GRCm39",
        SourceType = "HDF5",  ## check AnnotationHubData::getValidSourceTypes()
        SourceUrl = "https://maayanlab.cloud/archs4/",
        SourceVersion = "ARCHS4 v2.1.2",  ## ARCHS4 version
        Species = "Mus musculus",
        TaxonomyId = "10090",
        Coordinate_1_based = NA,
        DataProvider = "Ma'ayan Laboratory (https://labs.icahn.mssm.edu/maayanlab/)",
        Maintainer = "Panagiotis Papasaikas <panagiotis.papasaikas@fmi.ch>",
        RDataClass = "keras.engine.functional.Functional",
        DispatchClass = "keras::load_model_hdf5; require keras", ## AnnotationHub::DispatchClassList()
        Location_Prefix = "",
        RDataPath = "TRAINED_MODELS/DeltaEncoder_FT_ARCHS4_v212_Mouse"
    ),
    ## Human Context Encoder
    data.frame(
        Title = "DeltaEncoder_Human",
        Description = "Contrast (gene lcpm deltas i.e LFCs) encoder of human contrast experiments compiled from ARCHS4",
        BiocVersion = "3.17",
        Genome = "GRCh38",
        SourceType = "HDF5",  ## check AnnotationHubData::getValidSourceTypes()
        SourceUrl = "https://maayanlab.cloud/archs4/",
        SourceVersion = "ARCHS4 v2.1.2",  ## ARCHS4 version
        Species = "Homo sapiens",
        TaxonomyId = "9606",
        Coordinate_1_based = NA,
        DataProvider = "Ma'ayan Laboratory (https://labs.icahn.mssm.edu/maayanlab/)",
        Maintainer = "Panagiotis Papasaikas <panagiotis.papasaikas@fmi.ch>",
        RDataClass = "keras.engine.functional.Functional",
        DispatchClass = "keras::load_model_hdf5; require keras", ## AnnotationHub::DispatchClassList()
        Location_Prefix = "",
        RDataPath = "TRAINED_MODELS/DeltaEncoder_FT_ARCHS4_v212_Human"
    ),
    ## Mouse Contrast Decoder
    data.frame(
        Title = "DeltaDecoder_Mouse",
        Description = "Contrast (gene lcpm deltas i.e LFCs) decoder of mouse contrast experiments compiled from ARCHS4",
        BiocVersion = "3.17",
        Genome = "GRCm39",
        SourceType = "HDF5",  ## check AnnotationHubData::getValidSourceTypes()
        SourceUrl = "https://maayanlab.cloud/archs4/",
        SourceVersion = "ARCHS4 v2.1.2",  ## ARCHS4 version
        Species = "Mus musculus",
        TaxonomyId = "10090",
        Coordinate_1_based = NA,
        DataProvider = "Ma'ayan Laboratory (https://labs.icahn.mssm.edu/maayanlab/)",
        Maintainer = "Panagiotis Papasaikas <panagiotis.papasaikas@fmi.ch>",
        RDataClass = "keras.engine.functional.Functional",
        DispatchClass = "keras::load_model_hdf5; require keras", ## AnnotationHub::DispatchClassList()
        Location_Prefix = "",
        RDataPath = "TRAINED_MODELS/DeltaDecoder_FT_ARCHS4_v212_Mouse"
    ),
    ## Human Context Encoder
    data.frame(
        Title = "DeltaDecoder_Human",
        Description = "Contrast (gene lcpm deltas i.e LFCs) decoder of human contrast experiments compiled from ARCHS4",
        BiocVersion = "3.17",
        Genome = "GRCh38",
        SourceType = "HDF5",  ## check AnnotationHubData::getValidSourceTypes()
        SourceUrl = "https://maayanlab.cloud/archs4/",
        SourceVersion = "ARCHS4 v2.1.2",  ## ARCHS4 version
        Species = "Homo sapiens",
        TaxonomyId = "9606",
        Coordinate_1_based = NA,
        DataProvider = "Ma'ayan Laboratory (https://labs.icahn.mssm.edu/maayanlab/)",
        Maintainer = "Panagiotis Papasaikas <panagiotis.papasaikas@fmi.ch>",
        RDataClass = "keras.engine.functional.Functional",
        DispatchClass = "keras::load_model_hdf5; require keras", ## AnnotationHub::DispatchClassList()
        Location_Prefix = "",
        RDataPath = "TRAINED_MODELS/DeltaDecoder_FT_ARCHS4_v212_Human"
    ),
    
    ## Mouse ContrastDB Summarized Experiment
    data.frame(
        Title = "decomposed_contrasts_mouse",
        Description = "HDF5-based SummarizedExperiment with full and decomposed mouse contrast experiments compiled from ARCHS4 and associated metadata",
        BiocVersion = "3.17",
        Genome = "GRCm39",
        SourceType = "tar",  ## check AnnotationHubData::getValidSourceTypes()
        SourceUrl = "https://maayanlab.cloud/archs4/",
        SourceVersion = "ARCHS4 v2.1.2",  ## ARCHS4 version
        Species = "Mus musculus",
        TaxonomyId = "10090",
        Coordinate_1_based = NA,
        DataProvider = "Ma'ayan Laboratory (https://labs.icahn.mssm.edu/maayanlab/)",
        Maintainer = "Panagiotis Papasaikas <panagiotis.papasaikas@fmi.ch>",
        RDataClass = "keras.engine.functional.Functional",
        DispatchClass = "HDF5Array::loadHDF5SummarizedExperiment; require HDF5Array", ## AnnotationHub::DispatchClassList()
        Location_Prefix = "",
        RDataPath = "DECOMPOSED_CONTRASTS/mouse_v212_c100se.tar"
    ),
    ## Human ContrastDB Summarized Experiment
    data.frame(
        Title = "decomposed_contrasts_human",
        Description = "HDF5-based SummarizedExperiment with full and decomposed human contrast experiments compiled from ARCHS4 and associated metadata",
        BiocVersion = "3.17",
        Genome = "GRCh38",
        SourceType = "tar",  ## check AnnotationHubData::getValidSourceTypes()
        SourceUrl = "https://maayanlab.cloud/archs4/",
        SourceVersion = "ARCHS4 v2.1.2",  ## ARCHS4 version
        Species = "Homo sapiens",
        TaxonomyId = "9606",
        Coordinate_1_based = NA,
        DataProvider = "Ma'ayan Laboratory (https://labs.icahn.mssm.edu/maayanlab/)",
        Maintainer = "Panagiotis Papasaikas <panagiotis.papasaikas@fmi.ch>",
        RDataClass = "keras.engine.functional.Functional",
        DispatchClass = "HDF5Array::loadHDF5SummarizedExperiment; require HDF5Array", ## AnnotationHub::DispatchClassList()
        Location_Prefix = "",
        RDataPath = "DECOMPOSED_CONTRASTS/human_v212_c100se.tar"
    )
)

write.csv(
    df,
    file = "../extdata/metadata.csv",
    row.names = FALSE
)
