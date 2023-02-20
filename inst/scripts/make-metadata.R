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
    ## Mouse Context VAE weights
    data.frame(
        Title = "Weights_ContextVAE_Mouse",
        Description = "Weights for the full context (gene lcpms) VAE of mouse contrast experiments compiled from ARCHS4",
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
        RDataClass = "character",
        DispatchClass = "kerasHDF5ModelWeights", ## AnnotationHub::DispatchClassList()
        Location_Prefix = "https://zenodo.org/record/7554915/files/",
        RDataPath = "ContextVAE_ARCHS4_v212_Mouse.weights.hdf5"
    ),
    ## Ηuman Context VAE weights
    data.frame(
        Title = "Weights_ContextVAE_Human",
        Description = "Weights for the full context (gene lcpms) VAE of human contrast experiments compiled from ARCHS4",
        BiocVersion = "3.17",
        Genome = "GRCh38",
        SourceType = "HDF5",  ## check AnnotationHubData::getValidSourceTypes()
        SourceUrl = "https://maayanlab.cloud/archs4/",
        SourceVersion = "ARCHS4 v2.1.2",  ## ARCHS4 version
        Species = "Homo sapiens",
        TaxonomyId = "9609",
        Coordinate_1_based = NA,
        DataProvider = "Ma'ayan Laboratory (https://labs.icahn.mssm.edu/maayanlab/)",
        Maintainer = "Panagiotis Papasaikas <panagiotis.papasaikas@fmi.ch>",
        RDataClass = "character",
        DispatchClass = "kerasHDF5ModelWeights", ## AnnotationHub::DispatchClassList()
        Location_Prefix = "https://zenodo.org/record/7554915/files/",
        RDataPath = "ContextVAE_ARCHS4_v212_Human.weights.hdf5"
    ),
    ## Mouse Delta cVAE weights
    data.frame(
        Title = "Weights_DeltaCVAE_Mouse",
        Description = "Weights for the full contrast (gene lcpm deltas i.e LFCs) conditional VAE of mouse contrast experiments compiled from ARCHS4",
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
        RDataClass = "character",
        DispatchClass = "kerasHDF5ModelWeights", ## AnnotationHub::DispatchClassList()
        Location_Prefix = "https://zenodo.org/record/7554915/files/",
        RDataPath = "DeltaCVAE_FT_ARCHS4_v212_Mouse.weights.hdf5"
    ),
    ## Ηuman Delta cVAE weights
    data.frame(
        Title = "Weights_DeltaCVAE_Human",
        Description = "Weights for the full contrast (gene lcpm deltas i.e LFCs) conditional VAE of human contrast experiments compiled from ARCHS4",
        BiocVersion = "3.17",
        Genome = "GRCh38",
        SourceType = "HDF5",  ## check AnnotationHubData::getValidSourceTypes()
        SourceUrl = "https://maayanlab.cloud/archs4/",
        SourceVersion = "ARCHS4 v2.1.2",  ## ARCHS4 version
        Species = "Homo sapiens",
        TaxonomyId = "9609",
        Coordinate_1_based = NA,
        DataProvider = "Ma'ayan Laboratory (https://labs.icahn.mssm.edu/maayanlab/)",
        Maintainer = "Panagiotis Papasaikas <panagiotis.papasaikas@fmi.ch>",
        RDataClass = "character",
        DispatchClass = "kerasHDF5ModelWeights", ## AnnotationHub::DispatchClassList()
        Location_Prefix = "https://zenodo.org/record/7554915/files/",
        RDataPath = "DeltaCVAE_FT_ARCHS4_v212_Human.weights.hdf5"
    ),
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
        DispatchClass = "kerasHDF5Model", ## AnnotationHub::DispatchClassList()
        Location_Prefix = "https://zenodo.org/record/7554915/files/",
        RDataPath = "ContextEncoder_ARCHS4_v212_Mouse.hdf5"
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
        DispatchClass = "kerasHDF5Model", ## AnnotationHub::DispatchClassList()
        Location_Prefix = "https://zenodo.org/record/7554915/files/",
        RDataPath = "ContextEncoder_ARCHS4_v212_Human.hdf5"
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
        DispatchClass = "kerasHDF5Model", ## AnnotationHub::DispatchClassList()
        Location_Prefix = "https://zenodo.org/record/7554915/files/",
        RDataPath = "DeltaEncoder_FT_ARCHS4_v212_Mouse.hdf5"
    ),
    ## Human Contrast Encoder
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
        DispatchClass = "kerasHDF5Model", ## AnnotationHub::DispatchClassList()
        Location_Prefix = "https://zenodo.org/record/7554915/files/",
        RDataPath = "DeltaEncoder_FT_ARCHS4_v212_Human.hdf5"
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
        DispatchClass = "kerasHDF5Model", ## AnnotationHub::DispatchClassList()
        Location_Prefix = "https://zenodo.org/record/7554915/files/",
        RDataPath = "DeltaDecoder_FT_ARCHS4_v212_Mouse.hdf5"
    ),
    ## Human Contrast Decoder
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
        DispatchClass = "kerasHDF5Model", ## AnnotationHub::DispatchClassList()
        Location_Prefix = "https://zenodo.org/record/7554915/files/",
        RDataPath = "DeltaDecoder_FT_ARCHS4_v212_Human.hdf5"
    ),
    
    ## Full Mouse ContrastDB Summarized Experiment hdf5 representation
    data.frame(
        Title = "decomposed_contrasts_mouse_hdf5",
        Description = "HDF5 file with assays for the full and decomposed mouse contrast experiments compiled from ARCHS4. Ngenes x 100 column blocks.",
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
        RDataClass = "character",
        DispatchClass = "H5File", ## AnnotationHub::DispatchClassList()
        Location_Prefix = "https://zenodo.org/record/7554915/files/",
        RDataPath = "mouse_v212_NDF_c100assays.h5"
    ),
    data.frame(
        Title = "decomposed_contrasts_mouse_rds",
        Description = "Serialized version of decomposed_contrasts_mouse_hdf5 with row and column metadata for the assays",
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
        RDataClass = "character",
        DispatchClass = "Rds", ## AnnotationHub::DispatchClassList()
        Location_Prefix = "https://zenodo.org/record/7554915/files/",
        RDataPath = "mouse_v212_NDF_c100se.rds"
    ),
    
    ## Fill Human ContrastDB Summarized Experiment hdf5 representation
    data.frame(
        Title = "decomposed_contrasts_human_hdf5",
        Description = "HDF5 file with assays for the full and decomposed human contrast experiments compiled from ARCHS4. Ngenes x 100 column blocks.",
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
        RDataClass = "character",
        DispatchClass = "H5File", ## AnnotationHub::DispatchClassList()
        Location_Prefix = "https://zenodo.org/record/7554915/files/",
        RDataPath = "human_v212_NDF_c100assays.h5"
    ),
    data.frame(
        Title = "decomposed_contrasts_human_rds",
        Description = "Serialized version of tdecomposed_contrasts_human_hdf5 with row and column metadata for the assays",
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
        RDataClass = "character",
        DispatchClass = "Rds", ## AnnotationHub::DispatchClassList()
        Location_Prefix = "https://zenodo.org/record/7554915/files/",
        RDataPath = "human_v212_NDF_c100se.rds"
    ),
    
    ## "Toy" Mouse ContrastDB Summarized Experiment hdf5 representation
    data.frame(
        Title = "demo_decomposed_contrasts_mouse_hdf5",
        Description = "A heavily subsampled version of decomposed_contrasts_mouse_hdf5 used for demo purposes",
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
        RDataClass = "character",
        DispatchClass = "H5File", ## AnnotationHub::DispatchClassList()
        Location_Prefix = "https://zenodo.org/record/7554915/files/",
        RDataPath = "mouse_v212_NDF_c100_DEMOassays.h5"
    ),
    data.frame(
        Title = "demo_decomposed_contrasts_mouse_rds",
        Description = "Serialized version of demo_decomposed_contrasts_mouse_hdf5 with row and column metadata",
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
        RDataClass = "character",
        DispatchClass = "Rds", ## AnnotationHub::DispatchClassList()
        Location_Prefix = "https://zenodo.org/record/7554915/files/",
        RDataPath = "mouse_v212_NDF_c100_DEMOse.rds"
    ),
    
    ## "Toy" Human ContrastDB Summarized Experiment hdf5 representation
    data.frame(
        Title = "demo_decomposed_contrasts_human_hdf5",
        Description = "A heavily subsampled version of decomposed_contrasts_human_hdf5 used for demo purposes",
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
        RDataClass = "character",
        DispatchClass = "H5File", ## AnnotationHub::DispatchClassList()
        Location_Prefix = "https://zenodo.org/record/7554915/files/",
        RDataPath = "human_v212_NDF_c100_DEMOassays.h5"
    ),
    data.frame(
        Title = "demo_decomposed_contrasts_human_rds",
        Description = "Serialized version of demo_decomposed_contrasts_human_hdf5 with row and column metadata",
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
        RDataClass = "character",
        DispatchClass = "Rds", ## AnnotationHub::DispatchClassList()
        Location_Prefix = "https://zenodo.org/record/7554915/files/",
        RDataPath = "human_v212_NDF_c100_DEMOse.rds"
    )
)

write.csv(
    df,
    file = "../extdata/metadata.csv",
    row.names = FALSE
)

