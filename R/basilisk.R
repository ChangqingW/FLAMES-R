#' @importFrom basilisk BasiliskEnvironment
flames_env <- BasiliskEnvironment(
    envname = "flames_env", pkgname = "FLAMES",
    pip = c("fast-edit-distance==1.2.2", "blaze2==2.5.1", "matplotlib==3.9.1"),
    packages = c(
        "python==3.11.9",
        "numpy==2.1.1",
        "scipy==1.14.1",
        "pysam==0.22.1",
        "cutadapt==4.9",
        "tqdm==4.66.5",
        "pandas==2.2.3"
    ),
    channels = c("conda-forge", "bioconda")
)
