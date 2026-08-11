# Describe the provenance of this FLAMES build for logs / bug reports.
# Always returns a single string to display, most specific first: the exact
# commit for container and GitHub installs, otherwise the version plus the
# source (Bioconductor, or "unknown source" for a plain source install).
flames_build_revision <- function(pkgname = "FLAMES",
                                  desc = utils::packageDescription(pkgname)) {
  # Container -- see Dockerfile
  sha <- Sys.getenv("FLAMES_IMAGE_REVISION")
  if (nzchar(sha)) {
    return(paste0("container revision: ", sha))
  }
  # Installs from GitHub 
  # e.g. BiocManager::install("mritchielab/FLAMES") or remotes::install_github 
  remote_sha <- desc$RemoteSha
  if (is.null(remote_sha)) {
    remote_sha <- desc$GithubSHA1
  }
  if (!is.null(remote_sha) && !is.na(remote_sha) && nzchar(remote_sha)) {
    return(paste0("revision: ", remote_sha))
  }
  # Bioconductor
  repo <- desc$Repository
  if (!is.null(repo) && !is.na(repo) && grepl("^Bioconductor", repo)) {
    return(paste0("version ", desc$Version, " (", repo, ")"))
  }
  # Unknown (local build)
  ver <- desc$Version
  if (is.null(ver) || is.na(ver)) {
    ver <- "?"
  }
  paste0("version ", ver, " (unknown source)")
}

.onAttach <- function(libname, pkgname) {
  # when library()/require()
  rev <- flames_build_revision(pkgname)
  if (!is.null(rev)) {
    packageStartupMessage("FLAMES ", rev)
  }
}

# Log the FLAMES build revision
flames_log_revision <- function() {
  rev <- flames_build_revision()
  if (!is.null(rev)) {
    message("FLAMES ", rev)
  }
  invisible(NULL)
}
