cache_dir <- function() {
  inst_path <- tools::R_user_dir("FLAMES", "cache")
  pkg.v <- as.character(packageVersion("FLAMES"))
  file.path(path.expand(inst_path), pkg.v)
}

download_oarfish <- function(folder) {
  # get OS and architecture to determine the correct binary
  os <- Sys.info()[["sysname"]]
  arch <- Sys.info()[["machine"]]

  # https://github.com/COMBINE-lab/oarfish/releases
  binary_urls <- list(
    "Linux" = list(
      "x86_64" = c(
        url = "https://github.com/COMBINE-lab/oarfish/releases/download/v0.8.1/oarfish-x86_64-unknown-linux-gnu.tar.xz",
        sha256 = "b5bc3a577a2dfc53aa4f20d4ea34ce9a2ec42ec9509f747131a28a116a660e59"
      ),
      "arm64" = c(
        url = "https://github.com/COMBINE-lab/oarfish/releases/download/v0.8.1/oarfish-aarch64-unknown-linux-gnu.tar.xz",
        sha256 = "220a67b02af5d91f6de34111295a8de0412e8ccf40e1bb8526c67b5e3f938967"
      )
    ),
    "Darwin" = list(
      "x86_64" = c(
        url = "https://github.com/COMBINE-lab/oarfish/releases/download/v0.8.1/oarfish-x86_64-apple-darwin.tar.xz",
        sha256 = "654712ffcddd96f54a28dd2f7dec1a435f21f0fc09f27a58f3470d4935b1390f"
      ),
      "arm64" = c(
        url = "https://github.com/COMBINE-lab/oarfish/releases/download/v0.8.1/oarfish-aarch64-apple-darwin.tar.xz",
        sha256 = "28875820038bb29f20311d755d0354f2f5de1c8835a01583619bd7d0cc1ea6c0"
      )
    )
  )

  # check if the OS and architecture are supported
  if (!os %in% names(binary_urls)) {
    stop("Unsupported OS: ", os)
  }
  if (!arch %in% names(binary_urls[[os]])) {
    stop("Unsupported architecture: ", arch)
  }

  # download the binary
  url <- binary_urls[[os]][[arch]][1]
  sha256 <- binary_urls[[os]][[arch]][2]
  archive <- tempfile(fileext = ".tar.xz")
  if (!dir.exists(folder)) {
    dir.create(folder, recursive = TRUE)
  }
  download.file(url, archive, mode = "wb")

  # check the SHA256 checksum
  if (cli::hash_file_sha256(archive) != sha256) {
    stop("SHA256 checksum does not match for ", archive)
  }

  # extract the binary
  unarchive <- tempfile()
  dir.create(unarchive)
  if (grepl("\\.tar\\.xz$", archive)) {
    system(paste("tar -xf", shQuote(archive), "-C", shQuote(unarchive)))
  } else if (grepl("\\.zip$", archive)) {
    unzip(archive, exdir = unarchive)
  } else {
    stop("Unsupported file format: ", archive)
  }

  # find binary and move it to the target folder
  binary <- list.files(unarchive, pattern = "oarfish$", full.names = TRUE, recursive = TRUE)
  if (length(binary) != 1) {
    stop("Could not find oarfish binary in ", unarchive)
  }
  target_binary <- file.path(folder, "oarfish")
  file.copy(binary, target_binary)
  Sys.chmod(target_binary, mode = "0755")

  message("oarfish binary downloaded to ", target_binary)
  return(target_binary)
}

#' Find path to a binary
#' Wrapper for Sys.which to find path to a binary
#' @importFrom withr with_path
#' @importFrom basilisk obtainEnvironmentPath
#' @description
#' This function is a wrapper for \code{base::Sys.which} to find the path
#' to a command. It also searches within the \code{FLAMES} basilisk conda
#' environment. This function also replaces "" with \code{NA} in the
#' output of \code{base::Sys.which} to make it easier to check if the
#' binary is found.
#' @param command character, the command to search for
#' @return character, the path to the command or \code{NA}
#' @examples
#' find_bin("minimap2")
#' @export
find_bin <- function(command) {
  flames_bins <- cache_dir()
  which_command <- withr::with_path(
    new = c(flames_bins, system.file(package = "FLAMES", "bin")),
    action = "suffix",
    code = Sys.which(command)
  )
  # replace "" with NA
  which_command[which_command == ""] <- NA

  if ("oarfish" %in% command && is.na(which_command["oarfish"])) {
    if (interactive()) {
      message("oarfish not found. Do you want to download it? (Y/n)")
      answer <- readline(prompt = "> ")
      if (answer == "" || tolower(answer) == "y") {
        oarfish <- download_oarfish(flames_bins)
      } else {
        stop("oarfish not found and not downloaded.")
      }
    } else {
      message("oarfish not found. Downloading it now.")
      oarfish <- download_oarfish(flames_bins)
    }
    which_command["oarfish"] <- oarfish
  }

  return(which_command)
}
