# ---------- helpers ----------
.save_managed_upload <- function(upload, role, root) {

  if (is.null(upload) || nrow(upload) == 0) {
    return(NULL)
  }

  managed_name <- paste0(
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    "_",
    role,
    "_",
    basename(upload$name)
  )

  managed_path <- file.path(root, managed_name)

  ok <- file.copy(upload$datapath, managed_path, overwrite = TRUE)

  if (!ok) {
    stop("Failed to copy uploaded file.")
  }

  list(
    role = role,
    original_name = upload$name,
    managed_name = managed_name,
    managed_path = managed_path,
    size = upload$size
  )
}

.remove_managed_upload <- function(meta) {
  if (!is.null(meta) &&
      !is.null(meta$managed_path) &&
      file.exists(meta$managed_path)) {
    unlink(meta$managed_path, force = TRUE)
  }
  invisible(NULL)
}

.remove_dir_safely <- function(path) {
  if (!is.null(path) && dir.exists(path)) {
    unlink(path, recursive = TRUE, force = TRUE)
  }
  invisible(NULL)
}
