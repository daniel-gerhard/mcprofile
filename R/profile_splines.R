## Interpolating splines of the profiles in srdp, one per contrast (row of CM).
## Where a profile cannot be interpolated (e.g. fewer than 4 profile points),
## the element is a "try-error" and a single warning names the contrasts.
profile_splines <- function(srdp, CM){
  spl <- lapply(srdp, function(x){
    x <- na.omit(x)
    try(interpSpline(x[, 1], x[, 2]), silent = TRUE)
  })
  failed <- vapply(spl, inherits, logical(1), what = "try-error")
  if (any(failed)) {
    msg <- vapply(spl[failed], function(e) conditionMessage(attr(e, "condition")), character(1))
    warning("Could not interpolate the profile; results are NA for:\n",
            paste0("  ", rownames(CM)[failed], ": ", msg, collapse = "\n"),
            call. = FALSE)
  }
  spl
}
