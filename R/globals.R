# Global symbol declarations for dynamically resolved native entry points.
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "_OptiLCMS_PerformLipidDDADeco",
    "lipids_build_dot_buffer_cpp",
    "_OptiLCMS_lipids_build_dot_buffer_cpp",
    "lipidomics_db_search_precursor_cpp",
    "_OptiLCMS_lipidomics_db_search_precursor_cpp",
    "lipidomics_db_search_precursors_blob_cpp",
    "_OptiLCMS_lipidomics_db_search_precursors_blob_cpp"
  ))
}
