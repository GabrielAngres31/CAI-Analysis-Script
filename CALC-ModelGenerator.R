CALC.ModelGenerator <- function(df_in, fill_models_df) {
  fits_frame <- data.frame(accession = SWITCHBOARD.ACCESSIONLIST)
  for(model in fill_models_df$models) {
    new_models <- df_in %>%
      group_by(accession) %>%
      do(model = lm(model, data = .)) %>%
      setNames(c("accession", model)) %>%
      column_to_rownames("accession")
    fits_frame <- cbind(fits_frame, new_models[model])
  }
  return(fits_frame)
}