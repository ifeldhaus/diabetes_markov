setwd("/Users/isabellefeldhaus/Dropbox/Harvard/Dissertation/diabetes_modeling/diabetes_markov/")

OUTPUT_DIR = 'output'

filename_for_saving <- function(title, nlines, extension)
  sprintf('%s/%s_%dL_%s.%s', OUTPUT_DIR, title, nlines, format(Sys.time(), '%Y-%m-%d_%H%M'), extension)

save_df_to_csv <- function(df, title) {
  filename <- filename_for_saving(title, nlines = nrow(df), extension = 'csv')
  readr::write_csv(df, filename)
}

save_object_to_rdata <- function(object){
  title <- deparse(substitute(object))
  filename <- filename_for_saving(title, nlines = dim(object)[1], extension = 'RData')
  save(object, file = filename)
}

load_object_from_rdata <- function(filename){
  get(load(filename))
}
