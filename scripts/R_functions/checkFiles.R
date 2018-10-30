# check if all the files exists and if not, warn about the missing one and return just those that exists

checkFiles <- function(.files, .type){
    files_exists <- file.exists(.files)
    if( !all(files_exists) ){
        warning(paste("missing ", .type, "files:\n", paste('\t\t\t\t', .files[!files_exists], collapse = '\n')))
        .files <- .files[files_exists]
    }
    return(.files)
}