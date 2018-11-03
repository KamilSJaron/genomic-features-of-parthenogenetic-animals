make_data_frame <- function (variables) {
    asm_template <- matrix(ncol = length(variables))
    colnames(asm_template) <- variables
    asm_template
}