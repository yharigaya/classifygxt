#' Get the names of the eight models
#'
#' This function returns an ordered vector of character strings
#' corresponding to the names of the eight models.
#'
#' @export
#' @return A vector containing the eight model names.
#'
get_model_names <- function() {
    model.mat <- get_model_mat()
    get_names_from_matrix(model.mat)
}

#' Get the names of the 27 models
#'
#' This function returns an ordered vector of character strings
#' corresponding to the names of the 27 models accounting for the
#' sign of effect sizes.
#'
#' @export
#' @return A vector containing the eight model names.
#'
get_sign_names <- function() {
    sign.mat <- get_sign_mat()
    get_names_from_matrix(sign.mat)
}

#' Format input
#'
#' This function takes input and output files of
#' \href{https://github.com/broadinstitute/tensorqtl}{TensorQTL} and
#' returns a list object that can be used as input for
#' \code{\link{do_bms}}. It currenly only handles completely unpaired data,
#' where there is no repeated measurement per donor (subject). This
#' includes typical molecular data from clinical trials.
#'
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom stats residuals
#' @importFrom utils read.delim
#' @param qtl A character string specfying the name of the output file
#'     from interacton QTL mapping from TensorQTL.
#' @param pheno A character string specifynig the name of the file
#'     containing molecular phenotes. This must be the same as the
#'     input file used for TensorQTL.

#' @param geno A character string specifying the name of the genotype
#'     file in \href{https://www.cog-genomics.org/plink/2.0}{PLINK2}
#'     \code{.traw} format.
#' @param covar A character string specifyng the name of the file
#'     containing covariates. This must be the same as the input file
#'     used for TensorQTL.
#' @param int A character string specifying the name of the file
#'     containing the condition variable of interest, coded as {1,
#'     2}. This must be the same as the input file used for TensorQTL.
#' @param pheno.col A character string specifying the name of the
#'     column of \code{qtl} containing the feature IDs.
#' @param variant.col A character string specifying the name of the
#'     column of \code{qtl} containing the feature IDs.
#' @param feat.col A character string specifying the name of the
#'     column of \code{pheno} containing the feature IDs.
#' @param snp.col A character string specifying the name of the
#'     column of \code{geno} containing the SNP IDs.
#' @param pheno.rm A vector of integers corresponding to the column
#'     numbers for non-phenotype entries in \code{pheno}.
#' @param geno.rm A vector of integers corresponding to the column
#'     numbers for non-phenotype entries in \code{geno}.
#'
#' @return A list of lists containing:
#' \itemize{
#' \item{\code{y} - A vector of phenotypes.}
#' \item{\code{g} - A vector of genotypes.}
#' \item{\code{t} - A vector of treatment indicators.}
#' \item{\code{subject} - A vector of subject.}
#' \item{\code{feat.id} - A character string spcifyng the feature ID.}
#' \item{\code{snp.id} - A character string specifying the SNP ID.}
#' }
#'
format_input <- function(qtl, pheno, geno, covar, int,
                         pheno.col="phenotype_id",
                         variant.col="variant_id",
                         feat.col="gene_id",
                         snp.col="SNP",
                         pheno.rm=1:4,
                         geno.rm=1:6) {

    # read in and check the qtl file
    qtl <- qtl %>% read.delim

    qtl.check <- try(qtl[, pheno.col], silent=TRUE)
    if (class(qtl.check) != "character") {
        stop("check the `qtl` and `pheno.col` arguments")
    }

    # read in and check the covar file
    covar <- covar %>% read.delim(row.names=1)

    # convert the data frame to a matrix
    covar <- covar %>%
        as.matrix %>%
        t

    if (!all(is.numeric(covar))) {
        stop("check the `covar` argument")
    }

    # read in and check the interaction file
    int <- int %>% read.delim(header=FALSE)

    # get and format the condition variable
    t <- int[, 2]
    tmp <- t
    tmp[t == 1] <- 0
    tmp[t == 2] <- 1
    t <- tmp

    if (!all(t %in% c(0, 1))) {
        stop("check the `int` argument")
    }

    # read in the geno file
    geno <- geno %>% read.delim

    # read in the pheno file
    pheno <- pheno %>% read.delim

    # extract id columns (feat.vec, snp.vec)
    feat.vec <- try(pheno[, feat.col], silent=TRUE)
    snp.vec <- try(geno[, snp.col], silent=TRUE)

    if (class(feat.vec) != "character") {
        stop("check the `pheno` and `feat.col` arguments")
    }

    if (class(snp.vec) != "character") {
        stop("check the `geno` and `snp.col` arguments")
    }

    # remove unnecessary columns
    pheno <- pheno[, - pheno.rm]
    geno <- geno[, - geno.rm]

    if (!all(apply(pheno, 2, is.numeric))) {
        stop("check the `pheno` and `pheno.rm` argument")
    }

    if (!all(apply(geno, 2, is.numeric))) {
        stop("check the `geno` and `geno.rm` argument")
    }

    # get the subject names
    subject <- colnames(pheno)

    if (all(rownames(covar) != subject)) {
        stop("sample names in the covariate file do not match those in the phenotype file")
    }

    output <- lapply(
        seq_len(nrow(qtl)), format_each,
        qtl=qtl, pheno=pheno, geno=geno,
        pheno.col=pheno.col, variant.col=variant.col,
        covar=covar, t=t, subject=subject,
        feat.vec=feat.vec, snp.vec=snp.vec)

    output
}

format_each <- function(index, qtl, pheno, geno,
                        pheno.col, variant.col,
                        covar, t, subject,
                        feat.vec, snp.vec) {

    feat.id <- qtl[index, pheno.col]
    snp.id <- qtl[index, variant.col]

    # get the phenotype
    pheno.sel <- feat.vec == feat.id

    if (sum(pheno.sel) != 1) {
        stop("there is no match or more than one match in features")
    }

    y <- pheno[pheno.sel, ] %>%
        as.numeric %>%
        unname

    # get the genotype
    geno.sel <- snp.vec == snp.id

    if (sum(geno.sel) != 1) {
        stop("there is no match or more than one match in snps")
    }

    g <- geno[geno.sel, ] %>%
        as.numeric %>%
        unname

    # regress out the covariate
    d <- cbind(y, covar) %>%
        as.data.frame
    frml <- "y ~ ."
    y <- lm(frml, data=d) %>%
        residuals %>%
        unname

    output <- list(
        y=y, g=g, t=t, subject=subject,
        feat.id=feat.id, snp.id=snp.id)

    output
}
