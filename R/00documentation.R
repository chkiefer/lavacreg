#' A first example to illustrate 
#'
#' A dataset containing ... .
#'
#' @format A data frame with 871 rows and 9 variables:
#' \describe{
#'   \item{dv}{Count of correctly-answered items (dependent variable)}
#'   \item{treat}{Treatment group variable, where \code{0} is control and \code{2} is treatment}
#'   \item{z11}{First indicator of internal LoC}
#'   \item{z12}{Second indicator of internal LoC}
#'   \item{z21}{First indicator of external LoC}
#'   \item{z22}{Second indicator of external LoC}
#'   \item{z41}{First indicator of depression}
#'   \item{z42}{Second indicator of depression}
#'   \item{z43}{Third indicator of depression}
#' }
"example01"


#' lavacreg
#' 
#' Latent Variable Count Regression Models
#' 
#' @docType package
#' @author Christoph Kiefer <christoph.kiefer@uni-bielefeld.de>
#' @import Rcpp fastGHQuad pracma
#' @importFrom Rcpp evalCpp
#' @useDynLib lavacreg
#' @name lavacreg
NULL  