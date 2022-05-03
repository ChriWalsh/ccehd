#' Example data set to illustrate the estimator
#'
#' A simulated data set based on model used in scenario B of paper with
#' \eqn{N = 50} cross-sectional units, and \eqn{T = 50} time periods,
#' \eqn{K = 3} factors, and \eqn{p = 900} regressors.
#'
#' In contrast to the paper the additional irrelevant regressors are not
#' grouped by type. Instead they are resorted so that their types alternate.
#' Specifically, they alternate by type. That is for each \eqn{k = 0, ..., 298}:
#'       \itemize{
#'       \item{The \eqn{(4 + k * 3)}-th regressors only depend on the
#'       first factor.}
#'       \item{The \eqn{(5 + k * 3)}-th regressors only depend on the
#'       second factor.}
#'       \item{The \eqn{(6 + k * 3)}-th regressors only depend on the
#'       third factor.}
#'       }
#'
#'
#'
#' @format A list with :
#' \describe{
#' \item{data_example$y}{vector containing dependent variables for the
#'       \eqn{ N * T = 2500} observations.
#'       Sorted as \eqn{y_{11}, ... y_{1T},y_{21},...y_{2T},...,y_{NT}}.}
#' \item{data_example$x}{matrix containg \eqn{p = 900} regressors.
#'       Each column corresponds to one regressor sorted as the dependent variables.
#'       The first three regressors correspond to the relevant ones as in the
#'       paper.The remaining \eqn{3*299} regressors have been rearranged in
#'       comparison to the design in the paper.
#'       }
#' \item{data_example$f}{\eqn{T \times K} matrix containing the values for the
#'       \eqn{K=3} unobserved that were used to construct the dependent
#'       variables and the regressors.}
#' }
"data_example"
