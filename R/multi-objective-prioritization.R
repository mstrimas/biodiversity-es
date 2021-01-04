#' Multiobjective prioritization
#'
#' Generate a prioritisation using hierarchical multi-objective optimization.
#'
#' @param rij \code{matrix} or \code{\link[Matrix]{dgcMatrix}} object
#'   containing the amount of each feature in each planning unit. Columns
#'   correspond to planning units and rows correspond to features.
#'
#' @param obj \code{matrix} or \code{numeric} vector object containing
#'   the objective values for the prioritization. If the argument to \code{obj}
#'   is a vector, then these values correspond to the cost of each planning
#'   unit. If the argument to \code{obj} is a matrix, then these values
#'   correspond to the cost of each planning unit under different objectives.
#'   Specifically, each row corresponds to a different cost metric
#'   and each column corresponds to a different planning unit. Note that
#'   planning units with higher values are always considered more costly.
#'   The argument to \code{obj} should have cost data in the same order
#'   as the argument to \code{rij}. The argument to \code{obj} should have
#'   the objectives (rows) ordered according to their importance, such that
#'   the first row has the most important objective and the last row has
#'   the least important objective. No missing values (\code{NA}) are
#'   permitted.
#'
#' @param pu_locked_in \code{logical} vector object indicating which
#'   planning units should be locked into the solution (with \code{TRUE})
#'   values and which planning units are not locked into the solution
#'   (\code{FALSE} values). The argument to \code{pu_locked_in} should have
#'   one element for each planning unit, and in the same order as the
#'   argument to \code{rij}. No missing values (\code{NA}) are
#'   permitted.
#'
#' @param relative_target \code{numeric} vector object containing the
#'   the percentage representation targets for each feature. The
#'   argument to \code{relative_target} should have one element for each
#'   feature, and in the same order as the argument to \code{rij}.
#'   No missing values (\code{NA}) are
#'
#' @param gap \code{numeric} number of vector object containing the
#'   optimality gap(s).
#'   
#' @param flip_priority \code{logical} logical indicating whether the 
#'   the priority order of the objectives should be reversed.
#'
#' @details
#' This function requires the following R packages: assertthat, guorobi, and
#' Matrix.
#'
#' @return \code{list} object. This list has the following elements:
#'   \describe{
#'   \item{solution}{\code{numeric} vector indicating which planning units are
#'     selected in the prioritization with zeros and ones. These planning
#'     units are in the same order as the argument to \code{rij}.}
#'   \item{objval}{\code{numeric} vector containing the objective value(s)
#'     for each objective in the problem.}
#'   }
#'
#'
#' @examples
#' # create data for prioritisation,
#' # this prioritisation involves 2 features, 4 planning units, and 2 objectives
#' ## objective value (cost) data
#' cost <- rbind(matrix(c(2, 1, 1, 100), nrow = 1),
#'               matrix(c(1, 2.1, 2, 0), nrow = 1))
#' print(cost)
#'
#' ## specify which planning units should be (TRUE) locked in or (FALSE) not
#' locked_in <- c(FALSE, TRUE, FALSE, FALSE)
#' print(locked_in)
#'
#' ## create feature (rij) data,
#' ## these values indicate how much of each feature is contaned within
#' ## each planning unit (e.g. number of individuals, hectares of habitats, etc)
#' features <- rbind(matrix(c(2, 1, 1, 100), nrow = 1),
#'                   matrix(c(10, 10, 10, 100), nrow = 1))
#' print(features)
#'
#' ## set targets for features,
#' ## basically, this is saying we want 3 units of feature 1 and 12 units
#' ## of feature 2. We are expressing this as a percentage because the
#' ## input format for the function requires targets as percentages
#' targets <- c(3, 12) / rowSums(features)
#' print(targets)
#'
#' # find solution,
#' # the solution needs to be within 10% of optimality for the first objective,
#' # and then given the remaining pool of possible solutions, minimize
#' # the second objective (i.e. 0% of optimality)
#' r <- multiobjective_prioritization(rij = features, obj = cost,
#'                                    pu_locked_in = locked_in,
#'                                    relative_target = targets,
#'                                    gap = c(0.1, 0))
#'
#' # print solution,
#' # the solution should be: [1, 1, 0, 0]
#' print(r$solution)
#'
#' @export
multiobjective_prioritization <- function(rij, obj, pu_locked_in,
                                          relative_target, gap, flip_priority = FALSE,
                                          threads = 1) {
  # assert arguments are valid
  assertthat::assert_that(
    inherits(rij, "dgCMatrix") || is.matrix(rij), nrow(rij) > 0, ncol(rij) > 0,
    is.numeric(relative_target), assertthat::noNA(relative_target),
    length(relative_target) == nrow(rij),
    min(relative_target) >= 0, max(relative_target) <= 1,
    is.numeric(obj) || is.matrix(obj), assertthat::noNA(c(obj)),
    is.numeric(gap), assertthat::noNA(gap), all(gap >= 0),
    is.logical(pu_locked_in), length(pu_locked_in) == ncol(rij),
    assertthat::noNA(pu_locked_in))
  if (inherits(rij, "dgCMatrix")) {
    assertthat::assert_that(assertthat::noNA(rij@x))
  } else {
    assertthat::assert_that(assertthat::noNA(c(rij)))
  }
  if (is.matrix(obj))
    assertthat::assert_that(length(gap) == nrow(obj), ncol(obj) == ncol(rij))
  if (!is.matrix(obj))
    assertthat::assert_that(length(gap) == 1, length(obj) == ncol(rij))
  # preliminary processing
  ## convert relative targets to absolute targets
  absolute_target <- rowSums(rij) * relative_target
  ## set gurobi parameters
  params <- list(Presolve = 2,
                 Threads = threads,
                 MIPGap = gap[1])
  ## build gurobi model object
  model <- list()
  model$modelsense <- "min"
  model$lb <- as.numeric(pu_locked_in)
  model$ub <- rep(1, ncol(rij))
  model$vtype <- rep("B", ncol(rij))
  model$sense <- rep(">=", nrow(rij))
  model$rhs <- absolute_target
  model$A <- rij
  if (length(gap) == 1) {
    ## if single objective problem then just add objective
    model$obj <- obj
  } else {
    ## if multi-objective problem then we need to add data for each objective
    model$multiobj <- lapply(seq_along(gap), function(i) {
      list(objn = obj[i, ], priority = length(gap) + 1 + ifelse(flip_priority, i, i * -1), reltol = gap[i])
    })
  }
  # main processing
  o <- gurobi::gurobi(model = model, params = params)
  # exports
  list(solution = o$x[seq_len(ncol(rij))], objval = o$objval)
}

# code tests
library(testthat)
test_that("single objective (none locked)", {
  # data
  cost <- matrix(c(1, 2, 2), nrow = 1)
  locked_in <- rep(FALSE, 3)
  features <- rbind(matrix(c(2, 1, 1), nrow = 1),
                    matrix(c(10, 10, 10), nrow = 1))
  targets <- c(2, 10) / rowSums(features)
  # find solution
  r <- multiobjective_prioritization(rij = features, obj = cost,
                                     pu_locked_in = locked_in,
                                     relative_target = targets,
                                     gap = 0)
  # tests
  expect_equal(r$solution, c(1, 0, 0))
  expect_equal(r$objval, 1)
})

test_that("single objective (locked in)", {
  # data
  cost <- matrix(c(1, 2, 2), nrow = 1)
  locked_in <- c(FALSE, TRUE, FALSE)
  features <- rbind(matrix(c(2, 1, 1), nrow = 1),
                    matrix(c(10, 10, 10), nrow = 1))
  targets <- c(2, 10) / rowSums(features)
  # find solution
  r <- multiobjective_prioritization(rij = features, obj = cost,
                                     pu_locked_in = locked_in,
                                     relative_target = targets,
                                     gap = 0)
  # tests
  expect_equal(r$solution, c(1, 1, 0))
  expect_equal(r$objval, 3)
})

test_that("multiple objectives (none locked)", {
  # data
  cost <- rbind(matrix(c(2, 1, 1, 100), nrow = 1),
                matrix(c(1, 2.1, 2, 0), nrow = 1))
  locked_in <- rep(FALSE, 4)
  features <- rbind(matrix(c(2, 1, 1, 100), nrow = 1),
                    matrix(c(10, 10, 10, 100), nrow = 1))
  targets <- c(3, 12) / rowSums(features)
  # find solution
  r <- multiobjective_prioritization(rij = features, obj = cost,
                                     pu_locked_in = locked_in,
                                     relative_target = targets,
                                     gap = c(0.1, 0))
  # tests
  expect_equal(r$solution, c(1, 0, 1, 0))
  expect_equal(r$objval, c(3, 3))
})

test_that("multiple objectives (locked in)", {
  # data
  cost <- rbind(matrix(c(2, 1, 1, 100), nrow = 1),
                matrix(c(1, 2.1, 2, 0), nrow = 1))
  locked_in <- c(FALSE, TRUE, FALSE, FALSE)
  features <- rbind(matrix(c(2, 1, 1, 100), nrow = 1),
                    matrix(c(10, 10, 10, 100), nrow = 1))
  targets <- c(3, 12) / rowSums(features)
  # find solution
  r <- multiobjective_prioritization(rij = features, obj = cost,
                                     pu_locked_in = locked_in,
                                     relative_target = targets,
                                     gap = c(0.1, 0))
  # tests
  expect_equal(r$solution, c(1, 1, 0, 0))
  expect_equal(r$objval, c(3, 3.1))
})

test_that("multiple objectives (gap[1] = 0, second objective ignored)", {
  # data
  cost <- rbind(matrix(c(2, 1, 1.1, 100), nrow = 1),
                matrix(c(1, 500, 0, 0), nrow = 1))
  locked_in <- c(FALSE, TRUE, FALSE, FALSE)
  features <- rbind(matrix(c(2, 1, 1, 100), nrow = 1),
                    matrix(c(10, 10, 10, 100), nrow = 1))
  targets <- c(3, 12) / rowSums(features)
  # find solution
  r <- multiobjective_prioritization(rij = features, obj = cost,
                                     pu_locked_in = locked_in,
                                     relative_target = targets,
                                     gap = c(0, 0))
  # tests
  expect_equal(r$solution, c(1, 1, 0, 0))
  expect_equal(r$objval, c(3, 501))
})
