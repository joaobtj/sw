#' Normal Soil Water Balance
#'
#' Calculates the Normal (climatological) soil water balance for a 12 month period.
#'
#' @param month A string that identifies the column of the month
#' @param ET0 A string that identifies the column of the reference evapotranspiration
#' @param P A string that identifies the column of the precipitation
#' @param data A data.frame with (at least) three columns: month, reference evapotranspiration (et0) and precipitation (p)
#' @param awc Available water capacity of the soil
#' @param start_month Month of beginning of the soil water balance
#'
#' @export
#'
#' @examples
#' \dontrun{
#' swb_franca_sp <- tibble::tribble(
#'   ~m, ~et, ~p,
#'   "Jan", 117, 275,
#'   "Feb", 102, 218,
#'   "Mar", 104, 180,
#'   "Apr", 79, 60,
#'   "May", 60, 25,
#'   "Jun", 49, 20,
#'   "Jul", 54, 15,
#'   "Aug", 74, 12,
#'   "Sep", 93, 48,
#'   "Oct", 107, 113,
#'   "Nov", 108, 180,
#'   "Dec", 117, 245
#' )
#' normal_swb(month = "m", ET0 = "et", P = "p", data = swb_franca_sp, awc = 125, start_month = "Apr")
#' }
#'
normal_swb <- function(month, ET0, P, data, awc, start_month) {

  ## Reorder data.frame from star_month
  x <- which(data[[month]] == start_month)
  order_month <- data[[month]][c(
    x:12,
    seq_len(x - 1)
  )]
  data <- data[match(order_month, data[[month]]), ]

  ## Initialize data.frame
  data <- data.frame(data,
    pmet0 = rep(0, 12),
    L = rep(0, 12),
    A = rep(awc, 12),
    DA = rep(0, 12),
    ETa = rep(0, 12),
    D = rep(0, 12),
    E = rep(0, 12)
  )

  ## Calculate pmet0
  data$pmet0 <- data[[P]] - data[[ET0]]

  repeat{
    for (i in 1:12) {
      ## Calculate L and A
      if (data$pmet0[i] < 0) {
        ifelse(i == 1,
          data$L[i] <- data$L[12] - data$pmet0[i],
          data$L[i] <- data$L[i - 1] - data$pmet0[i]
        )
        data$A[i] <- awc * exp(-data$L[i] / awc)
        if (data$A[i] > awc) {
          data$A[i] <- awc
        }
      }
      if (data$pmet0[i] >= 0) {
        ifelse(i == 1,
          data$A[i] <- data$A[12] + data$pmet0[i],
          data$A[i] <- data$A[i - 1] + data$pmet0[i]
        )
        if (data$A[i] > awc) {
          data$A[i] <- awc
        }
        data$L[i] <- -awc * log(data$A[i] / awc)
      }
      ## Calculate DA
      ifelse(i == 1,
        data$DA[i] <- data$A[i] - data$A[12],
        data$DA[i] <- data$A[i] - data$A[i - 1]
      )
      ## Calculate ETa
      ifelse(data$pmet0[i] < 0,
        data$ETa[i] <- data[[P]][i] + abs(data$DA[i]),
        data$ETa[i] <- data[[ET0]][i]
      )
      ## Calculate D and E
      ifelse(data$pmet0[i] < 0,
        data$D[i] <- data[[ET0]][i] - data$ETa[i],
        data$E[i] <- data[[P]][i] - data[[ET0]][i] - data$DA[i]
      )
    }

    ## Check the conditions
    ## if satisfied, stop execution
    v1 <- sum(data[[ET0]]) - sum(data$ETa) - sum(data$D)
    v2 <- sum(data[[P]]) - sum(data$ETa) - sum(data$E)

    if (all(abs(v1) < 1e-10, abs(v2) < 1e-10)) break
  }

  return(data)
}
