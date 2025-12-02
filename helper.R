compute_expected2 <- function(counts,
                             exclude = NULL,
                             include.trend = TRUE,
                             trend.knots.per.year = 1/7,
                             extrapolate = TRUE,
                             harmonics = 2,
                             frequency = NULL,
                             weekday.effect = FALSE,
                             keep.components = TRUE,
                             verbose = TRUE){
  
  ## check column names
  if (!all(c("date", "outcome", "population") %in% names(counts))) stop("counts must have columns named date, outcome, and poulation.")
  
  if (any(table(counts$date)) > 1) stop("Each date can appear at most once.")
  
  if (any(c("expected", "excluded") %in% names(counts))) warning("expected and excluded columns will be overwritten.")
  
  if (any(is.na(counts$date)) | any(is.na(counts$outcome)) | any(is.na(counts$population)))
    stop("No NAs permited in date, outcome, or population columns.")
  
  if (!lubridate::is.Date(counts$date)) stop("date column must be class Date.")
  
  if (!is.numeric(counts$outcome)) stop("outome column must be counts.")
  
  if (!is.numeric(counts$population)) stop("population column must be numeric.")
  
  if (is.null(exclude)) {
    warning("No dates excluded. We recommend excluding at least the dates surrounding the event of interest.")
  }
  
  if (!identical(counts$date, dplyr::arrange(counts, date)$date)) stop("counts must be ordered by date.")
  
  if (include.trend & as.numeric(diff(range(counts$date)) / 365) < 5) warning("Including a trend in the model is not recommended with less than five years of data. Consider setting include.trend = FALSE.")
  
  ## helper function
  fourier_trend <- function(x, k = 3){
    H <- lapply(1:k, function(k){
      cbind(sin(2*pi*k/365*x), cos(2*pi*k/365*x))
    })
    res <- do.call(cbind, H)
    colnames(res) <- paste(rep(c("sin", "cos"), k), rep(1:k, each = 2), sep = "_")
    res
  }
  
  ## number of observations per year
  if (is.null(frequency)) {
    frequency <- round(365 / (as.numeric(diff(range(counts$date)))/nrow(counts)))
    if (verbose) message("No frequency provided, determined to be ", frequency, " measurements per year.")
  }
  
  if (frequency < 365 & weekday.effect) {
    warning("Modeling day effects is not recommended when frequency < 365. Consider setting weekday.effect = FALSE")
  }
  
  if (weekday.effect == TRUE) {
    if (frequency == 12) {
      message("Monthly data detected,  ingnoring weekday.effect argument.")
      weekday.effect <- FALSE
    }
    
    if (floor(frequency) == 52) {
      message("Weekly data detected,  ingnoring weekday.effect argument.")
      weekday.effect <- FALSE
    }
  }
  
  if (verbose) message("Overall death rate is ", signif(sum(counts$outcome, na.rm = TRUE)/sum(counts$population, na.rm = TRUE)*frequency*1000, 3), ".")
  
  if (mean(counts$outcome, na.rm = TRUE) < 1)
    warning("Average counts per unit of time is below 1. Confidence intervals are based on asymptotic approxmations that might not hold.")
  
  ## build design matrix
  # convert dates to time
  
  tt <- as.numeric(counts$date)
  index <- !(counts$date %in% exclude)
  last_index <- max(which(index))
  first_index <- min(which(index))
  
  if (!extrapolate) {
    ## index to later keep constant because we are not extrapolating
    if (last_index < length(index)) end_extrapolate_index <- seq(last_index + 1, length(index)) else end_extrapolate_index <- NULL
    if (first_index > 1) start_extrapolate_index <- seq(1, first_index - 1) else start_extrapolate_index <- NULL
  }
  
  # make trend basis (includes intercept)
  if (include.trend) {
    
    # compute knots
    # we use linear extrapolation outside the range of count$date
    years <- (max(tt[index]) - min(tt[index])) / 365
    nknots <- floor(years*trend.knots.per.year) + 1
    knots <- seq(min(tt[index]), max(tt[index]), length = nknots)
    
    if (nknots > 2) {
      knots <- knots[-c(1, length(knots))]
      x_t <- splines::ns(tt, knots = knots, intercept = TRUE, Boundary.knots = c(min(tt[index]), max(tt[index])))
    } else{
      knots <- c()
      x_t <- model.matrix(~tt)
    } 
  } else{ ## just an intercept
    x_t <- matrix(1, nrow = length(tt))
    trend.knots.per.year <- NA
    knots <- NA
  }
  
  # trend indices 
  i_t <- 1:ncol(x_t)
  
  #for harmonic model
  yd <- excessmort:::noleap_yday(counts$date)
  x_h <- fourier_trend(yd, k = harmonics)
  
  i_h <- ncol(x_t) + 1:ncol(x_h)
  
  ## build desing matrix
  if (weekday.effect) {
    ## weekday effects
    w <- factor(lubridate::wday(counts$date))
    contrasts(w) <- contr.sum(length(levels(w)), contrasts = TRUE)
    
    x_w <- model.matrix(~w)[, -1] ## intercept already in spline
    i_w <- ncol(x_t) + ncol(x_h) + 1:ncol(x_w)
    x <- cbind(x_t, x_h, x_w)
  } else{
    x <- cbind(x_t, x_h)
  }
  
  y <- counts$outcome
  n <- counts$population
  
  ## fit model
  
  fit <- glm(y[index] ~ x[index,]-1, offset = log(n[index]), family = "quasipoisson")
  dispersion <- pmax(1, summary(fit)$dispersion)
  
  # prepare stuff to return
  if (!extrapolate) {
    for (j in i_t) {
      x[end_extrapolate_index,j] <- x[last_index,j]  
      x[start_extrapolate_index,j] <- x[first_index,j]
    }
  }  
  expected <- exp(x %*% fit$coefficients) * n
  cova <- summary(fit)$cov.unscaled * dispersion
  log_expected_se <- apply(X = x, MARGIN = 1, function(xi){sqrt(t(xi) %*% cova %*% xi)})
  counts$log_expected_se = as.numeric(log_expected_se)
  
  counts <- dplyr::mutate(counts, expected = as.numeric(expected), log_expected_se = as.numeric(log_expected_se), excluded = !index)
  attr(counts, "dispersion") <- dispersion
  attr(counts, "trend.knots.per.year") <- trend.knots.per.year
  attr(counts, "knots") <- knots
  attr(counts, "harmonics") <- harmonics
  attr(counts, "frequency") <- frequency
  attr(counts, "weekday.effect") <- weekday.effect
  
  if (keep.components) {
    
    seasonal <- data.frame(day = seq(1, 365, length = frequency),
                           s = exp(fourier_trend(seq(1, 365, length = frequency), k = harmonics) %*% fit$coefficients[i_h]) - 1)
    
    if (!extrapolate) {
      for (j in 1:ncol(x_t)) {
        x_t[end_extrapolate_index,j] <- x_t[last_index,j]  
        x_t[start_extrapolate_index,j] <- x_t[first_index,j]
      }
    }  
    
    trend <- as.numeric(exp(x_t %*% fit$coefficients[i_t])*frequency*1000)
    
    if (weekday.effect) {
      w <- factor(1:7)
      contrasts(w) <- contr.sum(length(levels(w)), contrasts = TRUE)
      weekday <- data.frame(weekday = 1:7,
                            effect = exp(model.matrix(~w)[, -1] %*% fit$coefficients[i_w]) - 1)
    } else{
      weekday <- NULL
    }
    
    attr(counts, "components") <- list(trend = trend,
                                       log_trend_se = sqrt(diag(x_t %*% cova[i_t, i_t, drop = FALSE] %*% t(x_t))),
                                       seasonal = seasonal,
                                       weekday = weekday)
  }
  attr(counts, "class") <-  append("compute_expected", class(counts))
  attr(counts, "keep.components") <- keep.components
  
  return(counts)
}