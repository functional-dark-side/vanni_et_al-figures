nclass.all <- function(x, fun = median)
{
  fun(c(
    nclass.Sturges(x),
    nclass.scott(x),
    nclass.FD(x)
  ))
}

calc_bin_width <- function(x, ...)
{
  rangex <- range(x, na.rm = TRUE)
  (rangex[2] - rangex[1]) / nclass.all(x, ...)
}

get_ratio <- function(x = x, y = y, display = 4/3){
  ratio_display <- display
  ratio_values <- ((max(x) - min(x)))/((max(y) - min(y)))
  ratio_values/ratio_display
}

get_ratio_log10 <- function(x = x, y = y, display = 4/3){
  ratio_display <- display
  ratio_values <- (log10(max(x) - min(x)))/(log10(max(y) - min(y)))
  ratio_values/ratio_display
}
