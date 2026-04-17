############
#+ Functions I use alot. I put them here to save space in others scripts
############

# ****** Function: to get memory usage of objects in environment, ordered by size ******
memory_usage_summary <- function() {
  # Get names and sizes of all objects in the global environment
  obj_info <- sapply(ls(envir = .GlobalEnv), function(obj_name) {
    object.size(get(obj_name, envir = .GlobalEnv))
  })

  # Convert to data frame, sort by size in descending order
  obj_df <- data.frame(
    Object = names(obj_info),
    Size = obj_info
  ) %>%
    arrange(desc(Size)) %>%
    mutate(Size_MB = round(Size / (1024 ^ 2), 2),
           Size_GB = round(Size / 1e9, 3)) %>%
    filter(Size_MB>0)

  # Print the summary table
  print(obj_df, row.names = FALSE)
}


# ****** Function: calculate Euclidean distance ******
euclidean_distance <- function(x1, y1, x2, y2) {
  sqrt((x2 - x1)^2 + (y2 - y1)^2)
}


# ****** FUNCTION: ping banner notification with message ******
notify.me <- function(message = "R completed"){

  message <- paste0("'", message, "'")

  command <- paste0( "terminal-notifier -message ", message, "' -title 'R Notification' -subtitle 'Completion Notice' -sound /System/Library/Sounds/Ping.aiff'")
  return(system(command))
}

# ****** FUNCTION :: Cummulative max that ignores NAs ******
cummax_na <- function(x) {
  # Replace NA with the smallest possible number
  x[is.na(x)] <- -Inf
  cmx <- cummax(x)

  # Replace -Inf back to NA
  cmx[cmx == -Inf] <- NA
  return(cmx)
}

# ****** FUNCTION: clean dirty db download data formats ******
clean.it <- function(a_df, rm_bad_image_codes = TRUE){
  require(dplyr)
  require(tidyr)
  # clean psuedo missing values; fix values with additional characters; fix names with additional characters
  values_to_NA <- c( "null", "NULL", "Null", "na" , "Na","", " ", "\t", "\n", "  ")
  if(rm_bad_image_codes ) values_to_NA <-  c(values_to_NA,  "777", "7777", "888", "8888", "999", "9999", "7","8","9", "removed")

  names(a_df) <- gsub("X..|\\.$", "", names(a_df))

  a_df <- a_df %>%
    mutate_all(~case_when(
      . %in% values_to_NA ~ NA,
      TRUE ~ .
    )) %>%
    mutate_all(~gsub('=\"|\"$', "", .))

  return(a_df)
}

# x <- df[,3]
# ****** FUNCTION: convert columns of db to numeric if no error ******
convert_to_numeric_if_possible <- function(df) {
  df[] <- lapply(df, function(x) {
    if (is.character(x)) {
      converted <- suppressWarnings(as.numeric(x))

      # Check if all elements were successfully converted to numeric
      if (all(!is.na(converted) | is.na(x))) {
        return(converted)
      } else {
        return(x)  # Return as character if conversion isn't fully possible
      }
    } else {
      return(x)  # Return non-character columns as is
    }
  })
  df <- as.data.frame(df)
  return(df)
}

# ****** FUNCTION: convert vector to numeric if no error ******
as_numeric_if_possible <- function(x) {
  if (is.character(x)) {
    converted <- suppressWarnings(as.numeric(x))

    # Check if all elements were successfully converted to numeric
    if (all(!is.na(converted) | is.na(x))) {
      return(converted)
    } else {
      return(x)  # Return as character if conversion isn't fully possible
    }
  } else {
    return(x)  # Return non-character columns as is
  }
}

# ****** FUNCTION to return how long between function calls ******
how.long <- function(clear = FALSE){
  if(clear & ("how.long.time" %in% ls(envir = .GlobalEnv))){
    rm("how.long.time", envir = .GlobalEnv)
    how.long.time = Sys.time()
    return(assign("how.long.time", how.long.time, envir = .GlobalEnv))
  }else{
    if("how.long.time" %in% ls(envir = .GlobalEnv)){
      print(difftime(Sys.time() , how.long.time, units = "auto" ))
      rm("how.long.time", envir = .GlobalEnv)
    }else{
      how.long.time = Sys.time()
      return(assign("how.long.time", how.long.time, envir = .GlobalEnv))
    }
  }
}


# ****** FUNCTION:  plot >0 matrix w/out clustering ******
my.map <- function(mat,title = NULL, centered_on_zero = FALSE) {
  require(ComplexHeatmap)
  require(viridis)
  # if goes nevative, use default colors
  if(centered_on_zero) colors = NULL else colors = plasma(102)
  # plot
  Heatmap(mat, cluster_rows = F, cluster_columns = F, col = colors, row_names_side = "left", column_names_side = "top", column_title = title)
}


# ****** FUNCTION : add lines to base plot; on 0s and diag
add_lines <- function(...,log_transform_sig =FALSE, col=  "black", lty = 3){
  # Set default value to "xyd" if no argument is provided
  arg <- if (length(list(...)) == 0) "xyd" else list(...)[[1]]

  # Check for each component using grepl
  if (grepl("x", arg)) {
    abline(v = 0, lty = lty, col = col)
  }
  if (grepl("y", arg)) {
    abline(h = 0, lty = lty, col = col)
  }
  if (grepl("d", arg)) {
    abline(c(0,1), lty = 2, col = col)
  }
  if(grepl("sigX", arg)){
    abline(v = ifelse(log_transform_sig, -log10(.05),.05),
           lty = 4, col = 'darkgrey')
  }
  if(grepl("sigY", arg)){
    abline(h =  ifelse(log_transform_sig, -log10(.05),.05),
           lty = 4, col = 'darkgrey')
  }
  if(grepl("sigXY|sigYX", arg)){
    abline(v =  ifelse(log_transform_sig, -log10(.05),.05),
           lty = 4, col = 'darkgrey')
    abline(h =  ifelse(log_transform_sig, -log10(.05),.05),
           lty = 4, col = 'darkgrey')
  }
}



# MODELING FUNCTIONS----
# # ****** FUNCTION: to construct formula from terms ******
make.formula <- function(a_response, a_predictor, covariates) {
  # Check if the response variable starts with "I" and wrap in backticks if necessary
  if (grepl("^I", a_response)) {
    response_transformed <- paste0("`", a_response, "`")
  } else {
    response_transformed <- a_response
  }
  if(exists("covariates")){
    # Create the formula string
    formula_string <- paste(
      response_transformed, "~", a_predictor, "+",
      paste(covariates, collapse = " + ")
    )
  }else{
    formula_string <- paste( response_transformed, "~", a_predictor)
  }
  # Convert to formula
  as.formula(formula_string)
}

# ****** FUNCTION: to extract results  ******
# summaries = numeric_mod_summaries
extract_results <- function( summaries = list()){
  # empty df
  results_names <- colnames(coef( summaries[[1]] ))
  stat_cols <- matrix(NA, dimnames = list(1, results_names), ncol = length(results_names))
  results_df <- data.frame(Model = character(), Response = character(), Predictor = character,
                           # DF = integer(),
                           Term = character(),stringsAsFactors = FALSE) %>%
    cbind(stat_cols[-1,])


  # model names must be in format response~predictor or response ~ predictor...
  model_names <- names(summaries)
  response_names <-  model_names %>% str_split(., "~| ~ ") %>% sapply(., `[`, 1)
  predictor_names <- model_names %>% str_split(., "~| ~ ") %>% sapply(., `[`, 2)

  # fill df with model values
  for (i in seq(model_names)) {

    model_name <- model_names[i]
    response_name <-  response_names[i]
    predictor_name <- predictor_names[i]

    # # add degrees of freedom
    # if("df.null"%in% names(summaries[[model_name]])){
    #   model_df <- summaries[[model_name]]$df.null
    #   } else if("df"%in% names(summaries[[model_name]])){
    #     model_df <- summaries[[model_name]]$df[2]
    #    }else if("devcomp" %in% names(summaries[[model_name]])){
    #        model_df <- as.integer(summaries[[model_name]]$devcomp$dims["N"])
    #     } else{
    #       model_df <-NA
    #     }

    model_coef <- as.matrix(coef( summaries[[model_name]] ))

    a_result_df <-data.frame(model_name,
                             response_name,
                             predictor_name,
                             # model_df,
                             rownames(model_coef),
                             model_coef)

    names(a_result_df) <- names(results_df)

    # a_result_df <- data.frame(Model = model_name,
    #                        Term = rownames(model_coef),
    #                        Estimate = model_coef[, 1],
    #                        Std.Error = model_coef[, 2],
    #                        t.value = model_coef[, 3],
    #                        Pr = model_coef[, 4])

    results_df <- rbind(results_df, a_result_df)
  }

  # calculate adjusted p value
  results_df$q_value <- p.adjust(results_df[,grep("Pr\\(", names(results_df))], method = "BH")

  # # filter control variables from results
  # results_df <- results_df %>%
  # filter(!Term %in% c("(Intercept)", "SEX", "month", "patID_eye")) %>%
  # filter(!Model %in% c("patID_eye"))
  #
  # # filter pr of 0
  # results_df <- results_df %>%
  #   filter(Pr>0)

  return(results_df)
}

# ****** FUNCTION: take vector and convert to significance level asterisks ******
sig_asterisks <- function(vec, none = NA_character_){
  asterisks <- case_when(
    vec < .001 ~ "***",
    vec < .01 ~"**",
    vec < .05 ~ "*",
    vec < .1 ~ ".",
    TRUE ~ none)
  return(asterisks)
}


# ****** FUNCTION: ROTATE A N X 2 MATRIX A CERTAIN SPECIFIED DEGREES******
rotate_points <- function(points, angle_degrees) {

  if(angle_degrees<0) angle_degrees <- 360 + angle_degrees

  points <- as.matrix(points)
  # Convert angle from degrees to radians
  angle_radians <- angle_degrees * (pi / 180)

  # Create the rotation matrix
  rotation_matrix <- matrix(c(cos(angle_radians), -sin(angle_radians),
                              sin(angle_radians), cos(angle_radians)),
                            nrow = 2, ncol = 2)

  # Rotate the points
  rotated_points <- points %*% rotation_matrix
  return(rotated_points)
}


# ****** FUNCTION: mirror across x and axis ******
mirror <- function(points, over = NA, silence = FALSE){
  new_points <- as.matrix(points)
  if(ncol(new_points) == 2){
    if(over %in% c(1,"x")){
      new_points[, 2] <- new_points[,2]*-1
      return(new_points)
    }else if(over %in% c(2,"y")){
      new_points[, 1] <- new_points[,1]*-1
      return(new_points)
    }else{
      if(!silence) warning("Over argument not matched")
    }
    return(points)
  } else{
    if(!silence) warning("Points do not have 2 columns")
  }
}
# # test
# n <- 1000
# x <- 1:n
# points <- matrix(c(rnorm(n), rnorm(n)+ log(x)), ncol = 2, nrow = n)
# points <- scale(points)
# plot(points, col = "grey")
# points(mirror(points, "x"), col = rgb(1,.2,.2,.8))


# ****** FUNCTION : Filter for values between whiskers of boxplot ******
filter_outs <- function(vec, type = "identity"){
  bx <- boxplot(vec, plot = FALSE)
  # filter for values between whiskers
  filtered_vec <- vec[vec<bx$stats[5] & vec>bx$stats[1]]

  if(type %in% c("identity", "i", "I")){
    return(filtered_vec)
  }
  if(type %in% c("logical", "l", "L") ){
    return(vec %in%filtered_vec )
  }
  if(type %in% c("index", "ind", "idx")){
    return(which(vec %in%filtered_vec ))
  }
}



# plot reduced dim color by coloring for km_nn_clusters
cluster_cols <- c("1" = "#3A6FA9",
                  "2" = "#B2C6E5",
                  "3" = "#EE8635",
                  "4" = "#F2BF85",
                  "5" = "#509E3E",
                  "6" = "#A8DD93",
                  "7" = "#A12D27",
                  "8" = "#F09D99",
                  "9" = "#77599B",
                  "10" = "#C1B0D3",
                  "11" = "#6E4940",
                  "12" = "#C3A59D",
                  "13" = "#D47CBF",
                  "14" = "#EDB9D1",
                  "15" = "#818281",
                  "16" = "#C7C7C7"
)

# # ***** FUNCTION: assign non canonical task to a single core ******
# assign_task <- function(a_function) {
#   require(future)
#   # Plan for parallel execution
#   plan(multisession, workers = 1)
#   # Execute the function asynchronously
#   my_future <- future({
#     a_function()
#   })
#   # Return the future object for later retrieval
#   return(my_future)
# }
#
# # Define a task (e.g., saving a large image)
# save_large_image <- function() {
#   Sys.sleep(5)  # Simulate time taken to save a large image
#   return("Image saved successfully!")
# }
# how.long(T)
# # Assign the task to a single core
# my_future <- assign_task(save_large_image)
# how.long()
#
# # Run other code while the task runs
# print("Other code is running in the main session...")
# Sys.sleep(2)  # Simulate other work
# print("Still working...")
#
# # Retrieve the result when ready
# result <- value(my_future)
# print(result)

# ****** FUNCTION: grab the items of lookup where match found (similar to grep value = TRUE)******
grab_matches <- function(keys, lookup) {
  lookup[match(keys, names(lookup))]
}

# ***** FUNCTION: prints letter of word every s seconds until resolution ******
waiting <- function(expression_to_resolve, word = "WAITING", sep = " ", end = "\n", time = 1, msg = NULL) {
  # Capture the expression and evaluate it dynamically
  expr <- substitute(expression_to_resolve)

  lett <- strsplit(word, "")[[1]] # Split word into individual letters
  i <- 1

  while (!eval(expr, parent.frame())) {
    if (i > length(lett)) {
      i <- 1
      cat(end) # Print a new line after a complete cycle of letters
    }
    cat(sep, lett[i], sep = "") # Print the current letter with separator
    i <- i + 1
    Sys.sleep(time)
  }

  cat("\nDone!\n") # Print a final message once resolved
  if(!is.null(msg)){
    notify.me(msg)
  }
}

# ****** FUNCTION: find moving average across specified window ******
moving_average <- function(x, n) {
  if (n <= 0 | n > length(x)) stop("Invalid window size")
  filter <- rep(1 / n, n)
  stats::filter(x, filter, sides = 2)
}


# ****** FUNCTION: Distinct colors of length n ******----
# source("~/Documents/personal/R stuff/FIT FILE/distinct_palette.R")

distinct_palette <- function(n = NA, pal = "brewerPlus", add = "lightgrey") {
  stopifnot(rlang::is_string(pal))
  stopifnot(rlang::is_scalar_integerish(n) || identical(n, NA))
  stopifnot(rlang::is_na(add) || is.character(add) || is.numeric(add))

  # define valid palettes matched to retrieval functions
  palList <- list(
    brewerPlus = palBrewerPlus,
    kelly = palKelly,
    greenArmytage = palGreenArmytage
  )

  # match palette request
  pal <- rlang::arg_match0(arg = pal, values = names(palList))

  # get full palette
  palFun <- palList[[pal]]
  palCols <- palFun()

  # get n colors
  if (!rlang::is_na(n)) {
    if (n > length(palCols)) {
      stop("Palette '", pal, "' has ", length(palCols), " colors, not ", n)
    }
    palCols <- palCols[seq_len(n)]
  }

  # add last colour e.g. lightgrey default if requested
  if (!identical(add, NA)) {
    grDevices::col2rgb(add)
    palCols <- c(palCols, add)
  }
  return(palCols)
}

palBrewerPlus <- function() {
  c(
    # first 12 colours generated with:
    # RColorBrewer::brewer.pal(n = 12, name = "Paired")
    "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
    "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928",
    # vivid interlude
    "#1ff8ff", # a bright blue
    # "#FDFF00", # lemon (clashes with #FFFF99 on some screens)
    # "#00FF00", # lime (indistinguishable from bright blue on some screens)
    # next 8 colours generated with:
    # RColorBrewer::brewer.pal(n = 8, "Dark2")
    "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
    "#66A61E", "#E6AB02", "#A6761D", "#666666",
    # list below generated with iwanthue: all colours soft kmeans 20
    # with a couple of arbitrary tweaks by me
    "#4b6a53",
    "#b249d5",
    "#7edc45",
    "#5c47b8",
    "#cfd251",
    "#ff69b4", # hotpink
    "#69c86c",
    "#cd3e50",
    "#83d5af",
    "#da6130",
    "#5e79b2",
    "#c29545",
    "#532a5a",
    "#5f7b35",
    "#c497cf",
    "#773a27",
    "#7cb9cb",
    "#594e50",
    "#d3c4a8",
    "#c17e7f"
  )
}


palKelly <- function() {
  c(
    # "#f2f3f4", "#222222", # white and black removed
    "#f3c300", "#875692", "#f38400", "#a1caf1", "#be0032", "#c2b280",
    "#848482", "#008856", "#e68fac", "#0067a5", "#f99379", "#604e97",
    "#f6a600", "#b3446c", "#dcd300", "#882d17", "#8db600", "#654522",
    "#e25822", "#2b3d26"
  )
}

palGreenArmytage <- function() {
  c(
    "#F0A3FF", "#0075DC", "#993F00", "#4C005C", # "#191919", # black removed
    "#005C31", "#2BCE48", "#FFCC99", "#808080", "#94FFB5", "#8F7C00",
    "#9DCC00", "#C20088", "#003380", "#19A405", "#FFA8BB", "#426600",
    "#FF0010", "#5EF1F2", "#00998F", "#E0FF66", "#100AFF", "#990000",
    "#FFFF80", "#FFE100", "#FF5000"
  )
}

palScripps <- function(n=NULL){
  # pal <- c("#FAEE18","#FEC754","#F68F54","#F0614C","#A33244","#660B3B","#4B1949","#2b325D","#15476A","#1b6174","#12767e","#68aeb5","#bde4eb","#DDFBFA"   )
  # pal <- c("#F5F1A7","#fdc652","#f57d50","#79163a","#4a1846","#1c5570","#5ba4ac","#bde4eb","#DDFBFA")
  # pal <-   c(
  #   "#F5F1A7",
  #   "#FAEE09",
  #   "#FEC754",
  #   "#fdc652",
  #   "#F68F54",
  #   "#f57d50",
  #   "#F0614C",
  #   "#A33244",
  #   "#79163a",
  #   "#660B3B",
  #   "#4a1846",
  #   "#4B1949",
  #   "#2b325D",
  #   "#15476A",
  #   "#1c5570",
  #   "#1b6174",
  #   "#12767e",
  #   "#5ba4ac",
  #   "#68aeb5",
  #   "#bde4eb",
  #   "#DDFBFA"
  #   )
  pal =  c(
    "#FEF7AF",
    # "#FAEE09",
    # "#F5F1A7",
    # "#FEC754",
    "#fdc652",
    "#F68F54",
    "#f57d50",
    "#F0614C",
    "#A33244",
    "#79163a",
    "#660B3B",
    # "#4B1949",
    "#4a1846",
    "#2b325D",
    "#15476A",
    "#1c5570",
    "#1b6174",
    "#12767e",
    # "#5ba4ac",
    "#68aeb5",
    "#bde4eb",
    "#DDFBFA"
  )

  if(is.null(n)){
    barplot(rep(1, 100), col = palScripps(100), border = NA, space = 0)
  }else{
    return(colorRampPalette(pal)(n))
  }
}
# palScripps()


# ****** Function to reorder cols to show id cols first ******
reorder_id_cols <- function(df){
  # lookup<-"patID_eye|patID|patientId|eye|visitId|DATE"
  lookup<- c("patID_eye","pat","eye","visitI[Dd]","DATE_coalsced","DATE")
  vars <- sapply(lookup, function(l)grep(l, names(df), value = TRUE) %>% grep("HRDEX", ., value = TRUE, invert = TRUE)) %>%
    # sapply(`[`,1)
    unlist() %>%
    as.character()
  select(df, all_of(vars), everything())
}


# ****** FUNCTION to purrr::reduce() by breaking in to groups and running parralelly *********
parr_reduce <- function(...){
  require(purrr)
  require(parallel)
  args <- list(...)
  names(args)[sapply(args, function(x) is.list(x)|is.vector(x) )] <- '.x'
  names(args)[sapply(args, is.function)] <- '.f'

  if(!'.x'%in%names(args)){
    stop("please provide a list or vector to reduce")
  }else if(!'.f'%in% names(args)){
    stop("please provide a function")
  }else{

    inds <- cut(seq_along(args$.x), detectCores() - 1, labels = FALSE)
    groups <- split(seq_along(args$.x), inds)

    semi_reduced <- mclapply(groups, function(group_inds) {
      purrr::reduce(.x = args$.x[group_inds], .f = args$.f)
    }, mc.cores = detectCores() - 1)

    full_reduced <- purrr::reduce(semi_reduced, args$.f)

    return(full_reduced)
  }
}

# ****** FUNCTION to print set diff and intersection in up to 3 sets
setdiff_intersect <- function(sets, labels = names(sets), values = FALSE) {
  if (length(sets) < 2 || length(sets) > 3) {
    stop("Only 2 or 3 sets are supported.")
  }

  sets <- lapply(sets, unique) %>%
    lapply(na.omit)

  labels <- if (is.null(labels)) paste0("Set", seq_along(sets)) else labels

  if (length(sets) == 2) {
    a <- sets[[1]]
    b <- sets[[2]]
    la <- labels[1]
    lb <- labels[2]

    only_b <- setdiff(b, a)
    only_a <- setdiff(a, b)
    inter_ab <- intersect(a, b)

    cat(paste("Only", la, ":", length(only_a), "\n"))
    if (values) cat(paste(only_a, collapse = ", "), "\n")

    cat(paste("Only", lb, ":", length(only_b), "\n"))
    if (values) cat(paste(only_b, collapse = ", "), "\n")

    cat(paste(la, "∩", lb, ":", length(inter_ab), "\n"))
    if (values) cat(paste(inter_ab, collapse = ", "), "\n")

  } else if (length(sets) == 3) {
    a <- sets[[1]]
    b <- sets[[2]]
    c <- sets[[3]]
    la <- labels[1]
    lb <- labels[2]
    lc <- labels[3]

    list_out <- list(
     setdiff(a, union(b, c)),
     setdiff(b, union(a, c)),
     setdiff(c, union(a, b)),
     setdiff(intersect(a, b), c),
     setdiff(intersect(a, c), b),
     setdiff(intersect(b, c), a),
     intersect(intersect(a, b), c)
    )
    names(list_out) <- c( paste("Only", la),
                          paste("Only", lb),
                          paste("Only", lc),
                          paste(la, "∩", lb),
                          paste(la, "∩", lc),
                          paste(lb, "∩", lc),
                          paste(la, "∩", lb, "∩", lc)
                          )

    for (label in names(list_out)) {
      vals <- list_out[[label]]
      cat(paste0(label, ": ", length(vals), "\n"))
      if (values) cat(paste(vals, collapse = ", "), "\n\n")
    }
  }
}
# setdiff_intersect(list(a = letters[1:10], b = letters[5:15], c = letters[15:20]), values = TRUE)

# ****** FUNCTION to plot veritcal histogram

VerticalHist <- function(x, xscale = NULL, xwidth, hist,
                         fillCol = "gray80", lineCol = "gray40", bin_scale = 1) {
  ## x (required) is the x position to draw the histogram
  ## xscale (optional) is the "height" of the tallest bar (horizontally),
  ##   it has sensible default behavior
  ## xwidth (required) is the horizontal spacing between histograms
  ## hist (required) is an object of type "histogram"
  ##    (or a list / df with $breaks and $density)
  ## fillCol and lineCol... exactly what you think.
  binWidth <- (hist$breaks[2] - hist$breaks[1])

  if (is.null(xscale)) xscale <- xwidth * 0.90 / max(hist$density)
  n <- length(hist$density)
  x.l <- rep(x, n)
  x.r <- x.l + hist$density * xscale
  y.b <- hist$breaks[1:n]
  y.t <- hist$breaks[2:(n + 1)]

  rect(xleft = x.l, ybottom = y.b, xright = x.r, ytop = y.t,
       col = fillCol, border = lineCol)
}

# ***** FUNCTION spongebob mock text
mock_text <- function(text){
  text_split <- strsplit(tolower(text), "")[[1]]
  mocked <- paste0(ifelse(seq_along(text_split) %% 2 == 0, toupper(text_split), text_split), collapse = "")
  return(mocked)
}

# ****** FUNCTION print spongebob
spongebob <- function(byline = F){
x <- '  ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⣀⣀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
  ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣼⣿⣿⢯⡙⠳⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
  ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⡁⠀⠀⠀⠙⣆⠹⡆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
  ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠸⣇⠀⣀⠀⠀⠘⣦⠹⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
  ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠹⡄⠸⣇⡀⠀⠘⢧⡹⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
  ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⣄⠙⣇⠀⡀⠈⠳⡽⣆⠀⠀⠀⠀⠀⢀⣀⡀⠀⠀⠀⣀⣀⣀⣀⣀⠀⣀⣠⣤⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
  ⠀⠀⠀⠀⠀⠀⠀⠀⠀⣴⠖⠲⣦⣤⣿⣶⠻⠗⠃⠀⠀⢹⣾⡟⠲⣶⣖⣞⠋⠉⠉⠛⠒⣻⣿⣷⣦⡉⠙⠛⠉⣠⠞⠛⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
  ⠀⠀⠀⠀⠀⢀⣀⣴⡾⠁⠘⣿⣿⣿⣿⣿⣿⣄⡀⣀⣠⣾⣿⣷⡄⠙⣿⣿⣷⣤⣤⣄⡀⠈⠙⠛⠛⠁⠀⣠⠞⠁⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
  ⠀⠀⠀⣀⣴⠟⠁⠀⠀⠀⠀⢸⣿⣿⣦⠀⠀⣸⣿⡿⠟⠛⢹⣿⠃⠐⢿⣿⣿⣿⣿⣿⠟⠀⠀⠀⠀⣠⠼⠃⠀⠀⠀⠀⢧⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
  ⠀⣴⡿⠿⢦⣀⣀⡤⠶⠷⠶⠬⣉⡉⣹⡶⠟⠋⠁⠀⢀⣴⣿⠏⠠⢤⣠⣤⠾⠟⠛⡶⠤⠤⠤⠴⠾⡅⠀⠀⢀⣴⣶⣦⠘⣧⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
  ⠂⠹⣧⡀⣶⣶⣶⡆⠀⠀⠀⠀⠀⠀⢿⣷⣤⣤⣴⣶⣿⠿⠋⠀⠀⠀⠀⠀⣀⣤⠾⢿⣷⣀⡀⠀⠀⣿⠀⢠⣿⣿⣿⣿⡇⢸⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣀⠀⠀⠀⠀⠀
  ⠀⠀⠈⢹⣿⡭⠙⣿⠀⢠⡀⠀⠀⠀⠀⠉⠙⠛⠋⠉⠁⠀⠀⠀⠀⠀⡴⠛⠉⣷⡀⢠⣧⠈⠙⠂⠀⣿⠀⠸⣿⣿⣿⣿⠃⢸⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣴⠛⢻⡇⠀⠀⠀⠀
  ⠀⠀⢠⡾⣿⣴⠾⠛⠓⠛⠿⣤⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⣿⡶⠚⠛⠛⠛⠛⠲⣤⡀⠀⢻⡄⠀⠘⠛⠛⠁⠀⠀⠳⣄⠀⠀⠀⢀⣀⡀⣀⡴⠛⣵⠀⣸⡇⠀⠀⠀⠀
  ⠀⠀⢻⡿⣫⣤⣤⣄⡀⠀⠀⠀⢻⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡼⠁⠀⠀⠀⠀⣤⣤⣤⠈⢿⠀⠀⠙⢦⠀⠀⢠⣴⣿⣇⠀⢹⣆⣠⠟⠉⠀⠉⠻⣤⠟⣽⠀⣿⠀⣠⣤⣀⠀
  ⠀⠀⢸⡃⣿⣿⣷⢌⡷⠀⠀⠀⠀⣷⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠁⠀⠀⠀⠀⢾⣿⣿⣮⣷⣼⠀⠀⠀⠸⡇⠀⠀⠸⡿⠃⠀⢸⣿⠁⠀⠀⠀⠀⠀⣿⠀⡿⠀⣿⣸⣿⢏⣿⡦
  ⠀⠀⠘⢧⡙⠻⠿⠛⠁⠀⠀⢀⡾⠁⣀⠀⡀⠀⠀⠀⠀⠀⠀⠀⠘⢧⡀⠀⠀⠀⠈⠻⣿⣿⣯⠟⠀⠀⠀⠀⡇⠀⠀⠀⠀⢀⣠⣾⠁⠀⠀⠀⠀⢀⣤⠋⢀⡇⢀⣿⣿⡟⣼⢁⣿
  ⠀⠀⠀⠈⠹⢦⣄⣀⣀⡤⠶⠋⣀⡴⠋⣸⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢷⣶⠶⠶⠶⢛⣉⣥⠆⠀⠀⢀⠀⠹⣄⠀⠀⢠⣾⣿⡇⠀⠀⢀⣀⣶⣿⣷⣶⣿⡇⠸⣋⡤⠘⠃⣼⠇
  ⠀⠀⠀⠀⠀⠀⠙⢻⣏⣀⣴⠞⠋⠀⢠⡟⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠉⠛⠓⠛⠛⠉⠀⠀⠀⠀⠀⠀⠀⠈⢳⡀⢸⣿⣿⣿⡶⠒⠛⠛⢿⣿⣿⣿⣿⣧⠀⠟⠁⣠⣾⠁⠀
  ⠀⠀⠀⠀⠀⠀⣴⠟⠋⠁⠀⠀⠀⣠⠟⠀⣤⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⡇⠈⠻⠿⠋⠀⠀⣤⣄⠀⢹⠏⠀⠸⣿⣿⣶⣿⣿⣿⡆⠀
  ⠀⠀⠀⠀⠀⠀⠻⣇⣀⣀⣠⣤⠞⠁⠀⠘⠛⠀⠀⠀⠀⠀⢀⣤⣴⣶⣿⣿⣿⣿⣶⣤⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠛⠦⣄⠀⠀⠀⠘⣿⣿⠀⠸⣆⠀⠀⢿⣿⣿⣿⣿⣿⣿⡆
  ⠀⠀⠀⠀⠀⠀⠀⠈⠉⠙⡿⠁⠀⠀⠀⠀⠀⠀⠀⢀⣴⣾⣿⣿⣿⣿⡿⠿⠟⠛⠋⠉⣀⠀⠀⠀⣠⣴⣶⣶⣶⣦⣤⠀⠀⠈⢧⠀⠀⠀⠀⠀⠀⠀⣿⠀⠀⠘⣿⣿⣿⣿⡿⠿⠉
  ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣇⠀⠀⠀⠀⠀⠀⢀⣴⢿⣿⡿⠟⠋⠉⣁⣀⣠⡤⠶⠶⠛⠉⠀⠀⠀⠻⠿⠿⠿⠿⠟⠋⢀⡀⠀⠘⢧⣠⣀⡀⠀⠀⠀⣿⠀⠀⠀⣿⡿⠛⠉⠀⠀⠀
  ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣿⠀⠀⠀⠀⠀⣰⢿⡿⠛⢉⣲⣶⣶⠛⠉⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⣿⣿⡿⠀⠀⠀⠈⠹⠶⡶⢶⡘⣦⡴⠖⠋⠀⠀⠀⠀⠀⠀
  ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⣧⠀⠀⠀⣰⠇⠈⠛⠛⠉⠁⠀⠈⠙⠓⣶⣤⣤⣀⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⢀⣀⣠⣤⣤⣤⣶⣾⠛⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀
  ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⢷⡀⢠⠏⠀⠀⠀⠀⠀⠀⠀⠀⢠⣾⣿⣿⣿⣿⣿⣿⣿⠟⠉⠙⠛⠛⠛⠛⠛⠛⠋⠉⠉⣿⣿⣿⣿⣿⣿⣿⣿⠇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
  ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⠟⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣽⠟⢋⣿⣿⣿⣿⠟⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠛⢿⡏⢹⣿⣿⣿⠟⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
  ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢰⣞⢁⣾⣟⠉⠉⠉⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⡎⢀⣿⠉⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
  ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⣿⣿⣿⣿⡛⠦⣄⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⣇⣸⣷⣦⣀⣀⣠⣶⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
  ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⣿⣿⣿⣿⣿⣿⣷⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠛⠿⣿⣀⠉⢻⣿⣿⣷⠀⠀⠀⠀⠀⠀⠀⠀⠀
  ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⢻⣿⣿⣿⣿⣿⠿⠿⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⣀⣈⣿⣿⣿⣿⡿⠀⠀⠀⠀⠀⠀⠀⠀⠀
  ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠉⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣿⣿⣿⣿⣿⣿⡿⠛⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
  ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⠻⠿⠿⠟⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀'

  # Split x into lines
  lines <- strsplit(x, "\n")[[1]]

  if (byline) {
  # Otherwise, print each line with 1-second delay
    for (l in lines) {
      cat(l, "\n")
      Sys.sleep(1)
  }
  }else{
   cat(x)
 }
}
