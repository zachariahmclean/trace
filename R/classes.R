# Fragments class ---------------------------------------------------------
#' fragments object
#'
#' @description
#' An R6 Class representing a fragments object.
#'
#' @details
#' This is the parent class of both fragments_trace and fragments_repeats object. The idea is that shared fields and methods are both inherited from this object, but it is not itself directly used.
#' 
fragments <- R6::R6Class("fragments",
  public = list(
    #' @field unique_id unique id of the sanmple usually the file name
    unique_id = NA_character_,

    #' @field metrics_group_id sample grouping for metrics calculations. Associated with `add_metdata()`.
    metrics_group_id = NA_character_,

    #' @field metrics_baseline_control logical to indicate if sample is the baseline control. Associated with `add_metdata()`.
    metrics_baseline_control = FALSE,

    #' @field batch_run_id fragment analsyis run. Associated with `add_metdata()`.
    batch_run_id = NA_character_,

    #' @field batch_sample_id An id for the sample used as size standard for repeat calculation. Associated with `add_metdata()`.
    batch_sample_id = NA_character_,

    #' @description
    #' initialization function that is not used since the child classes are the main object of this package.
    #' @param unique_id unique_id
    initialize = function(unique_id) {
      if (length(unique_id) != 1) stop("Fragments must have a single unique id", call. = FALSE)
      self$unique_id <- unique_id
    },
    #' @description
    #' A function to print informative information to the console
    print = function() {
      print_helper(self,
        sample_attrs = c("unique_id", "batch_run_id", "metrics_group_id", "metrics_baseline_control", "size_standard", "batch_sample_id", "size_standard_repeat_length")
      )
    },
    #' @description
    #' plot the trace data
    #' @param show_peaks A logical to say if the called peaks should be overlayed on top of the trace. Only valid for fragments_repeats objects.
    #' @param x_axis Either "size" or "repeats" to indicate what should be plotted on the x-axis.
    #' @param xlim numeric vector length two specifying the x axis limits
    #' @param ylim numeric vector length two specifying the y axis limits
    #' @param height_color_threshold A thershold value to colour the peaks relative to the tallest peak. 
    #' @param plot_title A character string for setting the plot title. Defaults to the unique id of the object
    #' @return A base R plot
    plot_trace = function(show_peaks = TRUE,
                          x_axis = NULL,
                          ylim = NULL,
                          xlim = NULL,
                          height_color_threshold = 0.05,
                          plot_title = NULL) {
      plot_trace_helper(
        fragments = self,
        show_peaks = show_peaks,
        x_axis = x_axis,
        ylim = ylim,
        xlim = xlim,
        height_color_threshold = height_color_threshold,
        plot_title = plot_title)

    }
  ),
  private = list(
    min_bp_size = NULL,
    max_bp_size = NULL
  )
)

#' fragments_trace object
#'
#' @description
#' An R6 Class representing a fragments_trace object.
#'
#' @details
#' The idea behind this class is to store data for processing of the continuous trace-level information from an fsa file towards peak level data.
#' 
#' It also contains methods for plotting the ladder and traces
#' 
fragments_trace <- R6::R6Class(
  "fragments_trace",
  inherit = fragments,
  public = list(
    #' @field unique_id unique id of the sanmple usually the file name
    unique_id = NULL,

    #' @field fsa The whole fsa file, output from seqinr::read.abif()
    fsa = NULL,

    #' @field raw_ladder The raw data from the ladder channel
    raw_ladder = NULL,

    #' @field raw_data The raw data from the sample channel
    raw_data = NULL,

    #' @field scan The scan number
    scan = NULL,

    #' @field off_scale_scans vector indicating which scales were too big and off scale. Note can be in any channel
    off_scale_scans = NULL,

    #' @field ladder_df A dataframe of the identified ladder from `find_ladders()`. Scan is the scan number of peak and size is the associated bp size. 
    ladder_df = NULL,

    #' @field trace_bp_df A dataframe of bp size for every scan from `find_ladders()`.
    trace_bp_df = NULL,

    #' @field local_southern_mod Output from `local_southern()` function (not exported). It is bascially just a list of lm() after breaking up the ladder into chunks of three.
    local_southern_mod = NULL,

    #' @description
    #' Create a new fragments_trace.
    #' @param unique_id usually the file name
    #' @param fsa_file output from seqinr::read.abif()
    #' @param ladder_channel the name of the channel in the fsa file that contains the ladder data
    #' @param signal_channel the name of the channel in the fsa file that contains the trace data
    #' @return A new `fragments_trace` object.
    initialize = function(
      unique_id, 
      fsa_file,
      ladder_channel,
      signal_channel) {
        if (length(unique_id) != 1) stop("Fragments must have a single unique id", call. = FALSE)
        self$unique_id <- unique_id
        self$fsa <- fsa_file
        self$raw_ladder <- self$fsa$Data[[ladder_channel]]
        self$raw_data <- self$fsa$Data[[signal_channel]]
        self$scan <- 0:(length(self$fsa$Data[[signal_channel]]) - 1)
        self$off_scale_scans <- self$fsa$Data$OfSc.1
    },
    #' @description
    #' plot the ladder data
    #' @param xlim numeric vector length two specifying the x axis limits
    #' @param ylim numeric vector length two specifying the y axis limits
    #' @param plot_title A character string for setting the plot title. Defaults to the unique id of the object
    #' @return A base R plot
    plot_ladder = function(xlim = NULL, ylim = NULL,
                           plot_title = NULL) {
      plot_ladder_helper(
        self, xlim = xlim, ylim = ylim,
        plot_title = plot_title)
    },
    #' @description
    #' plot the raw data channels in the fsa file. It identifies every channel that has "DATA" in its name.
    #' @return A base R plot
    plot_data_channels = function(){
      plot_data_channels_helper(self)
    }
  )
)



#' fragments_repeats object
#'
#' @description
#' An R6 Class representing a fragments_repeats object.
#'
#' @details
#' The idea behind this class is to store data for processing of the peak level data towards calculation of repeat instability metrics.
#' 
#' It contains important setters and getters for alleles and index peaks. It's very important that the exactly correct size and repeat value is set for the alleles and index peak. This is used for subsetting etc, so if it's not exactly correct many functions would break.
#' 
#' It also contains methods for plotting the ladder and traces (if available).
#' 
fragments_repeats <- R6::R6Class(
  "fragments_repeats",
  inherit = fragments,
  public = list(
    #' @field trace_bp_df A dataframe of bp size for every scan from `find_ladders()`.
    trace_bp_df = NULL,

    #' @field peak_table_df A dataframe containing the fragment peak level informtation.
    peak_table_df = NULL,

    #' @field repeat_table_df A dataframe containing the fragment peak level informtation with the repeat size added. May or may not be the same as peak_table_df depending on what options are chosen in `call_repeats`.
    repeat_table_df = NULL,

    #' @description
    #' This returns a list with the allele information for this object.
    get_alleles = function(){
      alleles <- list(
        allele_1_size = private$allele_1_size,
        allele_1_height = private$allele_1_height,
        allele_1_repeat = private$allele_1_repeat,
        allele_2_size = private$allele_2_size,
        allele_2_height = private$allele_2_height,
        allele_2_repeat = private$allele_2_repeat
      )
      return(alleles)
    },

    #' @description
    #' This sets a single allele size/repeat. It searches through the appropriate peak table and finds the closest peak to the value that's provided.
    #' @param allele Either `1` or `2`, indicating which allele information should be set. Allele 1 is the only one used for repeat instability metrics calculations. 
    #' @param unit Either "size" or "repeats" to indicate if the value you're providing is bp size or repeat length.
    #' @param value Numeric vector (length one) of the size/repeat length to set.
    set_allele = function(allele, unit, value){

      if(!is.na(value)){
        if(is.null(self$repeat_table_df)){
          if(unit != "size") stop("Only size can be used to set alleles if repeats have not been called", call. = FALSE )
          df <- self$peak_table_df
        } else{
          df <- self$repeat_table_df
        }

        size_diff <- df[[unit]]- value
        allele_df <- df[which.min(abs(size_diff)), , drop = FALSE]

        if(nrow(allele_df) > 1){
          stop("More than one peak was selected with the value provided", call. = FALSE)
        }

        # Ensure the allele is either 1 or 2
        if (!(allele %in% c(1, 2))) {
          stop("Invalid 'allele' input. Please select between 1 or 2", call. = FALSE)
        }
      }
      # Dynamically construct the variable names and assign values
      private[[paste0("allele_", allele, "_size")]] <- ifelse(!is.na(value), allele_df$size, NA_real_)
      private[[paste0("allele_", allele, "_height")]] <- ifelse(!is.na(value), allele_df$height, NA_real_)
      private[[paste0("allele_", allele, "_repeat")]] <- ifelse(!is.null(self$repeat_table_df) && !is.na(value), allele_df$repeats, NA_real_)      
      private$find_main_peaks_used <- TRUE

      invisible(self)
    },

    #' @description
    #' This returns a list with the index peak information for this object.
    get_index_peak = function(){
      index <- list(
        index_repeat = private$index_repeat,
        index_height = private$index_height
      )
      return(index)
    },

    #' @description
    #' This sets the index repeat length. It searches through the repeat table and finds the closest peak to the value that's provided.
    #' @param value Numeric vector (length one) of the repeat length to set as index peak.
    set_index_peak = function(value){
      if(is.null(self$repeat_table_df)){
        stop("Index assignment requires repeats to be called", call. = FALSE )
      }

      if(!is.na(value)){
        size_diff <- self$repeat_table_df$repeats- value
        index_df <- self$repeat_table_df[which.min(abs(size_diff)), , drop = FALSE]

        if(nrow(index_df) > 1){
          stop("More than one peak was selected with the value provided", call. = FALSE)
        }
      }
      private$index_repeat <- ifelse(!is.na(value), index_df$repeats, NA_real_)
      private$index_height <- ifelse(!is.na(value), index_df$height, NA_real_)
      private$assigned_index_peak_used <- TRUE
      
      invisible(self)
    },

    #' @description
    #' This plots the peak/repeat table as a histogram
    #' @param xlim numeric vector length two specifying the x axis limits
    #' @param ylim numeric vector length two specifying the y axis limits
    #' @param plot_title A character string for setting the plot title. Defaults to the unique id of the object
    plot_fragments = function(ylim = NULL,
                              xlim = NULL,
                              plot_title = NULL) {
      plot_fragments_helper(self,
                            ylim = ylim,
                            xlim = xlim,
                            plot_title = plot_title)

    }
  ),
  private = list(
    find_main_peaks_used = FALSE,
    peak_regions = NA_real_,
    allele_1_size = NA_real_,
    allele_1_repeat = NA_real_,
    allele_1_height = NA_real_,
    allele_2_size = NA_real_,
    allele_2_repeat = NA_real_,
    allele_2_height = NA_real_,
    index_repeat = NA_real_,
    index_height = NA_real_,
    batch_correction_factor = NA_real_,
    repeats_not_called_reason = NA_character_,
    controls_repeats_df = NULL,
    assigned_index_peak_used = FALSE,
    index_samples = NULL
  )
)
