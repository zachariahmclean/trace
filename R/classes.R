# Fragments class ---------------------------------------------------------
#' fragments object
#'
#' @description
#' An R6 Class representing a fragments object.
#'
#' @details
#' This is the parent class of both fragments and fragments object. The idea is that shared fields and methods are both inherited from this object, but it is not itself directly used.
#' 
fragments <- R6::R6Class("fragments",
  public = list(
    #' @field unique_id unique id of the sample usually the file name
    unique_id = NA_character_,

    #' @field input_method pathway used to import data (either 'fsa', 'size', or 'repeats')
    input_method = NA_character_,

    #' @field metrics_group_id sample grouping for metrics calculations. Associated with `add_metadata()`.
    metrics_group_id = NA_character_,

    #' @field metrics_baseline_control logical to indicate if sample is the baseline control. Associated with `add_metadata()`.
    metrics_baseline_control = FALSE,

    #' @field batch_run_id fragment analysis run. Associated with `add_metadata()`.
    batch_run_id = NA_character_,

    #' @field batch_sample_id An id for the sample used as size standard for repeat calculation. Associated with `add_metadata()`.
    batch_sample_id = NA_character_,

    #' @field batch_sample_modal_repeat Validated repeat length for the modal repeat repeat in that sample. Associated with `add_metadata()`.
    batch_sample_modal_repeat = NA_real_,

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

    #' @field peak_table_df A dataframe containing the fragment peak level information.
    peak_table_df = NULL,

    #' @field repeat_table_df A dataframe containing the fragment peak level information with the repeat size added. May or may not be the same as peak_table_df depending on what options are chosen in `call_repeats`.
    repeat_table_df = NULL,

    #' @description
    #' initialization function that is not used since the child classes are the main object of this package.
    #' @param unique_id unique_id
    #' @param input_method pathway used to import data (either 'fsa', 'size', or 'repeats')
    #' @param object object to be inserted into either 'fsa', 'peak_table_df', or 'repeat_table_df'
    initialize = function(unique_id, input_method, object) {
      self$unique_id <- unique_id
      self$input_method <- input_method
      switch(input_method,
        fsa = self$fsa <- object,
        size = self$peak_table_df <- object,
        repeats = self$repeat_table_df <- object
      )
    },
    #' @description
    #' A function to print informative information to the console
    print = function() {
      print_helper(self,
        sample_attrs = c("unique_id",  "metrics_group_id", "metrics_baseline_control","batch_run_id", "batch_sample_id", "batch_sample_modal_repeat")
      )
    },
    #' @description
    #' plot the trace data
    #' @param show_peaks A logical to say if the called peaks should be plotted on top of the trace. Only valid for fragments objects.
    #' @param x_axis Either "size" or "repeats" to indicate what should be plotted on the x-axis.
    #' @param xlim numeric vector length two specifying the x axis limits
    #' @param ylim numeric vector length two specifying the y axis limits
    #' @param signal_color_threshold A threshold value to colour the peaks relative to the tallest peak. 
    #' @param plot_title A character string for setting the plot title. Defaults to the unique id of the object
    #' @return A base R plot
    plot_trace = function(show_peaks = TRUE,
                          x_axis = NULL,
                          ylim = NULL,
                          xlim = NULL,
                          signal_color_threshold = 0.05,
                          plot_title = NULL) {
      plot_trace_helper(
        fragments = self,
        show_peaks = show_peaks,
        x_axis = x_axis,
        ylim = ylim,
        xlim = xlim,
        signal_color_threshold = signal_color_threshold,
        plot_title = plot_title)

    },
    #' @description
    #' plot the ladder data
    #' @param xlim numeric vector length two specifying the x axis limits
    #' @param ylim numeric vector length two specifying the y axis limits
    #' @param plot_title A character string for setting the plot title. Defaults to the unique id of the object
    #' @return A base R plot
    plot_ladder = function(xlim = NULL, ylim = NULL, plot_title = NULL) {
      plot_ladder_helper(
      self, xlim = xlim, ylim = ylim,
      plot_title = plot_title)
    },
    #' @description
    #' plot the raw data channels in the fsa file. It identifies every channel that has "DATA" in its name.
    #' @return A base R plot
    plot_data_channels = function(){
      plot_data_channels_helper(self)
    },
    #' @description
    #' This returns a list with the allele information for this object.
    get_allele_peak = function(){
      #these have allele_ prefix, because in R if you just call something repeat, it causes many issues
      alleles <- list(
        allele_size = private$allele_size,
        allele_signal = private$allele_signal,
        allele_repeat = private$allele_repeat,
        allele_2_size = private$allele_2_size,
        allele_2_signal = private$allele_2_signal,
        allele_2_repeat = private$allele_2_repeat
      )
      
      return(alleles)
    },

    #' @description
    #' This sets a single allele size/repeat. It searches through the appropriate peak table and finds the closest peak to the value that's provided.
    #' @param allele Either `1` or `2`, indicating which allele information should be set. Allele 1 is the only one used for repeat instability metrics calculations. 
    #' @param unit Either "size" or "repeats" to indicate if the value you're providing is bp size or repeat length.
    #' @param value Numeric vector (length one) of the size/repeat length to set.
    set_allele_peak = function(allele, unit, value){

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
        
      }
      }

      # Dynamically construct the variable names and assign values
      if(allele == 1){
        private$allele_size <- ifelse(!is.na(value) && !is.null(allele_df$size), allele_df$size, NA_real_)
        private$allele_signal <- ifelse(!is.na(value), allele_df$signal, NA_real_)
        private$allele_repeat <- ifelse(!is.null(self$repeat_table_df) && !is.na(value), allele_df$repeats, NA_real_)   
      } else if(allele == 2){
        private$allele_2_size <- ifelse(!is.na(value) && !is.null(allele_df$size), allele_df$size, NA_real_)
        private$allele_2_signal <- ifelse(!is.na(value), allele_df$signal, NA_real_)
        private$allele_2_repeat <- ifelse(!is.null(self$repeat_table_df) && !is.na(value), allele_df$repeats, NA_real_)   
      } else{
        stop("Invalid 'allele' input. Please select between 1 or 2", call. = FALSE)
      }
      


      private$find_main_peaks_used <- TRUE

      invisible(self)
    },

    #' @description
    #' This returns a list with the index peak information for this object.
    get_index_peak = function(){
      index <- list(
        index_repeat = private$index_repeat,
        index_signal = private$index_signal
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

      if(!is.na(value) && nrow(self$repeat_table_df) > 0){
        size_diff <- self$repeat_table_df$repeats - value
        index_df <- self$repeat_table_df[which.min(abs(size_diff)), , drop = FALSE]

        if(nrow(index_df) > 1){
          stop("More than one peak was selected with the value provided", call. = FALSE)
        }
      } else{ 
        # deal with cases where nrow repeat_table_df == 0
        value <- NA_real_
      }
      private$index_repeat <- ifelse(!is.na(value), index_df$repeats, NA_real_)
      private$index_signal <- ifelse(!is.na(value), index_df$signal, NA_real_)
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
    # Fragments settings
    min_bp_size = NULL,
    max_bp_size = NULL,

    # allele data
    allele_size = NA_real_,
    allele_repeat = NA_real_,
    allele_signal = NA_real_,
    allele_2_size = NA_real_,
    allele_2_repeat = NA_real_,
    allele_2_signal = NA_real_,
    find_main_peaks_used = FALSE,

    # call_repeats data
    batch_correction_factor = NA_real_,
    repeat_correction_mod = NULL,
    repeat_correction_factor = NA_real_,
    repeats_not_called_reason = NA_character_,
    repeat_size = NA_real_,
    assay_size_without_repeat = NA_real_,

    #assign_index_peak data
    index_repeat = NA_real_,
    index_signal = NA_real_,
    index_samples = NULL,
    assigned_index_peak_used = FALSE,
    assigned_index_peak_grouped = NULL,

    # metrics calculation data
    metrics_qc_message = NA_character_

  )
)
