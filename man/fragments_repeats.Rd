% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes.R
\name{fragments_repeats}
\alias{fragments_repeats}
\title{fragments_repeats object}
\description{
An R6 Class representing a fragments_repeats object.
}
\details{
The idea behind this class is to store data for processing of the peak level data towards calculation of repeat instability metrics.

It contains important setters and getters for alleles and index peaks. It's very important that the exactly correct size and repeat value is set for the alleles and index peak. This is used for subsetting etc, so if it's not exactly correct many functions would break.

It also contains methods for plotting the ladder and traces (if available).
}
\section{Super class}{
\code{\link[trace:fragments]{trace::fragments}} -> \code{fragments_repeats}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{trace_bp_df}}{A dataframe of bp size for every scan from \code{find_ladders()}.}

\item{\code{peak_table_df}}{A dataframe containing the fragment peak level information.}

\item{\code{repeat_table_df}}{A dataframe containing the fragment peak level information with the repeat size added. May or may not be the same as peak_table_df depending on what options are chosen in \code{call_repeats}.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-fragments_repeats-get_allele_peak}{\code{fragments_repeats$get_allele_peak()}}
\item \href{#method-fragments_repeats-set_allele_peak}{\code{fragments_repeats$set_allele_peak()}}
\item \href{#method-fragments_repeats-get_index_peak}{\code{fragments_repeats$get_index_peak()}}
\item \href{#method-fragments_repeats-set_index_peak}{\code{fragments_repeats$set_index_peak()}}
\item \href{#method-fragments_repeats-plot_fragments}{\code{fragments_repeats$plot_fragments()}}
\item \href{#method-fragments_repeats-clone}{\code{fragments_repeats$clone()}}
}
}
\if{html}{\out{
<details open><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="trace" data-topic="fragments" data-id="initialize"><a href='../../trace/html/fragments.html#method-fragments-initialize'><code>trace::fragments$initialize()</code></a></span></li>
<li><span class="pkg-link" data-pkg="trace" data-topic="fragments" data-id="plot_trace"><a href='../../trace/html/fragments.html#method-fragments-plot_trace'><code>trace::fragments$plot_trace()</code></a></span></li>
<li><span class="pkg-link" data-pkg="trace" data-topic="fragments" data-id="print"><a href='../../trace/html/fragments.html#method-fragments-print'><code>trace::fragments$print()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-fragments_repeats-get_allele_peak"></a>}}
\if{latex}{\out{\hypertarget{method-fragments_repeats-get_allele_peak}{}}}
\subsection{Method \code{get_allele_peak()}}{
This returns a list with the allele information for this object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{fragments_repeats$get_allele_peak()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-fragments_repeats-set_allele_peak"></a>}}
\if{latex}{\out{\hypertarget{method-fragments_repeats-set_allele_peak}{}}}
\subsection{Method \code{set_allele_peak()}}{
This sets a single allele size/repeat. It searches through the appropriate peak table and finds the closest peak to the value that's provided.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{fragments_repeats$set_allele_peak(unit, value)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{unit}}{Either "size" or "repeats" to indicate if the value you're providing is bp size or repeat length.}

\item{\code{value}}{Numeric vector (length one) of the size/repeat length to set.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-fragments_repeats-get_index_peak"></a>}}
\if{latex}{\out{\hypertarget{method-fragments_repeats-get_index_peak}{}}}
\subsection{Method \code{get_index_peak()}}{
This returns a list with the index peak information for this object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{fragments_repeats$get_index_peak()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-fragments_repeats-set_index_peak"></a>}}
\if{latex}{\out{\hypertarget{method-fragments_repeats-set_index_peak}{}}}
\subsection{Method \code{set_index_peak()}}{
This sets the index repeat length. It searches through the repeat table and finds the closest peak to the value that's provided.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{fragments_repeats$set_index_peak(value)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{value}}{Numeric vector (length one) of the repeat length to set as index peak.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-fragments_repeats-plot_fragments"></a>}}
\if{latex}{\out{\hypertarget{method-fragments_repeats-plot_fragments}{}}}
\subsection{Method \code{plot_fragments()}}{
This plots the peak/repeat table as a histogram
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{fragments_repeats$plot_fragments(ylim = NULL, xlim = NULL, plot_title = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{ylim}}{numeric vector length two specifying the y axis limits}

\item{\code{xlim}}{numeric vector length two specifying the x axis limits}

\item{\code{plot_title}}{A character string for setting the plot title. Defaults to the unique id of the object}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-fragments_repeats-clone"></a>}}
\if{latex}{\out{\hypertarget{method-fragments_repeats-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{fragments_repeats$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
