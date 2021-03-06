% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scoreComposite.R
\name{scoreComposite}
\alias{scoreComposite}
\title{Composite Scores with Intervals}
\usage{
scoreComposite(returner_list, alpha = 0.05, type = "inst")
}
\arguments{
\item{returner_list}{List of lists, each element should be the return value (a list) from fitBabyMonitor}

\item{alpha}{Constructs (1-alpha)"\%" posterior intervals for the composite score}

\item{type}{Type of fitting to be done. 'inst' for institution composite scores
'subset_baseline' or 'subset_nobaseline' for subset rankings
'inst_subset_nobaseline' or 'inst_subset_baseline' for subset within institution rankings}
}
\description{
Constructs composite scores as well
as (1-alpha)% confidence intervals from scores on multiple performance indicators.
Composite scores constructed by taking the average of each score
}
