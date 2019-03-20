# Moment-Generating Function method example

An example of using the moment-generating function method for health economic evaluation in R. Also includes discrete event simulation code for comparison.

## Moment-Generating Function method

The moment-generating function (MGF) method is described in

> Snowsill T. A new method for model-based health economic evaluation utilising and extending moment-generating functions. *Med Decis Making* (under review).

Briefly, it expresses quantities of interest (e.g., total discounted costs) as functions of random variables representing event times and calculates their expected values by using the moment-generating functions for those random variables.

## Getting Started

The example depends on some R packages, all of which are available on CRAN:

```R
> install.packages(c("flexsurv", "tidyr", "dplyr", "purrr"))
```

Once you have the necessary packages, set your working directory to the directory containing the R scripts and source `mgf.R` to run the MGF method.

```R
> setwd("/path/to/scripts")
> source("mgf.R")
```

You will then have a data frame (actually a tibble) called `results` with the results of running the MGF method. Substitute `DES.R` for `mgf.R` to run the discrete event simulation instead.

## What is it doing?

The example included is a three state model (typical for health economic evaluations of cancer). The three states are *Stable disease*, *Progressive disease* and *Dead*.

There are three events included in the model:

* Disease progression (from *Stable disease* to *Progressive disease*)
* Death from stable disease (from *Stable disease* to *Dead*)
* Post-progrssion mortality (from *Progressive disease* to *Dead*)

These events are assumed to have time-to-event distributions (from time of entering the state) that are Weibull, Gompertz and log-normal respectively.

There are two treatment arms: *Treatment* and *Control*. The treatment has a proportional hazards impact on disease progression.

Costs are incurred for disease progression and death (one-off costs) as well as being incurred at a constant rate while in the *Stable disease* and *Progressive disease* states.

QALY weights are calculated from a baseline function (quadratic in age) multiplied by state-specific utility multipliers.

## Contributing

There are no plans to update or enhance this example, but you are free to use or adapt it as you wish (see License).

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
