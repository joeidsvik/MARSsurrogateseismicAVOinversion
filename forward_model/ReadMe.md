# Instructions for the Code and Variables

This repository contains scripts and functions for working with traveltime data.  
Below is a guide to the main function, input/output variables, and plotting examples.

---

## Main Script

The main function is:

```r
main_fm

## Functions and Variables

### `match_travletimes_mean_trend`

**Input**

* `traveltimes`: Matrix of traveltimes (must be passed as a matrix).
  *This matches the traveltimes with the corresponding prior mean.*

**Output**

* `prior_mean`: Matrix of dimension `((n × m) × 4)`

  * Column 1: `x_g`
  * Column 2: `x_o`
  * Column 3: `x_cl`
  * Column 4: `depth`

---

### `sim_forward_model`

**Inputs**

* `x_g`: Vector of gas saturations (Gaussian)
* `x_o`: Vector of oil saturations (Gaussian)
* `x_cl`: Vector of clay content (Gaussian)
* `traveltimes`: Section of the traveltimes matrix (single location, matrix, or vector).
  *Whatever is passed will be transformed into the appropriate shape.*

**Output**

* Vector of dimension `((n × m) × 2)`
  Format:

  ```
  (R_0, G)_1, (R_0, G)_2, ...
  ```

---

## Plotting Example

To plot traveltimes vs depth:

1. Extract the preferred area from the `traveltimes` matrix.
2. Extract the corresponding area from the `inline` and `xline` matrices.
3. Plot using, for example:

```r
image.plot(...)
