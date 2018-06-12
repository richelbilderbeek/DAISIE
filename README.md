# DAISIE

Dynamic Assembly of Island biota through Speciation, Immigration and Extinction in R

This is a development version before the official release on CRAN.

## Installing DAISIE

The DAISIE package has a stable version on CRAN and
a development version on GitHub.

### From CRAN

From within R, do:

```
install.packages("DAISIE")
```

### From GitHub

Because the DAISIE package is located in the folder `DAISIE`, do:

```
devtools::install_github("rsetienne/DAISIE/DAISIE")
```

## Using DAISIE as a package dependency

### From CRAN

To your DESCRIPTION file, add `DAISIE` as any normal package.

If your package directly uses `DAISIE`:

```
Imports:
  DAISIE
```

If your package uses `DAISIE` in its perepherals (e.g. vignettes and tests):

```
Suggests:
  DAISIE
```

### From GitHub

Because the DAISIE package is located in the folder `DAISIE`, do:

```
Removes:
  rsetienne/DAISIE/DAISIE
```
