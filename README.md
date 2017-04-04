# Marseille

[![Travis-CI Build Status](https://travis-ci.org/SWotherspoon/Marseille.svg?branch=master)](https://travis-ci.org/SWotherspoon/Marseille)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/SWotherspoon/Marseille?branch=master&svg=true)](https://ci.appveyor.com/project/SWotherspoon/Marseille)

Marseille provides basic support for the RAATD project. Typical users
will be more interested in
[RWalc](https://github.com/SWotherspoon/RWalc).

## Installing

Installing Marseille requires 

* a recent version of R (> 3.3.0),
* windows users will require [Rtools](https://cran.r-project.org/bin/windows/Rtools/),
* the [TMB](https://cran.r-project.org/web/packages/TMB/index.html)
  package and its dependencies.


Marseille is easily installed from GitHub using the devtools package. 

```R
devtools::install_github("SWotherspoon/Marseille")
```

If you don't have `devtools` installed already, install it first. 

```R
install.packages("devtools")
```

Marseille otherwise does not need devtools for normal use.


