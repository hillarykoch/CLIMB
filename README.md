# CLIMB: Composite LIkelihood eMpirical Bayes

## Dependencies
1.  To use the full extent of this package, you need to download the C++ library LEMON graph library (version 1.3.1); [download it here.](https://lemon.cs.elte.hu/trac/lemon/wiki/Downloads)

If you do not have them already, you do not need all of LEMON's dependencies to run the CLIMB package. To ease this during configuring LEMON, you can go into the INSTALL file after downloading LEMON and add options

```
-DLEMON_ENABLE_GLPK=NO
-DLEMON_ENABLE_COIN=NO
-DLEMON_ENABLE_ILOG=NO
```

<!---
LEMON citation:
Balázs Dezső, Alpár Jüttner, Péter Kovács. LEMON – an Open Source C++ Graph Template Library. Electronic Notes in Theoretical Computer Science, 264:23-45, 2011. Proceedings of the Second Workshop on Generative Technologies (WGT) 2010.
-->

2.  You need a compiler that has support for C++11, such as
    *   GCC: [see here, for example](https://www.gnu.org/software/gcc/projects/cxx-status.html#cxx11)
    *   clang: [see here](http://clang.llvm.org/cxx_status.html)

3. You also need a more recent version of the Julia programming language, version >= 1.0 (**CLIMB was developed with Julia version 1.0.2, but has been tested with version 1.5.3**):
    *   Download Julia [here](https://julialang.org/downloads/)
    *   Or, if you have Homebrew, from the terminal you can type
    ```console
    brew install --cask julia
    ```
4. Julia is used to support a Metropolis-within-Gibbs sampler. You will need to install the companion Julia module to this R package, called cgibbs.jl, to implement this. Please see [my other git repo](https://github.com/hillarykoch/cgibbs.jl) for the brief installation instructions.

## Getting the package
This package is currently only maintained on GitHub. To download just the R package portion of the software, you can open R and do the following commands:

```{r}
library(devtools)
install_github("hillarykoch/CLIMB")
library(CLIMB)
```
## Using the package
There is a manual walking through how to prepare your data, implement the CLIMB methodology, and perform some downstream analyses [here](https://hillarykoch.github.io/climb_page/index.html).
