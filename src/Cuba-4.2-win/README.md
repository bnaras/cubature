# Cuba 4.2 changes for Windows


The contents of this directory are the changes I (Balasubramanian
Narasimhan) made to `Cuba-4.2` library source to get it to compile
under Windows using
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) toolchain.

The changes are quite minor.

This organization allows me to ensure that the `Cuba-4.2` tree remains
untouched.

## How to effect the changes

All the files in this tree with the exception of `README.md` should be
copied over the analogous files in the original tree replacing the
originals where they exist.

The `Makevars.win` file in package `src` directory file essentially
performs this task for build the R package.

