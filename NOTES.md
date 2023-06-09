# Build Notes

Personal notes for synchronizing to a new version of underlying C
library.

## For Cubature

1. Download tagged tar.gz release
2. Replace `clencurt.h` with version generated using `M=16` rather
  than the default `M=19`. Easiest to copy our version of `clencurt.h`
  over
3. Move `Makefile` to `Makefile.orig`
4. Move in `Makefile` from our previous version over
5. Copy over any includes in the new C library release
6. Rest as usual

## For Cuba

This is now moved to a separate
[repo](https://github.com/bnaras/Cuba), allowing me to track the
changes that Thomas Hahn makes.

The changes needed for an R package are in the `R_pkg` branch and we
use that branch as submodule here, with Unix-like systems being the
default. Any changes to Windows system are isolated in the `Cuba-win`
directory and copied over using `Makefile.win`. 


