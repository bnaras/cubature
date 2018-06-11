# Build Notes

Personal notes for synchronizing to a new version of underlying C
library.

1. Download tagged tar.gz release
2. Replace `clencurt.h` with version generated using `M=16` rather
  than the default `M=19`. Easiest to copy our version of `clencurt.h`
  over
3. Move `Makefile` to `Makefile.orig`
4. Move in `Makefile` from our previous version over
5. Rest as usual
