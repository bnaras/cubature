# Build Notes

Personal notes for synchronizing to a new version of underlying C
library.

1. Download tagged tar.gz release
2. Replace `clencurt.h` with version generated using `M=16` rather
  than the default `M=19`. Easiest to copy our version of `clencurt.h`
  over
3. Move `Makefile` to `Makefile.orig`
4. Move in `Makefile` from our previous version over
5. Copy over any includes in the new C library release
6. Rest as usual

## For Cuba

I've tried to keep the sources pristine and apply changes to via
copies of files. However, some files were changed in place, notably,
`src/common/Parallel.c` and `src/suave/Integrate.c`. The news file
should have details.

