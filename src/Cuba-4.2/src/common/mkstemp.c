#define MKSTEMP

#include "stddecl.h"

int mkstemp(char *template)
{
    char *temp;

    temp = _mktemp(template);
    if (!temp)
        return -1;

    return _open(temp, _O_CREAT | _O_TEMPORARY | _O_EXCL | _O_RDWR, _S_IREAD | _S_IWRITE);
}
