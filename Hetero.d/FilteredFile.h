#ifndef _FILTEREDFILE_H_
#define _FILTEREDFILE_H_

#include <stdio.h>

class FilteredFile {
    FILE *file;
    char buffer[200];
    char *epos;
  public:
    FilteredFile(char *filename);
    int findToken(char *token);
    char *getLine();
    char *getLineAfterToken();
};

#endif
