#ifndef UTILITIES_H
#define UTILITIES_H

unsigned int countLinesInFile(FILE* fp)
{
    unsigned int lineCount = 0;
    char ch;

    while(!feof(fp))
    {
        // get a character
        ch = fgetc(fp);

        // if it is a newline, then increment linecount
        if(ch == '\n')
            lineCount++;
    }

    return lineCount;
}


#endif
