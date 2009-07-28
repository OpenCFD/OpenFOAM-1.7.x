%{
/*---------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 1991-2009 OpenCFD Ltd.
    \\/      M anipulation   |
------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    wmkdep

Description
    A fast dependency list generator that emulates the behaviour and
    output of cpp -M. However, the output contains no duplications and
    is ~40% faster than cpp.

    The algorithm uses flex to scan for includes and searches the files
    found.  Each file is entered into a hash table so that files are scanned
    only once.  This is why this program is faster than cpp.

Usage
    wmkdep [ -Idirectory ... -Idirectory ] filename

\*---------------------------------------------------------------------------*/

#define FILE_STACK_SIZE 300
#define HASH_TABLE_SIZE 500

#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <dirent.h>

void nextFile(const char* fileName);
void importFile(const char* fileName);
void importDirectory(const char* dirName);

#   undef yywrap	/* sometimes a macro by default */

%}

%x CMNT CFNAME SCFNAME JFNAME FFNAME
%%

"//".*\n                  ;            /* remove c++ style one line comments */

"/*"                      BEGIN(CMNT);    /* start removing c style comment */
<CMNT>.|\n                ;
<CMNT>"*/"                BEGIN(INITIAL); /* end removing c style comment */

^[ \t]*#[ \t]*include[ \t]+\"   BEGIN(CFNAME);  /* c-file name */
<CFNAME>[^"\n ]*        { BEGIN(INITIAL); nextFile(yytext); } /*"*/

^[ \t]*import[ \t]+       BEGIN(JFNAME);  /* java-file name */
<JFNAME>java.*;           BEGIN(INITIAL);
<JFNAME>org.*;            BEGIN(INITIAL);
<JFNAME>com.*;            BEGIN(INITIAL);
<JFNAME>sunw.*;           BEGIN(INITIAL);
<JFNAME>sun.*;            BEGIN(INITIAL);
<JFNAME>launcher.*;       BEGIN(INITIAL);
<JFNAME>[^"\n*]*;       { BEGIN(INITIAL); importFile(yytext); } /*"*/
<JFNAME>[^"\n]*\*;      { BEGIN(INITIAL); importDirectory(yytext); } /*"*/

"      "include[ \t]+\'   BEGIN(FFNAME);  /* FORTRAN-file name */
<FFNAME>[^']*           { BEGIN(INITIAL); nextFile(yytext); } /*'*/

.|\t|\n                   ;

%%


int nDirectories;
char** directories;
char* sourceFile = NULL;
char* sourceExt = NULL;
char* objectFile = NULL;
char* classFile = NULL;
char* depFile = NULL;

int main(int argc, char* argv[])
{
    char *dotPtr, *slashPtr;
    int i;

    if (argc == 1)
    {
        fprintf(stderr,"input file not supplied\n");
        exit(1);
    }

    sourceFile = (char*)malloc(strlen(argv[argc-1]) + 1);
    strcpy(sourceFile, argv[argc-1]);
    fprintf(stderr, "Making dependency list for source file %s\n", sourceFile);


    /* Get list of -I directories. */
    nDirectories = 0;

    for (i = 1; i < argc; i++)
    {
        if (strncmp(argv[i], "-I", 2) == 0)
        {
            nDirectories++;
        }
    }

    directories = (char**)malloc(sizeof(char*)*nDirectories);

    nDirectories = 0;

    for (i = 1; i < argc; i++)
    {
        if (strncmp(argv[i], "-I", 2) == 0)
        {
            directories[nDirectories++] = strdup(argv[i] + 2);
        }
    }



    if ((dotPtr = strrchr(sourceFile, '.')) == NULL)
    {
        fprintf
        (
            stderr,
            "Cannot find extension in source file name %s\n",
            sourceFile
        );
        exit(1);
    }

    if ((slashPtr = strrchr(sourceFile, '/')) == NULL)
    {
        slashPtr = sourceFile;
    }
    else
    {
        slashPtr++;
    }

    sourceExt = (char*)malloc(strlen(sourceFile));
    sourceExt[0] = 0;
    strncat
    (
        sourceExt,
        dotPtr+1,
        (&sourceFile[strlen(sourceFile) - 1] - dotPtr)/sizeof(char)
    );

    /* initialise depFile to zero and use strncat rather than strncpy
       because there is a bug in the SGI strncat that if 0 preceeds the .
       it inserts a space */
    depFile = (char*)malloc(strlen(sourceFile) + 3);
    depFile[0] = 0;
    strncat(depFile, sourceFile, (dotPtr - sourceFile)/sizeof(char));
    strcat(depFile, ".dep");

    if (strcmp(sourceExt, "java") == 0)
    {
        classFile = (char*)malloc(strlen(sourceFile) + 17);
        strcpy(classFile, "$(CLASSES_DIR)/");
        strncat(classFile, sourceFile, (dotPtr - sourceFile)/sizeof(char));
        strcat(classFile, ".class");

        printf("%s: %s\n", classFile, depFile);
    }
    else
    {
        objectFile = (char*)malloc(strlen(sourceFile) + 16);
        strcpy(objectFile, "$(OBJECTS_DIR)/");
        strncat(objectFile, slashPtr, (dotPtr - slashPtr)/sizeof(char));
        strcat(objectFile, ".o");

        printf("%s: %s\n", objectFile, depFile);
    }

    nextFile(sourceFile);
    yylex();


    for (i = 0; i < nDirectories; i++)
    {
        free(directories[i]);
    }
    free(directories);

    free(sourceFile);
    free(sourceExt);
    free(objectFile);
    free(classFile);
    free(depFile);

    return 0;
}


int currentBuffer = 0;                     /* Buffer pointer stack counter */
YY_BUFFER_STATE buffers[FILE_STACK_SIZE];  /* Buffer pointer stack */

/* file name entry in hash table */

struct FileName
{
    char* name;
    struct FileName* next;
};

struct FileName* fileHashTable[HASH_TABLE_SIZE];   /* File hash table */
struct FileName* dirHashTable[HASH_TABLE_SIZE];    /* Directory hash table */


/* lookup name in hash table, if not found insert in table */

int lookUp(struct FileName** hashTable, const char* p)
{
    int ii = 0;
    struct FileName* n;
    struct FileName* nn;

    /* hash */
    const char* pp = p;
    while (*pp) ii = ii<<1 ^ *pp++;
    if (ii < 0) ii = -ii;
    ii %= HASH_TABLE_SIZE;

    /* search */
    for (n=hashTable[ii]; n; n=n->next)
    {
        if (strcmp(p, n->name) == 0)
        {
            /* entry found so return true */
            return 1;
        }
    }

    /* insert */
    nn = (struct FileName*)malloc(sizeof(struct FileName));
    nn->name = (char*)malloc(strlen(p)+1);
    strcpy(nn->name, p);
    nn->next = hashTable[ii];
    hashTable[ii] = nn;

    /* entry not found, and therefore added. return false */
    return 0;
}



/* Add a directory name to the file name */

char* addDirectoryName(const char* directoryName, const char* fileName)
{
    char* pathName;
    pathName = (char*)malloc(strlen(directoryName) + strlen(fileName) + 2);

    strcpy(pathName, directoryName);

    if (directoryName[strlen(directoryName)-1] != '/')
    {
        strcat(pathName, "/");
    }

    strcat(pathName, fileName);

    return pathName;
}


/* open a file and create buffer and put on stack stack */

void nextFile(const char* fileName)
{
    int d;
    char* pathName;

    if (lookUp(fileHashTable, fileName)) return;

    if (currentBuffer >= FILE_STACK_SIZE)
    {
        fprintf
        (
            stderr,
            "depth of file search exceeds stack size %d "
            "while opening %s for file %s\n",
            FILE_STACK_SIZE, fileName, sourceFile
        );
        exit(1);
    }

    /* Pointer to new file which is set if the file is successfully opened */
    FILE* newyyin = NULL;

    if (!(newyyin = fopen(fileName, "r")))
    {
        for (d=0; d<nDirectories; d++)
        {
            pathName = addDirectoryName(directories[d], fileName);

            if ((newyyin = fopen(pathName, "r")))
            {
                printf("%s: %s\n", depFile, pathName);
                buffers[currentBuffer++] = YY_CURRENT_BUFFER;
                yy_switch_to_buffer(yy_create_buffer(newyyin, YY_BUF_SIZE));

                free(pathName);

                return;
            }

            free(pathName);
        }

        fprintf
        (
            stderr,
            "could not open file %s for source file %s\n",
            fileName, sourceFile
        );

        fflush(stdout);
    }
    else
    {
        printf("%s: %s\n", depFile, fileName);
        fflush(stdout);

        buffers[currentBuffer++] = YY_CURRENT_BUFFER;
        yy_switch_to_buffer(yy_create_buffer(newyyin, YY_BUF_SIZE));
    }
}


void dotToSlash(char* fileName)
{
    int i, len;

    len = strlen(fileName);

    for (i=0; i<len; i++)
    {
        if (fileName[i] == '.') fileName[i] = '/';
    }
}


void importFile(const char* fileName)
{
    char* javaFileName;
    int fileNameLen;

    fileNameLen = strlen(fileName);
    javaFileName = (char*)malloc(fileNameLen + 6);
    javaFileName[0] = 0;
    strncat(javaFileName, fileName, fileNameLen - 1);
    dotToSlash(javaFileName);
    strcat(javaFileName, ".java");

    nextFile(javaFileName);

    free(javaFileName);
}


void importDirectory(const char* dirName)
{
    int dirNameLen;
    char *uDirName, *path;
    DIR *source;
    struct dirent *list;

    if (lookUp(dirHashTable, dirName)) return;

    dirNameLen = strlen(dirName);
    uDirName = strdup(dirName);
    dotToSlash(uDirName);
    uDirName[dirNameLen-2] = 0;

    if ((source = opendir(uDirName)) == NULL)
    {
        fprintf
        (
            stderr,
            "could not open directory %s\n",
            uDirName
        );
        free(uDirName);

        return;
    }
    else
    {
        /* Read and parse all the entries in the directory */
        while ((list = readdir(source)) != NULL)
        {
            if
            (
                strstr(list->d_name, ".java")
             && !strstr(list->d_name, ".java~")
            )
            {
                path = addDirectoryName(uDirName, list->d_name);
                nextFile(path);
                free(path);
            }
        }

        closedir(source);
        free(uDirName);
    }
}


/* The lexer calls yywrap to handle EOF conditions */

int yywrap()
{
    /* Close the file for the buffer which has just reached EOF */
    /* This causes strange problems
    fclose(yyin);
    yyin = 0;
    */

    /* Delete the buffer */
    yy_delete_buffer(YY_CURRENT_BUFFER);

    /* Set buffer counter to previous buffer */
    currentBuffer--;

    if (currentBuffer >= 0) /* if buffer counter refers to a valid file */
    {
        /* reset input buffer to the previous buffer on the stack */
        yy_switch_to_buffer(buffers[currentBuffer]);

        /* Return to the normal state for the previous buffer on the stack */
        BEGIN(INITIAL);

        /* return 0 to inform lex to continue reading */
        return 0;
    }
    else                 /* else there are no more buffers on the stack */
    {
        /* return 1 to inform lex finish now that all buffers have been read */
        return 1;
    }
}


/*****************************************************************************/
