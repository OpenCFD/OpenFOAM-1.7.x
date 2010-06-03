/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

%{

#undef yyFlexLexer

 /* ------------------------------------------------------------------------ *\
   ------ local definitions
 \* ------------------------------------------------------------------------ */

#include "IFstream.H"
#include "triSurface.H"
#include "STLpoint.H"
#include "OSspecific.H"

using namespace Foam;

// Dummy yyFlexLexer::yylex() to keep the linker happy. It is not called
//! @cond dummy
int yyFlexLexer::yylex()
{
    FatalErrorIn("yyFlexLexer::yylex()")
        << "Should not have called this function"
        << abort(FatalError);
    return 0;
}
//! @endcond dummy

// Dummy yywrap to keep yylex happy at compile time.
// It is called by yylex but is not used as the mechanism to change file.
// See <<EOF>>
//! @cond dummy
#if YY_FLEX_SUBMINOR_VERSION < 34
extern "C" int yywrap()
#else
int yyFlexLexer::yywrap()
#endif
{
    return 1;
}
//! @endcond dummy


class STLLexer
:
    public yyFlexLexer
{
    // Private data

        label nTriangles_;
        short region_;       // current region
        short maxRegion_;    // max region

        label lineNo_;
        word startError_;

        DynamicList<STLpoint> STLpoints_;
        //DynamicList<STLpoint> STLnormals_;
        DynamicList<label> STLlabels_;
        HashTable<label, word> STLsolidNames_;


public:

    // Constructors

        //- From input stream and the approximate number of vertices in the STL
        STLLexer(istream* is, const label approxNpoints);


    // Member Functions

        //- The lexer function itself
        int lex();

    // Access

        label nTriangles() const
        {
            return nTriangles_;
        }

        DynamicList<STLpoint>& STLpoints()
        {
            return STLpoints_;
        }

        //DynamicList<STLpoint>& STLnormals()
        //{
        //    return STLnormals_;
        //}

        DynamicList<label>& STLlabels()
        {
            return STLlabels_;
        }

        const HashTable<label, word>& STLsolidNames() const
        {
            return STLsolidNames_;
        }
};


STLLexer::STLLexer(istream* is, const label approxNpoints)
:
    yyFlexLexer(is),

    nTriangles_(0),
    region_(-1),
    maxRegion_(0),

    lineNo_(1),

    STLpoints_(approxNpoints),
    //STLnormals_(approxNpoints),
    STLlabels_(approxNpoints)
{}


 /* ------------------------------------------------------------------------ *\
   ------ cppLexer::yylex()
 \* ------------------------------------------------------------------------ */

#define YY_DECL int STLLexer::lex()

%}

one_space             [ \t\f\r]
space                 {one_space}*
some_space            {one_space}+
cspace                ","{space}

alpha                 [_A-Za-z]
digit                 [0-9]
dec_digit             [0-9]
octal_digit           [0-7]
hex_digit             [0-9a-fA-F]

identifier            {alpha}({alpha}|{digit})*
integer               {dec_digit}+
label                 [1-9]{dec_digit}*
zeroLabel             {digit}*
signedInteger         [-+]?{integer}

word                  ([[:alnum:]]|[[:punct:]])*
string                {word}({some_space}{word})*

exponent_part         [eE][-+]?{digit}+
fractional_constant   [-+]?(({digit}*"."{digit}+)|({digit}+"."?))

double                (({fractional_constant}{exponent_part}?)|({digit}+{exponent_part}))
float                 {double}

x                     {float}
y                     {float}
z                     {float}

solid                 {space}("solid"|"SOLID"){space}
color                 {space}("color"|"COLOR"){some_space}{float}{some_space}{float}{some_space}{float}{space}
facet                 {space}("facet"|"FACET"){space}
normal                {space}("normal"|"NORMAL"){space}
point                 {space}{x}{some_space}{y}{some_space}{z}{space}
outerloop             {space}("outer"{some_space}"loop")|("OUTER"{some_space}"LOOP"){space}
vertex                {space}("vertex"|"VERTEX"){space}
endloop               {space}("endloop"|"ENDLOOP"){space}
endfacet              {space}("endfacet"|"ENDFACET"){space}
endsolid              {space}("endsolid"|"ENDSOLID")({some_space}{word})*


 /* ------------------------------------------------------------------------ *\
                      -----  Exclusive start states -----
 \* ------------------------------------------------------------------------ */

%option stack

%x readSolidName
%x readFacet
%x readNormal
%x readVertices
%x readVertex
%x stlerror

%%

%{
    // End of read character pointer returned by strtof
    char* endPtr;

    STLpoint normal;
    STLpoint vertex;
    label cmpt = 0;   // component index used for reading vertex

    static const char* stateNames[7] =
    {
        "reading solid",
        "reading solid name",
        "reading facet",
        "reading normal",
        "reading vertices",
        "reading vertex",
        "error"
    };

    static const char* stateExpects[7] =
    {
        "'solid', 'color', 'facet' or 'endsolid'",
        "<string>",
        "'normal', 'outer loop' or 'endfacet'",
        "<float> <float> <float>",
        "'vertex' or 'endloop'",
        "<float> <float> <float>",
        ""
    };
%}


 /* ------------------------------------------------------------------------ *\
                            ------ Start Lexing ------
 \* ------------------------------------------------------------------------ */

 /*                      ------ Reading control header ------                */

{solid} {
        BEGIN(readSolidName);
    }

<readSolidName>{string} {
        word solidName(Foam::string::validate<word>(YYText()));
        if (STLsolidNames_.found(solidName))
        {
            region_ = STLsolidNames_[solidName];
        }
        else
        {
            region_ = maxRegion_++;
            STLsolidNames_.insert(solidName, region_);
        }
        BEGIN(INITIAL);
    }

<readSolidName>{space}\n {
        word solidName("solid");
        if (STLsolidNames_.found(solidName))
        {
            region_ = STLsolidNames_[solidName];
        }
        else
        {
            region_ = maxRegion_++;
            STLsolidNames_.insert(solidName, region_);
        }

        lineNo_++;

        BEGIN(INITIAL);
    }

{color} {
    }

{facet} {
        BEGIN(readFacet);
    }

<readFacet>{normal} {
        BEGIN(readNormal);
    }

<readNormal>{point} {
        /*
        normal.x() = strtof(YYText(), &endPtr);
        normal.y() = strtof(endPtr, &endPtr);
        normal.z() = strtof(endPtr, &endPtr);
        STLnormals_.append(normal);
        */
        BEGIN(readFacet);
    }

<readFacet>{outerloop} {
        BEGIN(readVertices);
    }

<readVertices>{vertex} {
        BEGIN(readVertex);
    }

<readVertex>{space}{signedInteger}{space} {
        vertex[cmpt++] = atol(YYText());

        if (cmpt == 3)
        {
            cmpt = 0;
            STLpoints_.append(vertex);
            BEGIN(readVertices);
        }
    }

<readVertex>{space}{float}{space} {
        vertex[cmpt++] = atof(YYText());

        if (cmpt == 3)
        {
            cmpt = 0;
            STLpoints_.append(vertex);
            BEGIN(readVertices);
        }
    }

<readVertices>{endloop} {
        BEGIN(readFacet);
    }

<readFacet>{endfacet} {
        nTriangles_++;
        STLlabels_.append(region_);
        BEGIN(INITIAL);
    }

{endsolid} {
    }


 /* ------------------ Ignore remaining space and \n s. -------------------- */

<*>{space} {}
<*>\n      { lineNo_++; }


 /* ------------------- Any other characters are errors -------------------- */

<*>. {
        startError_ = YYText();
        yy_push_state(stlerror);
    }


 /* ---------------------------- Error handler ----------------------------- */

<stlerror>.* {
        yy_pop_state();
        FatalErrorIn
        (
            "triSurface::readSTLASCII(const fileName& STLfileName)"
        )   << "while " << stateNames[YY_START] << " on line " << lineNo_ << nl
            << "    expected " << stateExpects[YY_START]
            << " but found '" << startError_.c_str() << YYText() << "'"
            << exit(FatalError);
    }


 /*  ------------------------ On EOF terminate ----------------------------  */

<<EOF>> {
            yyterminate();
    }
%%


#include <fstream>

bool triSurface::readSTLASCII(const fileName& STLfileName)
{
    IFstream STLstream(STLfileName);

    if (!STLstream)
    {
        FatalErrorIn
        (
            "triSurface::readSTLASCII(const fileName&)"
        )   << "file " << STLfileName << " not found"
            << exit(FatalError);
    }

    // Create the lexer obtaining the approximate number of vertices in the STL
    // from the file size
    STLLexer lexer(&STLstream.stdStream(), Foam::fileSize(STLfileName)/400);
    while(lexer.lex() != 0)
    {}

    DynamicList<STLpoint>& STLpoints = lexer.STLpoints();

    /*
    DynamicList<STLpoint>& STLnormals = lexer.STLnormals();

    if (STLpoints.size() != 3*STLnormals.size())
    {
        FatalErrorIn
        (
            "triSurface::readSTLASCII(const fileName& STLfileName)"
        )   << "in file " << STLfileName << endl
            << "Problem: read " << STLnormals.size() << " normals"
            << " but " << STLpoints.size() << " points"
            << exit(FatalError);
    }
    */

    pointField rawPoints(STLpoints.size());

    forAll(rawPoints, i)
    {
        rawPoints[i] = STLpoints[i];
    }

    STLpoints.clear();

    setSize(lexer.nTriangles());
    DynamicList<label>& STLlabels = lexer.STLlabels();

    label pointi = 0;
    forAll(*this, i)
    {
        operator[](i)[0] = pointi++;
        operator[](i)[1] = pointi++;
        operator[](i)[2] = pointi++;
        operator[](i).region() = STLlabels[i];
    }

    STLlabels.clear();

    // Stitch all points within SMALL meters.
    stitchTriangles(rawPoints);

    // Convert solidNames into regionNames
    patches_.setSize(lexer.STLsolidNames().size());

    for
    (
        HashTable<label, word>::const_iterator iter =
            lexer.STLsolidNames().begin();
        iter != lexer.STLsolidNames().end();
        ++iter
    )
    {
        patches_[iter()].name() = iter.key();
    }

    // Fill in the missing information in the patches
    setDefaultPatches();

    return true;
}


 /* ------------------------------------------------------------------------ *\
    ------ End of STLToFoam.L
 \* ------------------------------------------------------------------------ */
