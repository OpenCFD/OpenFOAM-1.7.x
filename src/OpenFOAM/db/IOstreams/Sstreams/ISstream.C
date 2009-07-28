/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

\*---------------------------------------------------------------------------*/

#include "ISstream.H"
#include "int.H"
#include "token.H"
#include <cctype>


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

char Foam::ISstream::nextValid()
{
    char c = 0;

    while (true)
    {
        // Get next non-whitespace character
        while (get(c) && isspace(c))
        {}

        // Return if stream is bad
        if (bad() || isspace(c))
        {
            return 0;
        }

        // Is this the start of a C/C++ comment?
        if (c == '/')
        {
            // If cannot get another character, return this one
            if (!get(c))
            {
                return '/';
            }

            if (c == '/')
            {
                // This is the start of a C++ style one-line comment
                while (get(c) && c != '\n')
                {}
            }
            else if (c == '*')
            {
                // This is the start of a C style comment
                while (true)
                {
                    if (get(c) && c == '*')
                    {
                        if (get(c) && c == '/')
                        {
                            break;
                        }
                        else
                        {
                            putback(c);
                        }
                    }

                    if (!good())
                    {
                        return 0;
                    }
                }
            }
            else  // A lone '/' so return it.
            {
                putback(c);
                return '/';
            }
        }
        else  // c is a valid character so return it
        {
            return c;
        }
    }
}


Foam::Istream& Foam::ISstream::read(token& t)
{
    static char numberBuffer[100];

    // Return the put back token if it exists
    if (Istream::getBack(t))
    {
        return *this;
    }

    // Assume that the streams supplied are in working order.
    // Lines are counted by '\n'

    // Get next 'valid character': i.e. proceed through any white space
    // and/or comments until a semantically valid character is hit upon.

    char c = nextValid();

    // Set the line number of this token to the current stream line number
    t.lineNumber() = lineNumber();

    // return on error
    if (!c)
    {
        t.setBad();
        return *this;
    }

    // Analyse input starting with this character.
    switch (c)
    {
        // First check for punctuation characters.

        case token::END_STATEMENT :
        case token::BEGIN_LIST :
        case token::END_LIST :
        case token::BEGIN_SQR :
        case token::END_SQR :
        case token::BEGIN_BLOCK :
        case token::END_BLOCK :
        case token::COLON :
        case token::COMMA :
        case token::ASSIGN :
        case token::ADD :
     // case token::SUBTRACT : // Handled later as the posible start of a number
        case token::MULTIPLY :
        case token::DIVIDE :
        {
            t = token::punctuationToken(c);
            return *this;
        }

        // Strings: enclosed by double quotes.
        case token::BEGIN_STRING :
        {
            putback(c);
            string* sPtr = new string;

            if (!read(*sPtr).bad())
            {
                t = sPtr;
            }
            else
            {
                delete sPtr;
                t.setBad();
            }
            return *this;
        }

        // Numbers: do not distinguish at this point between Types.
        case '-' :
        case '.' :
        case '0' : case '1' : case '2' : case '3' : case '4' :
        case '5' : case '6' : case '7' : case '8' : case '9' :
        {
            bool isScalar = false;

            if (c == '.')
            {
                isScalar = true;
            }

            int i=0;
            numberBuffer[i++] = c;

            while
            (
                is_.get(c)
             && (
                    isdigit(c)
                 || c == '.'
                 || c == 'e'
                 || c == 'E'
                 || c == '+'
                 || c == '-'
                )
            )
            {
                numberBuffer[i++] = c;

                if (!isdigit(c))
                {
                    isScalar = true;
                }
            }
            numberBuffer[i] = '\0';

            setState(is_.rdstate());

            if (!is_.bad())
            {
                is_.putback(c);

                if (i == 1 && numberBuffer[0] == '-')
                {
                    t = token::punctuationToken(token::SUBTRACT);
                }
                else if (isScalar)
                {
                    t = scalar(atof(numberBuffer));
                }
                else
                {
                    long lt = atol(numberBuffer);
                    t = label(lt);

                    // If the integer is too large to be represented as a label
                    // return it as a scalar
                    if (t.labelToken() != lt)
                    {
                        isScalar = true;
                        t = scalar(atof(numberBuffer));
                    }
                }
            }
            else
            {
                t.setBad();
            }

            return *this;
        }

        // Should be a word (which can be a single character)
        default:
        {
            putback(c);
            word* wPtr = new word;

            if (!read(*wPtr).bad())
            {
                if (token::compound::isCompound(*wPtr))
                {
                    t = token::compound::New(*wPtr, *this).ptr();
                    delete wPtr;
                }
                else
                {
                    t = wPtr;
                }
            }
            else
            {
                delete wPtr;
                t.setBad();
            }
            return *this;
        }
    }
}


Foam::Istream& Foam::ISstream::read(char& c)
{
    c = nextValid();
    return *this;
}


Foam::Istream& Foam::ISstream::read(word& str)
{
    static const int maxLen = 1024;
    static const int errLen = 80; // truncate error message for readability
    static char buf[maxLen];

    register int i = 0;
    register int bc = 0;
    char c;

    while (get(c) && word::valid(c))
    {
        if (fail())
        {
            if (i < maxLen-1)
            {
                buf[i] = '\0';
            }
            else
            {
                buf[maxLen-1] = '\0';
            }
            buf[errLen] = '\0';

            FatalIOErrorIn("ISstream::read(word&)", *this)
                << "problem while reading word '" << buf << "'\n"
                << exit(FatalIOError);

            return *this;
        }

        if (i >= maxLen)
        {
            buf[maxLen-1] = '\0';
            buf[errLen] = '\0';

            FatalIOErrorIn("ISstream::read(word&)", *this)
                << "word '" << buf << "' ...\n"
                << "    is too long (max. " << maxLen << " characters)"
                << exit(FatalIOError);

            return *this;
        }

        if (c == token::BEGIN_LIST)
        {
            bc++;
        }
        else if (c == token::END_LIST)
        {
            bc--;

            if (bc == -1)
            {
                break;
            }
        }

        buf[i++] = c;
    }

    if (i == 0)
    {
        FatalIOErrorIn("ISstream::read(word&)", *this)
            << "invalid first character found : " << c
            << exit(FatalIOError);
    }

    buf[i] = '\0';        // Terminator
    str = buf;
    putback(c);

    return *this;
}


Foam::Istream& Foam::ISstream::read(string& str)
{
    static const int maxLen = 1024;
    static const int errLen = 80; // truncate error message for readability
    static char buf[maxLen];

    char c;

    if (!get(c))
    {
        buf[0] = '\0';

        FatalIOErrorIn("ISstream::read(string&)", *this)
            << "cannot read start of string"
            << exit(FatalIOError);

        return *this;
    }

    char endTok = token::END_STRING;

    // Note, we could also handle single-quoted strings here (if desired)
    if (c != token::BEGIN_STRING)
    {
        buf[0] = '\0';

        FatalIOErrorIn("ISstream::read(string&)", *this)
            << "Incorrect start of string character"
            << exit(FatalIOError);

        return *this;
    }

    register int i = 0;
    bool escaped = false;

    while (get(c))
    {
        if (c == endTok)
        {
            if (escaped)
            {
                escaped = false;
                i--;    // overwrite backslash
            }
            else
            {
                // done reading string
                buf[i] = '\0';
                str = buf;
                return *this;
            }
        }
        else if (c == token::NL)
        {
            if (escaped)
            {
                escaped = false;
                i--;    // overwrite backslash
            }
            else
            {
                buf[i] = '\0';
                buf[errLen] = '\0';

                FatalIOErrorIn("ISstream::read(string&)", *this)
                    << "found '\\n' while reading string \""
                    << buf << "...\""
                    << exit(FatalIOError);

                return *this;
            }
        }
        else if (c == '\\')
        {
            escaped = !escaped;    // toggle state (retains backslashes)
        }
        else
        {
            escaped = false;
        }

        buf[i] = c;
        if (i++ == maxLen)
        {
            buf[maxLen-1] = '\0';
            buf[errLen] = '\0';

            FatalIOErrorIn("ISstream::read(string&)", *this)
                << "string \"" << buf << "...\"\n"
                << "    is too long (max. " << maxLen << " characters)"
                << exit(FatalIOError);

            return *this;
        }
    }


    // don't worry about a dangling backslash if string terminated prematurely
    buf[i] = '\0';
    buf[errLen] = '\0';

    FatalIOErrorIn("ISstream::read(string&)", *this)
        << "problem while reading string \"" << buf << "...\""
        << exit(FatalIOError);

    return *this;
}


Foam::Istream& Foam::ISstream::read(label& val)
{
    is_ >> val;
    setState(is_.rdstate());
    return *this;
}


Foam::Istream& Foam::ISstream::read(floatScalar& val)
{
    is_ >> val;
    setState(is_.rdstate());
    return *this;
}


Foam::Istream& Foam::ISstream::read(doubleScalar& val)
{
    is_ >> val;
    setState(is_.rdstate());
    return *this;
}


// read binary block
Foam::Istream& Foam::ISstream::read(char* buf, std::streamsize count)
{
    if (format() != BINARY)
    {
        FatalIOErrorIn("ISstream::read(char*, std::streamsize)", *this)
            << "stream format not binary"
            << exit(FatalIOError);
    }

    readBegin("binaryBlock");
    is_.read(buf, count);
    readEnd("binaryBlock");

    setState(is_.rdstate());

    return *this;
}


Foam::Istream& Foam::ISstream::rewind()
{
    stream().rdbuf()->pubseekpos(0);

    return *this;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


std::ios_base::fmtflags Foam::ISstream::flags() const
{
    return is_.flags();
}


std::ios_base::fmtflags Foam::ISstream::flags(const ios_base::fmtflags f)
{
    return is_.flags(f);
}


// ************************************************************************* //
