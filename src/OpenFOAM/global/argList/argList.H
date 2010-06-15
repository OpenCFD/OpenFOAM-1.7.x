/*---------------------------------------------------------------------------*\
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

Class
    Foam::argList

Description
    Extract command arguments and options from the supplied
    @a argc and @a argv parameters.

    Sequences with "(" ... ")" are transformed into a stringList.
    For example,
    @verbatim
        program -listFiles \( *.txt \)
    @endverbatim
    would create a stringList:
    @verbatim
        ( "file1.txt" "file2.txt" ... "fileN.txt" )
    @endverbatim
    The backslash-escaping has been used to avoid shell expansions.

    @par Default command-line options
    @param -case \<dir\> \n
        select an case directory instead of the current working directory
    @param -parallel \n
        specify case as a parallel job
    @param -doc \n
        display the documentation in browser
    @param -srcDoc \n
        display the source documentation in browser
    @param -help \n
       print the usage

    The environment variable @b FOAM_CASE is set to the path of the
    global case (same for serial and parallel jobs).
    The environment variable @b FOAM_CASENAME is set to the name of the
    global case.

Note
    - Adjustment of the valid (mandatory) arguments
      by directly manipulating the static member argList::validArgs.
    - Adjustment of the valid options
      by directly manipulating the static member argList::validOptions.

SourceFiles
    argList.C

\*---------------------------------------------------------------------------*/

#ifndef argList_H
#define argList_H

#include "stringList.H"
#include "SubList.H"
#include "SLList.H"
#include "HashTable.H"
#include "word.H"
#include "fileName.H"
#include "parRun.H"
#include "IStringStream.H"

#include "sigFpe.H"
#include "sigInt.H"
#include "sigQuit.H"
#include "sigSegv.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class argList
{
    // Private data
        static bool bannerEnabled;

        stringList args_;
        HashTable<string> options_;

        word executable_;
        fileName rootPath_;
        fileName globalCase_;
        fileName case_;

        ParRunControl parRunControl_;

        // Signal handlers
        sigFpe sigFpe_;
        sigInt sigInt_;
        sigQuit sigQuit_;
        sigSegv sigSegv_;


    // Private member functions

        void getRootCase();

        //- Transcribe argv into internal args_
        //  return true if any "(" ... ")" sequences were captured
        bool regroupArgv(int& argc, char**& argv);


public:

    // Static data members

        //- A list of valid (mandatory) arguments
        static SLList<string> validArgs;

        //- A list of valid options
        static HashTable<string> validOptions;

        //- A list of valid parallel options
        static HashTable<string> validParOptions;

        //! @cond ignoreDocumentation
        class initValidTables
        {
        public:

            initValidTables();
        };
        //! @endcond ignoreDocumentation


    // Constructors

        //- Construct from argc and argv
        //  checking the arguments and options as requested
        argList
        (
            int& argc,
            char**& argv,
            bool checkArgs=true,
            bool checkOpts=true
        );


    // Destructor

        virtual ~argList();


    // Member functions

        // Access

            //- Name of executable
            const word& executable() const
            {
                return executable_;
            }

            //- Return root path
            const fileName& rootPath() const
            {
                return rootPath_;
            }

            //- Return case name
            const fileName& globalCaseName() const
            {
                return globalCase_;
            }

            //- Return case name (parallel run) or global case (serial run)
            const fileName& caseName() const
            {
                return case_;
            }

            //- Return the path
            fileName path() const
            {
                return rootPath()/caseName();
            }

            //- Return arguments
            const stringList& args() const
            {
                return args_;
            }

            //- Return additionl arguments,
            //  i.e. those additional to the executable itself
            stringList::subList additionalArgs() const;

            //- Return options
            const Foam::HashTable<string>& options() const
            {
                return options_;
            }

            //- Return the argument string associated with the named option
            const string& option(const word& opt) const
            {
                return options_.operator[](opt);
            }

            //- Return true if the named option is found
            bool optionFound(const word& opt) const
            {
                return options_.found(opt);
            }

            //- Return an IStringStream to the named option
            IStringStream optionLookup(const word& opt) const
            {
                return IStringStream(option(opt));
            }

            //- Read a value from the named option
            template<class T>
            T optionRead(const word& opt) const
            {
                T val;
                optionLookup(opt)() >> val;
                return val;
            }

            //- Read a value from the named option if present.
            //  Return true if the named option was found.
            template<class T>
            bool optionReadIfPresent(const word& opt, T& val) const
            {
                if (optionFound(opt))
                {
                    optionLookup(opt)() >> val;
                    return true;
                }
                else
                {
                    return false;
                }
            }

            //- Read a List of values from the named option
            template<class T>
            List<T> optionReadList(const word& opt) const
            {
                return readList<T>(optionLookup(opt)());
            }


        // Edit

            //- Disable emitting the banner information
            static void noBanner();

            //- Remove the parallel options
            static void noParallel();


        // Print

            //- Print usage
            void printUsage() const;

            //- Display documentation in browser
            //  Optionally display the application source code
            void displayDoc(bool source=false) const;


        // Check

            //- Check argument list
            bool check(bool checkArgs=true, bool checkOpts=true) const;

            //- Check root path and case path
            bool checkRootCase() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
