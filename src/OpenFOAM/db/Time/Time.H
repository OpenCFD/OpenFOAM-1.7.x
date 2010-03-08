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

Class
    Foam::Time

Description
    Class to control time during OpenFOAM simulations that is also the
    top-level objectRegistry.

SourceFiles
    Time.C
    TimeIO.C
    findInstance.C

\*---------------------------------------------------------------------------*/

#ifndef Time_H
#define Time_H

#include "TimePaths.H"
#include "objectRegistry.H"
#include "IOdictionary.H"
#include "FIFOStack.H"
#include "clock.H"
#include "cpuTime.H"
#include "TimeState.H"
#include "Switch.H"
#include "instantList.H"
#include "NamedEnum.H"
#include "typeInfo.H"
#include "dlLibraryTable.H"
#include "functionObjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class Time Declaration
\*---------------------------------------------------------------------------*/

class Time
:
    public clock,
    public cpuTime,
    public TimePaths,
    public objectRegistry,
    public TimeState
{
    // Private data

        //- The controlDict
        IOdictionary controlDict_;

public:

        //- Write control options
        enum writeControls
        {
            wcTimeStep,
            wcRunTime,
            wcAdjustableRunTime,
            wcClockTime,
            wcCpuTime
        };

        //- Stop-run control options
        enum stopAtControls
        {
            saEndTime,
            saNoWriteNow,
            saWriteNow,
            saNextWrite
        };

        //- Suported time directory name formats
        enum fmtflags
        {
            general = 0,
            fixed = ios_base::fixed,
            scientific = ios_base::scientific
        };


protected:

    // Protected data

        label  startTimeIndex_;
        scalar startTime_;
        scalar endTime_;

        static const NamedEnum<stopAtControls, 4> stopAtControlNames_;
        stopAtControls stopAt_;

        static const NamedEnum<writeControls, 5> writeControlNames_;
        writeControls writeControl_;

        scalar writeInterval_;

        label  purgeWrite_;
        mutable FIFOStack<word> previousOutputTimes_;

        //- Is the time currently being sub-cycled?
        bool subCycling_;

        //- If time is being sub-cycled this is the previous TimeState
        autoPtr<TimeState> prevTimeState_;

        //- Time directory name format
        static fmtflags format_;

        //- Time directory name precision
        static int precision_;

        //- Adjust the time step so that writing occurs at the specified time
        void adjustDeltaT();

        //- Set the controls from the current controlDict
        void setControls();

        //- Read the control dictionary and set the write controls etc.
        virtual void readDict();


private:

        //- Default write option
        IOstream::streamFormat writeFormat_;

        //- Default output file format version number
        IOstream::versionNumber writeVersion_;

        //- Default output compression
        IOstream::compressionType writeCompression_;

        //- Default graph format
        word graphFormat_;

        //- Is runtime modification of dictionaries allowed?
        Switch runTimeModifiable_;

        //- Instantiate a dummy class to cause the reading of dynamic libraries
        dlLibraryTable::readDlLibrary readLibs_;

        //- Function objects executed at start and on ++, +=
        mutable functionObjectList functionObjects_;


public:

    TypeName("time");

    //- The default control dictionary name (normally "controlDict")
    static word controlDictName;


    // Constructors

        //- Construct given name, rootPath and casePath
        Time
        (
            const word& name,
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName = "system",
            const word& constantName = "constant"
        );

        //- Construct given dictionary, rootPath and casePath
        Time
        (
            const dictionary& dict,
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName = "system",
            const word& constantName = "constant"
        );

        //- Construct given endTime, rootPath and casePath
        Time
        (
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName = "system",
            const word& constantName = "constant"
        );


    // Destructor

        //- Virtual destructor
        virtual ~Time();


    // Member functions

        // Database functions

            //- Return root path
            const fileName& rootPath() const
            {
                return TimePaths::rootPath();
            }

            //- Return case name
            const fileName& caseName() const
            {
                return TimePaths::caseName();
            }

            //- Return path
            fileName path() const
            {
                return rootPath()/caseName();
            }

            const dictionary& controlDict() const
            {
                return controlDict_;
            }

            virtual const fileName& dbDir() const
            {
                return fileName::null;
            }

            //- Return current time path
            fileName timePath() const
            {
                return path()/timeName();
            }

            //- Default write format
            IOstream::streamFormat writeFormat() const
            {
                return writeFormat_;
            }

            //- Default write version number
            IOstream::versionNumber writeVersion() const
            {
                return writeVersion_;
            }

            //- Default write compression
            IOstream::compressionType writeCompression() const
            {
                return writeCompression_;
            }

            //- Default graph format
            const word& graphFormat() const
            {
                return graphFormat_;
            }

            //- Read control dictionary, update controls and time
            virtual bool read();

            //- Read the objects that have been modified
            void readModifiedObjects();

            //- Return the location of "dir" containing the file "name".
            //  (eg, used in reading mesh data)
            //  If name is null, search for the directory "dir" only
            word findInstance
            (
                const fileName& dir,
                const word& name = word::null,
                const IOobject::readOption rOpt = IOobject::MUST_READ
            ) const;

            //- Search the case for valid time directories
            instantList times() const;

            //- Search the case for the time directory path
            //  corresponding to the given instance
            word findInstancePath(const instant&) const;

            //- Search the case for the time closest to the given time
            instant findClosestTime(const scalar) const;

            //- Search instantList for the time index closest to the given time
            static label findClosestTimeIndex(const instantList&, const scalar);

            //- Write using given format, version and compression
            virtual bool writeObject
            (
                IOstream::streamFormat,
                IOstream::versionNumber,
                IOstream::compressionType
            ) const;

            //- Write the objects now and continue the run
            bool writeNow();

            //- Write the objects now and end the run
            bool writeAndEnd();


        // Access

            //- Return time name of given scalar time
            static word timeName(const scalar);

            //- Return current time name
            virtual word timeName() const;

            //- Search a given directory for valid time directories
            static instantList findTimes(const fileName&);

            //- Return start time index
            virtual label startTimeIndex() const;

            //- Return start time
            virtual dimensionedScalar startTime() const;

            //- Return end time
            virtual dimensionedScalar endTime() const;

            //- Return the list of function objects
            const functionObjectList& functionObjects() const
            {
                return functionObjects_;
            }

            //- Return true if time currently being sub-cycled, otherwise false
            bool subCycling() const
            {
                return subCycling_;
            }

            //- Return previous TimeState if time is being sub-cycled
            const TimeState& prevTimeState() const
            {
                return prevTimeState_();
            }


        // Check

            //- Return true if run should continue,
            //  also invokes the functionObjectList::end() method
            //  when the time goes out of range
            //  @note
            //  For correct behaviour, the following style of time-loop
            //  is recommended:
            //  @code
            //      while (runTime.run())
            //      {
            //          runTime++;
            //          solve;
            //          runTime.write();
            //      }
            //  @endcode
            virtual bool run() const;

            //- Return true if run should continue and if so increment time
            //  also invokes the functionObjectList::end() method
            //  when the time goes out of range
            //  @note
            //  For correct behaviour, the following style of time-loop
            //  is recommended:
            //  @code
            //      while (runTime.loop())
            //      {
            //          solve;
            //          runTime.write();
            //      }
            //  @endcode
            virtual bool loop();

            //- Return true if end of run,
            //  does not invoke any functionObject methods
            //  @note
            //      The rounding heuristics near endTime mean that
            //      @code run() @endcode and @code !end() @endcode may
            //      not yield the same result
            virtual bool end() const;


        // Edit

            //- Reset the time and time-index to those of the given time
            virtual void setTime(const Time&);

            //- Reset the time and time-index
            virtual void setTime(const instant&, const label newIndex);

            //- Reset the time and time-index
            virtual void setTime
            (
                const dimensionedScalar&,
                const label newIndex
            );

            //- Reset the time and time-index
            virtual void setTime(const scalar, const label newIndex);

            //- Reset end time
            virtual void setEndTime(const dimensionedScalar&);

            //- Reset end time
            virtual void setEndTime(const scalar);

            //- Reset time step
            virtual void setDeltaT(const dimensionedScalar&);

            //- Reset time step
            virtual void setDeltaT(const scalar);

            //- Set time to sub-cycle for the given number of steps
            virtual TimeState subCycle(const label nSubCycles);

            //- Reset time after sub-cycling back to previous TimeState
            virtual void endSubCycle();

            //- Return non-const access to the list of function objects
            functionObjectList& functionObjects()
            {
                return functionObjects_;
            }


    // Member operators

        //- Set deltaT to that specified and increment time via operator++()
        virtual Time& operator+=(const dimensionedScalar&);

        //- Set deltaT to that specified and increment time via operator++()
        virtual Time& operator+=(const scalar);

        //- Prefix increment,
        //  also invokes the functionObjectList::start() or
        //  functionObjectList::execute() method, depending on the time-index
        virtual Time& operator++();

        //- Postfix increment, this is identical to the prefix increment
        virtual Time& operator++(int);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
