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

\*---------------------------------------------------------------------------*/

#include "readerDatabase.H"
#include "demandDrivenData.H"
#include "fvMesh.H"
#include "fvMeshSubset.H"
#include "Time.H"
#include "fileName.H"
#include "instant.H"
#include "cellSet.H"
#include "cellModeller.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const bool Foam::readerDatabase::debug_ = Foam::env("readerDatabase");


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Gets cell numbers of all polyHedra
void Foam::readerDatabase::getPolyHedra()
{
    const cellModel& tet = *(cellModeller::lookup("tet"));
    const cellModel& pyr = *(cellModeller::lookup("pyr"));
    const cellModel& prism = *(cellModeller::lookup("prism"));
    const cellModel& wedge = *(cellModeller::lookup("wedge"));
    const cellModel& tetWedge = *(cellModeller::lookup("tetWedge"));
    const cellModel& hex = *(cellModeller::lookup("hex"));

    DynamicList<label> polys(mesh().nCells()/100 + 1);

    const cellShapeList& cellShapes = mesh().cellShapes();

    forAll(cellShapes, celli)
    {
        const cellShape& cellShape = cellShapes[celli];
        const cellModel& cellModel = cellShape.model();

        if
        (
            (cellModel != tet)
         && (cellModel != pyr)
         && (cellModel != prism)
         && (cellModel != wedge)
         && (cellModel != tetWedge)
         && (cellModel != hex)
        )
        {
            polys.append(celli);
        }
    }

    Info<< "Found " << polys.size() << " polyhedral cells " << endl;
    polys_.transfer(polys);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
Foam::readerDatabase::readerDatabase()
:
    fieldviewNames_(10),
    runTimePtr_(NULL),
    meshPtr_(NULL),
    setName_(""),
    polys_(),
    volScalarNames_(),
    volVectorNames_()
{
    // Initialize name mapping table. See note on static in header file.
    fieldviewNames_.insert("alpha", "aalpha");
    fieldviewNames_.insert("Alpha", "AAlpha");
    fieldviewNames_.insert("fsmach", "ffsmach");
    fieldviewNames_.insert("FSMach", "FFSMach");
    fieldviewNames_.insert("re", "rre");
    fieldviewNames_.insert("Re", "RRe");
    fieldviewNames_.insert("time", "ttime");
    fieldviewNames_.insert("Time", "TTime");
    fieldviewNames_.insert("pi", "ppi");
    fieldviewNames_.insert("PI", "PPI");
    fieldviewNames_.insert("x", "xx");
    fieldviewNames_.insert("X", "XX");
    fieldviewNames_.insert("y", "yy");
    fieldviewNames_.insert("Y", "YY");
    fieldviewNames_.insert("z", "zz");
    fieldviewNames_.insert("Z", "ZZ");
    fieldviewNames_.insert("rcyl", "rrcyl");
    fieldviewNames_.insert("Rcyl", "RRcyl");
    fieldviewNames_.insert("theta", "ttheta");
    fieldviewNames_.insert("Theta", "TTheta");
    fieldviewNames_.insert("rsphere", "rrsphere");
    fieldviewNames_.insert("Rsphere", "RRsphere");
    fieldviewNames_.insert("k", "kk");
    fieldviewNames_.insert("K", "KK");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::readerDatabase::~readerDatabase()
{
    deleteDemandDrivenData(meshPtr_);
    deleteDemandDrivenData(runTimePtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::Time& Foam::readerDatabase::runTime() const
{
    if (!runTimePtr_)
    {
        FatalErrorIn("readerDatabase::runTime()")
            << "No database set" << abort(FatalError);
    }
    return *runTimePtr_;
}


const Foam::fvMesh& Foam::readerDatabase::mesh() const
{
    if (!meshPtr_)
    {
        FatalErrorIn("readerDatabase::runTime()")
            << "No mesh set" << abort(FatalError);
    }

    if (setName_.empty())
    {
        return *meshPtr_;
    }
    else
    {
        return meshPtr_->subMesh();
    }
}


const Foam::labelList& Foam::readerDatabase::polys() const
{
    return polys_;
}


const Foam::wordList& Foam::readerDatabase::volScalarNames() const
{
    return volScalarNames_;
}


const Foam::wordList& Foam::readerDatabase::volVectorNames() const
{
    return volVectorNames_;
}


const Foam::word& Foam::readerDatabase::getFvName(const word& foamName) const
{
    if (fieldviewNames_.found(foamName))
    {
        return fieldviewNames_[foamName];
    }
    else
    {
        return foamName;
    }
}


bool Foam::readerDatabase::setRunTime
(
    const fileName& rootDir,
    const fileName& caseName,
    const word& setName
)
{
    bool newDatabase = false;

    if (runTimePtr_)
    {
        if
        (
            (runTimePtr_->caseName() != caseName)
         || (runTimePtr_->rootPath() != rootDir)
         || (setName_ != setName)
        )
        {
            if (debug_)
            {
                Info<< "Deleting old mesh since deleting old database" << endl;
            }

            deleteDemandDrivenData(meshPtr_);

            if (debug_)
            {
                Info<< "Deleting old database for " << runTimePtr_->caseName()
                    << endl;
            }

            deleteDemandDrivenData(runTimePtr_);
        }
    }

    setName_ = setName;

    if (!runTimePtr_)
    {
        if (debug_)
        {
            Info<< "Deleting old mesh since loading new Time" << endl;
        }

        deleteDemandDrivenData(meshPtr_);

        if (debug_)
        {
            Info<< "Creating database for " << caseName << endl;
        }

        runTimePtr_ = new Time(Time::controlDictName, rootDir, caseName);

        newDatabase = true;
    }

    return newDatabase;
}


void Foam::readerDatabase::loadMesh()
{
    deleteDemandDrivenData(meshPtr_);

    Info<< "Loading new mesh" << endl;

    meshPtr_ = new fvMeshSubset
    (
        *runTimePtr_,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (setName_.size())
    {
        Info<< "Subsetting mesh based on cellSet " << setName_ << endl;

        fvMeshSubset& mesh = *meshPtr_;

        cellSet currentSet(mesh, setName_);

        mesh.setCellSubset(currentSet);
    }
    getPolyHedra();
}


Foam::polyMesh::readUpdateState Foam::readerDatabase::setTime
(
    const instant& timeInstance,
    const label timeIndex
)
{
    runTime().setTime(timeInstance, timeIndex);

    polyMesh::readUpdateState meshChange;

    if (meshPtr_)
    {
        // Update loaded mesh
        meshChange = meshPtr_->readUpdate();

        if (setName_.size() && meshChange != polyMesh::UNCHANGED)
        {
            Info<< "Subsetting mesh based on " << setName_ << endl;

            fvMeshSubset& mesh = *meshPtr_;

            cellSet currentSet(mesh, setName_);

            mesh.setCellSubset(currentSet);
        }

        if
        (
            (meshChange == polyMesh::TOPO_CHANGE)
         || (meshChange == polyMesh::TOPO_PATCH_CHANGE)
        )
        {
            getPolyHedra();
        }
    }
    else
    {
        // Force new mesh to be loaded for current time
        loadMesh();
        meshChange = polyMesh::TOPO_CHANGE;
    }

    return meshChange;
}


void Foam::readerDatabase::setFieldNames
(
   const wordList& vsNames,
   const wordList& vvNames
)
{
    volScalarNames_ = vsNames;
    volVectorNames_ = vvNames;
}




// ************************************************************************* //
