/*=========================================================================

  Program:   ParaView
  Module:    $RCSfile: vtkPVFoamServerSelectTimeSet.cxx,v $

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkPVFoamServerSelectTimeSet.h"

#include "vtkClientServerInterpreter.h"
#include "vtkObjectFactory.h"
#include "vtkPVProcessModule.h"
#include "vtkFoamReader.h"
#include "vtkDataArrayCollection.h"
#include "vtkDataArrayCollectionIterator.h"
#include "vtkClientServerStream.h"

#include <vtkstd/string>

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkPVFoamServerSelectTimeSet);
vtkCxxRevisionMacro(vtkPVFoamServerSelectTimeSet, "$Revision: 1.4 $");

//----------------------------------------------------------------------------
class vtkPVFoamServerSelectTimeSetInternals
{
public:
  vtkClientServerStream Result;
};

//----------------------------------------------------------------------------
vtkPVFoamServerSelectTimeSet::vtkPVFoamServerSelectTimeSet()
{
  this->Internal = new vtkPVFoamServerSelectTimeSetInternals;
}

//----------------------------------------------------------------------------
vtkPVFoamServerSelectTimeSet::~vtkPVFoamServerSelectTimeSet()
{
  delete this->Internal;
}

//----------------------------------------------------------------------------
void vtkPVFoamServerSelectTimeSet::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
const vtkClientServerStream&
vtkPVFoamServerSelectTimeSet::GetTimeSets(vtkFoamReader* reader)
{
  // Reset the stream for a new list of time sets.
  this->Internal->Result.Reset();

  // Get the time sets from the reader.
  vtkDataArrayCollection* timeSets = reader->GetTimeSets();

  // Iterate through the time sets.
  vtkDataArrayCollectionIterator* iter = vtkDataArrayCollectionIterator::New();
  iter->SetCollection(timeSets);
  for(iter->GoToFirstItem(); !iter->IsDoneWithTraversal();
      iter->GoToNextItem())
    {
    // Each time set is stored in one message.
    this->Internal->Result << vtkClientServerStream::Reply;
    vtkDataArray* da = iter->GetDataArray();
    for(int i=0; i < da->GetNumberOfTuples(); ++i)
      {
      this->Internal->Result << da->GetTuple1(i);
      }
    this->Internal->Result << vtkClientServerStream::End;
    }
  iter->Delete();

  // Return the stream.
  return this->Internal->Result;
}
