/*=========================================================================

  Program:   ParaView
  Module:    $RCSfile: vtkPVFoamServerSelectTimeSet.h,v $

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkPVFoamServerSelectTimeSet - Server-side helper for vtkPVFoamSelectTimeSet.
// .SECTION Description

#ifndef __vtkPVFoamServerSelectTimeSet_h
#define __vtkPVFoamServerSelectTimeSet_h

#include "vtkPVServerObject.h"

class vtkClientServerStream;
class vtkPVFoamServerSelectTimeSetInternals;
class vtkFoamReader;

class VTK_EXPORT vtkPVFoamServerSelectTimeSet : public vtkPVServerObject
{
public:
  static vtkPVFoamServerSelectTimeSet* New();
  vtkTypeRevisionMacro(vtkPVFoamServerSelectTimeSet, vtkPVServerObject);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Get a list the time sets provided by the given reader.
  const vtkClientServerStream& GetTimeSets(vtkFoamReader*);

protected:
  vtkPVFoamServerSelectTimeSet();
  ~vtkPVFoamServerSelectTimeSet();

  // Internal implementation details.
  vtkPVFoamServerSelectTimeSetInternals* Internal;
private:
  vtkPVFoamServerSelectTimeSet(const vtkPVFoamServerSelectTimeSet&); // Not implemented
  void operator=(const vtkPVFoamServerSelectTimeSet&); // Not implemented
};

#endif
