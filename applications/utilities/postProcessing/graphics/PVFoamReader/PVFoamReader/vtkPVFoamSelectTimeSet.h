/*=========================================================================

  Program:   ParaView
  Module:    $RCSfile: vtkPVFoamSelectTimeSet.h,v $

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkPVFoamSelectTimeSet - Special time selection widget used by PVFoamReaderModule
// .SECTION Description
// This is a PVWidget specially designed to be used with PVFoamReaderModule.
// It provides support for multiple sets. The time value selected by
// the user is passed to the Foam reader with a SetTimeValue() call.

#ifndef __vtkPVFoamSelectTimeSet_h
#define __vtkPVFoamSelectTimeSet_h

#include "vtkPVWidget.h"

class vtkKWLabel;
class vtkKWMenu;
class vtkKWLabeledFrame;
class vtkDataArrayCollection;
class vtkPVScalarListWidgetProperty;

class VTK_EXPORT vtkPVFoamSelectTimeSet : public vtkPVWidget
{
public:
  static vtkPVFoamSelectTimeSet* New();
  vtkTypeRevisionMacro(vtkPVFoamSelectTimeSet, vtkPVWidget);
  void PrintSelf(ostream& os, vtkIndent indent);

  virtual void Create(vtkKWApplication *pvApp);

  //BTX
  // Description:
  // Called when accept button is pushed.
  // Sets objects variable to the widgets value.
  // Adds a trace entry.  Side effect is to turn modified flag off.
  virtual void AcceptInternal(vtkClientServerID);
  //ETX

  // Description:
  // Called when the reset button is pushed.
  // Sets widget's value to the object-variable's value.
  // Side effect is to turn the modified flag off.
  virtual void ResetInternal();

  // Description:
  // Adds a script to the menu of the animation interface.
  virtual void AddAnimationScriptsToMenu(vtkKWMenu *menu, 
                                         vtkPVAnimationInterfaceEntry *ai);

  // Description:
  // Called whenthe animation method menu item is selected.
  // Needed for proper tracing.
  // It would be nice if the menu and cascade menus would trace
  // invokation of items (?relying of enumeration of menu items or label?)
  void AnimationMenuCallback(vtkPVAnimationInterfaceEntry *ai);

  // Description:
  // This is the labeled frame around the timeset tree.
  vtkGetObjectMacro(LabeledFrame, vtkKWLabeledFrame);

  // Description:
  // Label displayed on the labeled frame.
  void SetLabel(const char* label);
  const char* GetLabel();

  // Description:
  // Updates the time value label and the time ivar.
  void SetTimeValue(float time);
  vtkGetMacro(TimeValue, float);

  // Description:
  // Calls this->SetTimeValue () and Reader->SetTimeValue()
  // with currently selected time value.
  void SetTimeValueCallback(const char* item);

//BTX
  // Description:
  // Creates and returns a copy of this widget. It will create
  // a new instance of the same type as the current object
  // using NewInstance() and then copy some necessary state 
  // parameters.
  vtkPVFoamSelectTimeSet* ClonePrototype(vtkPVSource* pvSource,
                                 vtkArrayMap<vtkPVWidget*, vtkPVWidget*>* map);
//ETX

  // Description:
  // This serves a dual purpose.  For tracing and for saving state.
  virtual void Trace(ofstream *file);

  // Description:
  // Set/get the property to use with this widget.
  virtual void SetProperty(vtkPVWidgetProperty *prop);
  virtual vtkPVWidgetProperty* GetProperty();
  
  // Description:
  // Create the right property for use with this widget.
  virtual vtkPVWidgetProperty* CreateAppropriateProperty();
  
  // Description:
  // Set/get the command to pass the value to VTK.
  vtkSetStringMacro(SetCommand);
  vtkGetStringMacro(SetCommand);
  
  // Description:
  // Save this widget to a file.
  virtual void SaveInBatchScript(ofstream *file);
 
protected:
  vtkPVFoamSelectTimeSet();
  ~vtkPVFoamSelectTimeSet();

  vtkPVFoamSelectTimeSet(const vtkPVFoamSelectTimeSet&); // Not implemented
  void operator=(const vtkPVFoamSelectTimeSet&); // Not implemented

  vtkPVScalarListWidgetProperty *Property;

  char *SetCommand;
  
  vtkSetStringMacro(FrameLabel);
  vtkGetStringMacro(FrameLabel);

  vtkKWWidget* Tree;
  vtkKWWidget* TreeFrame;
  vtkKWLabel* TimeLabel;
  vtkKWLabeledFrame* LabeledFrame;

  void AddRootNode(const char* name, const char* text);
  void AddChildNode(const char* parent, const char* name, 
                    const char* text, const char* data);

  float TimeValue;
  char* FrameLabel;

  vtkDataArrayCollection* TimeSets;
  vtkClientServerID ServerSideID;

  // Fill the TimeSets collection with that from the actual reader.
  void SetTimeSetsFromReader();

//BTX
  virtual void CopyProperties(vtkPVWidget* clone, vtkPVSource* pvSource,
                              vtkArrayMap<vtkPVWidget*, vtkPVWidget*>* map);
//ETX
  
  int ReadXMLAttributes(vtkPVXMLElement* element,
                        vtkPVXMLPackageParser* parser);

  // Description:
  // An interface for saving a widget into a script.
  virtual void SaveInBatchScriptForPart(ofstream *file, vtkClientServerID);

};

#endif
