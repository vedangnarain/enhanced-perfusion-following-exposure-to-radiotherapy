/*=========================================================================

  Program:   VMTK
  Module:    $RCSfile: vtkvmtkDataSetItem.h,v $
  Language:  C++
  Date:      $Date: 2006/04/06 16:46:43 $
  Version:   $Revision: 1.4 $

  Copyright (c) Luca Antiga, David Steinman. All rights reserved.
  See LICENCE file for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm 
  for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
// .NAME vtkvmtkDataSetItem - ..
// .SECTION Description
// ..

#ifndef __vtkvmtkDataSetItem_h
#define __vtkvmtkDataSetItem_h

// JG: Add code to prevent deprecation warnings - June '17
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the vtk deprecated warning
#include "vtkObject.h"
#include "vtkvmtkItem.h"
#include "vtkDataSet.h"
//#include "vtkvmtkDifferentialGeometryWin32Header.h"
#include "vtkvmtkWin32Header.h"

class VTK_VMTK_DIFFERENTIAL_GEOMETRY_EXPORT vtkvmtkDataSetItem : public vtkvmtkItem 
{
public:

  vtkTypeMacro(vtkvmtkDataSetItem,vtkvmtkItem);

/*   vtkSetObjectMacro(DataSet,vtkDataSet); */
/*   vtkGetObjectMacro(DataSet,vtkDataSet); */
  void SetDataSet(vtkDataSet* dataSet) {this->DataSet = dataSet;};
  vtkDataSet* GetDataSet() {return this->DataSet;};

  vtkSetMacro(DataSetPointId,vtkIdType);
  vtkGetMacro(DataSetPointId,vtkIdType);

  // Description:
  // Build the item.
  virtual void Build() = 0;

  // Description:
  // Standard DeepCopy method.
  virtual void DeepCopy(vtkvmtkItem *src);

  vtkSetMacro(ReallocateOnBuild,int)
  vtkGetMacro(ReallocateOnBuild,int)
  vtkBooleanMacro(ReallocateOnBuild,int)

protected:
  vtkvmtkDataSetItem();
  ~vtkvmtkDataSetItem() {};

  vtkDataSet *DataSet;
  vtkIdType DataSetPointId;

  int ReallocateOnBuild;

private:
  vtkvmtkDataSetItem(const vtkvmtkDataSetItem&);  // Not implemented.
  void operator=(const vtkvmtkDataSetItem&);  // Not implemented.
};

#endif

