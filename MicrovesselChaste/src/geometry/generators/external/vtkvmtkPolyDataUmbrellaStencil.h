/*=========================================================================

  Program:   VMTK
  Module:    $RCSfile: vtkvmtkPolyDataUmbrellaStencil.h,v $
  Language:  C++
  Date:      $Date: 2006/04/06 16:46:44 $
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
// .NAME vtkvmtkPolyDataUmbrellaStencil - ..
// .SECTION Description
// ..

#ifndef __vtkvmtkPolyDataUmbrellaStencil_h
#define __vtkvmtkPolyDataUmbrellaStencil_h

// JG: Add code to prevent deprecation warnings - June '17
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the vtk deprecated warning
#include "vtkObject.h"
#include "vtkvmtkConstants.h"
#include "vtkvmtkPolyDataManifoldStencil.h"
#include "vtkvmtkWin32Header.h"

class VTK_VMTK_DIFFERENTIAL_GEOMETRY_EXPORT vtkvmtkPolyDataUmbrellaStencil : public vtkvmtkPolyDataManifoldStencil
{
public:

  static vtkvmtkPolyDataUmbrellaStencil *New();
  vtkTypeMacro(vtkvmtkPolyDataUmbrellaStencil,vtkvmtkPolyDataManifoldStencil);

  virtual vtkIdType GetItemType() {return VTK_VMTK_UMBRELLA_STENCIL;};

  void Build();

protected:
  vtkvmtkPolyDataUmbrellaStencil();
  ~vtkvmtkPolyDataUmbrellaStencil() {};

  void ScaleWithArea() {};

private:
  vtkvmtkPolyDataUmbrellaStencil(const vtkvmtkPolyDataUmbrellaStencil&);  // Not implemented.
  void operator=(const vtkvmtkPolyDataUmbrellaStencil&);  // Not implemented.
};

#endif

