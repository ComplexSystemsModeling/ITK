/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkShapedFloodFilledImageFunctionConditionalConstIterator.txx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkShapedFloodFilledImageFunctionConditionalConstIterator_txx
#define __itkShapedFloodFilledImageFunctionConditionalConstIterator_txx

#include "itkShapedFloodFilledImageFunctionConditionalConstIterator.h"

namespace itk
{
template< class TImage, class TFunction >
bool
ShapedFloodFilledImageFunctionConditionalConstIterator< TImage, TFunction >
::IsPixelIncluded(const IndexType & index) const
{
  return this->m_Function->EvaluateAtIndex(index);
}
} // end namespace itk

#endif
