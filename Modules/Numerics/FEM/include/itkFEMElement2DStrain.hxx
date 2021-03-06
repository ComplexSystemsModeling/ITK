/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#ifndef __itkFEMElement2DStrain_hxx
#define __itkFEMElement2DStrain_hxx

#include "itkFEMElement2DStrain.h"

namespace itk
{
namespace fem
{
template <class TBaseClass>
Element2DStrain<TBaseClass>
::Element2DStrain() : Superclass(), m_mat(0)
{
}

// ////////////////////////////////////////////////////////////////////////
/**
 * Methods related to the physics of the problem.
 */

template <class TBaseClass>
void
Element2DStrain<TBaseClass>
::GetStrainDisplacementMatrix(MatrixType & B, const MatrixType & shapeDgl) const
{
  unsigned int p;
  unsigned int Nn = this->GetNumberOfNodes();

  B.set_size(3, 2 * Nn);
  // Copy the shape function derivatives to the B matrix.
  for( unsigned int i = 0; i < Nn; i++ )
    {
    // Compute B index
    p = i << 1;

    // Compute B elements
    B[0][p]   = shapeDgl[0][i];
    B[0][p + 1] = 0;
    B[1][p]   = 0;
    B[1][p + 1] = shapeDgl[1][i];
    B[2][p]   = shapeDgl[1][i];
    B[2][p + 1] = shapeDgl[0][i];
    }
}

template <class TBaseClass>
void
Element2DStrain<TBaseClass>
::GetMassMatrix(MatrixType & Me) const
{
  // Call the parent's get matrix function
  Superclass::GetMassMatrix(Me);

  // Since parent class doesn't have the material properties,
  // we need to adjust Me matrix here for the density of the element.
  Me = Me * m_mat->GetDensityHeatProduct();
}

template <class TBaseClass>
void
Element2DStrain<TBaseClass>
::GetMaterialMatrix(MatrixType & D) const
{
  D.set_size(3, 3);

  /* Material properties matrix */
  Float fac = ( m_mat->GetThickness() * m_mat->GetYoungsModulus() )
    / ( ( 1 + m_mat->GetPoissonsRatio() ) * ( 1 - 2 * m_mat->GetPoissonsRatio() ) );
  D[0][0] = 1 - m_mat->GetPoissonsRatio();
  D[0][1] = m_mat->GetPoissonsRatio();
  D[0][2] = 0.0;

  D[1][0] = D[0][1];
  D[1][1] = D[0][0];
  D[1][2] = 0.0;

  D[2][0] = 0.0;
  D[2][1] = 0.0;
  D[2][2] = ( 1. - 2. * m_mat->GetPoissonsRatio() ) / 2.;

  D = D * fac;
}

template <class TBaseClass>
void
Element2DStrain<TBaseClass>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Materials: " << this->m_mat << std::endl;
}

}
}  // end namespace itk::fem

#endif // #ifndef __itkFEMElement2DStrain_hxx
