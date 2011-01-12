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
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkQuadEdgeMesh.h"
#include "itkCellInterface.h"
#include "itkPolylineCell.h"

#include <stdlib.h>


template < class TMesh >
int myTest()
{
  typedef TMesh                          MeshType;

  typedef typename MeshType::CellType             CellType;
  typedef typename CellType::CellAutoPointer      CellAutoPointer;

  typedef typename MeshType::PointType            PointType;
  typedef typename MeshType::PointIdentifier      PointIdentifier;
  typedef typename MeshType::CellIdentifier       CellIdentifier;
  typedef std::vector< PointIdentifier > PointIdList;

  typedef itk::PolylineCell< CellType >  PolylineCellType;

  // create the mesh
  typename MeshType::Pointer mesh = MeshType::New();

  // add the points to the mesh
  const int NumPoints = 7;

  PointType points[NumPoints];
  points[0][0] = 0.0;
  points[0][1] = 1.0;
  points[0][2] = 0.0;
  for( int i = 1; i < NumPoints; i++ )
    {
    points[i][0] = points[i-1][0] + 1.0;
    if( points[i-1][1] == 1.0 ) points[i][1] = 0.0;
    else points[i][1] = 1.0;
    points[i][2] = 0.0;
    }

  for(int i = 0; i < NumPoints; i++)
    {
    mesh->SetPoint( i, points[i] );
    }

  // create the polyline and add it to the mesh
  CellAutoPointer testCell;
  testCell.TakeOwnership( new PolylineCellType );

  for(int i = 0; i < NumPoints; i++)
    {
    testCell->SetPointId( i, i );
    }
  mesh->SetCell(0, testCell );

  if( mesh->GetNumberOfCells() != 1 )
    {
    std::cout << "Failed - The cell was not added to the container" << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;

}




int itkQuadEdgeMeshPolylineTest( int argc, char* argv[] )
{

  if( argc != 1 )
    {
    std::cout << "Does not require any argument." << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::Mesh< double, 3 >         MeshType1;
  typedef itk::QuadEdgeMesh< double, 3 > MeshType2;

  if( myTest< MeshType1 >() ) return EXIT_FAILURE;;
  if( myTest< MeshType2 >() ) return EXIT_FAILURE;;

  return EXIT_SUCCESS;

}
