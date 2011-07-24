#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshPoint.h"
#include "itkRegularSphereMeshSource.h"
#include "itkDefaultStaticMeshTraits.h"
#include "itkVTKPolyDataWriter.h"
#include "itkQuadEdgeMeshBoundaryEdgesMeshFunction.h"

//#include<conio.h>
#include<iostream>

int main( int argc, char ** argv )
{

  typedef itk::QuadEdgeMesh< float, 3>   MeshType;

  typedef MeshType::PointIdentifier      PointIdentifier;
  typedef MeshType::CellIdentifier       CellIdentifier;
  typedef MeshType::PointsContainer      PointsContainer;
  typedef MeshType::CellsContainer       CellsContainer;
  typedef MeshType::CellType             CellType;
  typedef MeshType::PointIdList          PointIdList;

  typedef CellType::PointIdConstIterator PointIdConstIterator;

  typedef PointsContainer::ConstIterator PointIterator;
  typedef CellsContainer::ConstIterator  CellIterator;

  typedef itk::RegularSphereMeshSource< MeshType >  SphereMeshSourceType;
  SphereMeshSourceType::Pointer  mySphereMeshSource = SphereMeshSourceType::New();

  typedef SphereMeshSourceType::PointType   PointType;
  typedef SphereMeshSourceType::VectorType  VectorType;

  typedef itk::QuadEdgeMeshBoundaryEdgesMeshFunction< MeshType > BoundaryLocatorType;

  PointType center;
  center.Fill( 0.0 );

  VectorType scale;
  scale.Fill( 1.0 );

  mySphereMeshSource->SetCenter( center );
  mySphereMeshSource->SetResolution( 3 );
  mySphereMeshSource->SetScale( scale );
  mySphereMeshSource->Modified();

 try
    {
    mySphereMeshSource->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Error during Sphere Mesh Update() " << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  MeshType::Pointer myMesh = mySphereMeshSource->GetOutput();
  MeshType::Pointer myDualMesh = MeshType::New();

  // take care of the primal cells => generate dual points
  const CellsContainer *cells = myMesh->GetCells();
  std::map< CellIdentifier ,PointIdentifier >   DualPointFromPrimalTriangleLUT;

  if( cells )
    {
    CellIterator cellIterator = cells->Begin();
    CellIterator cellEnd = cells->End();

    bool found = false;
    unsigned int numberOfDualPoints = 0;
    while( ( cellIterator != cellEnd ) && !found )
      {
      switch ( cellIterator.Value()->GetType() )
        {
        case 0: //VERTEX_CELL:
        case 1: //LINE_CELL:
        case 2: //TRIANGLE_CELL:
        case 3: //QUADRILATERAL_CELL:
          break;
        case 4: //POLYGON_CELL:
          if( cellIterator.Value()->GetNumberOfPoints() > 3 )
            {
            std::cout << "toto" << std::endl;
            found = true;
            }
          else
            {
            std::cout << "Found triangle." << std::endl;

            // compute dual point coordinate
            PointIdConstIterator current= cellIterator.Value()->PointIdsBegin();
            PointIdConstIterator end    = cellIterator.Value()->PointIdsEnd();
            PointType d_point;
            d_point[0] = 0.0;
            d_point[1] = 0.0;
            d_point[2] = 0.0;
            while( current != end )
              {
              PointType point = myMesh->GetPoint( *current );
              //std::cout << point << std::endl;
              for( unsigned int i =0; i < 3; i++ ) d_point[i] += point[i];
              current++;
              }
            for( unsigned int i =0; i < 3; i++ ) d_point[i] /= 3;

            // push dual point in dualPoints container
            myDualMesh->SetPoint( numberOfDualPoints, d_point );
            CellIdentifier cellIdentifier = cellIterator.Index();
            std::cout << "Triangle ID: " << cellIdentifier << std::endl;
            DualPointFromPrimalTriangleLUT[cellIdentifier] = numberOfDualPoints;
            numberOfDualPoints++;
            }
          break;
        case 7: //QUADRATIC_EDGE_CELL:
        case 8: //QUADRATIC_TRIANGLE_CELL:
          break;
        default:
          std::cerr << "Unhandled cell (volumic?)." << std::endl;
        }
      cellIterator++;
      }
    }

  // take care of the primal points => generate dual cells
  const PointsContainer *primalPoints = myMesh->GetPoints();
  const PointsContainer *dualPoints   = myDualMesh->GetPoints();

  if( primalPoints )
    {
    PointIterator pointIterator = primalPoints->Begin();
    PointIterator pointEnd      = primalPoints->End();

    while( pointIterator != pointEnd )
      {
      typedef MeshType::QEType QuadEdgeType;

      // grab the QEdge
      PointType point = pointIterator.Value();
      QuadEdgeType * start = point.GetEdge();
      QuadEdgeType * current = point.GetEdge();

      // create a point ID list to hold the dual point IDs while
      // we are iterating around a primal point to create the dual cell
      PointIdList pointidlist;

      // iterate around the o-ring
      do
        {
        // get the id of the face on the left
        QuadEdgeType::DualOriginRefType leftTriangle;
        leftTriangle = current->GetLeft();

        // push the dual point ID to the point ID list
        pointidlist.push_back( DualPointFromPrimalTriangleLUT[ leftTriangle ] );

        current = current->GetOnext();

        } while( current != start );

      // point list is complete, add the dual cell to the dual mesh;
      myDualMesh->AddFace( pointidlist );

      // next point
      pointIterator++;
      }
    }

  BoundaryLocatorType::Pointer boundaryEdges = BoundaryLocatorType::New();

  MeshType::EdgeListPointerType boundaryEdgesPointerList;
  boundaryEdgesPointerList = boundaryEdges->Evaluate( *myMesh );

  typedef itk::VTKPolyDataWriter<MeshType>   WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( myMesh );
  writer->SetFileName( "MyVeryFirstTriangularMesh.vtk" );
  writer->Write();

  writer->SetInput( myDualMesh );
  writer->SetFileName( "MyVeryFirstSimplexMesh.vtk" );
  writer->Write();
  // getch();
  return EXIT_SUCCESS;

}
