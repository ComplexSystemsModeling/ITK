#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshPoint.h"
#include "itkRegularSphereMeshSource.h"
#include "itkDefaultStaticMeshTraits.h"
#include "itkVTKPolyDataWriter.h"
#include "itkQuadEdgeMeshBoundaryEdgesMeshFunction.h"
#include "itkQuadEdgeMeshTopologyChecker.h"
#include "itkQuadEdgeMeshPolygonCell.h"

#include <iostream>

typedef itk::QuadEdgeMesh< float, 3>   MeshType;

template< class TMesh >
std::vector< typename TMesh::PointType > GeneratePointCoordinates( const unsigned int& iN );
template< class TMesh >
void CreateSquareQuadMesh( typename TMesh::Pointer mesh );
template< class TMesh >
void CreateSquareTriangularMesh( typename TMesh::Pointer mesh );
template< class TMesh >
void CreateTetraedronMesh( typename TMesh::Pointer mesh );
template< class TMesh >
void CreateSamosa( typename TMesh::Pointer mesh );

int main( int argc, char ** argv )
{

  const unsigned int dimension = 3;
  //  typedef itk::QuadEdgeMesh< float, 3>   MeshType;
  typedef MeshType::PointIdentifier      PointIdentifier;
  typedef MeshType::CellIdentifier       CellIdentifier;
  typedef MeshType::PointsContainer      PointsContainer;
  typedef MeshType::CellsContainer       CellsContainer;
  typedef MeshType::CellType             CellType;
  typedef MeshType::PointIdList          PointIdList;
  typedef MeshType::QEType               QuadEdgeType;

  typedef CellType::PointIdConstIterator PointIdConstIterator;

  typedef PointsContainer::ConstIterator PointIterator;
  typedef CellsContainer::ConstIterator  CellIterator;

  typedef itk::RegularSphereMeshSource< MeshType >  SphereMeshSourceType;
  SphereMeshSourceType::Pointer  mySphereMeshSource = SphereMeshSourceType::New();

  typedef SphereMeshSourceType::PointType   PointType;
  typedef SphereMeshSourceType::VectorType  VectorType;

  typedef itk::QuadEdgeMeshBoundaryEdgesMeshFunction< MeshType > BoundaryLocatorType;

  std::map< CellIdentifier ,PointIdentifier >   DualPointFromPrimalTriangleLUT;

  PointType center;
  center.Fill( 0.0 );

  VectorType scale;
  scale.Fill( 1.0 );

  mySphereMeshSource->SetCenter( center );
  mySphereMeshSource->SetResolution( dimension );
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

  MeshType::Pointer myDualMesh = MeshType::New();
  MeshType::Pointer myPrimalMesh = MeshType::New();
  //myPrimalMesh = MySphereMeshSource->GetOutput();
  CreateSquareTriangularMesh <MeshType> (myPrimalMesh);

   myPrimalMesh->LightWeightDeleteEdge( myPrimalMesh->FindEdge( 11, 12 ) );

  // take care of the primal cells => generate dual points
  const CellsContainer *primalCells = myPrimalMesh->GetCells();

  if( primalCells )
    {
    CellIterator cellIterator = primalCells->Begin();
    CellIterator cellEnd = primalCells->End();

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
            std::cout << "\n \n toto" << std::endl;
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
              PointType point = myPrimalMesh->GetPoint( *current );
              //std::cout << point << std::endl;
              for( unsigned int i =0; i < dimension; i++ ) d_point[i] += point[i];
              current++;
              }
            for( unsigned int i =0; i < dimension; i++ ) d_point[i] /= dimension;

            // push dual point in dualPoints container
            myDualMesh->SetPoint( numberOfDualPoints, d_point );
            CellIdentifier cellIdentifier = cellIterator.Index();
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
  const PointsContainer *primalPoints = myPrimalMesh->GetPoints();
  const PointsContainer *dualPoints   = myDualMesh->GetPoints();

  if( primalPoints )
    {
    PointIterator pointIterator = primalPoints->Begin();
    PointIterator pointEnd      = primalPoints->End();

    while( pointIterator != pointEnd )
      {
      // grab the QEdge
      PointType point = pointIterator.Value();
      QuadEdgeType * start = point.GetEdge();
      QuadEdgeType * current = point.GetEdge();

      // create a point ID list to hold the dual point IDs while
      // we are iterating around a primal point to create the dual cell
      PointIdList pointidlist;

     if( point.IsInternal() )
       {
       // iterate around the o-ring
       do
         {
         // get the id of the face on the left
         QuadEdgeType::DualOriginRefType leftTriangle = current->GetLeft();

         // push the dual point ID to the point ID list
         pointidlist.push_back( DualPointFromPrimalTriangleLUT[ leftTriangle ] );

         current = current->GetOnext();

         } while( current != start );

       // point list is complete, add the dual cell to the dual mesh;
       myDualMesh->AddFace( pointidlist );
       }
     // next point
     pointIterator++;
     }
   }

  BoundaryLocatorType::Pointer boundaryEdges = BoundaryLocatorType::New();
  MeshType::EdgeListPointerType boundaryEdgesPointerList = boundaryEdges->Evaluate( *myPrimalMesh );

  // for each boundary
  for( unsigned int i = 0; i < boundaryEdgesPointerList->size(); i++ )
    {
    // get the first edge (arbitrary) and remove it from the list
    QuadEdgeType* firstEdge = boundaryEdgesPointerList->front();
    QuadEdgeType* currentEdge = firstEdge;
    boundaryEdgesPointerList->pop_front();

    // circulate around the boundary and do what you have to do
    bool firstTime = true;
    PointIdentifier previousPointId;
    PointIdentifier currentPointId;
    PointIdentifier firstPointId = myDualMesh->GetNumberOfPoints();
    do
      {
     // HOMEWORK 1
     // create a new point in the middle of the edge
     PointIdentifier originPointId = currentEdge->GetOrigin();
     PointIdentifier destinationPointId = currentEdge->GetDestination();

     PointType originPoint = myPrimalMesh->GetPoint( originPointId );
     PointType destinationPoint = myPrimalMesh->GetPoint( destinationPointId );
     PointType currentPoint;
     for( unsigned int k =0; k<dimension; k++ )
       currentPoint[k] = (originPoint[k] + destinationPoint[k]) / 2 ;

    // add the new border point P1 to the dual point container
    currentPointId = myDualMesh->AddPoint(currentPoint);

    // find the dual point P2 associated with the face on the left
    QuadEdgeType::DualOriginRefType leftTriangle = currentEdge->GetRight();

    // add the dual edge P1-P2
    PointIdentifier leftDualPointId = DualPointFromPrimalTriangleLUT[ leftTriangle ];
    myDualMesh->AddEdge( currentPointId, leftDualPointId );

  // HOMEWORK 2

  // beware border cases
  // add the edge linking two new border dual points
  // either previous-current, or current-next
  // be carefull of first, respectingly last case
      if(firstTime == true)
 firstTime = false;
  else
    {
    // HOMEWORK 3

    // beware same border case
    // create and add cell to dual cell container
    myDualMesh->AddEdge( previousPointId, currentPointId );
    // grab the QEdge

    // create a point ID list to hold the dual point IDs while iterating to create the dual cell
    PointIdList pointidlist;
    pointidlist.push_back( previousPointId );
    QuadEdgeType *myEdge = currentEdge->GetOnext();
    do
      {
      QuadEdgeType::DualOriginRefType myleftTriangle = myEdge->GetLeft();
      PointIdentifier myleftDualPointId =  DualPointFromPrimalTriangleLUT[ myleftTriangle ];
      pointidlist.push_back( myleftDualPointId );
      myEdge = myEdge->GetOnext();
      } while( !myEdge->IsAtBorder() );

    pointidlist.push_back( currentPointId );

    // point list is complete, add the dual cell to the dual mesh;
    myDualMesh->AddFace( pointidlist );
    }

  previousPointId = currentPointId;

  // Update currentEdge
  currentEdge = currentEdge->GetLnext();

  } while (currentEdge != firstEdge);

PointIdList pointidlist;
pointidlist.push_back( previousPointId );
QuadEdgeType *myEdge = currentEdge->GetOnext();
do
  {
  QuadEdgeType::DualOriginRefType myleftTriangle = myEdge->GetLeft();
  PointIdentifier myleftDualPointId =  DualPointFromPrimalTriangleLUT[ myleftTriangle ];
 pointidlist.push_back( myleftDualPointId );
  myEdge = myEdge->GetOnext();
  } while( !myEdge->IsAtBorder() );

    pointidlist.push_back( firstPointId );

    // point list is complete, add the dual cell to the dual mesh;
    myDualMesh->AddFace( pointidlist );
    }

  typedef itk::VTKPolyDataWriter<MeshType>   WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( myPrimalMesh );
  writer->SetFileName( "E:\\TestSquareTriangularMesh.vtk" );
  writer->Write();

  writer->SetInput( myDualMesh );
  writer->SetFileName( "E:\\TestSquareTriangularSimplexMesh.vtk" );
  writer->Write();

  return EXIT_SUCCESS;

}

template< class TMesh >
std::vector< typename TMesh::PointType >
GeneratePointCoordinates( const unsigned int& iN )
{
  typedef typename TMesh::PointType        PointType;
  typedef typename PointType::CoordRepType CoordRepType;
  std::vector< PointType > oPt( iN * iN );

  for( unsigned int i = 0; i < iN; i++ )
  {
    for( unsigned int j = 0; j < iN; j++ )
    {
      oPt[ i * iN + j ][0] = static_cast< CoordRepType >( j );
      oPt[ i * iN + j ][1] = static_cast< CoordRepType >( i );
      oPt[ i * iN + j ][2] = static_cast< CoordRepType >( 0. );
    }
  }

  return oPt;
}

template< class TMesh >
void CreateSquareQuadMesh( typename TMesh::Pointer mesh )
{
  typedef TMesh                         MeshType;
  typedef typename MeshType::CellType   CellType;

  typedef itk::QuadEdgeMeshPolygonCell< CellType > QEPolygonCellType;

  if( mesh->GetNumberOfPoints( ) )
    {
    mesh->Clear( );
    mesh->ClearFreePointAndCellIndexesLists();
    }

  /////////////////////////////////////////////////////////////
  int expectedNumPts = 25;
  int expectedNumCells = 16;
  int simpleSquareCells[64] =
  {  0,  1,  6, 5,
     1,  2,  7, 6,
     2,  3,  8, 7,
     3,  4,  9, 8,
     5,  6, 11, 10,
     6,  7, 12, 11,
     7,  8, 13, 12,
     8,  9, 14, 13,
    10, 11, 16, 15,
    11, 12, 17, 16,
    12, 13, 18, 17,
    13, 14, 19, 18,
    15, 16, 21, 20,
    16, 17, 22, 21,
    17, 18, 23, 22,
    18, 19, 24, 23 };

  typedef typename MeshType::PointType PointType;

  std::vector< PointType > pts = GeneratePointCoordinates< MeshType >( 5 );

  for(int i=0; i<expectedNumPts; i++)
    {
    mesh->SetPoint( i, pts[i] );
    }

  typename CellType::CellAutoPointer cellpointer;
  QEPolygonCellType *poly;

  for(int i=0; i<expectedNumCells; i++)
    {
    poly = new QEPolygonCellType( 4 );
    cellpointer.TakeOwnership( poly );
    cellpointer->SetPointId( 0, simpleSquareCells[4*i] );
    cellpointer->SetPointId( 1, simpleSquareCells[4*i+1] );
    cellpointer->SetPointId( 2, simpleSquareCells[4*i+2] );
    cellpointer->SetPointId( 3, simpleSquareCells[4*i+3] );
    mesh->SetCell( i, cellpointer );
    }
}

template< class TMesh >
void CreateSquareTriangularMesh( typename TMesh::Pointer mesh )
{
  typedef TMesh                         MeshType;
  typedef typename MeshType::CellType   CellType;

  typedef itk::QuadEdgeMeshPolygonCell< CellType > QEPolygonCellType;

  if( mesh->GetNumberOfPoints( ) )
    {
    mesh->Clear( );
    mesh->ClearFreePointAndCellIndexesLists();
    }

  /////////////////////////////////////////////////////////////
  int expectedNumPts = 25;
  int expectedNumCells = 32;
  int simpleSquareCells[96] =
  {  0,  1,  6,
     0,  6,  5,
     1,  2,  7,
     1,  7,  6,
     2,  3,  8,
     2,  8,  7,
     3,  4,  9,
     3,  9,  8,
     5,  6, 11,
     5, 11, 10,
     6,  7, 12,
     6, 12, 11,
     7,  8, 13,
     7, 13, 12,
     8,  9, 14,
     8, 14, 13,
    10, 11, 16,
    10, 16, 15,
    11, 12, 17,
    11, 17, 16,
    12, 13, 18,
    12, 18, 17,
    13, 14, 19,
    13, 19, 18,
    15, 16, 21,
    15, 21, 20,
    16, 17, 22,
    16, 22, 21,
    17, 18, 23,
    17, 23, 22,
    18, 19, 24,
    18, 24, 23 };

  typedef typename TMesh::PointType PointType;
  std::vector< PointType > pts = GeneratePointCoordinates< TMesh >( 5 );

  for(int i=0; i<expectedNumPts; i++)
    {
    mesh->SetPoint( i, pts[i] );
    }

  typename CellType::CellAutoPointer cellpointer;
  QEPolygonCellType *poly;

  for(int i=0; i<expectedNumCells; i++)
    {
    poly = new QEPolygonCellType( 3 );
    cellpointer.TakeOwnership( poly );
    cellpointer->SetPointId( 0, simpleSquareCells[3*i] );
    cellpointer->SetPointId( 1, simpleSquareCells[3*i+1] );
    cellpointer->SetPointId( 2, simpleSquareCells[3*i+2] );
    mesh->SetCell( i, cellpointer );
    }
}

//----------------------------------------------------------------------------
template< class TMesh >
void CreateTetraedronMesh( typename TMesh::Pointer mesh )
{
  typedef TMesh                         MeshType;
  typedef typename MeshType::CellType   CellType;

  typedef itk::QuadEdgeMeshPolygonCell< CellType > QEPolygonCellType;

  if( mesh->GetNumberOfPoints( ) )
    {
    mesh->Clear( );
    mesh->ClearFreePointAndCellIndexesLists();
    }

  /////////////////////////////////////////////////////////////
  int expectedNumPts = 4;
  int expectedNumCells = 4;
  int simpleSquareCells[12] =
  {  0,  1,  2,
     1,  0,  3,
     1,  3,  2,
     2,  3,  0 };

  typedef typename TMesh::PointType PointType;
  std::vector< PointType > pts( 4 );
  int i(0);
  pts[i][0] = 0.; pts[i][1] = 1.; pts[i++][2] = 0.;
  pts[i][0] = 0.; pts[i][1] = -1.; pts[i++][2] = 0.;
  pts[i][0] = -1.; pts[i][1] = 0.; pts[i++][2] = 0.;
  pts[i][0] = 0.; pts[i][1] = 0.; pts[i++][2] = 1.;

  for( i=0; i<expectedNumPts; i++)
    {
    mesh->SetPoint( i, pts[i] );
    }

  typename CellType::CellAutoPointer cellpointer;
  QEPolygonCellType *poly;

  for( i=0; i<expectedNumCells; i++)
    {
    poly = new QEPolygonCellType( 3 );
    cellpointer.TakeOwnership( poly );
    cellpointer->SetPointId( 0, simpleSquareCells[3*i] );
    cellpointer->SetPointId( 1, simpleSquareCells[3*i+1] );
    cellpointer->SetPointId( 2, simpleSquareCells[3*i+2] );
    mesh->SetCell( i, cellpointer );
    }
}


//----------------------------------------------------------------------------
template< class TMesh >
void CreateSamosa( typename TMesh::Pointer mesh )
{
  typedef TMesh                         MeshType;
  typedef typename MeshType::CellType   CellType;

  typedef itk::QuadEdgeMeshPolygonCell< CellType > QEPolygonCellType;

  if( mesh->GetNumberOfPoints( ) )
    {
    mesh->Clear( );
    mesh->ClearFreePointAndCellIndexesLists();
    }

  /////////////////////////////////////////////////////////////
  int expectedNumPts = 3;
  int expectedNumCells = 2;
  int simpleSquareCells[6] =
  {  0,  1,  2,
     1,  0,  2 };

  typedef typename TMesh::PointType PointType;
  std::vector< PointType > pts( 3 );
  int i(0);
  pts[i][0] = 0.; pts[i][1] = 1.; pts[i++][2] = 0.;
  pts[i][0] = 0.; pts[i][1] = -1.; pts[i++][2] = 0.;
  pts[i][0] = -1.; pts[i][1] = 0.; pts[i++][2] = 0.;

  for( i=0; i<expectedNumPts; i++)
    {
    mesh->SetPoint( i, pts[i] );
    }

  typename CellType::CellAutoPointer cellpointer;
  QEPolygonCellType *poly;

  for( i=0; i<expectedNumCells; i++)
    {
    poly = new QEPolygonCellType( 3 );
    cellpointer.TakeOwnership( poly );
    cellpointer->SetPointId( 0, simpleSquareCells[3*i] );
    cellpointer->SetPointId( 1, simpleSquareCells[3*i+1] );
    cellpointer->SetPointId( 2, simpleSquareCells[3*i+2] );
    mesh->SetCell( i, cellpointer );
    }
}
