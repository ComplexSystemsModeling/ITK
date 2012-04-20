#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshWithDual.h"
#include "itkQuadEdgeMeshWithDualAdaptor.h"
#include "itkDefaultStaticMeshTraits.h"
#include "itkVTKPolyDataWriter.h"
#include "itkQuadEdgeMeshBoundaryEdgesMeshFunction.h"
#include "itkQuadEdgeMeshPolygonCell.h"
#include "itkQuadEdgeMeshEulerOperatorsTestHelper.h"

#include <iostream>

int main( int argc, char ** argv )
{
  //-----------------------------------------
  //  Define all the types we will work with
  //-----------------------------------------

  const unsigned int dimension = 3;

  typedef itk::QuadEdgeMeshWithDual< float, dimension > SimplexMeshType;

  typedef SimplexMeshType::PointIdentifier PointIdentifier;
  typedef SimplexMeshType::CellIdentifier  CellIdentifier;
  typedef SimplexMeshType::PointsContainer PointsContainer;
  typedef SimplexMeshType::CellsContainer  CellsContainer;
  typedef SimplexMeshType::CellType        CellType;
  typedef SimplexMeshType::PointType       PointType;
  typedef SimplexMeshType::PointIdList     PointIdList;
  typedef SimplexMeshType::QEType          QuadEdgeType;

  typedef CellType::PointIdConstIterator PointIdConstIterator;

  typedef PointsContainer::ConstIterator PointIterator;
  typedef CellsContainer::ConstIterator  CellIterator;

  // this will be use to find and track boundaries
  typedef itk::QuadEdgeMeshBoundaryEdgesMeshFunction< SimplexMeshType > BoundaryLocatorType;

  //--------------------------------------------------------------
  // Create the DataStructures that will hold the data in memory
  //--------------------------------------------------------------

  // LUT that will need to be integrated in the DataStructure
  std::map< CellIdentifier ,PointIdentifier >   DualPointFromPrimalTriangleLUT;

  // main DataStructure
  SimplexMeshType::Pointer myPrimalMesh = SimplexMeshType::New();

  //-----------------------------------------
  // Create an input mesh to toy around with
  //-----------------------------------------

  std::cout << "Create a square triangulated planar mesh." << std::endl;
  CreateSquareTriangularMesh< SimplexMeshType >( myPrimalMesh );

  std::cout << "Poke a hole in the mesh to have two boundaries." << std::endl;
  myPrimalMesh->LightWeightDeleteEdge( myPrimalMesh->FindEdge( 11, 12 ) );

  //-------------------------------------------------------
  // First pass: dual points for 2D cells (triangles here)
  //-------------------------------------------------------

  std::cout << "Generation of dual points: barycenter of primal cells" << std::endl;

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
            std::cout << "We found a polygon, this is not handled right now." << std::endl;
            found = true;
            }
          else
            {
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
              for( unsigned int i =0; i < dimension; i++ )
                {
                d_point[i] += point[i];
                }
              current++;
              }
            for( unsigned int i =0; i < dimension; i++ )
              {
              d_point[i] /= dimension;
              }
            // push dual point in dualPoints container
            myPrimalMesh->SetDualPoint( numberOfDualPoints, d_point );
            CellIdentifier cellIdentifier = cellIterator.Index();
            DualPointFromPrimalTriangleLUT[cellIdentifier] = numberOfDualPoints;

 			// step 1 - loop on all edges of cell
			// step 2 - Get left face of each edge (dualpoint) rot or GetLeft
			// setp 3 - assign cellIdentifier and pointIdentifier to edge (QuadEdgeMeshTraits::GeometricalQuadEdge::OriginRefType)

			SimplexMeshType::QEType::GeometricalQuadEdge::DualOriginRefType dualOriginRef;
			dualOriginRef.first = cellIdentifier;
			dualOriginRef.second = numberOfDualPoints;

			SimplexMeshType::CellAutoPointer cellPointer;
			myPrimalMesh->GetCell( cellIdentifier, cellPointer );
			//QuadEdgeType *currentEdge = cellPointer->GetEdgeRingEntry();
			PointIdConstIterator cellPointsItertor = cellPointer->GetPointIds();
			
			PointType primalPoint;
			myPrimalMesh->GetPoint( cellPointsItertor.value, primalPoint);
			QuadEdgeType * currentEdge = primalPoint.GetEdge();
			QuadEdgeType *firstEdge = currentEdge;
			do
			  {
			  currentEdge->SetLeft( dualOriginRef );
              currentEdge = currentEdge->GetLnext();
			  }while( currentEdge != firstEdge );

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

  //-------------------------------------------------------
  // Second pass: dual cells (polygons) for primal points
  //-------------------------------------------------------

  std::cout << "Generate Dual Cells from the primal points" << std::endl;

  const PointsContainer *primalPoints = myPrimalMesh->GetPoints();
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
          pointidlist.push_back( DualPointFromPrimalTriangleLUT[ leftTriangle.first ] );
		  pointidlist.push_back( leftTriangle.second );

          current = current->GetOnext();

          } while( current != start );

        // point list is complete, add the dual cell to the dual mesh;
        myPrimalMesh->AddDualFace( pointidlist );
        }
      // next point
      pointIterator++;
      }
    }

  //-----------------------------------------
  // last pass: Treat the borders as 1D cells
  //-----------------------------------------

  BoundaryLocatorType::Pointer boundaryEdges = BoundaryLocatorType::New();
  SimplexMeshType::EdgeListPointerType boundaryEdgesPointerList = boundaryEdges->Evaluate( *myPrimalMesh );

  // for each boundary
  unsigned int i = 0;
  while( !boundaryEdgesPointerList->empty() )
    {

    std::cout << "Boundary #" << i << std::endl;
    i++;

    // get the first edge and remove it from the list
    QuadEdgeType* firstEdge = boundaryEdgesPointerList->front();
    QuadEdgeType* currentEdge = firstEdge;
    boundaryEdgesPointerList->pop_front();

    // circulate around the boundary and do what you have to do
    bool firstTime = true;
    PointIdentifier previousPointId;
    PointIdentifier currentPointId;
    PointIdentifier firstPointId = myPrimalMesh->GetNumberOfDualPoints();
    do
      {
      // create a new point in the middle of the edge
      // this is a dual point. Cells are sample in 2d, boundaries are sample in 1D
      PointIdentifier originPointId = currentEdge->GetOrigin().first;
      PointIdentifier destinationPointId = currentEdge->GetDestination().first;

      PointType originPoint = myPrimalMesh->GetPoint( originPointId );
      PointType destinationPoint = myPrimalMesh->GetPoint( destinationPointId );
      PointType currentPoint;
      for( unsigned int k =0; k<dimension; k++ )
        {
        currentPoint[k] = (originPoint[k] + destinationPoint[k]) / 2 ;
        }

      // add the new border point P1 to the dual point container
      currentPointId = myPrimalMesh->AddDualPoint( currentPoint );

      // find the dual point P2 associated with the face on the left
      QuadEdgeType::DualOriginRefType leftTriangle = currentEdge->GetRight();

      // add the dual edge P1-P2
      PointIdentifier leftDualPointId = DualPointFromPrimalTriangleLUT[ leftTriangle.first ];
      myPrimalMesh->AddDualEdge( currentPointId, leftDualPointId );

      // add the edge linking two new border dual points
      // either previous-current, or current-next
      if( firstTime == true )
        firstTime = false;
      else
        {
        myPrimalMesh->AddDualEdge( previousPointId, currentPointId );

        // create a point ID list to hold the dual point IDs while iterating to create the dual cell
        PointIdList pointidlist;
        pointidlist.push_back( previousPointId );
        QuadEdgeType *myEdge = currentEdge->GetOnext();
        do
          {
          QuadEdgeType::DualOriginRefType myleftTriangle = myEdge->GetLeft();
          PointIdentifier myleftDualPointId =  DualPointFromPrimalTriangleLUT[ myleftTriangle.first ];
          pointidlist.push_back( myleftDualPointId );
          myEdge = myEdge->GetOnext();
          }
        while( !myEdge->IsAtBorder() );

        pointidlist.push_back( currentPointId );

        // point list is complete, add the dual cell to the dual mesh;
        myPrimalMesh->AddDualFace( pointidlist );
        }
      previousPointId = currentPointId;

      // Update currentEdge
      currentEdge = currentEdge->GetLnext();
      }
    while (currentEdge != firstEdge);

    PointIdList pointidlist;
    pointidlist.push_back( previousPointId );
    QuadEdgeType *myEdge = currentEdge->GetOnext();
    do
      {
      QuadEdgeType::DualOriginRefType myleftTriangle = myEdge->GetLeft();
      PointIdentifier myleftDualPointId =  DualPointFromPrimalTriangleLUT[ myleftTriangle.first ];
      pointidlist.push_back( myleftDualPointId );
      myEdge = myEdge->GetOnext();
      }
    while( !myEdge->IsAtBorder() );

    pointidlist.push_back( firstPointId );

    // point list is complete, add the dual cell to the dual mesh;
    myPrimalMesh->AddDualFace( pointidlist );
    }

  //-----------------------------------------------------
  // Write the original mesh, and the computed dual mesh
  //-----------------------------------------------------

  typedef itk::VTKPolyDataWriter<SimplexMeshType> MeshWriterType;
  MeshWriterType::Pointer writer1 = MeshWriterType::New();
  writer1->SetInput( myPrimalMesh );
  writer1->SetFileName( "TestSquareTriangularMesh.vtk" );
  writer1->Write();

  typedef itk::QuadEdgeMeshWithDualAdaptor< SimplexMeshType >  AdaptorType;
  AdaptorType* adaptor = new AdaptorType();
  adaptor->SetInput( myPrimalMesh );

  typedef itk::VTKPolyDataWriter< AdaptorType > DualMeshWriterType;
  DualMeshWriterType::Pointer writer2 = DualMeshWriterType::New();
  writer2->SetInput( adaptor );
  writer2->SetFileName( "TestSquareTriangularSimplexMesh.vtk" );
  writer2->Write();

  return EXIT_SUCCESS;
}
