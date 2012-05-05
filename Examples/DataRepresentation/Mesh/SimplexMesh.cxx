#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshWithDual.h"
#include "itkQuadEdgeMeshWithDualAdaptor.h"
#include "itkDefaultStaticMeshTraits.h"
#include "itkVTKPolyDataWriter.h"
#include "itkQuadEdgeMeshBoundaryEdgesMeshFunction.h"
#include "itkQuadEdgeMeshPolygonCell.h"
#include "itkQuadEdgeMeshEulerOperatorsTestHelper.h"
#include "itkQuadEdgeMeshToQuadEdgeMeshFilter.h"


#include <iostream>

namespace itk
{
/**
 * \class QuadEdgeMeshToQuadEdgeMeshWithDualFilter
 * \brief TODO
 * \ingroup ITKQuadEdgeMeshFiltering
 */
template< class TInputMesh, class TOutputMesh=TInputMesh >
class ITK_EXPORT QuadEdgeMeshToQuadEdgeMeshWithDualFilter:
  public QuadEdgeMeshToQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
{
public:
  typedef QuadEdgeMeshToQuadEdgeMeshWithDualFilter                    Self;
  typedef SmartPointer< Self >                                        Pointer;
  typedef SmartPointer< const Self >                                  ConstPointer;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter< TInputMesh, TOutputMesh > Superclass;

  /** Run-time type information (and related methods).   */
  itkTypeMacro(QuadEdgeMeshToQuadEdgeMeshWithDualFilter, QuadEdgeMeshToQuadEdgeMeshFilter);

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro(Self);

protected:
  QuadEdgeMeshToQuadEdgeMeshWithDualFilter() {};

  virtual ~QuadEdgeMeshToQuadEdgeMeshWithDualFilter() {};

  void GenerateData()
  {
  // use the superclass method
  this->CopyInputMeshToOutputMesh();
  };

  void PrintSelf(std::ostream & os, Indent indent) const
  {
    Superclass::PrintSelf(os, indent);
    //os << indent << "AbsoluteTolerance2: " << m_AbsoluteTolerance2 << std::endl;
  }

private:
  QuadEdgeMeshToQuadEdgeMeshWithDualFilter(const Self &);
  void operator=(const Self &);
};
}

template< typename MeshType >
bool
ComputeDualPointFromFaceAndSet(
  MeshType* myPrimalMesh,
  typename MeshType::CellsContainer::ConstIterator cellIterator
  );

template< typename MeshType >
bool
ComputeDualPointsForAllPolygons(
  MeshType* myPrimalMesh
  );

template< typename MeshType >
bool
ComputeDualPolygonsForAllPoints(
  MeshType* myPrimalMesh
  );

template< typename MeshType >
bool
CreateDualCellOfBorderPoint(
  MeshType*                          myPrimalMesh,
  typename MeshType::PointIdentifier firstPointId,
  typename MeshType::PointIdentifier previousPointId,
  typename MeshType::QEType*         currentEdge
  );

template< typename MeshType >
bool
CreateDualOfBorderPointsAndEdges(
  MeshType*                          myPrimalMesh,
  typename MeshType::QEType*         currentEdge
  );

int main( int, char ** )
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
  typedef SimplexMeshType::PolygonCellType PolygonCellType;

  typedef SimplexMeshType::EdgeListPointerType EdgeListPointerType;

  typedef CellType::PointIdConstIterator PointIdConstIterator;

  typedef PointsContainer::ConstIterator PointIterator;
  typedef CellsContainer::ConstIterator  CellIterator;

  // this will be use to find and track boundaries
  typedef itk::QuadEdgeMeshBoundaryEdgesMeshFunction< SimplexMeshType > BoundaryLocatorType;

  //--------------------------------------------------------------
  // Create the DataStructures that will hold the data in memory
  //--------------------------------------------------------------

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
  ComputeDualPointsForAllPolygons< SimplexMeshType >( myPrimalMesh );

  //-------------------------------------------------------
  // Second pass: dual cells (polygons) for primal points
  //-------------------------------------------------------

  std::cout << "Generate Dual Cells from the primal points" << std::endl;
  ComputeDualPolygonsForAllPoints< SimplexMeshType >( myPrimalMesh );

  //-----------------------------------------
  // last pass: Treat the borders as 1D cells
  //-----------------------------------------

  BoundaryLocatorType::Pointer boundaryEdges = BoundaryLocatorType::New();
  EdgeListPointerType boundaryEdgesPointerList = boundaryEdges->Evaluate( *myPrimalMesh );

  // for each boundary
  unsigned int i = 0;
  while( !boundaryEdgesPointerList->empty() )
    {
    std::cout << "Boundary #" << i++ << std::endl;


    CreateDualOfBorderPointsAndEdges< SimplexMeshType >(
      myPrimalMesh, boundaryEdgesPointerList->front() );
    boundaryEdgesPointerList->pop_front();
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
//--------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------
template< typename MeshType >
bool
ComputeDualPointFromFaceAndSet(
  MeshType* myPrimalMesh,
  typename MeshType::CellsContainer::ConstIterator cellIterator
  )
{
  // NOTE ALEX: to extract from MeshType
  const unsigned int dimension = 3;

  typedef typename MeshType::CellIdentifier  CellIdentifier;
  typedef typename MeshType::PointIdentifier PointIdentifier;
  typedef typename MeshType::CellType        CellType;
  typedef typename MeshType::PointType       PointType;
  typedef typename MeshType::QEType          QuadEdgeType;
  typedef typename MeshType::PolygonCellType PolygonCellType;

  typedef typename CellType::PointIdConstIterator PointIdConstIterator;

  typedef typename QuadEdgeType::DualOriginRefType DOrgRefType;
  typedef typename MeshType::CellAutoPointer CellAutoPointer;

  // 1. compute dual point coordinate and push it to the container
  PointIdConstIterator current= cellIterator.Value()->PointIdsBegin();
  PointIdConstIterator end    = cellIterator.Value()->PointIdsEnd();
  PointType d_point;
  for( unsigned int i = 0; i < 3; i++ ) // dimension; i++ )
    {
    d_point[i] = 0.0;
    }
  while( current != end )
    {
    PointType point = myPrimalMesh->GetPoint( *current );
    for( unsigned int i = 0; i < 3; i++ ) // dimension; i++ )
      {
      d_point[i] += point[i];
      }
    current++;
    }
  for( unsigned int i =0; i < 3; i++ ) // dimension; i++ )
    {
    d_point[i] /= dimension;
    }
  PointIdentifier d_point_id = myPrimalMesh->AddDualPoint( d_point );

  // 2. Compute the new OriginRefType and set all the QEdges

  CellIdentifier cellIdentifier = cellIterator.Index();
  DOrgRefType dualOriginRef = DOrgRefType( cellIdentifier, d_point_id );


  // NOTE ALEX: isn't there a method in QE to do that?
  CellAutoPointer cellPointer;
  myPrimalMesh->GetCell( cellIdentifier, cellPointer );
  PolygonCellType* myCell = dynamic_cast< PolygonCellType* >( cellPointer.GetPointer() );
  QuadEdgeType *currentEdge = myCell->GetEdgeRingEntry();
  QuadEdgeType *firstEdge = currentEdge;
  do
    {
    currentEdge->SetLeft( dualOriginRef );
    currentEdge = currentEdge->GetLnext();
    }
  while( currentEdge != firstEdge );

  return true;
}
//--------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------
template< typename MeshType >
bool
ComputeDualPointsForAllPolygons(
  MeshType* myPrimalMesh
  )
{
  typedef typename MeshType::CellsContainer       CellsContainer;
  typedef typename CellsContainer::ConstIterator  CellIterator;

  const CellsContainer *primalCells = myPrimalMesh->GetCells();
  if( primalCells )
    {
    CellIterator cellIterator = primalCells->Begin();
    CellIterator cellEnd = primalCells->End();
    while( cellIterator != cellEnd )
      {
      ComputeDualPointFromFaceAndSet< MeshType >( myPrimalMesh, cellIterator );
      cellIterator++;
      }
    }
  return true;
}
//--------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------
template< typename MeshType >
bool
ComputeDualPolygonsForAllPoints(
  MeshType* myPrimalMesh
  )
{
  typedef typename MeshType::PointsContainer       PointsContainer;
  typedef typename PointsContainer::ConstIterator  PointIterator;
  typedef typename MeshType::PointIdList           PointIdList;
  typedef typename MeshType::PointType             PointType;
  typedef typename MeshType::QEType                QuadEdgeType;
  typedef typename MeshType::CellIdentifier        CellIdentifier;
  typedef typename MeshType::PointIdentifier       PointIdentifier;
  typedef typename QuadEdgeType::OriginRefType     OriginRefType;

  // Get hold of the point container
  const PointsContainer *primalPoints = myPrimalMesh->GetPoints();

  // if it s not empty, proceed
  if( primalPoints )
    {

    // for all points
    PointIterator pointIterator = primalPoints->Begin();
    PointIterator pointEnd      = primalPoints->End();
    while( pointIterator != pointEnd )
      {

      // borders are treated separately
      PointType point = pointIterator.Value();
      if( point.IsInternal() )
        {

        // grab the QEdge
        QuadEdgeType * start   = point.GetEdge();
        QuadEdgeType * current = start;

        // create a point ID list to hold the dual point IDs while
        // we are iterating around a primal point to create the dual cell
        PointIdList pointidlist;

        // iterate around the o-ring and record dual point Ids
        do
          {
          pointidlist.push_back( current->GetLeft().second );
          current = current->GetOnext();
          }
        while( current != start );

        // point list is complete, add the dual cell to the dual mesh;
        CellIdentifier DualFaceID = myPrimalMesh->AddDualFace( pointidlist );

        // compute new origin
        OriginRefType newORF = OriginRefType( current->GetOrigin().first, DualFaceID );
        // iterate around the o-ring and set the new origin
        do
          {
          current->SetOrigin( newORF );
          current = current->GetOnext();
          }
        while( current != start );
        }

      // next point
      pointIterator++;
      }

    } // endof if( primalPoints )

  return true;
}
//--------------------------------------------------------------------------------------------------


//--------------------------------------------------------------------------------------------------
template< typename MeshType >
bool
CreateDualCellOfBorderPoint(
  MeshType*                          myPrimalMesh,
  typename MeshType::PointIdentifier firstPointId,
  typename MeshType::PointIdentifier previousPointId,
  typename MeshType::QEType*         currentEdge
  )
{
  typename MeshType::PointIdList pointidlist;
  pointidlist.push_back( previousPointId );
  typename MeshType::QEType *myEdge = currentEdge->GetOnext();
  do
    {
    pointidlist.push_back( myEdge->GetLeft().second );
    myEdge = myEdge->GetOnext();
    }
  while( !myEdge->IsAtBorder() );
  pointidlist.push_back( firstPointId );

  // point list is complete, add the dual cell to the dual mesh;
  myPrimalMesh->AddDualFace( pointidlist );

  // NOTE ALEX: TODO here we have to reset the OriginRefType
  // loop around the onext, border or not

  return true;
}
//--------------------------------------------------------------------------------------------------


//--------------------------------------------------------------------------------------------------
template< typename MeshType >
bool
CreateDualOfBorderPointsAndEdges(
  MeshType*                          myPrimalMesh,
  typename MeshType::QEType*         currentEdge
  )
{
  typedef typename MeshType::PointIdentifier     PointIdentifier;
  typedef typename MeshType::QEType              QEType;
  typedef typename MeshType::PointType           PointType;

  QEType*         firstEdge    = currentEdge;
  bool            firstTime    = true;
  PointIdentifier firstPointId = myPrimalMesh->GetNumberOfDualPoints();
  PointIdentifier previousPointId;
  PointIdentifier currentPointId;
  do
    {
    // 1. always create a new point in the middle of the edge

    // this is a dual point. Cells are sampled in 2d, boundaries are sampled in 1D
    PointIdentifier originPointId      = currentEdge->GetOrigin().first;
    PointIdentifier destinationPointId = currentEdge->GetDestination().first;

    // border primal points
    PointType originPoint = myPrimalMesh->GetPoint( originPointId );
    PointType destinationPoint = myPrimalMesh->GetPoint( destinationPointId );

    PointType currentPoint;
    // NOTE ALEX: TODO extract dimension from MeshTye
    for( unsigned int k =0; k < 3; k++ ) // dimension; k++ )
      {
      currentPoint[k] = (originPoint[k] + destinationPoint[k]) / 2 ;
      }
    // add the new border dual point P1 to the dual point container
    currentPointId = myPrimalMesh->AddDualPoint( currentPoint );

    // 2. always add the edge between this (1D) point, and previous (2D) point

    // add the dual edge P1-P2
    // NOTE ALEX: do we know on which side the hole is. Is it stable?
    myPrimalMesh->AddDualEdge( currentPointId, currentEdge->GetRight().second );

    // 3. Almost always add the dual edge along the border,
    // in which case we also create the dual cell.

    // add the edge linking two new border dual points
    if( firstTime == true )
      {
      firstTime = false;
      }
    else
      {
      // NOTE ALEX: how to deal with OriginRefType in this case?
      // use the EdgeCellContainer ID?
      myPrimalMesh->AddDualEdge( previousPointId, currentPointId );

      CreateDualCellOfBorderPoint< MeshType >(
        myPrimalMesh, currentPointId, previousPointId, currentEdge );
      }
    previousPointId = currentPointId;
    currentEdge = currentEdge->GetLnext();
    }
  while( currentEdge != firstEdge );

  // NOTE ALEX: are we missing a dual edge here?
  CreateDualCellOfBorderPoint< MeshType >(
    myPrimalMesh, firstPointId, previousPointId, currentEdge );

  return true;
}
