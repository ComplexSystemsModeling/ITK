#include "itkQuadEdgeMeshBoundaryEdgesMeshFunction.h"
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

  // our types

  // NOTE ALEX: to extract from MeshType
  static const unsigned int dimension = 3;

  typedef typename TInputMesh::Pointer             InputMeshPointerType;

  typedef typename TOutputMesh::Pointer            OutputMeshPointerType;
  typedef typename TOutputMesh::PointIdList        PointIdList;
  typedef typename TOutputMesh::CellIdentifier     CellIdentifier;
  typedef typename TOutputMesh::PointIdentifier    PointIdentifier;
  typedef typename TOutputMesh::PointsContainer    PointsContainer;
  typedef typename TOutputMesh::CellsContainer     CellsContainer;
  typedef typename TOutputMesh::CellType           CellType;
  typedef typename TOutputMesh::PointType          PointType;
  typedef typename TOutputMesh::QEType             QuadEdgeType;
  typedef typename TOutputMesh::PolygonCellType    PolygonCellType;
  typedef typename TOutputMesh::CellAutoPointer    CellAutoPointer;
  typedef typename CellType::PointIdConstIterator  PointIdConstIterator;
  typedef typename QuadEdgeType::OriginRefType     OriginRefType;
  typedef typename QuadEdgeType::DualOriginRefType DOrgRefType;
  typedef typename CellsContainer::ConstIterator   CellIterator;
  typedef typename PointsContainer::ConstIterator  PointIterator;

  // this will be use internally to find and track boundaries
  typedef typename itk::QuadEdgeMeshBoundaryEdgesMeshFunction< TOutputMesh > BoundaryLocatorType;
  typedef typename TOutputMesh::EdgeListPointerType                          EdgeListPointerType;

protected:
  QuadEdgeMeshToQuadEdgeMeshWithDualFilter() {};

  virtual ~QuadEdgeMeshToQuadEdgeMeshWithDualFilter() {};

  void GenerateData()
  {
    // use the superclass method to copy input over output
    // this ensure the const correctness of the inputs
    std::cout << "Filter: Copy input to output." << std::endl;
    this->CopyInputMeshToOutputMesh();

    //-------------------------------------------------------
    // First pass: dual points for 2D cells (triangles here)
    //-------------------------------------------------------
    std::cout << "Filter: Generation of dual points: barycenter of primal cells" << std::endl;
    ComputeDualPointsForAllPolygons( );

    //-------------------------------------------------------
    // Second pass: dual cells (polygons) for primal points
    //-------------------------------------------------------
    std::cout << "Filter: Generate Dual Cells from the primal points" << std::endl;
    ComputeDualPolygonsForAllPoints( );

    //-----------------------------------------
    // last pass: Treat the borders as 1D cells
    //-----------------------------------------
    std::cout << "Filter: Handle Borders." << std::endl;
    typename BoundaryLocatorType::Pointer boundaryEdges = BoundaryLocatorType::New();
    OutputMeshPointerType myPrimalMesh = this->GetOutput();
    EdgeListPointerType boundaryEdgesPointerList = boundaryEdges->Evaluate( *myPrimalMesh );

    // for each boundary
    unsigned int i = 0;
    while( !boundaryEdgesPointerList->empty() )
      {
      std::cout << "Filter: Boundary #" << i++ << std::endl;
      CreateDualOfBorderPointsAndEdges( boundaryEdgesPointerList->front() );
      boundaryEdgesPointerList->pop_front();
      }
  };

  void PrintSelf(std::ostream & os, Indent indent) const
  {
    Superclass::PrintSelf(os, indent);
    //os << indent << "AbsoluteTolerance2: " << m_AbsoluteTolerance2 << std::endl;
  }

private:
  bool ComputeDualPointsForAllPolygons( )
  {
    OutputMeshPointerType myPrimalMesh = this->GetOutput();

    const CellsContainer *primalCells = myPrimalMesh->GetCells();
    if( primalCells )
      {
      CellIterator cellIterator = primalCells->Begin();
      CellIterator cellEnd = primalCells->End();
      while( cellIterator != cellEnd )
        {
        ComputeDualPointFromFaceAndSet( cellIterator );
        cellIterator++;
        }
      }
    return true;
  }

  bool CreateDualOfBorderPointsAndEdges( QuadEdgeType* currentEdge )
  {
    OutputMeshPointerType myPrimalMesh = this->GetOutput();

    QuadEdgeType*   firstEdge    = currentEdge;
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

        CreateDualCellOfBorderPoint( currentPointId, previousPointId, currentEdge );
        }
      previousPointId = currentPointId;
      currentEdge = currentEdge->GetLnext();
      }
    while( currentEdge != firstEdge );

    // NOTE ALEX: are we missing a dual edge here?
    CreateDualCellOfBorderPoint( firstPointId, previousPointId, currentEdge );

    return true;
  }

  bool ComputeDualPolygonsForAllPoints()
  {
    OutputMeshPointerType myPrimalMesh = this->GetOutput();

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

  bool CreateDualCellOfBorderPoint(
    PointIdentifier firstPointId,
    PointIdentifier previousPointId,
    QuadEdgeType*   currentEdge
    )
  {
    OutputMeshPointerType myPrimalMesh = this->GetOutput();

    PointIdList pointidlist;
    pointidlist.push_back( previousPointId );
    QuadEdgeType *myEdge = currentEdge->GetOnext();
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

  bool ComputeDualPointFromFaceAndSet( CellIterator cellIterator )
  {
    OutputMeshPointerType myPrimalMesh = this->GetOutput();

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


  // usual ITK stuff
  QuadEdgeMeshToQuadEdgeMeshWithDualFilter(const Self &);
  void operator=(const Self &);
};
}
