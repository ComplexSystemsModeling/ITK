#include "itkSimplexMeshGeometry.h"
#include "itkFixedArray.h"

#include <algorithm>
#include <iostream>
#include <vector>
#include <set>

namespace itk
{

template< class TMesh >
class QuadEdgeMeshWithDualToSimplexMeshAdaptor
{
public:
  /** Standard typedefs. */
  typedef QuadEdgeMeshWithDualToSimplexMeshAdaptor      Self;
  typedef Self*                        Pointer;
  typedef const Self*                  ConstPointer;

  /** to foul the itkSuperclassTraitMacro */
  typedef TMesh Superclass;

  /** Convenient constants obtained from MeshTraits. */
  itkStaticConstMacro(PointDimension, unsigned int,
                      Superclass::PointDimension);
  itkStaticConstMacro(MaxTopologicalDimension, unsigned int,
                      Superclass::MaxTopologicalDimension);

  /** those types need to be defined for the itkDeformableSimplexMesh3DFilter */
  typedef typename Superclass::PixelType               PixelType;
  typedef typename Superclass::PointType               PointType;
  typedef typename Superclass::CellType                CellType;
  typedef typename Superclass::PointIdentifier         PointIdentifier;
  typedef typename Superclass::PointIdIterator         PointIdIterator;
  typedef typename Superclass::CellTraits              CellTraits;
  typedef typename Superclass::PointsContainer         PointsContainer;
  typedef typename Superclass::CellsContainer          CellsContainer;
  typedef typename Superclass::PointsContainerPointer  PointsContainerPointer;
  typedef typename Superclass::CellsContainerPointer   CellsContainerPointer;

  /** API functions from itkSimplexMesh */

  /** */
  typedef CovariantVector< typename VectorType::ValueType, 3 > CovariantVectorType;

  /** definition for a set of neighbor indices */
  typedef std::vector< SizeValueType >        NeighborListType;
  typedef std::set< SizeValueType >           NeighborSetType;
  typedef typename NeighborSetType::iterator  NeighborSetIterator;

  /** map containing a SimplexMeshGeomtry data object for each mesh */
  typedef typename SimplexMeshGeometry::IndexArray                  IndexArray;
  typedef itk::MapContainer< SizeValueType, SimplexMeshGeometry * > GeometryMapType;
  typedef typename GeometryMapType::Pointer	                        GeometryMapPointer;
  typedef typename GeometryMapType::Iterator                        GeometryMapIterator;
  typedef typename GeometryMapType::ConstIterator                   GeometryMapConstIterator;


  /** map containing a SimplexMeshGeometry data object for each mesh point */
  //typedef itk::MapContainer< SizeValueType, SimplexMeshGeometry * > GeometryMapType;
  typedef Self::PointsContainer                    GeometryMapType;
  typedef typename GeometryMapType::Pointer        GeometryMapPointer;
  typedef typename GeometryMapType::Iterator       GeometryMapIterator;
  typedef typename GeometryMapType::ConstIterator  GeometryMapConstIterator;


  /** set the map of geometrydata to the new pointer */
  itkSetMacro(GeometryData, GeometryMapPointer);

  /** returns the current map of geometrydata */
  itkGetConstReferenceMacro(GeometryData, GeometryMapPointer);

  /**
   * Get the geometry data for a specified point
   */
  //GeometryMapPointer & GetGeometryData( void );
  /* ****** Alex! should we have Macro or function definition ??? */

  /**
   * Get the three direct neighbors of a point
   */
  IndexArray GetNeighbors(PointIdentifier pointId) const
    {
    IndexArray neighborArray;
    PointType point;
    int i = 0;
    if( m_Mesh->GetPoint( pointId, pointt )
      {
      if( point.IsInternal() )
        {
        QuadEdgeType *start = point.GetEdge();
        QuadEdgeType *current = start;
        do
          {
          //neighborArray[i++] = current->GetSym().first;
          neighborArray[i++] = current->GetDestination().first;
          current = current.GetOnext();
          }
        while( current != start );
        }
      }
    return neighborArray;
    }
  /**
   * Get all neighbor points with a specified radius
   */
  NeighborListType* GetNeighbors( PointIdentifier idx, unsigned int radius, NeighborListType *list)
    {
    if( list == NULL )
      {
      list = new NeighborListType();
      IndexArray neighborArray = GetNeighbors( idx );
      list->push_back( neighborArray[0] );
      list->push_back( neighborArray[1] );
      list->push_back( neighborArray[2] );

      if( radius > 0 )
        {
        list = GetNeighbors( neighborArray[0], radius - 1, list );
        list = GetNeighbors( neighborArray[1], radius - 1, list );
        list = GetNeighbors( neighborArray[2], radius - 1, list );
        }
      NeighborListType::Iterator it = std::find( list->begin(), list->end(), idx );
      if( it != list->end() )
        list->erase( it );

      return list;
      }
    else
      {
      IndexArray neighborArray = GetNeighbors( idx );
      NeighborListType::Iterator foundIt1 = std::find( list->begin(), list->end(), neighborArray[0] );
      NeighborListType::Iterator foundIt2 = std::find( list->begin(), list->end(), neighborArray[1] );
      NeighborListType::Iterator foundIt3 = std::find( list->begin(), list->end(), neighborArray[2] );

      bool found1 = false, found2 = false, found3 = false;
      if( foundIt1 != endIt )
        found1 = true;
      if( foundIt2 != endIt )
        found2 = true;
      if( foundIt3 != endIt )
        found3 = true;

      if( != found1 )
        list->push_back( neighborArray[0] );
      if( != found2 )
        list->push_back( neighborArray[1] );
      if( != found3 )
        list->push_back( neighborArray[2] );

      if( radius == 0 )
        return list;
      else
        {
        list = GetNeighbors( neighborArray[0], radius - 1, list );
        list = GetNeighbors( neighborArray[0], radius - 1, list );
        list = GetNeighbors( neighborArray[0], radius - 1, list );
        return list;
        }
      }
    }

  /** compute the normal vector in the specified mesh point */
  CovariantVectorType ComputeNormal(PointIdentifier idx) const;


  /** API that will be used by itkDeformableSimplexMesh3DFilter */
  const PointIdentifier  GetNumberOfPoints() const { return m_Mesh->GetDualPoints()->size(); };
  const CellsContainerPointer  GetCells()    const { return m_Mesh->GetDualCells();  };
  const PointsContainerPointer GetPoints()   const { return m_Mesh->GetDualPoints(); };
  void SetInput( TMesh* mesh ) { m_Mesh = mesh; };

private:
  TMesh* m_Mesh;

  GeometryMapPointer m_GeometryData;

};

}
