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
  typedef std::vector< SizeValueType >					NeighborListType;
  typedef std::set< SizeValueType >						NeighborSetType;
  typedef typename NeighborSetType::iterator			NeighborSetIterator;

  /** definition for array of indices. */
  typedef typename SimplexMeshGeometry::IndexArray		IndexArray;

  /** map containing a SimplexMeshGeometry data object for each mesh point */
  typedef itk::MapContainer< SizeValueType, SimplexMeshGeometry * > GeometryMapType;

  /** smartpointer def for the geometry map */
  typedef typename GeometryMapType::Pointer				GeometryMapPointer;

  /** iterator definition for iterating over a geometry map */
  typedef typename GeometryMapType::Iterator			GeometryMapIterator;
  typedef typename GeometryMapType::ConstIterator		GeometryMapConstIterator;


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
  IndexArray GetNeighbors(PointIdentifier pointId) const;

  /**
   * Get all neighbor points with a specified radius
   */
  NeighborListType * GetNeighbors(PointIdentifier pointId, unsigned int radius, NeighborListType *list = NULL) const;

  /** compute the normal vector in the specified mesh point */
  CovariantVectorType ComputeNormal(PointIdentifier idx) const;


  /** API that will be used by itkDeformableSimplexMesh3DFilter */
  const PointIdentifier  GetNumberOfPoints() const { return GetPoints()->size(); };
  const CellsContainerPointer  GetCells()    const { return m_Mesh->GetDualCells();  };
  const PointsContainerPointer GetPoints()   const { return m_Mesh->GetDualPoints(); };
  void SetInput( TMesh* mesh ) { m_Mesh = mesh; };

private:
  TMesh* m_Mesh;

};

}
