#include <iostream>

namespace itk
{

template< class TMesh >
class QuadEdgeMeshWithDualAdaptor
{
public:
  /** Standard typedefs. */
  typedef QuadEdgeMeshWithDualAdaptor  Self;
  typedef Self*                        Pointer;
  typedef const Self*                  ConstPointer;

  /** to foul the itkSuperclassTraitMacro */
  typedef TMesh Superclass;

  /** Convenient constants obtained from MeshTraits. */
  itkStaticConstMacro(PointDimension, unsigned int,
                      Superclass::PointDimension);
  itkStaticConstMacro(MaxTopologicalDimension, unsigned int,
                      Superclass::MaxTopologicalDimension);

  /** those types need to be defined for the itkVTKPolydataWriter */
  itkSuperclassTraitMacro( PixelType );
  itkSuperclassTraitMacro( PointType );
  itkSuperclassTraitMacro( CellType  );
  itkSuperclassTraitMacro( PointIdentifier );
  itkSuperclassTraitMacro( PointIdIterator );
  itkSuperclassTraitMacro( CellTraits      );
  itkSuperclassTraitMacro( PointsContainer );
  itkSuperclassTraitMacro( CellsContainer  );
  itkSuperclassTraitMacro( PointsContainerPointer );
  itkSuperclassTraitMacro( CellsContainerPointer  );

  /** API that will be used by itkVTKPolydataWriter */
  const PointIdentifier  GetNumberOfPoints() const { return GetPoints()->size(); };
  const CellsContainerPointer  GetCells()    const { return m_Mesh->GetDualCells();  };
  const PointsContainerPointer GetPoints()   const { return m_Mesh->GetDualPoints(); };
  void SetInput( TMesh* mesh ) { m_Mesh = mesh; };

private:
  TMesh* m_Mesh;

};

}
