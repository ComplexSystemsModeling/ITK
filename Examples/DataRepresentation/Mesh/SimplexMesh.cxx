#include "itkQuadEdgeMeshWithDual.h"
#include "itkQuadEdgeMeshToQuadEdgeMeshWithDualFilter.h"
#include "itkQuadEdgeMeshWithDualAdaptor.h"
#include "itkVTKPolyDataWriter.h"
#include "itkQuadEdgeMeshEulerOperatorsTestHelper.h"

#include <iostream>

int main( int, char ** )
{
  //-----------------------------------------
  //  Define all the types we will work with
  //-----------------------------------------

  // two base tyes for meshes
  const unsigned int dimension = 3;
  typedef float PixelType;

  // the mesh type
  typedef itk::QuadEdgeMeshWithDual< PixelType, dimension > SimplexMeshType;

  // the primal to primal+dual filter
  typedef itk::QuadEdgeMeshToQuadEdgeMeshWithDualFilter< SimplexMeshType > FillDualFilterType;

  // all the filters to write the result
  typedef itk::VTKPolyDataWriter<SimplexMeshType> MeshWriterType;
  typedef itk::QuadEdgeMeshWithDualAdaptor< SimplexMeshType >  AdaptorType;
  typedef itk::VTKPolyDataWriter< AdaptorType > DualMeshWriterType;

  //--------------------------------------------------------------
  // Create the DataStructures that will hold the data in memory
  //--------------------------------------------------------------

  // main DataStructure
  SimplexMeshType::Pointer myPrimalMesh = SimplexMeshType::New();

  //-----------------------------------------
  // Create an input mesh to toy around with
  //-----------------------------------------

  std::cout << "Main: Create a square triangulated planar mesh." << std::endl;
  CreateSquareTriangularMesh< SimplexMeshType >( myPrimalMesh );

  std::cout << "Main: Poke a hole in the mesh to have two boundaries." << std::endl;
  myPrimalMesh->LightWeightDeleteEdge( myPrimalMesh->FindEdge( 11, 12 ) );

  //------------
  // Do the job
  //------------

  std::cout << "Main: Apply filter." << std::endl;
  FillDualFilterType::Pointer fillDual = FillDualFilterType::New();
  fillDual->SetInput( myPrimalMesh );
  try
    {
    fillDual->Update( );
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Main: Exception thrown while processing the input." << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  //-----------------------------------------------------
  // Write the original mesh, and the computed dual mesh
  //-----------------------------------------------------

  std::cout << "Main: Write primal mesh." << std::endl;
  MeshWriterType::Pointer writer1 = MeshWriterType::New();
  writer1->SetInput( fillDual->GetOutput() );
  writer1->SetFileName( "TestSquareTriangularMesh.vtk" );
  writer1->Write();

  AdaptorType* adaptor = new AdaptorType();
  adaptor->SetInput( fillDual->GetOutput() );

  std::cout << "Main: Write dual mesh." << std::endl;
  DualMeshWriterType::Pointer writer2 = DualMeshWriterType::New();
  writer2->SetInput( adaptor );
  writer2->SetFileName( "TestSquareTriangularSimplexMesh.vtk" );
  writer2->Write();

  //-----------------------------------------------------
  // and ... we're outta here.
  //-----------------------------------------------------

  return EXIT_SUCCESS;
}
//--------------------------------------------------------------------------------------------------
