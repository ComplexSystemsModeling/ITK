#include "itkQuadEdgeMeshWithDual.h"
#include "itkQuadEdgeMeshToQuadEdgeMeshWithDualFilter.h"
#include "itkQuadEdgeMeshWithDualAdaptor.h"
#include "itkVTKPolyDataWriter.h"
#include "itkVTKPolyDataReader.h"
#include "itkQuadEdgeMeshEulerOperatorsTestHelper.h"

#include <iostream>

int main( int argc, char **argv )
{
  //-----------------------------------------
  //  Check inputs
  //-----------------------------------------

  bool UseDefaultMesh = true;

  if( argc == 1 )
    {
    std::cout << "No input mesh file provided, using default planar mesh.";
    std::cout << std::endl;
    }

  if( argc == 2 )
    {
    // TODO: check that files exists
    UseDefaultMesh = false;
    }

   if( argc > 2 )
    {
    // TODO Print usage
    exit(0);
    }

  //-----------------------------------------
  //  Define all the types we will work with
  //-----------------------------------------

  // two base tyes for meshes
  const unsigned int dimension = 3;
  typedef float PixelType;

  // the mesh type
  typedef itk::QuadEdgeMeshWithDual< PixelType, dimension > MeshType;

  // the primal to primal+dual filter
  typedef itk::QuadEdgeMeshToQuadEdgeMeshWithDualFilter< MeshType > FillDualFilterType;

  // all the filters to write the result
  typedef itk::VTKPolyDataWriter< MeshType >            MeshWriterType;
  typedef itk::QuadEdgeMeshWithDualAdaptor< MeshType >  AdaptorType;
  typedef itk::VTKPolyDataWriter< AdaptorType >         DualMeshWriterType;

  //--------------------------------------------------------------
  // Create the DataStructures that will hold the data in memory
  //--------------------------------------------------------------

  // main DataStructure
  MeshType::Pointer myMesh = MeshType::New();

  //-----------------------------------------
  // Create an input mesh to toy around with
  //-----------------------------------------

  if( UseDefaultMesh )
    {
    // NOTE ALEX: nothing prevents us from proposing a nonplanar option
    // by using a command line argument to switch to Regular Sphere e.g.
    std::cout << "Main: Create a square triangulated planar mesh." << std::endl;
    CreateSquareTriangularMesh< MeshType >( myMesh );

    // NOTE ALEX: this is not random. We should be carefull here and write
    // something more generic by deleting the first edge in the edge container
    // for example. This code could then be reused for any mesh.
    std::cout << "Main: Poke a hole in the mesh to have two boundaries." << std::endl;
    myMesh->LightWeightDeleteEdge( myMesh->FindEdge( 11, 12 ) );
    }
  else
    {
    // NOTE ALEX: we suppose the file exist
    // we suppose the extension is provided
    // we suppose it's legacy vtk file.

    typedef itk::VTKPolyDataReader< MeshType > MeshReaderType;
    MeshReaderType::Pointer reader = MeshReaderType::New();
    reader->SetFileName( argv[1] );
    reader->Update();
    // NOTE ALEX: do we need a delete here?
    myMesh = reader->GetOutput();
    }

  //------------
  // Do the job
  //------------

  std::cout << "Main: Apply filter." << std::endl;
  FillDualFilterType::Pointer fillDual = FillDualFilterType::New();
  fillDual->SetInput( myMesh );
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

  // By default, the new datastructure is seen by ITK filters as the primal mesh
  // you can use it in a regular pipeline. Most of the existing filters will be compatible.
  std::cout << "Main: Write primal mesh." << std::endl;
  MeshWriterType::Pointer writer1 = MeshWriterType::New();
  writer1->SetInput( fillDual->GetOutput() );
  writer1->SetFileName( "TestSquareTriangularMesh.vtk" );
  writer1->Write();

  // if you want to access the dual with existing filters
  // you can just use the adaptor
  AdaptorType* adaptor = new AdaptorType();
  adaptor->SetInput( fillDual->GetOutput() );

  // then use a regular pipeline
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
