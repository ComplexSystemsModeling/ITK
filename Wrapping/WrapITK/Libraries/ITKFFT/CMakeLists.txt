include(FindFFTW.cmake)

WRAP_LIBRARY(ITKFFT)

set(WRAPPER_LIBRARY_DEPENDS ITKCommon  ITKStatistics ITKTransform ITKImageFunction ITKImageFilterBase ITKImageGrid ITKFiniteDifference ITKCurvatureFlow ITKAnisotropicSmoothing ITKMesh ITKIOBase ITKSpatialObjects ITKImageCompose ITKImageStatistics ITKPath ITKImageIntensity ITKSmoothing ITKImageGradient ITKImageFeature ITKThresholding ITKConnectedComponents ITKIO-PNG ITKMathematicalMorphology ITKNarrowBand ITKDistanceMap ITKSignedDistanceFunction ITKOptimizers ITKImageCompare ITKQuadEdgeMesh ITKFastMarching ITKLevelSets ITKAntiAlias ITKPolynomials ITKBiasCorrection ITKBioCell ITKIO-RAW ITKBlox ITKClassifiers ITKDeformableMesh ITKDeformationField ITKDiffusionTensorImage ITKEigen ITKRegistrationCommon ITKFEM ITKPDEDeformableRegistration ITKFEMRegistration)

# find needed libs
set(libs ${ITK_LIBRARIES})

if(WRAP_complex_float AND WRAP_float)
  if(FFTWF_LIB)
    set(libs ${libs} ${FFTWF_LIB})
  endif(FFTWF_LIB)
endif(WRAP_complex_float AND WRAP_float)

if(WRAP_complex_double AND WRAP_double)
  if(FFTWD_LIB)
    set(libs ${libs} ${FFTWD_LIB})
  endif(FFTWD_LIB)
endif(WRAP_complex_double AND WRAP_double)

set(WRAPPER_LIBRARY_LINK_LIBRARIES ${ITK_LIBRARIES})

AUTO_INCLUDE_MODULES()
END_WRAP_LIBRARY()
