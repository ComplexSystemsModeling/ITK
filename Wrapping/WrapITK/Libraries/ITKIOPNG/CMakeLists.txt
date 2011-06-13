WRAP_LIBRARY(ITKIOPNG)

set(WRAPPER_LIBRARY_GROUPS
  ITKIOBase
)

set(WRAPPER_LIBRARY_DEPENDS ITKCommon  ITKStatistics ITKTransform ITKImageFunction ITKImageFilterBase ITKImageGrid ITKFiniteDifference ITKCurvatureFlow ITKAnisotropicSmoothing ITKMesh ITKIO-Base ITKSpatialObjects ITKImageCompose ITKImageStatistics ITKPath ITKImageIntensity ITKSmoothing ITKImageGradient ITKImageFeature ITKThresholding ITKConnectedComponents)

set(WRAPPER_LIBRARY_LINK_LIBRARIES ${ITK_LIBRARIES})

AUTO_INCLUDE_MODULES()
END_WRAP_LIBRARY()
