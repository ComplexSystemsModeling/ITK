WRAP_LIBRARY(ITKPath)

set(WRAPPER_LIBRARY_GROUPS
  itkPathBase
  itkParametricPath
  itkPolyLineParametricPath
)

set(WRAPPER_LIBRARY_DEPENDS ITKCommon  ITKStatistics ITKTransform ITKImageFunction ITKImageFilterBase ITKImageGrid ITKFiniteDifference ITKCurvatureFlow ITKAnisotropicSmoothing ITKMesh ITKIOBase ITKSpatialObjects ITKImageCompose ITKImageStatistics)

set(WRAPPER_LIBRARY_LINK_LIBRARIES ${ITK_LIBRARIES})

AUTO_INCLUDE_MODULES()

END_WRAP_LIBRARY()
