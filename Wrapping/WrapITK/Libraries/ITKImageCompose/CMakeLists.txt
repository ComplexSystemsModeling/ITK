WRAP_LIBRARY(ITKImageCompose)

set(WRAPPER_LIBRARY_DEPENDS ITKCommon  ITKStatistics ITKTransform ITKImageFunction ITKImageFilterBase ITKImageGrid ITKFiniteDifference ITKCurvatureFlow ITKAnisotropicSmoothing ITKMesh ITKIOBase ITKSpatialObjects)

set(WRAPPER_LIBRARY_LINK_LIBRARIES ${ITK_LIBRARIES})

AUTO_INCLUDE_MODULES()
END_WRAP_LIBRARY()
