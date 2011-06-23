WRAP_LIBRARY(ITKSmoothing)

set(WRAPPER_LIBRARY_GROUPS
)

set(WRAPPER_LIBRARY_DEPENDS ITKCommon  ITKStatistics ITKTransform ITKImageFunction ITKImageFilterBase ITKImageGrid ITKFiniteDifference ITKCurvatureFlow ITKAnisotropicSmoothing ITKMesh ITKIOBase ITKSpatialObjects ITKImageCompose ITKImageStatistics ITKPath ITKImageIntensity)

set(WRAPPER_LIBRARY_LINK_LIBRARIES ${ITK_LIBRARIES})

AUTO_INCLUDE_MODULES()

END_WRAP_LIBRARY()
