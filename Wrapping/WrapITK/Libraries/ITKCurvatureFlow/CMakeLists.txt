WRAP_LIBRARY(ITKCurvatureFlow)

set(WRAPPER_LIBRARY_GROUPS
  itkCurvatureFlowImageFilter
  itkMinMaxCurvatureFlowImageFilter
)

set(WRAPPER_LIBRARY_DEPENDS ITKCommon  ITKStatistics ITKTransform ITKImageFunction ITKImageFilterBase ITKImageGrid ITKFiniteDifference)

set(WRAPPER_LIBRARY_LINK_LIBRARIES ${ITK_LIBRARIES})

AUTO_INCLUDE_MODULES()
END_WRAP_LIBRARY()
