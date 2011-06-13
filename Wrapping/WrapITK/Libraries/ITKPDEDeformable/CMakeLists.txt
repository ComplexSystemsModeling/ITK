WRAP_LIBRARY(itkPDEDeformable)

set(WRAPPER_LIBRARY_DEPENDS
  VXLNumerics
  ITKCommon
  Numerics
  Transforms
  SpatialObject
  Resize
  Calculators
)

set(WRAPPER_LIBRARY_LINK_LIBRARIES ${ITK_LIBRARIES})

AUTO_INCLUDE_MODULES()
END_WRAP_LIBRARY()
