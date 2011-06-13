WRAP_LIBRARY(ITKDeformationField)

set(WRAPPER_LIBRARY_DEPENDS
  VXLNumerics
  ITKCommon
  Numerics
  itkPDEDeformable
  Transforms
  SpatialObject
  Resize
  Calculators
)

set(WRAPPER_LIBRARY_LINK_LIBRARIES ${ITK_LIBRARIES})

AUTO_INCLUDE_MODULES()
END_WRAP_LIBRARY()
