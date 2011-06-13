WRAP_LIBRARY(ITKNonunitReview)

set(WRAPPER_LIBRARY_GROUPS
  itkImageToPathFilter
  itkLabelObjectLine
  itkLabelObject
  itkShapeLabelObject
  itkStatisticsLabelObject
  itkLabelMap
  ITKLabelMapBase
  itkLabelMapFilter
  itkInPlaceLabelMapFilter
  itkChangeRegionLabelMapFilter
)

set(WRAPPER_LIBRARY_DEPENDS
  VXLNumerics
  ITKCommon
  itkPDEDeformable
  ITKTransform
  ITKIOSpatialObjects
)

set(WRAPPER_LIBRARY_LINK_LIBRARIES ${ITK_LIBRARIES})

AUTO_INCLUDE_MODULES()
END_WRAP_LIBRARY()
