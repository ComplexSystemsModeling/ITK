/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkDisplacementFieldToBSplineImageFilter_h
#define __itkDisplacementFieldToBSplineImageFilter_h

#include "itkImageToImageFilter.h"

#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkPointSet.h"
#include "itkVector.h"

namespace itk
{

/**
 * \class DisplacementFieldToBSplineImageFilter
 * \brief Class which takes a displacement field image and smooths it
 * using B-splines.  The inverse can also be estimated.
 *
 * \author Nick Tustison
 *
 * \ingroup ITKDisplacementField
 */

template <class TInputImage, class TOutputImage = TInputImage>
class DisplacementFieldToBSplineImageFilter
  : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  typedef DisplacementFieldToBSplineImageFilter            Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage>    Superclass;
  typedef SmartPointer<Self>                               Pointer;
  typedef SmartPointer<const Self>                         ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Extract dimension from input image. */
  itkStaticConstMacro( ImageDimension, unsigned int, TInputImage::ImageDimension );

  typedef TInputImage                          InputFieldType;
  typedef TOutputImage                         OutputFieldType;

  typedef InputFieldType                       DisplacementFieldType;
  typedef OutputFieldType                      InverseDisplacementFieldType;

  /** Image typedef support. */
  typedef typename OutputFieldType::PixelType     PixelType;
  typedef typename OutputFieldType::PixelType     VectorType;
  typedef typename OutputFieldType::RegionType    RegionType;
  typedef typename OutputFieldType::IndexType     IndexType;

  typedef typename OutputFieldType::PointType     PointType;
  typedef typename OutputFieldType::SpacingType   SpacingType;
  typedef typename OutputFieldType::PointType     OriginType;
  typedef typename OutputFieldType::SizeType      SizeType;
  typedef typename OutputFieldType::DirectionType DirectionType;

  /** B-spline smoothing filter argument typedefs */
  typedef PointSet<VectorType, ImageDimension>    PointSetType;

  /** B-sline filter typedefs */
  typedef BSplineScatteredDataPointSetToImageFilter<
    PointSetType, OutputFieldType>                          BSplineFilterType;
  typedef typename BSplineFilterType::WeightsContainerType  WeightsContainerType;
  typedef typename BSplineFilterType::PointDataImageType    DisplacementFieldControlPointLatticeType;
  typedef typename BSplineFilterType::ArrayType             ArrayType;

  /** Set the displacement field */
  void SetDisplacementField( const InputFieldType * field )
    {
    this->SetInput( field );
    }

  /**
   * Get the deformation field.
   */
  const InputFieldType* GetDisplacementField() const
    {
    return this->GetInput( 0 );
    }

  /**
   * Get the displacement field control point lattice.
   */
  const DisplacementFieldControlPointLatticeType * GetDisplacementFieldControlPointLattice() const
    {
    return this->m_DisplacementFieldControlPointLattice;
    }

  /**
   * Set the spline order defining the bias field estimate.  Default = 3.
   */
  itkSetMacro( SplineOrder, unsigned int );

  /**
   * Get the spline order defining the bias field estimate.  Default = 3.
   */
  itkGetConstMacro( SplineOrder, unsigned int );

  /**
   * Set the control point grid size definining the B-spline estimate of the
   * scalar bias field.  In each dimension, the B-spline mesh size is equal
   * to the number of control points in that dimension minus the spline order.
   * Default = 4 control points in each dimension for a mesh size of 1 in each
   * dimension.
   */
  itkSetMacro( NumberOfControlPoints, ArrayType );

  /**
   * Get the control point grid size definining the B-spline estimate of the
   * scalar bias field.  In each dimension, the B-spline mesh size is equal
   * to the number of control points in that dimension minus the spline order.
   * Default = 4 control points in each dimension for a mesh size of 1 in each
   * dimension.
   */
  itkGetConstMacro( NumberOfControlPoints, ArrayType );

  /**
   * Set the number of fitting levels.  One of the contributions of N4 is the
   * introduction of a multi-scale approach to fitting. This allows one to
   * specify a B-spline mesh size for initial fitting followed by a doubling of
   * the mesh resolution for each subsequent fitting level.  Default = 1 level.
   */
  itkSetMacro( NumberOfFittingLevels, ArrayType );

  /**
   * Set the number of fitting levels.  One of the contributions of N4 is the
   * introduction of a multi-scale approach to fitting. This allows one to
   * specify a B-spline mesh size for initial fitting followed by a doubling of
   * the mesh resolution for each subsequent fitting level.  Default = 1 level.
   */
  void SetNumberOfFittingLevels( unsigned int n )
    {
    ArrayType nlevels;

    nlevels.Fill( n );
    this->SetNumberOfFittingLevels( nlevels );
    }

  /**
   * Get the number of fitting levels.  One of the contributions of N4 is the
   * introduction of a multi-scale approach to fitting. This allows one to
   * specify a B-spline mesh size for initial fitting followed by a doubling of
   * the mesh resolution for each subsequent fitting level.  Default = 1 level.
   */
  itkGetConstMacro( NumberOfFittingLevels, ArrayType );

  /**
   * Estimate the inverse field instead of the forward field.  Default = false.
   */
  itkBooleanMacro( EstimateInverse );
  itkSetMacro( EstimateInverse, bool );
  itkGetConstMacro( EstimateInverse, bool );

  /**
   * Enforce stationary boundary conditions.  Default = false.
   */
  itkBooleanMacro( EnforceStationaryBoundary );
  itkSetMacro( EnforceStationaryBoundary, bool );
  itkGetConstMacro( EnforceStationaryBoundary, bool );

protected:

  /** Constructor */
  DisplacementFieldToBSplineImageFilter();

  /** Deconstructor */
  virtual ~DisplacementFieldToBSplineImageFilter();

  /** Standard print self function **/
  void PrintSelf( std::ostream& os, Indent indent ) const;

  /** preprocessing function */
  void GenerateData();

private:
  DisplacementFieldToBSplineImageFilter( const Self& ); //purposely not implemented
  void operator=( const Self& );                 //purposely not implemented

  bool                                         m_EstimateInverse;
  bool                                         m_EnforceStationaryBoundary;
  unsigned int                                 m_SplineOrder;
  ArrayType                                    m_NumberOfControlPoints;
  ArrayType                                    m_NumberOfFittingLevels;

  typename DisplacementFieldControlPointLatticeType::Pointer    m_DisplacementFieldControlPointLattice;


};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDisplacementFieldToBSplineImageFilter.hxx"
#endif

#endif
