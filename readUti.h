#pragma once
#include "itkImage.h"
#include "itkMetaImageIO.h"
#include <iostream>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBinaryMask3DMeshSource.h"
#include "itkImage.h"
#include "itkMesh.h"
#include <stdlib.h>  
#include "GL/glut.h"
#include <math.h>  
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkImagePCAShapeModelEstimator.h"
#include "itkResampleImageFilter.h" 
#include "itkTranslationTransform.h"
#include "itkNumericTraits.h"
#include "itkCovariantVector.h"
#include "itkIdentityTransform.h"
#include "itkBSplineInterpolateImageFunction.h"
#include <itkImageToVTKImageFilter.h>
#include <vtkSmartPointer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkImageData.h>
#include <vtkImageMapper3D.h>
#include <vtkImageCast.h>
#include <vtkMetaImageWriter.h>
#include <vtkMetaImageReader.h>
#include <vtkImageMandelbrotSource.h>
#include <vtkImageActor.h>
#include <vtkCamera.h>  
#include <vtkVolumeRayCastCompositeFunction.h>
#include <vtkVolumeRayCastMapper.h>
#include<vtkPiecewiseFunction.h> 
#include<vtkVolumeProperty.h>
#include <vtkColorTransferFunction.h>
#include <vtkImageData.h>
#include <vtkMarchingCubes.h>
#include <vtkVectorNorm.h>
#include <vtkDataSetMapper.h>
#include <exception>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkCleanPolyData.h>
#include "vtkImageResample.h"
#include <vtkDelaunay3D.h>
#include <vtkSmartPointer.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkIterativeClosestPointTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkLandmarkTransform.h>
#include <vtkMath.h>
#include <vtkMatrix4x4.h>
#include <fstream>
#include "vtkTransform.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyDataReader.h"
#include "vtkXMLPolyDataWriter.h"
#include <vtkTransformPolyDataFilter.h>
#include <vtkHomogeneousTransform.h>
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkCleanPolyData.h"
#include "vtkCellLocator.h"
#include <vector>
#include "itkImageDuplicator.h"
#include <itkAddImageFilter.h>
#include "itkMultiplyImageFilter.h"
#include "itkCenteredAffineTransform.h"
#include "itkScalarImageToCooccurrenceMatrixFilter.h"
using std::vector;
using std::exception;
using std::string;
using std::string;
const unsigned int InputDimension = 3;
//分割结果图像-类型定义
typedef signed char SegedPixelType;
typedef itk::Image<SegedPixelType,InputDimension> SegedImageType;
typedef itk::ImageFileReader< SegedImageType > SegedReaderType;
typedef itk::ImageFileWriter< SegedImageType > SegedWriterType;
typedef itk::ImageRegionIterator< SegedImageType > IteratorTypeSeged;
//原始图像-类型定义
typedef signed short OriginPixelType;
typedef itk::Image<OriginPixelType,InputDimension> OriginImageType;
typedef itk::ImageFileReader<OriginImageType> OriginReaderType;
typedef itk::ImageFileWriter<OriginImageType> OriginWriterType;
typedef itk::ImageRegionIterator<OriginImageType> IteratorTypeOrigin;
//浮点数图像定义
typedef float FloatPixelType;
typedef itk::Image<FloatPixelType,InputDimension> FloatImageType;
typedef itk::ImageFileReader<FloatImageType> FloatReaderType;
typedef itk::ImageFileWriter<FloatImageType> FloatWriterType;
typedef itk::ImageRegionIterator<FloatImageType> IteratorTypeFloat;
//PCA类型定义
typedef itk::ImagePCAShapeModelEstimator<FloatImageType,   FloatImageType >  PCAestimatorType;
//重采样类型定义
typedef itk::ResampleImageFilter<FloatImageType, FloatImageType> ResampleImageFilterType;
typedef itk::IdentityTransform<double, 3>
	IdentityTransformType;
typedef itk::BSplineInterpolateImageFunction<FloatImageType, double, double>
	BSplineInterpolatorType;
 typedef itk::LinearInterpolateImageFunction<
                       FloatImageType, double >  LinearInterpolatorType;

//transform  类型定义
typedef itk::CenteredAffineTransform< double, InputDimension >  TransformType;
//纹理模型 类型定义
typedef itk::Statistics::ScalarImageToCooccurrenceMatrixFilter<OriginImageType>
    Image2CoOccuranceType;

SegedImageType::Pointer readSegedImage(const char* fileName);
bool writeSegedImage(const char* fileName,SegedImageType::Pointer image);

OriginImageType::Pointer readOriginImage(const char* fileName);
bool writeOriginImage(const char* fileName,OriginImageType::Pointer image);

FloatImageType::Pointer readFloatImage(const char* fileName);
bool writeFloatImage(const char* fileName,FloatImageType::Pointer image);



string getNumOfFile(string pre,string post,int num);
OriginImageType* FloatImageToOrigin(FloatImageType* floatImage);
OriginImageType* levelsetFloatImageToOrigin(FloatImageType* floatImage);

FloatImageType* cloneEmptyImageFromFile();
//vtk
vtkRenderer* getRendererLoadFile(string filePath,bool isSmooth);

void writePolyData(vtkPolyData* poly,string pathName);
vtkPolyData* readPolyData(string pathName);
vector<int> getCorArray(vtkPolyData* src,vtkPolyData* des,vector<int> srcCor);
vtkPolyData* getPolyDataLoadImage(vtkImageData* imageData);

vtkPolyData* getPolyDataLoadFile(string filePath);
vtkRenderer* getRenderer(vtkPolyData* polyData);
void icpPoly(vtkPolyData*& polySrc,vtkPolyData*& polyDes);

//vtk and itk

vtkImageData* itkFloatImageToVtk(FloatImageType* imageItk);
vtkImageData* itkOriginImageToVtk(OriginImageType* imageItk);

//pca
FloatImageType* getImageMode(FloatImageType* mean,FloatImageType* princ,double para);
//const int princNum = 10;
FloatImageType* getImageModePca(double* para);
FloatImageType* transform(FloatImageType* image,double *angle,double scale,TransformType::CenterType center);
double transformForFitness(double* para,double *angle,int* pos,double scale);
double	 transformForFitnessWithSave(double* para,double *angle,int* pos,double scale);
class  PcaData
{
public:
	static const int princNum =10;
	static const int constNum = 233;
	FloatImageType* princModeArray[princNum];
	FloatImageType* mean;
	static PcaData* getInstance()
	{
		if(!pcaData)
			pcaData = new PcaData();
		return pcaData;
	}
	static PcaData* pcaData;
	static FloatImageType* ratioImage;
	static bool loadRatioImage(string ratioName);
	static bool deleteModel();
private:
	
	PcaData();


};
//const int PcaData::princNum;
//PcaData *PcaData::pcaData = NULL;
