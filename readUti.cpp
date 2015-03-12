#include "readUti.h"


PcaData *PcaData::pcaData = NULL;
FloatImageType *PcaData::ratioImage = NULL;
PcaData::PcaData()
{
	std::string pre="D:/liver segmentation/data/distance3/";
	std::string post=".mhd";

	this->mean = readFloatImage((pre+"mean.mhd").c_str());

	for(int index = 1;index<=princNum;++index)
	{
		this->princModeArray[index-1] = readFloatImage(getNumOfFile(pre,post,index).c_str());
	}

}
bool PcaData::loadRatioImage(string ratioName)
{
	PcaData::ratioImage = readFloatImage(ratioName.c_str());
	return true;
}
bool PcaData::deleteModel()
{
	pcaData->mean->Delete();
	for(int index = 0;index<princNum;++index)
		pcaData->princModeArray[index]->Delete();
	pcaData=NULL;
	return 1;
}
string getNumOfFile(string pre,string post,int num)
{
	char numStr[10];
	itoa(num,numStr,10);
	return pre+numStr+post;
}
SegedImageType::Pointer readSegedImage(const char* fileName)
{
	SegedReaderType::Pointer reader = SegedReaderType::New();
	reader->SetFileName( fileName);

	typedef itk::MetaImageIO ImageIOType;
	ImageIOType::Pointer mhdImageIO = ImageIOType::New();
	reader->SetImageIO( mhdImageIO);


	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "exception in file reader " << std::endl;
		std::cerr << e << std::endl;
		return NULL;
	}
	return reader->GetOutput();

}
bool writeSegedImage(const char* fileName,SegedImageType::Pointer image)
{
	typedef itk::MetaImageIO ImageIOType;
	ImageIOType::Pointer mhdImageIO = ImageIOType::New();

	SegedWriterType::Pointer writer = SegedWriterType::New();
	writer->SetFileName(fileName);
	writer->SetImageIO(mhdImageIO);
	writer->SetInput(image);
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "exception in file reader " << std::endl;
		std::cerr << e << std::endl;
		return false;
	}
	return true;
}

OriginImageType::Pointer readOriginImage(const char* fileName)
{
	OriginReaderType::Pointer reader = OriginReaderType::New();
	reader->SetFileName( fileName);

	typedef itk::MetaImageIO ImageIOType;
	ImageIOType::Pointer mhdImageIO = ImageIOType::New();
	reader->SetImageIO( mhdImageIO);


	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "exception in file reader " << std::endl;
		std::cerr << e << std::endl;
		return NULL;
	}
	return reader->GetOutput();

}
bool writeOriginImage(const char* fileName,OriginImageType::Pointer image)
{
	typedef itk::MetaImageIO ImageIOType;
	ImageIOType::Pointer mhdImageIO = ImageIOType::New();

	OriginWriterType::Pointer writer = OriginWriterType::New();
	writer->SetFileName(fileName);
	writer->SetImageIO(mhdImageIO);
	writer->SetInput(image);
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "exception in file reader " << std::endl;
		std::cerr << e << std::endl;
		return false;
	}
	return true;
}



FloatImageType::Pointer readFloatImage(const char* fileName)
{
	FloatReaderType::Pointer reader = FloatReaderType::New();
	reader->SetFileName( fileName);

	typedef itk::MetaImageIO ImageIOType;
	ImageIOType::Pointer mhdImageIO = ImageIOType::New();
	reader->SetImageIO( mhdImageIO);


	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "exception in file reader " << std::endl;
		std::cerr << e << std::endl;
		return NULL;
	}
	auto re = reader->GetOutput();
	reader->Register();
	mhdImageIO->Register();

	re->Update();
	return re;
}
bool writeFloatImage(const char* fileName,FloatImageType::Pointer image)
{
	typedef itk::MetaImageIO ImageIOType;
	ImageIOType::Pointer mhdImageIO = ImageIOType::New();

	FloatWriterType::Pointer writer = FloatWriterType::New();
	writer->SetFileName(fileName);
	writer->SetImageIO(mhdImageIO);
	writer->SetInput(image);
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "exception in file reader " << std::endl;
		std::cerr << e << std::endl;
		return false;
	}
	return true;
}

OriginImageType* FloatImageToOrigin(FloatImageType* floatImage)
{
	auto origin = OriginImageType::New();
	origin->SetRegions(floatImage->GetRequestedRegion());
	origin->SetSpacing(floatImage->GetSpacing());
	origin->Allocate();
	IteratorTypeOrigin iterOri(origin,origin->GetRequestedRegion());
	IteratorTypeFloat iterFloat(floatImage,floatImage->GetRequestedRegion());
	iterFloat.GoToBegin();
	iterOri.GoToBegin();
	while (!iterFloat.IsAtEnd())
	{
		if(iterFloat.Get()>0)
			iterOri.Set(100);
		else
		{
			iterOri.Set(0);
		}
		++iterFloat;
		++iterOri;
	}
	origin->Register();
	return origin;
}
OriginImageType* levelsetFloatImageToOrigin(FloatImageType* floatImage)
{
	auto origin = OriginImageType::New();
	origin->SetRegions(floatImage->GetRequestedRegion());
	origin->SetSpacing(floatImage->GetSpacing());
	origin->Allocate();
	IteratorTypeOrigin iterOri(origin,origin->GetRequestedRegion());
	IteratorTypeFloat iterFloat(floatImage,floatImage->GetRequestedRegion());
	iterFloat.GoToBegin();
	iterOri.GoToBegin();
	while (!iterFloat.IsAtEnd())
	{
		if(iterFloat.Get()<0)
			iterOri.Set(100);
		else
		{
			iterOri.Set(0);
		}
		++iterFloat;
		++iterOri;
	}
	origin->Register();
	//floatImage->Delete();
	return origin;
}

FloatImageType* cloneEmptyImageFromFile()
{

	auto origin = FloatImageType::New();
	FloatImageType::RegionType region;
	region.SetSize(0,512);
	region.SetSize(1,512);
	region.SetSize(2,178);
	origin->SetRegions(region);
	FloatImageType::SpacingType spacing;
	spacing.SetElement(0,0.74219);
	spacing.SetElement(1,0.74219);
	spacing.SetElement(2,1.5);
	origin->SetSpacing(spacing);
	origin->Allocate(1);

	origin->Register();
	return origin;
}










//vtk




void writePolyData(vtkPolyData* poly,string pathName)
{
	auto polyWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
	polyWriter->SetFileName(pathName.c_str());
	polyWriter->SetInputData(poly);
	polyWriter->Update();

}

vtkPolyData* readPolyData(string pathName)
{
	auto polyReader = vtkPolyDataReader::New();
	polyReader->SetFileName(pathName.c_str());

	polyReader->Update();
	auto poly = polyReader->GetOutput();
	auto num = poly->GetNumberOfPoints();
	return poly;
}
vector<int> getCorArray(vtkPolyData* src,vtkPolyData* des,vector<int> srcCor)
{
	vector<int> curCor;
	auto locator = vtkSmartPointer<vtkCellLocator>::New();
	locator->SetDataSet(des);
	locator->SetNumberOfCellsPerBucket(1);
	locator->BuildLocator();
	double outPoint[3];
	int cell_id;
	int sub_id;
	double dist2;
	for(vector<int>::iterator iter = srcCor.begin();iter != srcCor.end();++iter)
	{
		int index = *iter;
		locator->FindClosestPoint(src->GetPoint(index),
			outPoint,
			cell_id,
			sub_id,
			dist2);
		curCor.push_back(sub_id);
	}
	return curCor;

}



vtkPolyData* getPolyDataLoadImage(vtkImageData* imageData)
{
	auto resample = vtkSmartPointer<vtkImageResample>::New();
	int reductionFactor = 0.25;
	resample->SetAxisOutputSpacing(0,imageData->GetSpacing()[0]*4);
	resample->SetAxisOutputSpacing(1,imageData->GetSpacing()[1]*4);
	resample->SetAxisOutputSpacing(2,imageData->GetSpacing()[2]*4);

	resample->SetInputData( imageData );
	resample->SetAxisMagnificationFactor(0, reductionFactor);
	resample->SetAxisMagnificationFactor(1, reductionFactor);
	resample->SetAxisMagnificationFactor(2, reductionFactor);

	resample->Update();
	auto imageDataPost = resample->GetOutput();
	//cout<<imageData->GetDimensions()[0]<<" "<<imageData->GetDimensions()[1]<<" "<<imageData->GetDimensions()[2]<<endl;
	//cout<<imageDataPost->GetDimensions()[0]<<" "<<imageDataPost->GetDimensions()[1]<<" "<<imageDataPost->GetDimensions()[2]<<endl;
	vtkSmartPointer<vtkImageCast> castFilter = 
		vtkSmartPointer<vtkImageCast>::New();
	/*用filter处理把有符号的short转换为无符号的short*/
	castFilter->SetOutputScalarTypeToUnsignedInt();
	castFilter->SetInputData(imageDataPost);
	castFilter->Update();
	imageData= castFilter->GetOutput();
	auto image = imageData;
	//auto iso = vtkSmartPointer<vtkMarchingCubes>::New();  
	auto iso = vtkMarchingCubes::New();
	iso->SetInputData(image);  
	iso->SetNumberOfContours(1);  
	iso->SetValue(0,50);  
	iso->ComputeGradientsOn();  
	iso->ComputeNormalsOff();
	//iso->ComputeNormalsOn();
	//iso->ComputeScalarsOff();  
	iso->ComputeScalarsOn();
	iso->Update();

	auto polyData = iso->GetOutput();
	return polyData;
}
void icpPoly(vtkPolyData*& polySrc,vtkPolyData*& polyDes)
{
	vtkSmartPointer<vtkIterativeClosestPointTransform> icp = 
		vtkSmartPointer<vtkIterativeClosestPointTransform>::New();
	icp->SetSource(polySrc);
	icp->SetTarget(polyDes);
	icp->GetLandmarkTransform()->SetModeToRigidBody();
	icp->SetMaximumNumberOfLandmarks(20000);

	icp->SetMaximumNumberOfIterations(300);

	icp->StartByMatchingCentroidsOn();
	icp->Modified();
	icp->Update();

	//get the resulting transformation matrix (this matrix takes the source points to the target points)
	vtkSmartPointer<vtkMatrix4x4> M = icp->GetMatrix();

	//std::cout << "The resulting matrix is: " << *M << std::endl;

	auto trans = vtkTransform::New();
	trans->SetMatrix(M);
	trans->Update();
	//std::cout << "The trans is: " << *trans << std::endl;
	auto transFilter = vtkTransformPolyDataFilter::New();
	transFilter->SetInputData(polySrc);
	transFilter->SetTransform(trans);
	transFilter->Update();
	polySrc = transFilter->GetOutput();
}

vtkPolyData* getPolyDataLoadFile(string filePath)
{
	vtkSmartPointer<vtkMetaImageReader> reader = 
		vtkSmartPointer<vtkMetaImageReader>::New();
	reader->SetFileName(filePath.c_str());
	reader->Update();
	auto imageData = reader->GetOutput();
	//resample
	auto resample = vtkSmartPointer<vtkImageResample>::New();
	int reductionFactor = 0.25;
	resample->SetAxisOutputSpacing(0,imageData->GetSpacing()[0]*4);
	resample->SetAxisOutputSpacing(1,imageData->GetSpacing()[1]*4);
	resample->SetAxisOutputSpacing(2,imageData->GetSpacing()[2]*4);

	resample->SetInputData( imageData );
	resample->SetAxisMagnificationFactor(0, reductionFactor);
	resample->SetAxisMagnificationFactor(1, reductionFactor);
	resample->SetAxisMagnificationFactor(2, reductionFactor);

	resample->Update();
	auto imageDataPost = resample->GetOutput();
	//cout<<imageData->GetDimensions()[0]<<" "<<imageData->GetDimensions()[1]<<" "<<imageData->GetDimensions()[2]<<endl;
	//cout<<imageDataPost->GetDimensions()[0]<<" "<<imageDataPost->GetDimensions()[1]<<" "<<imageDataPost->GetDimensions()[2]<<endl;
	vtkSmartPointer<vtkImageCast> castFilter = 
		vtkSmartPointer<vtkImageCast>::New();
	/*用filter处理把有符号的short转换为无符号的short*/
	castFilter->SetOutputScalarTypeToUnsignedInt();
	castFilter->SetInputData(imageDataPost);
	castFilter->Update();
	imageData= castFilter->GetOutput();
	auto image = imageData;
	//auto iso = vtkSmartPointer<vtkMarchingCubes>::New();  
	auto iso = vtkMarchingCubes::New();
	iso->SetInputData(image);  
	iso->SetNumberOfContours(1);  
	iso->SetValue(0,50);  
	iso->ComputeGradientsOn();  
	iso->ComputeNormalsOff();
	//iso->ComputeNormalsOn();
	//iso->ComputeScalarsOff();  
	iso->ComputeScalarsOn();
	iso->Update();

	auto polyData = iso->GetOutput();
	return polyData;
}
vtkRenderer* getRenderer(vtkPolyData* polyData)
{
	auto cubeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();  

	cubeMapper->SetInputData(polyData); 
	auto cubeActor = vtkSmartPointer<vtkActor>::New();  
	cubeActor->SetMapper(cubeMapper);  
	auto camera = vtkSmartPointer<vtkCamera>::New();  
	camera->SetPosition(1,1,1);  
	camera->SetFocalPoint(0,0,0);  
	//auto renderer = vtkSmartPointer<vtkRenderer>::New();
	vtkRenderer *renderer = vtkRenderer::New();
	renderer->AddActor(cubeActor);  
	renderer->SetActiveCamera(camera);  
	renderer->ResetCamera();  
	renderer->SetBackground(1,1,1);  
	return renderer;
}
//vtk and itk
vtkImageData* itkFloatImageToVtk(FloatImageType* imageItk)
{
	auto filter = itk::ImageToVTKImageFilter<FloatImageType>::New();
	//filter->Register();
	filter->SetInput(imageItk);
	filter->Update();
	auto image=filter->GetOutput();
	image->Register(image);
	return filter->GetOutput();
}

vtkImageData* itkOriginImageToVtk(OriginImageType* imageItk)
{
	auto filter = itk::ImageToVTKImageFilter<OriginImageType>::New();
	//filter->Register();
	filter->SetInput(imageItk);
	filter->Update();
	auto image=filter->GetOutput();
	image->Register(image);
	return filter->GetOutput();
}

//pca
FloatImageType* getImageMode(FloatImageType* mean,FloatImageType* princ,double para)
{
	mean = PcaData::getInstance()->mean;
	princ = PcaData::getInstance()->princModeArray[0];



	auto addFilter = itk::AddImageFilter<FloatImageType,FloatImageType,FloatImageType>::New();
	auto multiplyFilter = itk::MultiplyImageFilter<FloatImageType,FloatImageType,FloatImageType>::New();
	multiplyFilter->SetInput1(princ);
	multiplyFilter->SetConstant2(para);
	multiplyFilter->Update();
	addFilter->SetInput1(mean);
	addFilter->SetInput2(multiplyFilter->GetOutput());
	addFilter->Update();
	auto result = addFilter->GetOutput();
	result->Register();
	return result;

}

FloatImageType* getImageModePca(double* para)
{

	auto mean = PcaData::getInstance()->mean;
	auto princ = PcaData::getInstance()->princModeArray;
	auto pricNum = PcaData::getInstance()->princNum;


	for(int index = 0;index <pricNum;++index)
	{
		auto addFilter = itk::AddImageFilter<FloatImageType,FloatImageType,FloatImageType>::New();
		auto multiplyFilter = itk::MultiplyImageFilter<FloatImageType,FloatImageType,FloatImageType>::New();
		multiplyFilter->SetInput1(princ[index]);
		multiplyFilter->SetConstant2(para[index]);
		multiplyFilter->Update();
		addFilter->SetInput1(mean);
		addFilter->SetInput2(multiplyFilter->GetOutput());
		addFilter->Update();
		auto temp = mean;
		mean = addFilter->GetOutput();
		mean->Register();
		if(index>0)
			temp->Delete();

	}

	return mean;

}
FloatImageType* transformRotate(FloatImageType* image,double *angle)
{

	auto filter = ResampleImageFilterType::New();
	auto interpolator = LinearInterpolatorType::New();
	filter->SetInterpolator(interpolator);
	filter->SetDefaultPixelValue(1);
	auto trans = TransformType::New();
	auto spacing = image->GetSpacing();
	auto origin  = image->GetOrigin();
	auto size =
		image->GetLargestPossibleRegion().GetSize();
	auto
		direction  = image->GetDirection();
	filter->SetOutputOrigin( origin );
	filter->SetOutputSpacing( spacing );
	filter->SetOutputDirection( direction );
	filter->SetSize( size );
	TransformType::OutputVectorType translation1;
	TransformType::OutputVectorType translation2;
	TransformType::CenterType center;
	center[0] = origin[0] + spacing[0] * size[0] / 2.0;
	center[1] = origin[1] + spacing[1] * size[1] / 2.0;
	center[2] = origin[2] + spacing[2] * size[2] / 2.0;
	trans->SetCenter(center);
	TransformType::OutputVectorType axieX;
	axieX[0]=1;
	axieX[1]=0;
	axieX[2]=0;
	trans->Rotate3D(axieX,angle[0],false);
	axieX[0]=0;
	axieX[1]=1;
	axieX[2]=0;
	trans->Rotate3D(axieX,angle[1],false);
	axieX[0]=0;
	axieX[1]=0;
	axieX[2]=1;
	trans->Rotate3D(axieX,angle[2],false);
	filter->SetTransform(trans);
	filter->SetInput(image);
	filter->Update();
	auto re = filter->GetOutput();
	re->Register();
	return re;
}

FloatImageType* transformScaling(FloatImageType* image,double scale,TransformType::CenterType center)
{

	auto filter = ResampleImageFilterType::New();
	auto interpolator = LinearInterpolatorType::New();
	filter->SetInterpolator(interpolator);
	filter->SetDefaultPixelValue(1);
	auto trans = TransformType::New();
	auto spacing = image->GetSpacing();
	auto origin  = image->GetOrigin();
	auto size =
		image->GetLargestPossibleRegion().GetSize();
	auto
		direction  = image->GetDirection();
	filter->SetOutputOrigin( origin );
	filter->SetOutputSpacing( spacing );
	filter->SetOutputDirection( direction );
	filter->SetSize( size );



	trans->SetCenter(center);

	trans->Scale(scale,false);


	filter->SetTransform(trans);
	filter->SetInput(image);
	filter->Update();
	auto re = filter->GetOutput();
	re->Register();
	return re;
}
FloatImageType* transform(FloatImageType* image,double *angle,double scale,TransformType::CenterType center)
{
	auto filter = ResampleImageFilterType::New();
	auto interpolator = LinearInterpolatorType::New();
	filter->SetInterpolator(interpolator);
	filter->SetDefaultPixelValue(1);
	auto trans = TransformType::New();
	auto spacing = image->GetSpacing();
	auto origin  = image->GetOrigin();
	auto size =
		image->GetLargestPossibleRegion().GetSize();
	auto
		direction  = image->GetDirection();
	filter->SetOutputOrigin( origin );
	filter->SetOutputSpacing( spacing );
	filter->SetOutputDirection( direction );
	filter->SetSize( size );



	trans->SetCenter(center);

	trans->Scale(scale,false);
	TransformType::OutputVectorType axieX;
	axieX[0]=1;
	axieX[1]=0;
	axieX[2]=0;
	trans->Rotate3D(axieX,angle[0],false);
	axieX[0]=0;
	axieX[1]=1;
	axieX[2]=0;
	trans->Rotate3D(axieX,angle[1],false);
	axieX[0]=0;
	axieX[1]=0;
	axieX[2]=1;
	trans->Rotate3D(axieX,angle[2],false);

	filter->SetTransform(trans);
	filter->SetInput(image);
	filter->Update();
	auto re = filter->GetOutput();
	re->Register();
	return re;
}
double	 transformForFitness(double* para,double *angle,int* pos,double scale)
{
	auto image = getImageModePca(para);
	auto spacing = image->GetSpacing();
	auto origin  = image->GetOrigin();
	auto size =
		image->GetLargestPossibleRegion().GetSize();
	TransformType::CenterType center;
	center[0] = origin[0] + spacing[0] * size[0] / 2.0;
	center[1] = origin[1] + spacing[1] * size[1] / 2.0;
	center[2] = origin[2] + spacing[2] * size[2] / 2.0;
	center[0] +=pos[0];
	center[1] +=pos[1];
	center[2] +=pos[2];

	PcaData::deleteModel();
	//提前分配空间给几个大图像
	auto emptyImage = cloneEmptyImageFromFile();

	FloatImageType::SizeType size0 = image->GetLargestPossibleRegion().GetSize();
	auto size1 = emptyImage->GetLargestPossibleRegion().GetSize();
	for(int index0 =0;index0<size0.GetElement(0);index0++)
	{
		for(int index1=0;index1<size0.GetElement(1);index1++)
		{
			for(int index2=0;index2<size0.GetElement(2);index2++)
			{
				FloatImageType::IndexType index ;
				index[0] = index0;
				index[1] = index1;
				index[2] = index2;
				auto value = image->GetPixel(index);
				index[0]+=pos[0];
				index[1]+=pos[1];
				index[2]+=pos[2];
				if(index[0]>= size1.GetElement(0) || index[1]>= size1.GetElement(1) || index[2]>= size1.GetElement(2))
					continue;
				emptyImage->SetPixel(index,value);
			}
		}
	}


	//放缩
	auto result =transform(emptyImage,angle,scale,center);
	//result->Register();
	if(!PcaData::ratioImage)
	{
		cout<<"ratioImage is empty!";
		return 0;
	}
	auto ratio = PcaData::ratioImage;
	IteratorTypeFloat ratioIter(ratio,ratio->GetRequestedRegion());
	IteratorTypeFloat resultIter(result,ratio->GetRequestedRegion());
	ratioIter.GoToBegin();
	resultIter.GoToBegin();
	double fitness= 0;
	while (!ratioIter.IsAtEnd())
	{
		if(resultIter.Get()<0)
			fitness+=ratioIter.Get();
		++ratioIter;
		++resultIter;
	}
	writeFloatImage("D:/liver segmentation/data/segResultPcaFloat.mhd",result);
	auto originImage = levelsetFloatImageToOrigin(result);
	writeOriginImage("D:/liver segmentation/data/segResultPca.mhd",originImage);

	image->Delete();
	emptyImage->Delete();
	originImage->Delete();
	result->Delete();

	return fitness;



}
double	 transformForFitnessWithSave(double* para,double *angle,int* pos,double scale)
{
	auto image = getImageModePca(para);





	auto spacing = image->GetSpacing();
	auto origin  = image->GetOrigin();
	auto size =
		image->GetLargestPossibleRegion().GetSize();
	TransformType::CenterType center;
	center[0] = origin[0] + spacing[0] * size[0] / 2.0;
	center[1] = origin[1] + spacing[1] * size[1] / 2.0;
	center[2] = origin[2] + spacing[2] * size[2] / 2.0;

	//把旋转后的图像放入要分割图像大小的空图像中


	auto re = transformRotate(image,angle);
	//image->Delete();
	//return re;
	//释放模型占据的内存
	PcaData::deleteModel();
	//提前分配空间给几个大图像
	auto emptyImage = cloneEmptyImageFromFile();
	FloatImageType::SizeType size0 = re->GetLargestPossibleRegion().GetSize();
	auto size1 = emptyImage->GetLargestPossibleRegion().GetSize();
	for(int index0 =0;index0<size0.GetElement(0);index0++)
		for(int index1=0;index1<size0.GetElement(1);index1++)
			for(int index2=0;index2<size0.GetElement(2);index2++)
			{
				FloatImageType::IndexType index ;
				index[0] = index0;
				index[1] = index1;
				index[2] = index2;
				auto value = re->GetPixel(index);
				index[0]+=pos[0];
				index[1]+=pos[1];
				index[2]+=pos[2];
				if(index[0]>= size1.GetElement(0) || index[1]>= size1.GetElement(1) || index[2]>= size1.GetElement(2))
					continue;
				emptyImage->SetPixel(index,value);
			}
			//emptyImage->Register();

			center[0] +=pos[0];
			center[1] +=pos[1];
			center[2] +=pos[2];
			//放缩
			auto result =transformScaling(emptyImage,scale,center);
			//result->Register();
			if(!PcaData::ratioImage)
			{
				cout<<"ratioImage is empty!";
				return 0;
			}
			auto ratio = PcaData::ratioImage;
			IteratorTypeFloat ratioIter(ratio,ratio->GetRequestedRegion());
			IteratorTypeFloat resultIter(result,ratio->GetRequestedRegion());
			ratioIter.GoToBegin();
			resultIter.GoToBegin();
			double fitness= 0;
			while (!ratioIter.IsAtEnd())
			{
				if(resultIter.Get()<0)
					fitness+=ratioIter.Get();
				++ratioIter;
				++resultIter;
			}
			writeFloatImage("D:/liver segmentation/data/segResultPcaFloat.mhd",result);
			auto originImage = levelsetFloatImageToOrigin(result);
			writeOriginImage("D:/liver segmentation/data/segResultPca.mhd",originImage);
			
			image->Delete();
			emptyImage->Delete();
			originImage->Delete();
			result->Delete();
			re->Delete();
			return fitness;



}
