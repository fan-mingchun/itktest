#include "../itk_test/readUti.h"
int main(int argc,char ** argv)
{
	std::string pre="D:/liver segmentation/data/distance/";
	std::string preWrite="D:/liver segmentation/data/distanceExample/";
	std::string post=".mhd";
	int trainNum = 20;
	int fileNum =1;
	auto resampleFilter =   ResampleImageFilterType::New();
	FloatImageType::SizeType outputSize ;
	outputSize[0] = 50;
	outputSize[1] = 50;
	outputSize[2] = 50;
	resampleFilter->SetSize(outputSize);
	
	auto transform =
		IdentityTransformType::New();
	auto interpolator = BSplineInterpolatorType::New();
	resampleFilter->SetTransform(transform);
	//resampleFilter->SetInterpolator(interpolator);
	while (fileNum <= trainNum)
	{
		char num[10];
		itoa(fileNum,num,10);
		std::string name= pre+num+post;
		auto moveImage = readFloatImage(name.c_str());
		//filter->SetInput(fileNum-1,moveImage);
		
		resampleFilter->SetInput(moveImage);
		resampleFilter->SetOutputSpacing(moveImage->GetSpacing());
		resampleFilter->Update();
		auto writeImage = resampleFilter->GetOutput();

		std::string namew= preWrite+num+post;
		writeFloatImage(namew.c_str(),writeImage);
		fileNum++;

	}
	//filter->Update();
	//auto image = filter->GetOutput(1);
	//writeFloatImage("D:/liver segmentation/data/distance/pcaTest.mhd",image);
	return 0;
}