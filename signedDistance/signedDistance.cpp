#include "../readUti.h"
int main(int argc,char ** argv)
{
	std::string pre="D:/liver segmentation/data/ssmsegLevel2/";
	std::string post=".mhd";
	int fileNum = 20;
	int fileIndex = 15;
	while (fileIndex <= fileNum)
	{
		auto origin = readOriginImage(getNumOfFile(pre,post,fileIndex).c_str());
		auto sign=itk::SignedMaurerDistanceMapImageFilter<OriginImageType,FloatImageType>::New();
		sign->SetInput(origin);
		sign->Update();
		std::string pre="D:/liver segmentation/data/distance2/";
		std::string post=".mhd";
		auto output= sign->GetOutput();
		writeFloatImage(getNumOfFile(pre,post,fileIndex).c_str(),output);
		++fileIndex;
	}
	
	
	return 0;
}