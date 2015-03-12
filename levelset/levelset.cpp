#include "../readUti.h"


int main(int argc,char ** argv)
{
	std::string pre="D:/liver segmentation/data/distance2/";
	std::string post=".mhd";
	auto inputImage = readFloatImage("D:/liver segmentation/data/distance2/1.mhd");
	int fileNum =2;
	IteratorTypeFloat iter(inputImage,inputImage->GetRequestedRegion());
	iter.GoToBegin();
	while (fileNum < 21)
	{
		char num[10];
		itoa(fileNum,num,10);
		std::string namew= pre+num+post;
		auto moveImage = readFloatImage(namew.c_str());
		IteratorTypeFloat iterMove(moveImage,moveImage->GetRequestedRegion());
		iterMove.GoToBegin();
		while (!iter.IsAtEnd())
		{
			iter.Set(iter.Get()+iterMove.Get());
			++iter;
			++iterMove;
		}

		fileNum++;
		std::cout<<fileNum<<std::endl;
		iter.GoToBegin();
	}
	iter.GoToBegin();
	while (!iter.IsAtEnd())
	{
		iter.Set(iter.Get()/20);
		++iter;

	}
	writeFloatImage("D:/liver segmentation/data/distance2/mean.mhd",inputImage);
	iter.GoToBegin();
	while (!iter.IsAtEnd())
	{
		if(iter.Get() < 0)
		{
			iter.Set(1);
		}
		else
		{
			iter.Set(0);
		}
		
		++iter;

	}
	writeFloatImage("D:/liver segmentation/data/distance3/meanImage.mhd",inputImage);
	return 0;
}