#include "../itk_test/readUti.h"
int main(int argc,char ** argv)
{
	auto image=readFloatImage("D:/liver segmentation/data/distance/1.mhd");
	IteratorTypeFloat iter(image,image->GetRequestedRegion());
	iter.GoToBegin();
	int a = 0;
	while (a++ < 5)
	{
		std::cout<<iter.Get()<<std::endl;
		iter++;
	}
	return 0;
}