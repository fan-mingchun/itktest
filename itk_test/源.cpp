#include "../readUti.h"
int main(int argc,char ** argv)
{
	//read image

	//auto imageFixed = readSegedImage("D:/liver segmentation/data/newTraining/liver-seg001.mhd");
	//IteratorTypeSeged iterFixed(imageFixed,imageFixed->GetRequestedRegion());

	//auto tempFixed = SegedImageType::New();
	//tempFixed ->CopyInformation(imageFixed);
	//tempFixed->SetRegions(imageFixed->GetRequestedRegion());
	//tempFixed->Allocate(1);
	//IteratorTypeSeged iterTempFixed(tempFixed,tempFixed->GetRequestedRegion());
	//iterFixed.GoToBegin();
	//iterTempFixed.GoToBegin();
	//while(!iterFixed.IsAtEnd())
	//{
	//	if(iterFixed.Get()>0)
	//	{
	//		iterFixed.Set(1);
	//		iterTempFixed.Set(1);
	//	}
	//		
	//	++iterFixed;
	//	++iterTempFixed;
	//}
	//
	std::string namePre = "D:/liver segmentation/data/testtrans/t";
	std::string namePost = "/result.0.mhd";
	int fileNum = 1;
	while (fileNum <2)
	{
		char num[10];
		itoa(fileNum,num,10);
		std::string name = namePre+num+namePost;
		
		auto imageMove = readOriginImage(name.c_str());
		auto imageSave = OriginImageType::New();
		imageSave->SetRegions(imageMove->GetRequestedRegion());
		imageSave->SetSpacing(imageMove->GetSpacing());
		imageSave->Allocate(0);
		IteratorTypeOrigin iterMove(imageMove,imageMove->GetRequestedRegion());
		IteratorTypeOrigin iterSave(imageSave,imageSave->GetRequestedRegion());
		iterMove.GoToBegin();
		iterSave.GoToBegin();
		while(!iterMove.IsAtEnd())
		{
			
			auto movePixel = iterMove.Get();
			//if(movePixel >0)
				iterMove.Set(movePixel);
			/*else
			{
				iterMove.Set(0);
			}*/
			//iterSave.Set(movePixel);
			
			++iterMove;
			++iterSave;
		}
		std::string pre="D:/liver segmentation/data/test3.4/";
		std::string post=".mhd";
		std::string namew= pre+num+post;
		writeOriginImage(namew.c_str(),imageMove);
		//imageMove->Delete();
		++fileNum;
		std::cout<<fileNum<<std::endl;
	}

	int a;
	std::cin>>a;

}