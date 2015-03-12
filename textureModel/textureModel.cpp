#include "../readUti.h"

int main()
{
	
	typedef itk::Offset<3> offsetType;
	offsetType::OffsetValueType offsetArray[3]={1,1,1};
	itk::Offset<3> offset;
	offset.SetOffset(offsetArray);

	auto filter = Image2CoOccuranceType::New();
	filter->SetOffset(offset);
	return 0;
}