#include "../readUti.h"

int main(int argc,char ** argv)
{
	std::string pre=
		"D:/liver segmentation/data/distance3/";
	std::string post=".mhd";
	auto mean = readFloatImage((pre+"mean.mhd").c_str());
	auto	princ = readFloatImage
		(getNumOfFile(pre,post,1).c_str());
	const int size = 10;
	int displayNum = size*2+1;
	FloatImageType** resultImageArray = 
		new FloatImageType*[displayNum];
	double constant = 233;
	int index = 0;
	double para[PcaData::princNum];
	double pos[3];
	double angle[3];
	angle[0]=1.57;
	angle[1]=0;
	angle[2]=0;
	pos[0] = 50;
	pos[1] = 50;
	pos[2] = 50;
	double scale = 1.8;
	srand(time(0));

	PcaData::loadRatioImage("D:/liver segmentation/data/knnresult/ratio1.mhd");
	for(int i = -size;i <=size;++i)
	{
		for(int index1 = 0;index1 <PcaData::princNum;++index1)
			para[index1] = (((rand()%200)-100)/100.0)*constant;
		//resultImageArray[index++] = getImageModePca(para);
		//resultImageArray[index++] = transformForFitness(para,angle,pos,scale);
		cout<<transformForFitness(para,angle,pos,scale)<<endl;

	}
	/*vtkRenderWindow *renWin = vtkRenderWindow::New();
	for(int i = 0;i < displayNum;++i)
	{
		auto originImage = levelsetFloatImageToOrigin(resultImageArray[i]);
		auto vtkImage = itkOriginImageToVtk(originImage);
		auto poly = getPolyDataLoadImage(vtkImage);
		auto renderer = getRenderer(poly);
		renderer->SetViewport(i*(1.0/displayNum),0,(i+1)*(1.0/displayNum),1);
		renWin->AddRenderer(renderer);
	}
	vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();  
	iren->SetRenderWindow(renWin);  
	renWin->SetSize(2400,600);  
	renWin->Render();  
	iren->Start();   
	renWin->Delete();  
	iren->Delete(); */
	return 0;
}
