
#include "../readUti.h"


int main()
{
	string pre = "D:/liver segmentation/data/ssmsegLevel/";
	string post = ".mhd";
	auto iktImage = readFloatImage("D:/liver segmentation/data/distance/meanImage.mhd");
	auto origin = FloatImageToOrigin(iktImage);
	auto vtkImage = itkOriginImageToVtk(origin);
	vtkRenderWindow *renWin = vtkRenderWindow::New(); 

	auto poly =getPolyDataLoadImage(vtkImage);
	auto renderer= getRenderer(poly);
	renWin->AddRenderer(renderer);

	vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();  
	iren->SetRenderWindow(renWin);  
	renWin->SetSize(800,800);  
	renWin->Render();  
	iren->Start();   
	renWin->Delete();  
	iren->Delete();  
	return 0;
}

void simpliedPolyData(vtkPolyData* polydata)
{

}