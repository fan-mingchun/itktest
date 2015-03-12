#include "../readUti.h"


int main()
{
	string pre = "D:/liver segmentation/data/ssmsegLevel/";
	string newPre = "D:/liver segmentation/data/ssmsegLevel2/";
	string post = ".raw";
	int xmax,ymax,zmax,xmin,ymin,zmin;
	int indexMax,indexMin;
	xmin = ymin=zmin = 10000;
	xmax = ymax = zmax = -1;
	indexMax = -1;
	indexMin = 1000;
	int filenum = 20;
	int fileIndex = 1;
	int wid = 512;
	int hei = 512;
	int dep = 178;
	int size = wid*hei*dep;
	short *data = new short[size];
	while (fileIndex <=filenum)
	{
		string name = getNumOfFile(pre,post,fileIndex);
		ifstream input(name.c_str(),ios::binary);
		input.read((char*)data,sizeof(short)*size);

		for(int x = 0;x < wid;++x)
			for(int y = 0;y <hei;++y)
				for(int z = 0;z <dep ;++z)
				{
					if(data[z*wid*hei+y*wid+x])
					{
						if(x < xmin)
							xmin=x;
						if(y < ymin)
							ymin = y;
						if(z < zmin)
							zmin = z;
						if(x > xmax)
							xmax=x;
						if(y > ymax)
							ymax = y;
						if(z > zmax)
							zmax = z;
					}
				}
				fileIndex++;
	}
	int oldwid = wid;
	int oldhei = hei;
	int olddep = dep;
	wid = xmax-xmin+1;
	hei = ymax-ymin+1;
	dep = zmax-zmin+1;
	
	int oldSize = size;
	size = wid*hei*dep;
	indexMin = zmin*wid*hei+ymin*hei+xmin;

	fileIndex = 1;
	while (fileIndex <=filenum)
	{
		string name = getNumOfFile(pre,post,fileIndex);
		ifstream input(name.c_str(),ios::binary);

		input.read((char*)data,sizeof(short)*oldSize);
		short* newData = new short[size];
		//memset(newData,0,sizeof(short)*size);
		for(int x = xmin;x <= xmax;++x)
			for(int y = ymin;y <= ymax;++y)
				for(int z = zmin;z <= zmax ;++z)
				{
					newData[(z-zmin)*wid*hei+(y-ymin)*wid+(x-xmin)]=data[z*oldwid*oldhei+y*oldwid+x];
					
				}
				ofstream output(getNumOfFile(newPre,post,fileIndex).c_str(),ios::binary);
				output.write((char*)newData,sizeof(short)*size);
				fileIndex++;
		delete [] newData;
	}
	cout<<wid<<endl;
	cout<<hei<<endl;
	cout<<dep<<endl;
	return 0;
}