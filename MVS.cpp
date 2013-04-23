// MVS.cpp : Defines the entry point for the console application.
#include "stdafx.h"
#include "Settings.h"

#include <iostream>
#include <sstream>

using namespace std;

int _tmain(int argc, char* argv[])
{
	//mainCalib("in_VID5.xml");
	string path;
	Settings settings;
	if (argc >= 2) {
		path = argv[1];
	} else {
		path = "xml\\default.xml";
	}
	const string constpath = path;
    FileStorage fs(constpath, FileStorage::READ); // Read the settings
    if (!fs.isOpened())
    {
        cout << "Could not open the configuration file: \"" << path << "\"" << endl;
        return -1;
    }
    fs["Settings"] >> settings;
    fs.release(); 
	
	cout<<path<<"\n"<<settings.imageList.size()<<"\n";
	return 0;
}

