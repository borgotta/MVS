#pragma once

#include "stdafx.h"

#include <iostream>
#include <sstream>

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>


using namespace cv;
using namespace std;

class Settings
{
public:
	Settings(void): goodInput(false) {}
	~Settings(void);
	public:

    void write(FileStorage& fs) const;                      //Write serialization for this class

    void read(const FileNode& node);                          //Read serialization for this class

    void interprate();

    Mat nextImage();

	//double * nextMat();

    static bool readStringLists( const string& filename, vector<string>& images, vector<string>& txts);

public:

    string outputFileName;      // The name of the file where to write
    string input;               // The input ->

    vector<string> imageList;
	vector<string> txtList;
    bool goodInput;
    int flag;

private:
    string patternToUse;

};
void read(const FileNode& node, Settings& x, const Settings& default_value);
void write(FileStorage& fs, const std::string&, const Settings& x);