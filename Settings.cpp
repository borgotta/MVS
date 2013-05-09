#include "StdAfx.h"
#include "Settings.h"

namespace MVS {
	//void Settings::init(const Reconstructor& r) {
	//	rec = r;
	//}
	void Settings:: write(FileStorage& fs) const                        //Write serialization for this class
	{
		fs << "{" 
			<< "Write_outputFileName"  << outputFileName
			<< "Input" << input
			<< "CPU" << CPU
			<< "}";
	}
	void Settings:: read(const FileNode& node)                          //Read serialization for this class
	{
		node["Write_outputFileName"] >> outputFileName;
		node["Input"] >> input;
		node["CPU"] >> CPU;

		interprate();
	}

	void write(FileStorage& fs, const std::string&, const Settings& x)
	{
		x.write(fs);
	}
	void read(const FileNode& node, Settings& x, const Settings& default_value = Settings())
	{
		if(node.empty())
			x = default_value;
		else
			x.read(node);
	}
	bool Settings::readStringLists( const string& filename, vector<string>& images, vector<string>& txts) {
		try {
			images.clear();
			txts.clear();
			FileStorage fs(filename, FileStorage::READ);
			if( !fs.isOpened() )
				return false;
			FileNode n = fs["images"];
			if( n.type() != FileNode::SEQ )
				return false;
			FileNodeIterator it = n.begin(), it_end = n.end();
			for( ; it != it_end; ++it ) {
				images.push_back((string)*it);
			}
			n = fs["txt"];
			if( n.type() != FileNode::SEQ )
				return false;
			it = n.begin(), it_end = n.end();
			for( ; it != it_end; ++it ) {
				txts.push_back((string)*it);
			}
		}
		catch (Exception e) {
			cout << "Caught exception number:  " << e.code << endl;
			return false;
		}
		return true;
	}
	void Settings::interprate() {
		if (!readStringLists(input,imageList,txtList)) {
			goodInput = false;
			return;
		}
		//if (!loadPhotos(rec.ps)) {
		//	goodInput = false;
		//	return;
		//}
		goodInput = true; //TODO implement, mock for now
	}
	bool Settings::loadPhotos(PhotoSet &ps) {
		if (imageList.empty() || txtList.empty()) {
			return false;
		}
		
		ps.init(imageList,txtList);
		return true;

	}
	Settings::~Settings(void)
	{
	}
};