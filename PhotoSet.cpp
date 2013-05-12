#include "stdafx.h"
#include "PhotoSet.h"

using namespace MVS;

Photo PhotoSet::getPhoto(int i) {
	return photoList[i];
}
Mat PhotoSet::getImage(int i) {
	return photoList[i].image;
}
void PhotoSet::init(const vector<string> &iL, const vector<string> &tL) {
	int size = 0;
	if (iL.size() > tL.size()) {
		size = tL.size();
	} else size = iL.size();
	//main filling
	photoList = vector<Photo>(size);
	for (int i = 0; i < size; i++) {
		Photo temp;
		temp.init(iL[i], tL[i], 0);
		photoList[i] = temp;
	}
}
int PhotoSet::size(void) {
	return photoList.size();
}
int PhotoSet::getWidth(const int index){
	return photoList[index].image.cols;
}
int PhotoSet::getHeight(const int index) {
	return photoList[index].image.rows;
}