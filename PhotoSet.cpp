#include "stdafx.h"
#include "PhotoSet.h"

using namespace MVS;

Photo PhotoSet::getPhoto(int i) {
	return photoList[i];
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