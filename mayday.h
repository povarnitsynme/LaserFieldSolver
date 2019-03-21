#ifndef _mayday_h
#define _mayday_h
#include <string>
#include <fstream>
#include <stdlib.h>

using std::string;
using std::ofstream;
using std::endl;

extern ofstream pofferr;
void mayDayMessage(const string msg) {

	if(pofferr.is_open())
		pofferr << msg << endl;


}


void mayDayAbort(const string msg) {

	ofstream timedata;
	do {
		timedata.open("time.dat");
	} while (!timedata.is_open());

	if(timedata.is_open()) {
		timedata << msg << endl;
       		timedata.close();
	}
	exit(666);

}



#endif


