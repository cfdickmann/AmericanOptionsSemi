#include <stdio.h>
#include "AmericanOption.h"
#include <cstring>
#include <string.h>
#include <iostream>

using namespace std;

int main( int argc, char* args[]) {
	int runden=1;
	AmericanOption AMO;

	bool wieder=false;
	for(int i=0;i<argc;++i)
	{
		string arg=args[i];
		bool geaendert=false;
		if(! arg.compare("-w"))                     {geaendert=true;wieder=true;}
		if(! arg.compare("-10"))                    {geaendert=true;runden=10;}
		if(! arg.compare("-20"))                    {geaendert=true;runden=20;}
		if(! arg.compare("-50"))                    {geaendert=true;runden=50;}
		if(! arg.compare("-100"))                   {geaendert=true;runden=100;}
		if(! arg.compare("-1000"))                  {geaendert=true;runden=1000;}
		if(! arg.compare("-verfaelscht"))           {geaendert=true;AMO.verfaelscht=true;}
		if(! arg.compare("-verbose"))               {geaendert=true;AMO.verbose=true;}
		if(! arg.compare("-loadAlphas"))            {geaendert=true;AMO.loadAlphas=true;}
		if(! arg.compare("-zehnmal"))               {geaendert=true;runden=10;};
		if(! arg.compare("-fuenfzigmal"))           {geaendert=true;runden=50;};
		if(! arg.compare("-10T")) 		            {geaendert=true;AMO.zehnT=true;};
		if(! arg.compare("-4T")) 		            {geaendert=true;AMO.zehnT=true;};
		if(! arg.compare("-semi"))                  {geaendert=true;AMO.semiinf=true;}
		if(i>0 && !geaendert){printf("Unverständliche Parameter!\n");return 0;}
	}
	if (wieder) {
		printf("Anzahl der Wiederholungen\n");
		cin >> runden;
	}
	for(int i=0;i<runden;++i){
		if (AMO.semiinf)AMO.semi();
	}
}
