/* File: vfs.cpp
 * Author: CRE
 * Last Edited: Sun Apr  2 16:06:10 2017
 */

#include "crelib/crelib.h"
#include <stdio.h>
#include <map>
#include <string>
#include <sstream>
#include <math.h>
using namespace std;
using namespace cre;

#define BUFFER_SIZE 10240000

map<string, uint> FreqS;

void getFreqInFile(map<string, uint> &FreqS, const char * FileName)
{
	FILE * InFile=fopen(FileName, "r");
	if (InFile==NULL) die("no such file %s.", FileName);

	char * Buffer=myalloc(BUFFER_SIZE, char);

	while (fgets(Buffer, BUFFER_SIZE, InFile))
	{
		if(Buffer[0]=='#') continue;//ignore header lines
		if(Buffer[0]=='\n') continue;//in case for blank lines
		stringstream Ss;
		Ss.str(Buffer);
		string Item;
		string Freq;
		for (uint i=0;i<13;++i)
		{
			getline(Ss, Item, '\t');
			if (Item==string("Chromesome")) break;
			if (i==9) Freq=Item;
			if (i==10)
			{
				stringstream TSs;
				TSs.str(Item);
				uint Number;
				TSs>>Number;
				FreqS[Freq]+= Number;
			}
		}


		if (feof(InFile)) break;
	}


	fclose(InFile);
	free (Buffer);
}

double getC(double p, double MAF, uint N, map<string, uint> &FreqS, uint64 T=0xFFFFFFFFFFFFFFFF)
{
	double C=0;
	double r=0;
	double F=0;
	if (T==0xFFFFFFFFFFFFFFFF)
	{
		T=0;
		for (map<string,uint>::iterator i=FreqS.begin();i!=FreqS.end();++i)
		{
			if (i->first=="") continue;
			if (i->first[0]<'0'||i->first[0]>'9') continue;
			T+=i->second;
		}
	}
	if (T==0) die("T==0?");
	for (map<string,uint>::iterator i=FreqS.begin();i!=FreqS.end();++i)
	{
		if (i->first=="") continue;
		if (i->first[0]<'0'||i->first[0]>'9') continue;
		stringstream Ss;
		Ss.str(i->first);
		Ss>>r;
		if(r>=MAF)
		{
			C+=(1.0-pow(1-r*p, 2*N))*((double)(i->second))/((double)T);
			F+=((double)(i->second))/((double)T);
		}
	}
	C/=F;
	return C;
}

void usage()
{
	die("Usage: vfs a.tsv [b.tsv c.tsv...] > statistics.txt");
}

int main (int argc, char ** argv)
{
	if (argc==1) usage();
	for (int i=1;i<argc;++i)
	{
		getFreqInFile(FreqS, argv[i]);
	}
	uint64 Total=0;
	for (map<string,uint>::iterator i=FreqS.begin();i!=FreqS.end();++i)
	{
		if (i->first=="") continue;
		if (i->first[0]<'0'||i->first[0]>'9') continue;
		if (i!=FreqS.begin()) printf("\n");
		printf("%s:%u", i->first.c_str(), i->second);
		Total+=i->second;
	}
	printf("\nTotal:%llu",Total);
	printf ("\n1070 samples with MAF=0.001, variant find ratio=0.8, coverage=%lf", getC(0.8,0.001,1070,FreqS,Total));
	printf ("\n1000 samples with MAF=0.001, variant find ratio=0.8, coverage=%lf", getC(0.8,0.001,1000,FreqS,Total));
	return 0;
}
