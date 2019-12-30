#pragma once
#include "head.h"
class Cordic
{
public:
	Cordic();
	int cordic_sin(int theta);
	int cordic_cos(int theta);
	~Cordic();
private:
	int Pi_div2 = ((int)(Pi * 32768) >> 1); // pi/2  536870912=2^29
	int Pi_2 = (Pi * 32768);
	const short atanget[16] = { 25736,15193,8027,4075,2045,1024,512,256,128,64,32,16,8,4,2,1 }; // *32768  2^14
	int x[16], y[16];
	
};

