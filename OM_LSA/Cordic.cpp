#include "Cordic.h"



Cordic::Cordic()
{
}
int Cordic::cordic_sin(int theta){
	int negtive_flag = 1; 
	int sin_value, d;
	x[0] = 652032943;  // 0.607253*1073741824  30
	y[0] = 0;
	while (theta < -Pi_2) {
		theta += Pi_2; 
		negtive_flag *= -1; 
	}
	while (theta > Pi_2) {
		theta -= Pi_2; 
		negtive_flag *= -1; 
	}
	while (theta > Pi_div2) {
		theta = Pi_2 - theta; 
	}
	while (theta < -Pi_div2) {
		theta = -Pi_2 - theta; 
	}

	for (int i = 0; i < 15; i++) {
		if (theta > 0)d = 1;
		else d = -1;
		x[i + 1] = x[i] - d * (y[i] >> i);
		y[i + 1] = y[i] + d * (x[i] >> i);
		theta = theta - d * atanget[i];
	}
	sin_value = negtive_flag * y[15];
	return sin_value;
}

int Cordic::cordic_cos(int theta)
{
	int negtive_flag = 1;
	int cos_value, d;
	x[0] = 652032943;  // 0.607253*256
	y[0] = 0;
	while (theta < -Pi_2) {
		theta += Pi_2;
		negtive_flag *= -1;
	}
	while (theta > Pi_2) {
		theta -= Pi_2;
		negtive_flag *= -1;
	}
	while (theta > Pi_div2) {
		theta = Pi_2 - theta;
		negtive_flag *= -1;
	}
	while (theta < -Pi_div2) {
		theta = -Pi_2 - theta;
		negtive_flag *= -1;
	}

	for (int i = 0; i < 15; i++) {
		if (theta > 0)d = 1;
		else d = -1;
		x[i + 1] = x[i] - d * (y[i] >> i);
		y[i + 1] = y[i] + d * (x[i] >> i);
		theta = theta - d * atanget[i];
	}
	cos_value = negtive_flag * x[15];
	return cos_value;
}

Cordic::~Cordic()
{
}
