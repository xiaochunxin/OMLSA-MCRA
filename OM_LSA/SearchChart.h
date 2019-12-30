#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include<iostream>
#include "G_calculate.h"
using namespace std;

template<class T>
T* G_calculate::file_read(const char* Filename) {
	//ff:要打开的文件  
	FILE* ff;
	fopen_s(&ff,Filename, "rb");
	if (NULL == ff) {
		printf("open out file err4! \n");
		system("pause");
	}
	//查表积分法文件大小计算   
	fseek(ff, 0, 2);		 //使文件指针指向end     // ftell(fpIn);   
	int DataLength = ftell(ff);  //ftell 用于得到文件位置指针当前位置相对于文件首的偏移字节数
	fseek(ff, 0, 0);
	T* T_int_value = new T[DataLength / sizeof(T)];
	int filecount=fread(T_int_value, sizeof(T), DataLength / sizeof(T), ff);
	fclose(ff);
	return T_int_value;
}

int G_calculate::expintpow_solution(int v_subscript) {
	int vec = 0;
	int g = 0;

	vec = ((__int64)v_subscript * 100 ) >>24;  // / 0.0001;
	vec = vec < 1 ? 1 : vec;
	 
	//g = (m_int_value[vec - 1]) ;
	g= (m_int_value1[vec - 1]);
	return g;
}

int G_calculate::subexp_solution(int v_subscript) {

	int vec = 0;
	int g = 0;

	vec = ((__int64)v_subscript * 100) >> 24;  // / 0.0001;
	vec = vec < 1 ? 1 : vec;
	//g = (m_expsub_value[vec - 1]);
	g = (m_expsub_value1[vec - 1]);
	return g;
}

int G_calculate::Gvalue_solution(int Gh1_subscript,int pp_subscript) {
	int veci = 0,vecj=0;    // j:m_pp  i:m_Gh1
	int g = 0;				
	veci = min(Gh1_subscript * 100 >> 14,6999);
	vecj = max((pp_subscript * 100 >> 14)-1,0);

	//g= dp1[veci][vecj];
	//cout << sizeof(m_G_value) << sizeof(m_G_value) / sizeof(m_G_value[0]);
	g = m_G_value[vecj* 7000 + veci];
	return g;
}
