#define _CRT_SECURE_NO_WARNINGS
#include<iostream>
#include"Datablock_Read.h"
#include<vector>
#include"head.h"
#include"time.h"

using namespace std;
void procese(const char *szFileIn, const char *szFileOut,int sample_rate)
{
	FILE *fpIn=fopen(szFileIn, "rb");
	if (NULL == fpIn)
	{
		printf("open src file err1 \n");
	}
	fseek(fpIn, 44, 0);
	FILE *fpOut=fopen(szFileOut, "wb");  //出错
	if (NULL == fpOut)
	{
		printf("open out file err2! \n");
	}

	short pInBuffer[section_max];
	memset(pInBuffer, 0, section_max*sizeof(short));
	short pOutBuffer[section_max+frame_max];
	memset(pOutBuffer, 0, (section_max + frame_max) * sizeof(short));

	int Out_Length = 0;		//输出数据长度
	short voice_channel = 1;
	short abnormal_flag=0;
	int MaxDataLen = 10000;
	Datablock_Read dtr(sample_rate, voice_channel, MaxDataLen);
	
	if (abnormal_flag < 0){}
	else { 
		int current_length = 0;

		clock_t start, finish;
		double duration;
		start = clock(); 
		while (!feof(fpIn))
		{
			int read_length = 4000; // 896;  
			int nread = fread(pInBuffer, sizeof(short), read_length, fpIn);

			if (nread > 0)
			{
				abnormal_flag = dtr.Data_procese(pInBuffer, pOutBuffer, read_length, Out_Length);

				if (abnormal_flag < 0) break;
				if (Out_Length > 1) {
					fwrite(pOutBuffer, sizeof(short), Out_Length, fpOut);
				}
				memset(pInBuffer, 0, MaxDataLen * sizeof(short));
				memset(pOutBuffer, 0, (MaxDataLen + frame_max) * sizeof(short));
				current_length += read_length;
			}
		}
		finish = clock();
		duration = (double)(finish - start) / CLOCKS_PER_SEC;
		printf("%f seconds", duration);

	}
	system("pause");
    fclose(fpIn);
	fclose(fpOut);
}
//MY_B4_FFT* MY_B4_FFT::NSingleton = new MY_B4_FFT; //饿汉式  单例模式 

int main(const int argc, const char * argv[]) {
	//当vs读取wav文件时，会将wav的44字节头文件也读取进来，影响了部分程序的读写正确性
	procese("cleanadf10db.wav", "out.pcm",32000);
	for (int i = 0; i < argc; i++) {
	//printf("%d ", nParams);
	printf("%s\n", argv[i]);
	}
	printf("降噪完成");

	system("pause");
	return 0;
}