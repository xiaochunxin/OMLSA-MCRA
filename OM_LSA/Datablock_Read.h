#pragma once
#include"base4_fft.h"
#include<iostream>
#include"LSA_denoise.h"
#include"head.h"
// pInBuffer:输入数据暂存，空间需要接口外部开辟，大小为一次读取数据所占空间的最大值
// pOutBuffer:输出数据暂存，空间需要接口外部开辟，大小为pInBuffer+一帧的数据，可取一帧数据的最大值4096
class Datablock_Read{
public:
	Datablock_Read(int sample_rate, short channels,int MaxDataLen);
	//降噪程序将会返回一值，若该值小于0则程序出现异常，需要跳过降噪程序回到主循环
	short Data_procese(short* pInBuffer, short* pOutBuffer,int read_length, int& out_length);
	~Datablock_Read();
private:
	short m_derr_code;  //错误代码返回值
	//出错范围：-1到-14 base4_fft   -15到-18 LSA_denoise
	//-19到-49 G_calculate  -50到-58 Datablock_Read
	int m_maxdata;
	short m_inc, m_wlen,m_blockInd,m_inc_move;
	short m_channels,m_wlen15,m_inc2;
	int m_sample_rate,m_data_rest_length ;
	short* m_data_in;
	short* m_data_storage;
	int* m_process_storage;
	short* m_buffer;
	short* m_data_resize;
	short* m_DoubDataBuffer;
	short* m_data_out;
	LSA_denoise LSA;
	short Initial(int MaxDataLen);

};

