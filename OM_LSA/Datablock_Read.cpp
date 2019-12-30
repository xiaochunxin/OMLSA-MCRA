#define _CRT_SECURE_NO_WARNINGS
#include "Datablock_Read.h"
#include "head.h"
using namespace std;

Datablock_Read::Datablock_Read(int sample_rate,short channels,int MaxDataLen)
{
	m_channels = channels;
	m_sample_rate = sample_rate;
	m_derr_code=Initial(MaxDataLen);
}

short Datablock_Read::Initial(int MaxDataLen) {
	m_maxdata = MaxDataLen;
	if (m_sample_rate < 23000)	//采样率小于23K
	{
		m_wlen = 256;
		m_inc = 128;  //帧长为4的次方
	}
	else {
		m_wlen = 1024;
		m_inc = 512;
	}
	m_wlen15 = m_wlen +(m_wlen>>1);
	m_inc2 = m_inc <<1;

	m_inc_move = log(m_inc) / log(2);

	m_DoubDataBuffer = new short[MaxDataLen + m_wlen];
	if (!m_DoubDataBuffer) {
		return -52;
	}
	m_data_in = new short[m_wlen15];
	if (!m_data_in) {
		return -53;
	}
	m_data_storage = new short[m_wlen15];
	if (!m_data_storage) {
		return -54;
	}
	m_process_storage=new int[m_wlen15];
	if (!m_process_storage) {
		return -55;
	}
	m_buffer = new short[m_wlen15];
	if (!m_buffer) {
		return -56;
	}
	m_data_out = new short[m_wlen15];
	if (!m_data_out) {
		return -57;
	}
	m_data_resize = new short[MaxDataLen + m_wlen];
	if (!m_data_resize) {
		return -58;
	}

	memset(m_process_storage, 0, m_wlen15 * sizeof(int));
	memset(m_DoubDataBuffer, 0, (MaxDataLen + m_wlen) * sizeof(short));
	m_data_rest_length = 0;
	m_blockInd = 0;
	short m_ierr_code = LSA.Initialize(m_wlen15); //15--49
	if (m_ierr_code < 0)return m_ierr_code;
	//My_fft->initial(m_wlen); //	FFT初始化  单例转换
	return 0;
}

short Datablock_Read::Data_procese(short* pInBuffer, short* pOutBuffer,int read_length,int& Out_Length){
	if (m_derr_code < 0) return m_derr_code;
	if (pInBuffer == NULL)return -50;
	if (pOutBuffer == NULL)return -51;

	if (m_channels == 2)
		read_length = (read_length >> 1);
	int current_total_data = read_length + m_data_rest_length;  // 重组后的数据个数

	if (current_total_data < m_wlen15) {	//  数据不满一帧
		if(m_channels==1)
		memcpy(m_data_storage + m_data_rest_length, pInBuffer, (read_length) * sizeof(short));
		else if (m_channels == 2) {
			for (int i = 0; i < read_length; i++) 
				m_data_storage[i + m_data_rest_length] = (pInBuffer[i << 1] + pInBuffer[(i << 1) + 1]) >> 1;
		}

		Out_Length = 0;
		m_data_rest_length = current_total_data;
	}
	else {
	// 数据整合  考虑第一帧数据的存储
	// 先将保留区数据导入
	memcpy(m_data_resize, m_data_storage, (m_data_rest_length) * sizeof(short));
	// 再将读取数据导入
	if (m_channels == 1){
		memcpy(m_data_resize + m_data_rest_length, pInBuffer, read_length * sizeof(short));
	}
	else if (m_channels == 2) {
		for (int i = 0; i < read_length; i++) {
			m_data_resize[i + m_data_rest_length] = (pInBuffer[i << 1] + pInBuffer[(i << 1) + 1]) >> 1;   //2^15
		}
	}

 	int data_use = 0;				 // 计算读取的数据中已经使用过的数据-以帧头计数
	int current_frame = ((current_total_data - m_wlen) >> m_inc_move) + 1;  // 当前数据所能构成的帧数
	//对段总帧数取偶数  奇数部分存进暂存区
	current_frame = (current_frame | 1) - 1;

	m_data_rest_length = current_total_data -( current_frame*m_inc);  // 更新剩余数据 尾帧中之后
	Out_Length = current_frame * m_inc;	  // 输出到末段中

	for (int k = 0; k < current_frame-1; k+=2, m_blockInd+=2){
		memcpy(m_data_in, m_data_resize + data_use, m_wlen15 * sizeof(short)); //分帧
		memset(m_data_out, 0, sizeof(short)*m_wlen15);
		m_derr_code=LSA.Denoise_process(m_data_in, m_data_out, m_blockInd);

		if (m_derr_code < 0) return m_derr_code;
		//输出要求：重叠好的两帧
 		int block_inc = k * m_inc;  //段的当前位置
		if (current_frame == 2) {   
			for (int i = 0; i < m_inc; i++) {
				m_DoubDataBuffer[block_inc + i] += m_data_out[i];
				m_DoubDataBuffer[block_inc + i] += m_process_storage[i];
				m_DoubDataBuffer[block_inc + m_inc+i] += m_data_out[i+m_inc];
				m_process_storage[i] = m_data_out[i+m_wlen];
			}
		}
		else {
			if (k == 0) {                    
				for (int i = 0; i < m_inc; i++) {
					m_DoubDataBuffer[block_inc + i] += m_data_out[i];
					m_DoubDataBuffer[block_inc + i] += m_process_storage[i];
					m_DoubDataBuffer[block_inc + m_inc + i] += m_data_out[i + m_inc];
					m_DoubDataBuffer[block_inc + m_wlen + i ] += m_data_out[i + m_wlen];
				}
			}
			else if (k == current_frame-2 ){
				for (int i = 0; i < m_inc; i++) {   
					m_DoubDataBuffer[block_inc + i] += m_data_out[i];
					m_DoubDataBuffer[block_inc + m_inc + i] += m_data_out[i + m_inc];
					m_process_storage[i] = m_data_out[i + m_wlen];
				}
			}
			else {   //段中间帧
				for (int i = 0; i < m_wlen15; i++) {
					m_DoubDataBuffer[block_inc + i] += m_data_out[i];
				}
			}
		}
		if (m_channels == 1) {
			for (int i = 0; i < Out_Length; i++) {
				pOutBuffer[i] = m_DoubDataBuffer[i];
			}
		}
		else if (m_channels == 2){
			for (int i = 0; i < Out_Length; i++) {
				pOutBuffer[i << 1] = m_DoubDataBuffer[i];
				pOutBuffer[(i << 1) + 1] = m_DoubDataBuffer[i];  //双通道输出 
			}
		}
		data_use += m_wlen;
	}
	if (m_channels == 2)
		Out_Length = Out_Length << 1;  
	//存入残余数据   
	memcpy(m_data_storage, m_data_resize + data_use, (m_data_rest_length) * sizeof(short));
	memset(m_DoubDataBuffer, 0, (m_maxdata + m_wlen) * sizeof(short));
	}
	return 0;
}

	 
Datablock_Read::~Datablock_Read()
{
	if (!m_data_storage) {
		delete[] m_data_storage;
		m_data_storage = NULL;
	}
	if (!m_data_in) {
		delete[] m_data_in;
		m_data_in = NULL;
	}
	if (!m_data_resize) {
		delete[] m_data_resize;
		m_data_resize = NULL;
	}
	if (!m_buffer) {
		delete[] m_buffer;
		m_buffer = NULL;
	}
	if (!m_DoubDataBuffer) {
		delete[] m_DoubDataBuffer;
		m_DoubDataBuffer = NULL;
	}
	if (!m_data_out) {
		delete[] m_data_out;
		m_data_out = NULL;
	}
}



