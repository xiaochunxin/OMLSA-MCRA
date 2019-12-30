#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include "LSA_denoise.h"
#include<iostream>
#include"time.h"
#include<algorithm>

using namespace std;

LSA_denoise::LSA_denoise()
{
}

short LSA_denoise::Initialize(short wlen)
{
	m_lwlen = wlen;
	m_linc13 = m_lwlen/3;
	m_linc23 = m_linc13<<1;//标准帧长
	m_winData = new Complex_num[m_lwlen<<1];
	if (!m_winData) {
		return -15;
	}
	m_ns_hn = new int[m_lwlen];
	if (!m_ns_hn) {
		return -16;
	}
	m_lerr_code=Gc.Initialize(m_linc23);//将真正的wlen长度传入
	if (m_lerr_code < 0) return m_lerr_code;
	m_lerr_code=MyN_fft.initial(m_linc23);  // 1--14
	if (m_lerr_code <0) return m_lerr_code;  

	for (int i = 1; i < (m_linc23+1); ++i) {   //  2^30= 1073741824
		m_ns_hn[i - 1] = sqrt((0.5 - 0.5 * cos(2.0 * Pi*(i) / (m_linc23 + 1))) * 1073741824);// 1073741824;
	}

	return 0;
}


short  LSA_denoise::Denoise_process( short* data_in, short* data_out , int blockInd)
{
	if (data_in == NULL)return -17;
	if (data_out == NULL) return -18;
	
	for (int i = 0; i < m_linc23; i++){  //取前半帧和后半帧分为两个实序列
		m_winData[i].real = ((__int64)data_in[i] *m_ns_hn[i]) >> 9; //放大2^6 
		m_winData[i].imag = ((__int64)data_in[i + m_linc13] * m_ns_hn[i]) >> 9; //放大2^6 
	} 
 	m_lerr_code=MyN_fft.base4_fft(m_winData, 1);  //将两帧完全分开
	if (m_lerr_code < 0) return m_lerr_code;

	m_lerr_code=Gc.G_calculate_process(m_winData, blockInd);  //19--49
	m_lerr_code=Gc.G_calculate_process(m_winData + m_linc23, blockInd+1);
	if (m_lerr_code <0) return m_lerr_code;

	m_lerr_code=MyN_fft.base4_fft(m_winData, -1);
	m_lerr_code = MyN_fft.base4_fft(m_winData + m_linc23, -1);
	if (m_lerr_code < 0) return m_lerr_code;

	for (int i = 0; i < m_linc23; i++) {
		data_out[i]+= ((__int64)m_winData[i].real * m_ns_hn[i]) >> 21;
		data_out[i+m_linc13] += ((__int64)m_winData[i + m_linc23].real * m_ns_hn[i]) >> 21;
	}

	return 0;
}
LSA_denoise::~LSA_denoise()
{
	if (!m_ns_hn) {
		delete[] m_ns_hn;
		m_ns_hn = NULL;
	}
	if (!m_winData) {
		delete[] m_winData;
		m_winData = NULL;
	}

	//fout.close();
}

//abc[i] = m_G[i] / 1.6384;;
//bcd[i] = m_Gh1[i] / 1.6384;
//cde[i] = m_pr_SNR[i] / 1.6384;
//efg[i] = m_pp[i] / 1.6384;
//aa[i] = m_v[i] / 1.677;;
//bb[i] = m_q[i] / 1.6384;
//cc[i] = m_E_pr_SNR[i] / 1.6384;;
//kds[i] = m_post_SNR[i] / 1.6384;;
//dd[i] = m_abs_Y[i] >> 6;
//ee[i] = m_lamda_d[i] >> 6;
//ff[i] = m_M[i] >> 6;

//for (int k = 0; k < m_lwlen; k++) {
//	abc[k] = m_q[k] / 1.6384;
//	bcd[k] = m_plocal[k] / 1.6384;
//	cde[k] = m_pglobal[k] / 1.6384;
//	efg[k] = m_E_pr_SNR[k] / 1.6384;
//	aa[k] = m_cosen_local[k] / 1.6384;
//	bb[k] = m_cosen_global[k] / 1.6384;
//	cc[k] = m_cosen[k] / 1.6384;
//}

