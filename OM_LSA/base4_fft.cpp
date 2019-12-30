#pragma once
#include "base4_fft.h"
#include <math.h>
#include<iostream>

MY_B4_FFT::MY_B4_FFT() {
}

short MY_B4_FFT::initial(int N) {
	m_fwlen = N;  //真实帧长的初始化
	m_finc = m_fwlen>>1;
	m_M4 = log(m_fwlen) / log(4);
	int N_pow = pow(4, m_M4);
	if (N_pow - m_fwlen!=0  || m_fwlen<=0|| m_fwlen> 16384) {  //2^14 =16384   2^12=4096
		return -1;
	}  
	m_Buffer_cos = new int[m_fwlen], m_Buffer_sin = new int[m_fwlen];
	m_sort4_count = new int[m_fwlen];
	m_sort_temp_r = new int[m_fwlen];
	m_sort_temp_i = new int[m_fwlen];
	if (!m_Buffer_cos) {
		return -2;
	}
	if (!m_Buffer_sin) {
		return -3;
	}
	if (!m_sort4_count) {
		return -4;
	}
	if (!m_sort_temp_i) {
		return -5;
	}
	if (!m_sort_temp_r) {
		return -6;
	}
	m_e1_position = new int[m_fwlen];
	if (!m_e1_position) {
		return -7;
	}
	m_o1_position = new int[m_fwlen];
	if (!m_o1_position) {
		return -8;
	}
	m_e2_position = new int[m_fwlen];
	if (!m_e2_position) {
		return -9;
	}
	m_o2_position = new int[m_fwlen];
	if (!m_o2_position) {
		return -10;
	}
	m_yr = new int[m_fwlen];
	if (!m_yr) {
		return -11;
	}
	m_yi = new int[m_fwlen];
	if (!m_yi) {
		return -12;
	}
	m_ifft_move_bit = m_M4<<1;  // 等价于  /m_fwlen

	for (int i = 0; i < m_fwlen; i++) {  
		m_Buffer_cos[i] = (coc.cordic_cos((int)(-2 * Pi*i * 32768)>>m_ifft_move_bit));
		m_Buffer_sin[i] = (coc.cordic_sin((int)(-2 * Pi*i * 32768)>>m_ifft_move_bit));
	}
	Base4_Sort();
	return 0;
}

void MY_B4_FFT::Base4_Sort(){
	int  bit_rest, space, add;
	memset(m_sort4_count, 0, sizeof(int)*m_fwlen);
	for (int l = 0; l < m_M4; l++)
	{
		space = pow(4, l);
		add = pow(4, m_M4 - l - 1);
		for (int i = 0; i < m_fwlen; i++){
			bit_rest = (i / space) % 4;
			if (bit_rest != 0)			   
				m_sort4_count[i] += add * bit_rest;// 位余项，用于乘上逆序后的对应位
		}
	}
}

short MY_B4_FFT::base4_fft(Complex_num *x, int sign) {
	int i, j, k, tr, ti, wi2, wi3;
	if (x == NULL) {
		return -13;//
	}
	if (!abs(sign)) {
		return -14;//
	}
	//对于输入数据序列进行倒位序变换
	int ctemp;
	for (int i = 0; i < m_fwlen; i++){
		m_sort_temp_r[i] = x[i].real;
		m_sort_temp_i[i] = x[i].imag;
	}
	for (int i = 0; i < m_fwlen; i++){
		ctemp = m_sort4_count[i];
		x[i].real = m_sort_temp_r[ctemp];
		x[i].imag = m_sort_temp_i[ctemp];
	}

	int BlockLen, BlockNum, BlockLen2;  int wi, F1, F2; int F3, F0; Complex_num X0, X1, X2, X3;
	for (i = 1; i <= m_M4; i++) {
		BlockNum = pow(4, m_M4 - i);
		BlockLen = pow(4, i);		// 组间间隔+组内数据个数  interval_2
		BlockLen2 = BlockLen >> 2;	// 组内蝶形个数  interval_1
		for (j = 0; j < BlockNum; j++) {	  // 对每一个组的基4蝶形循环计算
			for (k = 0; k < BlockLen2; k++) { // 为组内的第k个蝶形
				F0 = k + j * BlockLen;		  // 每一蝶形的第一个数据
				F1 = F0 + BlockLen2;
				F2 = F1 + BlockLen2;
				F3 = F2 + BlockLen2;
				wi = BlockNum * k;			  // 旋转因子的指数
				wi2 = wi << 1;
				wi3 = (wi << 1) + wi;

				X0.real = x[F0].real;
				X0.imag = x[F0].imag;
				X1.real = (((__int64)x[F1].real*m_Buffer_cos[wi]) >> 30) - (((__int64)x[F1].imag*m_Buffer_sin[wi]) >> 30);
				X1.imag = (((__int64)x[F1].imag*m_Buffer_cos[wi]) >> 30) + (((__int64)x[F1].real*m_Buffer_sin[wi]) >> 30);
				X2.real = (((__int64)x[F2].real*m_Buffer_cos[wi2]) >> 30) - (((__int64)x[F2].imag*m_Buffer_sin[wi2]) >> 30);
				X2.imag = (((__int64)x[F2].imag*m_Buffer_cos[wi2]) >> 30) + (((__int64)x[F2].real*m_Buffer_sin[wi2]) >> 30);
				X3.real = (((__int64)x[F3].real*m_Buffer_cos[wi3]) >> 30) - (((__int64)x[F3].imag*m_Buffer_sin[wi3]) >> 30);
				X3.imag = (((__int64)x[F3].imag*m_Buffer_cos[wi3]) >> 30) + (((__int64)x[F3].real*m_Buffer_sin[wi3]) >> 30);

				x[F0].real = X0.real + X1.real + X2.real + X3.real;
				x[F0].imag = X0.imag + X1.imag + X2.imag + X3.imag;
				x[F1].real = X0.real + X1.imag - X2.real - X3.imag;
				x[F1].imag = X0.imag - X1.real - X2.imag + X3.real;
				x[F2].real = X0.real - X1.real + X2.real - X3.real;
				x[F2].imag = X0.imag - X1.imag + X2.imag - X3.imag;
				x[F3].real = X0.real - X1.imag - X2.real + X3.imag;
				x[F3].imag = X0.imag + X1.real - X2.imag - X3.real;
			}
		}
	}

	if (sign == -1)
	{
		x[0].real = x[0].real >> m_ifft_move_bit;
		x[0].imag = x[0].imag >> m_ifft_move_bit;
		for (i = 1; i <= m_finc; i++) {
			m_fly_tempr = (x[m_fwlen - i].real) >> m_ifft_move_bit;
			m_fly_tempi = (x[m_fwlen - i].imag) >> m_ifft_move_bit;
			x[m_fwlen - i].real = (x[i].real) >> m_ifft_move_bit;
			x[m_fwlen - i].imag = (x[i].imag) >> m_ifft_move_bit;
			x[i].real = m_fly_tempr;
			x[i].imag = m_fly_tempi;
		}
	}
	else if (sign == 1) {  
		//对于两实序列傅里叶变换的还原  
		for (int i = 0; i < m_fwlen; i++) {
			if (i == 0) {
				m_yr[i] = x[i].real; 
				m_yi[i] = x[i].imag;
			}
			else {
				m_yr[i] = x[m_fwlen - i].real; 
				m_yi[i] = x[m_fwlen - i].imag;
			}
			m_e1_position[i] = (x[i].real + m_yr[i]) >> 1;
			m_o1_position[i] = (x[i].real - m_yr[i]) >> 1;
			m_e2_position[i] = (x[i].imag + m_yi[i]) >> 1;
			m_o2_position[i] = (x[i].imag - m_yi[i]) >> 1;
		}
		for (int i = 0; i < m_fwlen; i++) {
			x[i].real = m_e1_position[i];
			x[i].imag = m_o2_position[i];
			x[i + m_fwlen].real = m_e2_position[i];
		 	x[i + m_fwlen].imag = -m_o1_position[i];
		}
	}
	return 0;
}

MY_B4_FFT::~MY_B4_FFT() {
	if (!m_Buffer_cos) {
		delete[] m_Buffer_cos;
		m_Buffer_cos = NULL;
	}
	if (!m_Buffer_sin) {
		delete[] m_Buffer_sin;
		m_Buffer_sin = NULL;
	}
	if (!m_sort4_count) {
		delete[] m_sort4_count;
		m_sort4_count = NULL;
	}
	if (!m_sort_temp_r) {
		delete[] m_sort_temp_r;
		m_sort_temp_r = NULL;
	}
	if (!m_sort_temp_i) {
		delete[] m_sort_temp_i;
		m_sort_temp_i = NULL;
	}
	if (!m_e1_position) {
		delete[] m_e1_position;
		m_e1_position = NULL;
	}
	if (!m_o1_position) {
		delete[] m_o1_position;
		m_o1_position = NULL;
	}
	if (!m_yr) {
		delete[] m_yr;
		m_yr = NULL;
	}
	if (!m_e2_position) {
		delete[] m_e2_position;
		m_e2_position = NULL;
	}
	if (!m_o2_position) {
		delete[] m_o2_position;
		m_o2_position = NULL;
	}
	if (!m_yi) {
		delete[] m_yi;
		m_yi = NULL;
	}
}

void MY_B4_FFT::fftfix(Complex_num *x, int sign) {

	//--------------------------------------------------------------------------
	//按时间抽取法的fft变换，输入反序，输出正序

	int i, j, k, u = 0, l = 0, wi = 0, n1, tr, ti, N = m_fwlen; //j第二层循环（子块中的每个蝶形的循环计数）
									 //k第一层循环（横向fft变换阶数，为log2（N）N为总采样点数
									 //u 蝶形上标x[upper],l 蝶形下标x[lower]，wi旋转因子下标wn[wi]
	int SubBlockNum, SubBlockStep = 1;
	//SubBlockNum当前k层子块数量，SubBlockStep当前k层不同子块的相同位置元素间间隔

	int tempr, tempi;
	N = N << 1;

	n1 = m_fwlen - 1;
	for (j = 0, i = 0; i < n1; i++)
	{              // 位反转运算
		if (i < j)
		{
			tr = x[i].real;
			ti = x[i].imag;

			x[i].real = x[j].real;
			x[i].imag = x[j].imag;

			x[j].real = tr;
			x[j].imag = ti;
		}
		k = N / 2;
		while (k < (j + 1))
		{
			j = j - k;
			k = k / 2;
		}
		j = j + k;
	}

	for (k = N; k > 1; k = (k >> 1)) {				//第一个循环，代表log2(k)阶的变换

		SubBlockNum = k >> 1;				//子块个数为所做点数的一半
		SubBlockStep = SubBlockStep << 1;	//子块间同等地位的元素间隔以2为倍数递增
		wi = 0;							//旋转因子初始化
		for (j = 0; j < SubBlockStep >> 1; j++) {	//第二层循环，更新j值，j表示各个子块的第j个蝶形。
			//因为每个子块的同地位蝶形具有相同的wn，所以用第二层循环控制wn
			for (u = j; u < N; u += SubBlockStep) {	//第三层循环，循环于各个子块间的第j个蝶形，计算所有蝶形。
				//直到下标u越界。(u>N)		
				l = u + (SubBlockStep >> 1);//下标l计算

				//同时计算两实序列的FFT  cos-isin
				tempr = (((__int64)x[l].real*m_Buffer_cos[wi]) >> 30) - (((__int64)x[l].imag*m_Buffer_sin[wi]) >> 30);
				//蝶形x[u]=x[u]+x[l]*Wn,x[l]=x[u]-x[l]*Wn的复数计算
				tempi = (((__int64)x[l].imag*m_Buffer_cos[wi]) >> 30) + (((__int64)x[l].real*m_Buffer_sin[wi]) >> 30);

				x[l].real = x[u].real - tempr;  //迭代，每次更新一次旋转因子带来的变化
				x[l].imag = x[u].imag - tempi;
				x[u].real = x[u].real + tempr;
				x[u].imag = x[u].imag + tempi;
			}
			wi += SubBlockNum;	//第二层循环更新wi值
		}
	}

	if (sign == -1)
	{
		for (i = 0; i < N; i++)
		{
			x[i].real /= N;
			x[i].imag /= N;
		}
	}
}