#pragma once
#include "head.h"
#include "Cordic.h"

typedef struct
{
	int real;   // 8192
	int imag;
} Complex_num;


class MY_B4_FFT {
public:
	MY_B4_FFT();
	~MY_B4_FFT();
	short initial(int N);
	short base4_fft(Complex_num *x, int sign);
	/*static MY_B4_FFT* getInstance() {
		return NSingleton;
	}*/
	void fftfix(Complex_num *x, int sign);
private:
	Cordic coc;
	int* m_Buffer_cos, *m_Buffer_sin;
	int* m_sort4_count;
	int *m_e1_position, *m_o1_position;
	int *m_e2_position, *m_o2_position;
	int* m_yr,*m_yi;
	int* m_sort_temp_r;
	int* m_sort_temp_i;
	int m_fly_tempr, m_fly_tempi;
	int m_fwlen, m_finc, m_M4;
	int m_value_bit;
	int m_ifft_move_bit ;
	int m_value_limit;
	void Base4_Sort();

	//static MY_B4_FFT* NSingleton;
};