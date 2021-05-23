#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>
#define M 400		//x column
#define N 200		//y string
#define L 200		//Монте-Карло (кол-во траекторий)
#define scatt 12	// акты рассеяния
#define c 1500.		//скорость звука
#define sigmatimes 5	//множитель
#define mu (0.018 * sigmatimes)	//коэф. затухания
#define h 0.4					//шаг
#define PI 3.1415
#define maximum (0.018 * sigmatimes)	


double m_max(double mass[N][M]) {

	double max = -100000000000;

	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < M; ++j) {
			if (mass[i][j] > max)
				max = mass[i][j];
		}
	}
	return max;
}

double m_min(double mass[N][M]) {

	double min = 100000000000;

	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < M; ++j) {
			if (mass[i][j] < min)
				min = mass[i][j];
		}
	}
	return min;
}

bool mass_print_int_bin(const char* f_name, double sigma_min, double sigma_max, double u_im[N][M]) {//print array to binary file with scale from 0 to 255

	FILE* f_out = fopen(f_name, "wb");
	if (!f_out)
		return 0;
	char s_out[M];
	int i = 0;
	int j = 0;
	char value = 0;
	double v = 0;
	double base = 42;
	for (i = 0; i < N; i++) {
		for (j = 0; j < M; j++) {
			if (u_im[i][j] <= sigma_min)
				value = 0;
			else if (u_im[i][j] >= sigma_max)
				value = 255;
			else value = char((u_im[i][j] - sigma_min) / (sigma_max - sigma_min) * 255.);
			/*else {
				v = (powl(base, ((u_im[i][j] - sigma_min) / (sigma_max - sigma_min))) - 1)*(sigma_max - sigma_min) / (base - 1) + sigma_min;
				value = char((v - sigma_min) / (sigma_max - sigma_min)*255.);
			}*/
			s_out[j] = value;
		}
		fwrite(&s_out, sizeof(char), M, f_out);
	}
	fclose(f_out);
	return true;
}