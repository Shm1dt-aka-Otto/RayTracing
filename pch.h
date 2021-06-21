#define _CRT_SECURE_NO_WARNINGS
#include "framework.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <amp_math.h>
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
