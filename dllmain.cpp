#include "pch.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <amp_math.h>

extern "C"

using namespace std;

ofstream fout("test.txt");

double masI[N][M];
double masSig[N][M];

double m[N][M];

#define TYPE double
#define VECTOR_SIZE 2

class rand_191_bit {			//класс 
	long long X[2][4];
	long long m_1;
	long long m_2;

	long double val();//return value in range [0..1]
public:
	rand_191_bit();
	long long mod(long long x, long long y);
	long double random_val(double a, double b);
};

rand_191_bit::rand_191_bit() {
	srand(clock());
	for (int i = 0; i < 1.0e6; i++);
	for (int i = 0; i < 2; i++)
		for (int j = 1; j < 4; j++)
			X[i][j] = clock() + rand();
	m_1 = 4294967087;
	m_2 = 4294944443;
};

long double rand_191_bit::val() {
	X[0][0] = mod((1403580 * X[0][2] - 810728 * X[0][3]), m_1);
	X[1][0] = mod((527612 * X[1][1] - 1370589 * X[1][3]), m_2);
	long long z = mod((X[0][0] - X[1][0]), 4294967087);
	long double u;
	u = double(z) / 4294967087.;
	for (int i = 0; i < 2; i++)
		for (int j = 3; j > 0; j--)
			X[i][j] = X[i][j - 1];
	return u;
}

long double rand_191_bit::random_val(double a, double b) {
	return (a + val() * (b - a));
}

long long rand_191_bit::mod(long long x, long long y) {
	if (y < 0)
		y = -y;
	if (x < 0)
		return x % y + y;
	else
		return x % y;
}

//класс для работы с векторами
class Vect {
	int size_v;
	double* data;

public:
	Vect() {
		data = new double[VECTOR_SIZE];
		data[0] = 0;
		data[1] = 0;
		size_v = VECTOR_SIZE;
	}

	Vect(int size) {
		size_v = size;
		data = new double[size_v];
		for (int i = 0; i < size_v; i++) {
			data[i] = 0;
		}
	}

	Vect(double a, double b) {
		size_v = 2;
		data = new double[size_v];
		data[0] = a;
		data[1] = b;
	}

	Vect(const Vect& a) {
		size_v = a.size_v;
		data = new double[size_v];
		for (int i = 0; i < size_v; i++) {
			data[i] = a.data[i];
		}
	}

	~Vect() {
		delete[] data;
	}

	double SizeVector() {
		return size_v;
	}

	TYPE GetElem(int i) {						// получение элемента вектора
		for (int j = 0; j < size_v; ++j)
			if (i == j)
				return data[i];
	}

	const Vect operator+(Vect a) {
		Vect x;
		for (int i = 0; i < VECTOR_SIZE; i++) {
			x.data[i] = data[i] + a.data[i];
		}
		return x;
	}

	const Vect& operator=(const Vect a) {

		for (int i = 0; i < VECTOR_SIZE; i++) {
			data[i] = a.data[i];
		}
		return *this;
	}


	Vect& AddScal(double s, int i) {
		data[i] = data[i] + s;
		return *this;
	}


	const Vect operator-(Vect a) {
		Vect x;
		for (int i = 0; i < VECTOR_SIZE; i++) {
			x.data[i] = data[i] - a.data[i];
		}
		return x;
	}

	const Vect operator*(double a) {				//умножение на число
		Vect x;
		for (int i = 0; i < VECTOR_SIZE; i++) {
			x.data[i] = data[i] * a;
		}
		return x;
	}

	const double operator*(Vect a) {
		double res = 0;
		for (int i = 0; i < VECTOR_SIZE; i++) {
			res = res + a.data[i] * data[i];
		}
		return res;
	}

	const Vect operator/(double a) {
		Vect x;
		for (int i = 0; i < VECTOR_SIZE; i++) {
			x.data[i] = data[i] / a;
		}
		return x;
	}

	const Vect operator/(Vect a) {
		Vect x;
		for (int i = 0; i < VECTOR_SIZE; i++) {
			x.data[i] = data[i] / a.data[i];
		}
		return x;
	}

	Vect& MulsVector(const Vect x, double s) {		//
		for (int i = 0; i < VECTOR_SIZE; i++) {
			data[i] = x.data[i] * s;
		}
		return *this;
	}

	Vect& DivVec(Vect a, double s) {
		for (int i = 0; i < VECTOR_SIZE; i++)
			data[i] = a.data[i] / s;
		return *this;
	}

	double DotVector(Vect a) {
		double res = 0;
		for (int i = 0; i < VECTOR_SIZE; i++) {
			res = res + a.data[i] * data[i];
		}
		return res;
	}

	double LengthVector() {
		double b = 0;
		for (int i = 0; i < VECTOR_SIZE; i++) {
			b = b + data[i] * data[i];
		}
		return sqrt(b);
	}
};



//координаты: (y,x)
double Sigma(Vect y) {

	Vect a(-40, -20);
	Vect f(20, 10);
	Vect b(40, -40);
	Vect e(-55, -40);
	Vect p(-45, -75);
	Vect q(60, -50);

	Vect onebubble(-70, -65);
	Vect twobubble(-60, -70);
	Vect threebubble(-40, -75);
	Vect fourbubble(5, -70);
	Vect fivebubble(20, -65);
	Vect sixbubble(55, -75);

	Vect toponebubble(-76, -5);
	Vect toptwobubble(-68, -2);
	Vect topthreebubble(-35, -6);
	Vect topfourbubble(5, -3);
	Vect topfivebubble(25, -9);
	Vect topsixbubble(72, -6);

	Vect mediumonebubble(-75, -40);
	Vect mediumtwobubble(-60, -50);
	Vect mediumthreebubble(-45, -45);

	double inclus1, inclus2, inclus3, inclus4, inclus5;
	double inc1, inc2, inc3, inc4, inc5, inc6;
	double topinc1, topinc2, topinc3, topinc4, topinc5, topinc6;
	double mediuminc1, mediuminc2, mediuminc3;

	//inclus1 = sqrt((pow((y.GetElem(0) + a.GetElem(0)), 2) + pow((y.GetElem(1) + a.GetElem(1)), 2)));
	inclus1 = sqrt(pow((y.GetElem(0) + a.GetElem(0)), 2) / pow(f.GetElem(0), 2) + pow((y.GetElem(1) + a.GetElem(1)), 2) / pow(f.GetElem(1), 2));
	inclus2 = sqrt((pow((y.GetElem(0) + b.GetElem(0)), 2) + pow((y.GetElem(1) + b.GetElem(1)), 2)));
	inclus3 = sqrt((pow((y.GetElem(0) + e.GetElem(0)), 2) + pow((y.GetElem(1) + e.GetElem(1)), 2)));
	inclus4 = sqrt((pow((y.GetElem(0) + p.GetElem(0)), 2) + pow((y.GetElem(1) + p.GetElem(1)), 2)));
	inclus5 = sqrt((pow((y.GetElem(0) + q.GetElem(0)), 2) + pow((y.GetElem(1) + q.GetElem(1)), 2)));

	inc1 = sqrt((pow((y.GetElem(0) + onebubble.GetElem(0)), 2) + pow((y.GetElem(1) + onebubble.GetElem(1)), 2)));
	inc2 = sqrt((pow((y.GetElem(0) + twobubble.GetElem(0)), 2) + pow((y.GetElem(1) + twobubble.GetElem(1)), 2)));
	inc3 = sqrt((pow((y.GetElem(0) + threebubble.GetElem(0)), 2) + pow((y.GetElem(1) + threebubble.GetElem(1)), 2)));
	inc4 = sqrt((pow((y.GetElem(0) + fourbubble.GetElem(0)), 2) + pow((y.GetElem(1) + fourbubble.GetElem(1)), 2)));
	inc5 = sqrt((pow((y.GetElem(0) + fivebubble.GetElem(0)), 2) + pow((y.GetElem(1) + fivebubble.GetElem(1)), 2)));
	inc6 = sqrt((pow((y.GetElem(0) + sixbubble.GetElem(0)), 2) + pow((y.GetElem(1) + sixbubble.GetElem(1)), 2)));

	topinc1 = sqrt((pow((y.GetElem(0) + toponebubble.GetElem(0)), 2) + pow((y.GetElem(1) + toponebubble.GetElem(1)), 2)));
	topinc2 = sqrt((pow((y.GetElem(0) + toptwobubble.GetElem(0)), 2) + pow((y.GetElem(1) + toptwobubble.GetElem(1)), 2)));
	topinc3 = sqrt((pow((y.GetElem(0) + topthreebubble.GetElem(0)), 2) + pow((y.GetElem(1) + topthreebubble.GetElem(1)), 2)));
	topinc4 = sqrt((pow((y.GetElem(0) + topfourbubble.GetElem(0)), 2) + pow((y.GetElem(1) + topfourbubble.GetElem(1)), 2)));
	topinc5 = sqrt((pow((y.GetElem(0) + topfivebubble.GetElem(0)), 2) + pow((y.GetElem(1) + topfivebubble.GetElem(1)), 2)));
	topinc6 = sqrt((pow((y.GetElem(0) + topsixbubble.GetElem(0)), 2) + pow((y.GetElem(1) + topsixbubble.GetElem(1)), 2)));

	mediuminc1 = sqrt((pow((y.GetElem(0) + mediumonebubble.GetElem(0)), 2) + pow((y.GetElem(1) + mediumonebubble.GetElem(1)), 2)));
	mediuminc2 = sqrt((pow((y.GetElem(0) + mediumtwobubble.GetElem(0)), 2) + pow((y.GetElem(1) + mediumtwobubble.GetElem(1)), 2)));
	mediuminc3 = sqrt((pow((y.GetElem(0) + mediumthreebubble.GetElem(0)), 2) + pow((y.GetElem(1) + mediumthreebubble.GetElem(1)), 2)));
	//Верхние пузырьки
	if (topinc1 < 3)
		return (mu * 0.1 + 0.004) * sigmatimes;
	if (topinc2 < 1.5)
		return (mu * 0.2 + 0.001) * sigmatimes;
	if (topinc3 < 4)
		return (mu * 0.2 - 0.0016) * sigmatimes;
	if (topinc4 < 2.6)
		return (mu * 0.1 + 0.0032) * sigmatimes;
	if (topinc5 < 4.5)
		return (mu * 0.1 + 0.0025) * sigmatimes;
	if (topinc6 < 2.3)
		return (mu * 0.2 + 0.0018) * sigmatimes;
	//маленькие пузырьки в среднем слое
	if (mediuminc1 < 2)
		return (mu * 0.1 + 0.003) * sigmatimes;
	if (mediuminc2 < 2.6)
		return (mu * 0.2 - 0.001) * sigmatimes;
	if (mediuminc3 < 3.2)
		return (mu * 0.1 + 0.0014) * sigmatimes;
	//Верхний слой
	//if (fabs(y.GetElem(1)) < 10)
	//return mu * 0.2;
	//Основные неоднородности
	if (inclus1 < 1)
		return mu;
	if (inclus2 < 4.5)
		return mu * 0;
	if (inclus3 < 3.5)
		return mu * 0;
	if (inclus4 < 3)
		return mu * 0;
	if (inclus5 < 3)
		return mu * 0;
	//Маленькие пузырьки в нижнем слое
	if (inc1 < 5)
		return (mu * 0.2 + 0.0022) * sigmatimes;
	if (inc2 < 3.5)
		return (mu * 0.2 - 0.001) * sigmatimes;
	if (inc3 < 2.5)
		return (mu * 0.1 + 0.002) * sigmatimes;
	if (inc4 < 4)
		return (mu * 0.1 + 0.0032) * sigmatimes;
	if (inc5 < 2.62)
		return (mu * 0.2 + 0.0019) * sigmatimes;
	if (inc6 < 3)
		return (mu * 0.1 + 0.0018) * sigmatimes;
	//Нижний слой
	//if (fabs(y.GetElem(1)) > 60)
	//return mu * 0.2;
	//Основной
	return mu * 0.1;
}

double TauF(Vect r, Vect k, double s) {

	double n = 0, d = 0;
	double w = 0;
	n = (c * s) * (c * s) - pow(r.LengthVector(), 2);
	d = (c * s) - r * k;
	w = n / (2 * d);
	return w;
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
			s_out[j] = value;
		}
		fwrite(&s_out, sizeof(char), M, f_out);
	}
	fclose(f_out);
	return true;
}

int __declspec(dllexport) main() {

	rand_191_bit value;
	int i, j, g;
	Vect tmp;	//вспомогательная
	Vect r;
	Vect a;		//вспомогательная
	double tau, teta = 1, t, t0, tprev;
	double I1 = 0;
	double sig = 0, fi, eps = h / 10;
	double result = 0, res = 0;
	double multiplier = 0;
	double denom = 0, numerator = 0;
	double m[N][M];

	double error[scatt + 1];		//среднеквадратичная
	double errorexact[scatt + 1];		//среднеквадратичная
	double errormax[scatt + 1];
	double reserrormax[scatt + 1];
	double scattering[scatt + 1];		//отдельные рассеяния
	double scatteringsum[scatt + 1];
	double sigma[scatt + 1], difference[scatt + 1];				//
	Vect k, kprev, pointprev;
	clock_t time;

	for (int i = 0; i <= scatt; i++) {
		error[i] = 0;
		errorexact[i] = 0;
		scattering[i] = 0;
		sigma[i] = 0;
		errormax[i] = 0;
		reserrormax[i] = 0;
		scatteringsum[i] = 0;
	}
	time = clock();

	for (int i = 1; i < N; i++) {		//string			
		r = Vect(-80, i * h);
		for (int j = M - 1; j >= 0; j--) {	//column		
			r.AddScal(h, 0);
			res = 0;
			for (int d = 1; d <= scatt; d++) {
				scattering[d] = 0;
				scatteringsum[d] = 0;
			}
			t0 = (2 * r.LengthVector()) / c;				//t0		
			multiplier = c * exp(-mu * c * t0) / (2 * PI);
			for (g = 0; g < L; g++) {			//по Монте-Карло
				t = t0;
				Vect point(0, 0);					//r0
				k = r / -r.LengthVector();			//k0
				teta = 1;							//teta1
				tau = TauF(point, k, t);
				a = point - k * tau;
				scattering[1] = Sigma(a) / (c * t - point * k);
				for (int s = 2; s <= scatt; s++) {		//по рассеяниям					
					tau = value.random_val(0, TauF(point, k, t));
					fi = value.random_val(0, 2 * PI);
					pointprev = point;
					point = point - k * tau;
					kprev = k;
					Vect k(cos(fi), sin(fi));
					tprev = t;
					t = t - tau / c;
					teta = teta * Sigma(point) * TauF(pointprev, kprev, tprev);

					a = point - k * TauF(point, k, t);
					numerator = Sigma(a);
					denom = c * t - point * k;
					scattering[s] = scattering[s] + teta * numerator / denom;					//отдельные результаты
				}
			}

			scattering[1] = scattering[1] * multiplier;
			for (int d = 2; d <= scatt; d++) {
				scattering[d] = scattering[d] * multiplier / L;
			}
			for (int d = 1; d <= scatt; d++)
				scatteringsum[d] = scatteringsum[d - 1] + scattering[d];

			for (int d = 1; d <= scatt; d++) {
				sigma[d] = (4 * PI * r.LengthVector() * scatteringsum[d]) / (exp(-2 * mu * r.LengthVector()) * c);
				if (sigma[d] > maximum)
					sigma[d] = maximum;
				difference[d] = fabs(sigma[d] - Sigma(r));
				error[d] += pow(Sigma(r) - sigma[d], 2);
				errorexact[d] += pow(Sigma(r), 2);
				errormax[d] = fabs(Sigma(r) - sigma[d]);// / maximum;
				if (errormax[d] > reserrormax[d])
					reserrormax[d] = errormax[d];
			}

			m[i][j] = sigma[scatt];
		}
		cout << i << endl;
	}

	for (int s = 1; s <= scatt; s++) {
		error[s] = sqrt(error[s] / errorexact[s]);
	}

	mass_print_int_bin("sigma_real_data.bin", 0, maximum, m);

	fout << "Standard deviation:	";
	for (int d = 1; d <= scatt; d++) {
		fout << error[d] << "	";
	}
	fout << endl;
	fout << "Max:	";
	for (int d = 1; d <= scatt; d++) {
		fout << reserrormax[d] << "	";
	}
	fout << endl;
	time = clock() - time;
	fout << "Time:	" << (double)time / CLOCKS_PER_SEC / 60 << endl;

	fout.close();
	cout << (double)time / CLOCKS_PER_SEC / 60 << endl;
	//system("bin_array_to_image_bmp_with_pal_log_dyn_name.exe");
	system("pause");
	return 0;
}

