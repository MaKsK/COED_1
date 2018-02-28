#include <iostream>

using namespace std;
const int M = 1e2;
template<typename T>
void matrixinput(int n, int m, T a[M][M])
{
	cout.precision(5);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
			cout << a[i][j] << " ";
		cout << endl;
	}
}
double t(double r, double n)
{
	return r*sqrt(n - 2) / sqrt(1 - r*r);
}
int main()
{
	//111-120
	int N, p;
	double Z[M][M] = { 0 };

	cout << "N = ";
	cin >> N;
	cout << "p = ";
	cin >> p;


	cout << "Z(" << N << " x " << p << " ) = " << endl;
	double zm[M] = { 0 };
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < p; j++)
		{
			cin >> Z[i][j];
			zm[j] += Z[i][j];
		}
	}
	for (int i = 0; i < p; i++) zm[i] /= N;
	double zd[M] = { 0 };
	for (int i = 0; i < p; i++)
	{
		for (int j = 0; j < N; j++)
		{
			zd[i] += pow(Z[j][i] - zm[j], 2);
		}
	}
	for (int i = 0; i < p; i++) zd[i] /= N;
	cout << "Среднее по столбцам:" << endl;
	for (int i = 0; i < p; i++) cout << zm[i] << " ";
	cout << endl;

	cout << "Дисперсия по столбцам:" << endl;
	for (int i = 0; i < p; i++) cout << zd[i] << " ";
	cout << endl;

	double sqrt_zd[M];

	for (int i = 0; i < p; i++) sqrt_zd[i] = sqrt(zd[i]);

	double X[M][M];
	for (int i = 0; i < N; i++)
		for (int j = 0; j < p; j++)
			X[i][j] = (Z[i][j] - zm[j]) / sqrt_zd[j];

	cout << "Стандартизованная матрица: " << endl;
	matrixinput(N, p, X);
	double covm[M][M];
	for (int i = 0; i < p; i++)
		for (int j = 0; j < p; j++)
		{
			for (int k = 0; k < N; k++)
				covm[i][j] += (Z[k][i] - zm[i])*(Z[k][j] - zm[j]);
			covm[i][j] /= N;
		}
	cout << "Ковариационная матрица: " << endl;
	matrixinput(N, p, covm);
	double R[M][M];
	for (int i = 0; i < p; i++)
		for (int j = 0; j < p; j++)
		{
			for (int k = 0; k < N; k++)
				R[i][j] += X[k][i]*X[k][j];
			R[i][j] /= N;
		}
	cout << "Корреляционная матрица: " << endl;
	matrixinput(N, p, R);

	string H[M][M];
	double ttable = 2.0095752;
	for (int i = 0; i < p; i++)
		for (int j = 0; j < p; j++)
		{
			if (t(R[i][j], N) < ttable) H[i][j] = "H0";
			else H[i][j] = "H1";
		}
	cout << "Проверка гипотез о значимости КК: " << endl;
	matrixinput(N, p, R);
	cin >> N;
	return 0;
}