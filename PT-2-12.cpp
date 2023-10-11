#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;

void copyArray(int* to, int* from, int size)
{
	for (size_t i = 0; i < size; i++)
	{
		to[i] = from[i];
	}
}

bool isNumber(string value) {
	for (size_t i = 0; i < value.size(); i++) {
		if (!isdigit(value[i])) {
			return false;
		}
	}
	return true;
}

class Matrix {
private:
	int** matrix;
	int** tempmatrix;
	int size; // Размер матрицы (количество строк или столбцов)
	int nnz;
	int* AD;
	int* AU;
	int* AL;
	int* LJ;
	int* LI;
public:
	Matrix(const char* filename) 
	{
		ifstream file(filename);

		if (!file.is_open()) {
			cout << "Не удается открыть файл." << endl;
			size = 0;
			matrix = nullptr;
			return;
		}

		file >> size;
		string value;

		// Выделяем память под матрицу
		matrix = new int* [size];
		for (int i = 0; i < size; ++i) {
			matrix[i] = new int[size];
		}

		// Считываем элементы матрицы из файла
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < size; ++j) {
				file >> value;
				if (isNumber(value)) {
					matrix[i][j] = atoi(value.c_str());
				}
				else {
					cout << "Введено неверное значение!" << endl;
					exit(-1);
				}
			}
		}

		file.close();

		countNonZeroElements();
		packMatrix();
	}

	Matrix() 
	{
		// Выделяем память под матрицу
		matrix = new int* [size];
		for (int i = 0; i < size; ++i) {
			matrix[i] = new int[size];
		}

		countNonZeroElements();
		packMatrix();
	}

	Matrix(const Matrix& right)
	{
		copy(right);
	}

	~Matrix() {
		if (matrix != nullptr) {
			for (int i = 0; i < size; ++i) {
				delete[] matrix[i];
			}
			delete[] matrix;
		}

		clearMemory();

	}

	void printMatrix() {
		if (matrix == nullptr) {
			cout << "Матрица не инициализирована." << endl;
			return;
		}

		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < size; ++j) {
				cout << matrix[i][j] << " ";
			}
			cout << endl;
		}

		cout << "Сжатый формат:\n";

		string outStr = "AD: ";
		for (size_t i = 0; i < size; i++)
		{
			outStr += to_string(AD[i]) + " ";
		}
		cout << outStr + "\n";

		outStr = "AU: ";
		for (size_t i = 0; i < nnz / 2; i++)
		{
			outStr += to_string(AU[i]) + " ";
		}
		cout << outStr + "\n";

		outStr = "AL: ";
		for (size_t i = 0; i < nnz / 2; i++)
		{
			outStr += to_string(AL[i]) + " ";
		}
		cout << outStr + "\n";

		outStr = "LJ: ";
		for (size_t i = 0; i < nnz / 2; i++)
		{
			outStr += to_string(LJ[i]) + " ";
		}
		cout << outStr + "\n";

		outStr = "LI: ";
		for (size_t i = 0; i < size; i++)
		{
			outStr += to_string(LI[i]) + " ";
		}
		cout << outStr + "\n";
	}

	void countNonZeroElements() {
		nnz = 0;
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < size; ++j) {
				if (matrix[i][j] != 0 && i != j) {
					nnz++;
				}
			}
		}
	}

	void packMatrix() {
		allocateMemory();

		int adIndex = 0;
		int auIndex = 0;
		int alIndex = 0;
		int ljIndex = 0;
		int liIndex = 0;

		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < size; ++j) {
				if (i == j) {
					AD[adIndex++] = matrix[i][j];
				}
				else if (i < j && matrix[i][j] != 0) {
					AU[auIndex++] = matrix[i][j];
				}
				else if (i > j && matrix[i][j] != 0) {
					AL[alIndex++] = matrix[i][j];
				}
			}
		}

		for (int j = 0; j < size; ++j) {
			LJ[j] = -1;
			for (int i = 0; i < size; ++i) {
				if (i < j && matrix[i][j] != 0) {
					LJ[ljIndex++] = j + 1;
					break;
				}
			}
		}

		int elementId = 1;
		for (int i = 0; i < size; i++) // Массив LI
		{
			int tempId = elementId;
			bool isNull = true;
			for (int j = i + 1; j < size; j++)
			{
				if (matrix[i][j] != 0)
				{
					++elementId;
					isNull = false;
				}
			}

			if (isNull)
				LI[i] = 0;
			else
				LI[i] = tempId;
		}

		int lastElementId = nnz / 2 + 1;
		for (int i = size - 1; i >= 0; i--)
		{
			if (LI[i] == 0)
				LI[i] = lastElementId;
			else
				lastElementId = LI[i];
		}
	}

	void gaussianElimination() {
		cout << "Before Gauss:\n";
		printMatrix();

		int numRows = size;
		int numColumns = size;

		for (int i = 0; i < numRows; i++) {
			// Поиск максимального элемента в текущем столбце
			double maxValue = abs(matrix[i][i]);
			int maxRow = i;
			for (int j = i + 1; j < numRows; j++) {
				if (abs(matrix[j][i]) > maxValue) {
					maxValue = abs(matrix[j][i]);
					maxRow = j;
				}
			}

			// Перестановка строк
			for (int k = i; k < numColumns + 1; k++) {
				double temp = matrix[maxRow][k];
				matrix[maxRow][k] = matrix[i][k];
				matrix[i][k] = temp;
			}

			// Приведение матрицы к треугольному виду
			for (int j = i + 1; j < numRows; j++) {
				double coefficient = -matrix[j][i] / matrix[i][i];
				for (int k = i; k < numColumns + 1; k++) {
					if (i == k) {
						matrix[j][k] = 0;
					}
					else {
						matrix[j][k] += coefficient * matrix[i][k];
					}
				}
			}
		}


		// Обратный ход метода Гаусса
		int solution = 1;
		//std::vector<double> solution(numColumns);
		for (size_t i = 0; i < size; i++)
		{
			solution *= matrix[i][i];
		}

		// Вывод решения
		std::cout << "Solution: " + to_string(solution);
		std::cout << std::endl;
		cout << "After Gauss:\n";
		printMatrix();
	}

private:
	void allocateMemory()
	{
		AD = new int[size];
		AU = new int[nnz / 2];
		AL = new int[nnz / 2];
		LJ = new int[nnz / 2];
		LI = new int[size];
	}
	void clearMemory()
	{
		delete[] AD;
		delete[] AU;
		delete[] AL;
		delete[] LI;
		//delete[] LJ;
	}

	void copy(const Matrix& right)
	{
		nnz = right.nnz;
		size = right.size;

		clearMemory();
		allocateMemory();
		
		copyArray(AD, right.AD, size);
		copyArray(AU, right.AU, nnz / 2);
		copyArray(AL, right.AL, nnz / 2);
		copyArray(LJ, right.LJ, nnz / 2);
		copyArray(LI, right.LI, size);
	}
};


int main() {
	setlocale(LC_ALL, "Ru");

	const char* filename = "input.txt";

	Matrix matrix(filename);

	matrix.printMatrix();

	matrix.gaussianElimination();
	//matrix.printPackedMatrix();


	return 0;

}


