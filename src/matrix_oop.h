#ifndef MATRIXPLUS_MATRIX_OOP_H_
#define MATRIXPLUS_MATRIX_OOP_H_
#include <cmath>
#include <exception>
#include <iostream>
class S21Matrix {
 public:
  bool EqMatrix(const S21Matrix& other) const;
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();
  S21Matrix();
  S21Matrix(int rows_, int cols_);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other);
  ~S21Matrix();
  S21Matrix operator+(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other);
  S21Matrix operator*(const S21Matrix& other);
  S21Matrix operator*(const double num);
  bool operator==(const S21Matrix& other);
  S21Matrix& operator=(const S21Matrix& other);
  S21Matrix& operator=(S21Matrix&& other);
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix& operator*=(const S21Matrix& other);
  S21Matrix& operator*=(const double num);
  double& operator()(int i, int j);
  const double& operator()(int i, int j) const;
  void SetRows(int rows);
  void SetCols(int cols);
  int GetRows() const { return rows_; }
  int GetCols() const { return cols_; }
  friend S21Matrix operator*(const double num, const S21Matrix& other);

  private:
  int rows_, cols_;
  double** matrix_;
  void Allocate();
  void RemoveMatrix();
  bool IsEqualSize(const S21Matrix& other) const;
  bool IsSquareMatrix() const;
  void GetMinor(int row, int col, S21Matrix& minor);
  void CopyMatrix(const S21Matrix& other);
};

#endif

