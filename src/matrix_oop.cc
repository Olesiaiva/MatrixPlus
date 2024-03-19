#include "matrix_oop.h"

bool S21Matrix::EqMatrix(const S21Matrix& other) const {
  if (!IsEqualSize(other)) {
    return false;
  }
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        if (fabs(matrix_[i][j] - other.matrix_[i][j]) > 1e-7) {
          return false;
        }
      }
    }
 
  return true;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (!IsEqualSize(other)) {
    throw std::invalid_argument("ERROR: Matrix dimensions do not match");
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (!IsEqualSize(other)) {
    throw std::invalid_argument("ERROR: Matrix dimensions do not match");
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (rows_ != other.cols_) {
    throw std::invalid_argument("ERROR: Matrix rows and columns do not match");
  }
  S21Matrix result(rows_, other.cols_);
  for (int i = 0; i < result.rows_; i++) {
    for (int j = 0; j < result.cols_; j++) {
      for (int k = 0; k < cols_; k++) {
        result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
      }
    }
  }

  *this = std::move(result);
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix temp(cols_, rows_);
  for (int i = 0; i < cols_; i++) {
    for (int j = 0; j < rows_; j++) {
      temp.matrix_[i][j] = matrix_[j][i];
    }
  }
  return temp;
}

S21Matrix S21Matrix::CalcComplements() {
  if (!IsSquareMatrix()) {
    throw std::invalid_argument("ERROR: Matrix should be square");
  }
  S21Matrix minor(rows_ - 1, cols_ - 1);
  S21Matrix result(rows_, cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      GetMinor(i, j, minor);
      result.matrix_[i][j] = pow(-1.0, i + j) * minor.Determinant();
    }
  }
  return result;
}

void S21Matrix::GetMinor(int row, int col, S21Matrix& minor) {
  for (int i = 0, m_i = 0; i < rows_; i++) {
    if (i != row) {
      for (int j = 0, m_j = 0; j < cols_; j++) {
        if (j != col) {
          minor.matrix_[m_i][m_j] = matrix_[i][j];
          m_j++;
        }
      }
      m_i++;
    }
  }
}
double S21Matrix::Determinant() {
  if (!IsSquareMatrix()) {
    throw std::invalid_argument("ERROR: Matrix should be square");
  }
  int sign = 1;
  double det = 1.0;
  for (int i = 0; i < rows_; i++) {
    int mainElementRow = i;
    for (int j = i + 1; j < rows_; j++) {
      if (std::fabs(matrix_[i][i]) < std::fabs(matrix_[j][mainElementRow])) {
        mainElementRow = j;
      }
    }
    if (mainElementRow != i) {
      std::swap(matrix_[i], matrix_[mainElementRow]);
      sign *= -1;
    }
    if (std::fabs(matrix_[i][i]) < 1e-10) {
      return 0;
    }
    double mainValue = matrix_[i][i];
    for (int j = i + 1; j < rows_; j++) {
      double factor = matrix_[j][i] / mainValue;
      for (int k = i; k < rows_; k++) {
        matrix_[j][k] -= factor * matrix_[i][k];
      }
    }
    det *= mainValue;
  }
  det *= sign;
  return det;
}

S21Matrix S21Matrix::InverseMatrix() {
  if (!IsSquareMatrix()) {
    throw std::invalid_argument("ERROR: Matrix should be square");
  }
  S21Matrix comp_temp = CalcComplements();
  S21Matrix tran_temp = comp_temp.Transpose();
  double det = this->Determinant();
  if (det == 0) {
    throw std::invalid_argument("ERROR: Matrix is not invertible");
  }
  tran_temp.MulNumber(1 / det);
  return tran_temp;
}

S21Matrix::S21Matrix() : S21Matrix(3, 3) {}

S21Matrix::S21Matrix(const S21Matrix& other): rows_(other.rows_), cols_(other.cols_) {
  Allocate();
  CopyMatrix(other);
}
S21Matrix::S21Matrix(S21Matrix&& other): rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.rows_ = other.cols_ = 0;
  other.matrix_ = nullptr;
}

void S21Matrix::Allocate() {
  matrix_ = new double*[rows_]();
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_]();
  }
}
S21Matrix::S21Matrix(int rows, int cols)
    : rows_(rows), cols_(cols), matrix_(nullptr) {
  if (rows_ <= 0 or cols_ <= 0) {
    throw std::invalid_argument("ERROR: Invalid matrix size");
  }
  Allocate();
}
void S21Matrix::SetRows(int rows) {
  if (rows <= 0) {
    throw std::invalid_argument("ERROR: Invalid matrix size");
  }
  if (rows_ != rows) {
    S21Matrix temp(*this);
    RemoveMatrix();
    rows_ = rows;
    cols_ = temp.cols_;
    Allocate();
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        if (i < temp.rows_) matrix_[i][j] = temp.matrix_[i][j];
      }
    }
  }
}
void S21Matrix::SetCols(int cols) {
  if (cols <= 0) {
    throw std::invalid_argument("ERROR: Invalid matrix size");
  }

  if (cols_ != cols) {
    S21Matrix temp(*this);
    RemoveMatrix();
    cols_ = cols;
    rows_ = temp.rows_;
    Allocate();
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        if (j < temp.cols_) matrix_[i][j] = temp.matrix_[i][j];
      }
    }
  }
}


const double& S21Matrix::operator()(int i, int j) const { 
  if (i < 0 or j < 0 or i >= rows_ or j >= cols_) {
    throw std::out_of_range("ERROR: Index out of range");
  }
  return matrix_[i][j]; }

double& S21Matrix::operator()(int i, int j)  { 
  if (i < 0 or j < 0 or i >= rows_ or j >= cols_) {
    throw std::out_of_range("ERROR: Index out of range");
  }
  return matrix_[i][j]; }

bool S21Matrix::IsEqualSize(const S21Matrix& other) const {
  return ((rows_ == other.rows_) and (cols_ == other.cols_));
}

bool S21Matrix::IsSquareMatrix()const { return rows_ == cols_; }

S21Matrix::~S21Matrix() { RemoveMatrix(); }
void S21Matrix::RemoveMatrix() {
    for (int i = 0; i < rows_; i++) {
      delete[] matrix_[i];
    }
    delete[] matrix_;
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  S21Matrix result(*this);
  result.SumMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  S21Matrix result(*this);
  result.SubMatrix(other);
  return result;
}
S21Matrix S21Matrix::operator*(const S21Matrix& other) {
  S21Matrix result(*this);
  result.MulMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(const double num) {
  S21Matrix result(*this);
  result.MulNumber(num);
  return result;
}
S21Matrix operator*(const double num, const S21Matrix& other) {
  S21Matrix result(std::move(other));
  result.MulNumber(num);
  return result;
}

bool S21Matrix::operator==(const S21Matrix& other) { return EqMatrix(other); }

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if(this != &other){
  RemoveMatrix();
  rows_ = other.rows_;
  cols_ = other.cols_;
  Allocate();
  CopyMatrix(other);
  }
  return *this;
}

S21Matrix& S21Matrix::operator=(S21Matrix&& other) {
  if(this != &other){
  RemoveMatrix();
  rows_ = other.rows_;
  cols_ = other.cols_;
  matrix_ = other.matrix_;
  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;
}
return *this;
}

void S21Matrix::CopyMatrix(const S21Matrix& other) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      if (i < other.rows_) matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) { 
this->SumMatrix(other);
return *this;
 }
S21Matrix& S21Matrix::operator-=(const S21Matrix& other) { 
this->SubMatrix(other);
return *this;
}
S21Matrix& S21Matrix::operator*=(const S21Matrix& other) { 
this->MulMatrix(other);
return *this;
}
S21Matrix& S21Matrix::operator*=(const double num) {
this->MulNumber(num);
return *this;
}