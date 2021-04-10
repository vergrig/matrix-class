#include <vector>
#include <ostream>
#include <algorithm>
#include <iostream>
#include "Polynomial.h"
#include "Rational.h"
using std::vector;

template<typename T>
class Matrix {
 private:
    std::vector<std::vector<T>> data;
    size_t rows = 0, columns = 0;

 public:
    Matrix(const T& elem) {
        data.resize(1);
        rows = 1;
        columns = 1;
        data[0].push_back(elem);
    }

    Matrix(const std::vector<std::vector<T>>& base) {
        data.resize(base.size());
        rows = base.size();
        if (rows != 0) {
            columns = base[0].size();
        }
        data = base;
    }

    Matrix& operator+=(const Matrix& other) {
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < columns; ++j) {
                data[i][j] += other.data[i][j];
            }
        }
        return *this;
    }

    friend Matrix operator+(Matrix first, const Matrix& second) {
        first += second;
        return first;
    }

    Matrix& operator-=(const Matrix& other) {
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < columns; ++j) {
                data[i][j] -= other.data[i][j];
            }
        }
        return *this;
    }

    friend Matrix operator-(Matrix first, const Matrix& second) {
        first -= second;
        return first;
    }

    template<typename N>
    Matrix& operator*=(const N& scalar) {
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < columns; ++j) {
                data[i][j] *= scalar;
            }
        }
        return *this;
    }

    template<typename N>
    friend Matrix operator*(Matrix first, const N& scalar) {
        first *= scalar;
        return first;
    }

    template<typename N>
    friend Matrix operator*(const N& scalar, Matrix first) {
        first *= scalar;
        return first;
    }

    Matrix& operator*=(const Matrix& other) {
        std::vector<std::vector<T>> new_data;
        new_data.assign(rows, std::vector<T>(other.columns, 0));
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < other.columns; ++j) {
                for (size_t k = 0; k < columns; ++k) {
                    new_data[i][j] += data[i][k] * other.data[k][j];
                }
            }
        }
        data = new_data;
        columns = other.columns;
        return *this;
    }

    friend Matrix operator*(Matrix first, const Matrix& second) {
        first *= second;
        return first;
    }

    Matrix& transpose() {
        std::vector<std::vector<T>> new_data;
        new_data.assign(columns, std::vector<T>(rows));
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < columns; ++j) {
                new_data[j][i] = data[i][j];
            }
        }
        size_t tmp = rows;
        rows = columns;
        columns = tmp;
        data = new_data;
        return *this;
    }

    Matrix transposed() const {
        Matrix tmp = *this;
        tmp.transpose();
        return tmp;
    }

    std::pair<size_t, size_t> size() const {
        return {rows, columns};
    }

    template<typename U>
    std::vector<U> solve(const std::vector<U>& b) {
        std::vector<std::vector<U>> new_data(rows);
        std::vector<U> res = b, ans(columns);
        std::vector<size_t> place(columns, rows);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < columns; ++j) {
                new_data[i].push_back(static_cast<U>(data[i][j]));
            }
        }
        size_t cur_row = 0;
        for (size_t cur_col = 0; cur_col < columns; ++cur_col) {
            size_t optimal = cur_row;
            for (size_t i = cur_row; i < rows; ++i) {
                if (abs(new_data[i][cur_col]) > abs(new_data[optimal][cur_col])) {
                    optimal = i;
                }
            }
            if (new_data[optimal][cur_col] != static_cast<U>(0)) {
                for (size_t i = cur_col; i < columns; ++i) {
                    std::swap(new_data[optimal][i], new_data[cur_row][i]);
                }
                std::swap(res[optimal], res[cur_row]);
                place[cur_col] = cur_row;
                for (size_t i = 0; i < rows; ++i) {
                    U mult;
                    if (i != cur_row) {
                        mult = new_data[i][cur_col] / new_data[cur_row][cur_col];
                        for (size_t j = cur_col; j < columns; ++j) {
                            new_data[i][j] -= new_data[cur_row][j] * mult;
                        }
                        res[i] -= res[cur_row] * mult;
                    }
                }
                ++cur_row;
                if (cur_row >= rows) {
                    break;
                }
            }
        }
        for (size_t i = 0; i < columns; ++i) {
            if (place[i] != rows) {
                ans[i] = res[place[i]] / new_data[place[i]][i];
            }
        }
        return ans;
    }

    Matrix& slice(size_t st, size_t sz) {
        vector<vector<T>> new_data(rows);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < columns; ++j) {
                if (j < st || j >= st + sz) {
                    new_data[i].push_back(data[i][j]);
                }
            }
        }
        data = new_data;
        columns = columns - sz;
        return *this;
    }

    Matrix& gauss() {
        std::vector<std::vector<T>> new_data(rows);
        std::vector<size_t> place(columns, rows);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < columns; ++j) {
                new_data[i].push_back(data[i][j]);
            }
        }
        size_t cur_row = 0;
        for (size_t cur_col = 0; cur_col < columns; ++cur_col) {
            size_t optimal = cur_row;
            for (size_t i = cur_row; i < rows; ++i) {
                if (new_data[i][cur_col].absolute() > new_data[optimal][cur_col].absolute()) {
                    optimal = i;
                }
            }
            if (new_data[optimal][cur_col] != static_cast<T>(0)) {
                for (size_t i = cur_col; i < columns; ++i) {
                    std::swap(new_data[optimal][i], new_data[cur_row][i]);
                }
                place[cur_col] = cur_row;
                for (size_t i = 0; i < rows; ++i) {
                    T mult;
                    if (i != cur_row) {
                        mult = new_data[i][cur_col] / new_data[cur_row][cur_col];
                        for (size_t j = cur_col; j < columns; ++j) {
                            new_data[i][j] -= new_data[cur_row][j] * mult;
                        }
                    }
                }
                ++cur_row;
                if (cur_row >= rows) {
                    break;
                }
            }
        }
        for (size_t i = 0; i < columns; ++i) {
            if (place[i] != rows) {
                T rem = new_data[place[i]][i];
                for (size_t j = 0; j < columns; ++j) {
                    new_data[place[i]][j] /= rem;
                }
            }
        }
        data = new_data;
        return *this;
    }

    T det() const {
        if (rows == 2) {
            return data[0][0] * data[1][1] - data[0][1] * data[1][0];
        }
        T summ = 0;
        for (size_t j = 0; j < columns; ++j) {
            if (j % 2 == 0) {
                summ += data[0][j] * (remove(*this, 0, j)).det();
            } else {
                summ -= data[0][j] * (remove(*this, 0, j)).det();
            }
        }
        return summ;
    }

    friend Matrix remove(Matrix m, size_t xi, size_t xj) {
        vector<vector<T>> newdata = vector<vector<T>>(m.rows - 1, vector<T>(m.columns - 1));
        size_t ci = 0, cj = 0;
        for (size_t i = 0; i < m.rows; ++i) {
            if (i != xi) {
                for (size_t j = 0; j < m.columns; ++j) {
                    if (j != xj) {
                        newdata[ci][cj++] = m.data[i][j];
                    }
                }
                ++ci;
                cj = 0;
            }
        }
        return Matrix(newdata);
    }

    friend std::ostream& operator<<(std::ostream& out, const Matrix<T>& m) {
        for (size_t i = 0; i < m.rows; ++i) {
            if (i != 0) {
                out << "\n";
            }
            for (size_t j = 0; j < m.columns; ++j) {
                if (j != 0) {
                    out << "\t";
                }
                out << m.data[i][j];
            }
        }
        return out;
    }
};

using namespace std;

int main() {
    freopen("input.txt", "r", stdin);

    //inverse A
    int n, m;
    cin >> m >> n;
    vector<vector<Rational>> tmp1(m, vector<Rational>(n, 0));
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            int x;
            cin >> x;
            tmp1[i][j] = Rational(x);
        }
    }
    Matrix<Rational> a(tmp1);
    Matrix<Rational> ainv = a.gauss();
    ainv.slice(0, n / 2);
    cout << "\t#1\n";
    cout << "Inversed A:\n" << ainv << "\n\n";

    //B * A^(-1)
    cin >> m >> n;
    vector<vector<Rational>> tmp2(m, vector<Rational>(n, 0));
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            int x;
            cin >> x;
            tmp2[i][j] = Rational(x);
        }
    }
    Matrix<Rational> b(tmp2);
    cout << "B:\n" << b << "\n\n";
    cout << "B * A^(-1) = linear map matrix O:\n";
    Matrix<Rational> o = b * ainv;
    cout << b * ainv << "\n\n";

    cout << "Row echelon form of O:\n";
    cout << o.gauss() << "\n\n";

    cout << "\t#2\n";
    cin >> m >> n;
    vector<vector<Rational>> tmp3(m, vector<Rational>(n, 0));
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            int x;
            cin >> x;
            tmp3[i][j] = Rational(x);
        }
    }
    Matrix<Rational> A(tmp3);
    cout << "Row echelon form of linear map matrix A:\n";
    cout << A.gauss() << "\n\n";

    cin >> m >> n;
    vector<vector<Rational>> tmp4(m, vector<Rational>(n, 0));
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            int x;
            cin >> x;
            tmp4[i][j] = Rational(x);
        }
    }
    Matrix<Rational> C1(tmp4);

    cin >> m >> n;
    vector<vector<Rational>> tmp5(m, vector<Rational>(n, 0));
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            int x;
            cin >> x;
            tmp5[i][j] = Rational(x);
        }
    }
    Matrix<Rational> D(tmp5);

    cin >> m >> n;
    vector<vector<Rational>> tmp6(m, vector<Rational>(n, 0));
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            int x;
            cin >> x;
            tmp6[i][j] = Rational(x);
        }
    }
    Matrix<Rational> C2(tmp6);

    cout << "Check C1 * D * C2^(-1) == A:\n";
    cout << C1 * D * (C2.gauss().slice(0, n/2)) * 2 << "\n\n";

    cout << "\t#3\n";
    cin >> m >> n;

    vector<vector<Rational>> tmp7(m, vector<Rational>(n, 0));
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            int x;
            cin >> x;
            tmp7[i][j] = Rational(x);
        }
    }
    Matrix<Rational> P(tmp7);
    cout << "Inverse P:\n";
    Matrix<Rational> Pinv = (P.gauss().slice(0, n/2));
    cout << Pinv << "\n\n";

    cin >> m >> n;
    vector<vector<Rational>> tmp8(m, vector<Rational>(n, 0));
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            int x;
            cin >> x;
            tmp8[i][j] = Rational(x);
        }
    }
    Matrix<Rational> Eps(tmp8);
    cout << "Inverse Eps:\n";
    Matrix<Rational> Epsinv = (Eps.gauss().slice(0, n/2));
    cout << (Epsinv) << "\n\n";

    cin >> m >> n;
    vector<vector<Rational>> tmp9(m, vector<Rational>(n, 0));
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            int x;
            cin >> x;
            tmp9[i][j] = Rational(x);
        }
    }
    cout <<"Functiona Alpha:\n";
    Matrix<Rational> alpha(tmp9);
    cout << alpha * Epsinv << "\n\n";

    cin >> m >> n;
    vector<vector<Rational>> tmp10(m, vector<Rational>(n, 0));
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            int x;
            cin >> x;
            tmp10[i][j] = Rational(x);
        }
    }
    cout <<"Vector H:\n";
    Matrix<Rational> h(tmp10);
    cout <<Pinv * h << "\n\n";

    cout << "Alpha(H):\n";
    cout << (alpha * Epsinv) * (Pinv * h) << "\n\n";

}