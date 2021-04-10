#include <vector>
#include <map>
#include <ostream>
#include <algorithm>
#include <iostream>

template <typename T>
class Polynomial {
 private:
    std::map<size_t, T> data;

 public:
    void format() {
        std::vector<size_t> del;
        for (auto[deg, coef] : data) {
            if (coef == static_cast<T>(0)) {
                del.push_back(deg);
            }
        }
        for (size_t i : del) {
            data.erase(i);
        }
    }

    Polynomial(const std::vector<T>& coefs) {
        for (size_t i = 0; i < coefs.size(); ++i) {
            data[i] = coefs[i];
        }
        format();
    }

    Polynomial(const T coef) {
        data[0] = coef;
        format();
    }

    Polynomial() {
        data[0];
        format();
    }

    bool operator==(const Polynomial other) const {
        return data == other.data;
    }

    bool operator!=(const Polynomial other) const {
        return !(*this == other);
    }

    Polynomial& operator*=(const Polynomial& other) {
        Polynomial tmp;
        for (auto[deg1, coef1] : data) {
            for (auto[deg2, coef2] : other.data) {
                tmp.data[deg1 + deg2] += coef1 * coef2;
            }
        }
        *this = tmp;
        format();
        return *this;
    }

    friend Polynomial operator*(Polynomial first, const Polynomial& second) {
        first *= second;
        return first;
    }

    Polynomial& operator+=(const Polynomial& other) {
        Polynomial tmp;
        for (auto[deg, coef] : data) {
            tmp.data[deg] += coef;
        }
        for (auto[deg, coef] : other.data) {
            tmp.data[deg] += coef;
        }
        *this = tmp;
        format();
        return *this;
    }

    friend Polynomial operator+(Polynomial first, const Polynomial& second) {
        first += second;
        return first;
    }

    Polynomial& operator-=(const Polynomial other) {
        Polynomial tmp;
        for (auto[deg, coef] : data) {
            tmp.data[deg] += coef;
        }
        for (auto[deg, coef] : other.data) {
            tmp.data[deg] -= coef;
        }
        *this = tmp;
        format();
        return *this;
    }

    friend Polynomial operator-(Polynomial first, const Polynomial& second) {
        first -= second;
        return first;
    }

    Polynomial& operator/=(const Polynomial other) {
        Polynomial div;
        if (other.Degree() == -1) {
            *this = div;
            format();
            return *this;
        }
        while (this->Degree() >= other.Degree()) {
            int old_deg = this->Degree();
            T coef = (data.rbegin())->second / (other.data.rbegin())->second;
            int deg = this->Degree() - other.Degree();
            Polynomial tmp;
            tmp.data[deg] = coef;
            div += tmp;
            (*this) -= other * tmp;
            if (this->Degree() == old_deg) {
                data.erase(data.rbegin()->first);
            }
        }
        *this = div;
        format();
        return *this;
    }

    friend Polynomial operator/(Polynomial first, const Polynomial& second) {
        first /= second;
        return first;
    }

    Polynomial& operator%=(const Polynomial other) {
        Polynomial div;
        if (other.Degree() == -1) {
            format();
            return *this;
        }
        while (this->Degree() >= other.Degree()) {
            int old_deg = this->Degree();
            T coef = (data.rbegin())->second / (other.data.rbegin())->second;
            int deg = this->Degree() - other.Degree();
            Polynomial tmp;
            tmp.data[deg] = coef;
            div += tmp;
            (*this) -= other * tmp;
            if (this->Degree() == old_deg) {
                data.erase(data.rbegin()->first);
            }
        }
        format();
        return *this;
    }

    friend Polynomial operator%(Polynomial first, const Polynomial& second) {
        first %= second;
        return first;
    }

    friend Polynomial operator,(Polynomial first, Polynomial second) {
        if (second.Degree() == -1) {
            if (first.Degree() == -1) {
                return first;
            }
            T mult = first.data.rbegin()->second;
            for (auto&[deg, coef] : first.data) {
                coef /= mult;
            }
            return first;
        }
        return (second, first % second);
    }

    T operator[](size_t i) const {
        if (data.count(i) != 0) {
            return data.at(i);
        } else {
            return static_cast<T>(0);
        }
    }

    T operator()(T x) const {
        T summ = static_cast<T>(0);
        T pow = static_cast<T>(1);
        size_t curdeg = 0;
        for (auto[deg, coef] : data) {
            while (curdeg < deg) {
                pow *= x;
                ++curdeg;
            }
            summ += coef * pow;
        }
        return summ;
    }

    friend Polynomial operator&(const Polynomial& first, const Polynomial& second) {;
        Polynomial tmp(static_cast<T>(0));
        Polynomial mult(static_cast<T>(1));
        size_t curdeg = 0;
        for (auto[deg, coef] : first.data) {
            while (curdeg < deg) {
                mult *= second;
                ++curdeg;
            }
            if (coef != static_cast<T>(0)) {
                tmp += (coef * mult);
            }
        }
        tmp.format();
        return tmp;
    }

    int Degree() const {
        if (data.size() == 0) {
            return -1;
        } else {
            return data.rbegin()->first;
        }
    }

    friend std::ostream& operator<<(std::ostream& out, const Polynomial& m) {
        if (m.Degree() == -1) {
            out << "0";
            return out;
        }
        for (auto i = m.data.crbegin(); i != m.data.crend(); ++i) {
            if (i->second == static_cast<T>(0)) {
                continue;
            }
            if (i->second > static_cast<T>(0) && i->first != (m.data.rbegin()->first)) {
                out << "+";
            } else if (i->second < static_cast<T>(0)) {
                if (i->second == static_cast<T>(-1) && i->first != 0) {
                    out << "-";
                }
            }
            if (i->first > 0 && i->second != static_cast<T>(1) && i->second != static_cast<T>(-1)) {
                out << i->second << "*";
            }
            if (i->first == 1) {
                out << "x";
            } else if (i->first > 1) {
                out << "x^" << i->first;
            } else if (i->first == 0) {
                out << i->second;
            }
        }
        return out;
    }

    const typename std::map<size_t, T>::const_iterator begin() const {
        return data.begin();
    }

    const typename std::map<size_t, T>::const_iterator end() const {
        return data.end();
    }
};