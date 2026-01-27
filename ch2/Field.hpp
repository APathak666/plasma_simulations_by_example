#ifndef _FIELD_H
#define _FIELD_H

#include <ostream>
#include <iostream>

template <typename T>
struct vec3 {
	vec3 (const T u, const T v, const T w) : d{u,v,w} {}
	vec3 (const T a[3]) : d{a[0],a[1],a[2]} {}
	vec3 (): d{ 0,0,0 } {}
	T& operator [] (int i) { return d[i]; }
	T operator () (int i) const { return d[i]; }
	vec3<T>& operator = (double s) { d[0] = s; d[1] = s; d[2] = s; return (*this);}
	vec3<T>& operator += (vec3<T> o) { d[0] += o[0]; d[1] += o[1]; d[2] += o[2]; return(*this); }
	vec3<T>& operator -= (vec3<T> o) { d[0] -= o[0]; d[1] -= o[1]; d[2] -= o[2]; return(*this); }

    protected:
        T d[3];
};

using double3 = vec3<double>;
using int3 = vec3<int>;

//vec3-vec3 operations
template<typename T>	//addition of two vec3s
vec3<T> operator+(const vec3<T>& a, const vec3<T>& b) {
	return vec3<T> (a(0)+b(0),a(1)+b(1),a(2)+b(2));	}
template<typename T>	//subtraction of two vec3s
vec3<T> operator-(const vec3<T>& a, const vec3<T>& b) {
	return vec3<T> (a(0)-b(0),a(1)-b(1),a(2)-b(2));	}
template<typename T>	//element-wise multiplication of two vec3s
vec3<T> operator*(const vec3<T>& a, const vec3<T>& b) {
	return vec3<T> (a(0)*b(0),a(1)*b(1),a(2)*b(2));	}
template<typename T>	//element wise division of two vec3s
vec3<T> operator/(const vec3<T>& a, const vec3<T>& b) {
	return vec3<T> (a(0)/b(0),a(1)/b(1),a(2)/b(2));	}

//vec3 - scalar operations
template<typename T>		//scalar multiplication
vec3<T> operator*(const vec3<T> &a, T s) {
	return vec3<T>(a(0)*s, a(1)*s, a(2)*s);}
template<typename T>		//scalar multiplication 2
vec3<T> operator*(T s,const vec3<T> &a) {
	return vec3<T>(a(0)*s, a(1)*s, a(2)*s);}

//output
template<typename T>	//ostream output
std::ostream& operator<<(std::ostream &out, vec3<T>& v) {
	out<<v[0]<<" "<<v[1]<<" "<<v[2];
	return out;
}

template <typename T>
class Field_ {
    public:
        Field_(int ni, int nj, int nk) : ni{ni}, nj{nj}, nk{nk} {
            data = new T **[ni];
            for (int i = 0; i < ni; i++) {
                data[i] = new T *[nj];
                for (int j = 0; j < nj; j++) {
                    data[i][j] = new T [nk];
                }
            }

            (*this) = 0;
        }

        Field_(const Field_ &other) :
        Field_{other.ni, other.nj, other.nk} {
            for (int i = 0; i < ni; i++)
                for (int j = 0; j < nj; j++)
                    for (int k = 0; k < nk; k++)
                        data[i][j][k] = other(i, j, k);
        }

        Field_(Field_ &&other) :
        ni{other.ni}, nj{other.nj}, nk{other.nk} {
            data = other.data;
            other.data = nullptr;
        }

        Field_& operator = (Field_ &&f) {
            if (data) {
                for (int i = 0; i < ni; i++) {
                    for (int j = 0; j < nj; j++) {
                        delete[] data[i][j];
                    }
                    delete[] data[i];
                }
                delete[] data;
            }
            data = f.data;
            f.data = nullptr;
            return *this;
        }

        ~Field_() {
            if (data == nullptr) return;
            for (int i = 0; i < ni; i++) {
                for (int j = 0; j < nj; j++) {
                    delete[] data[i][j];
                }
                delete[] data[i];
            }
            delete[] data;
            data = nullptr;
        }

        T ** operator [] (int i) { return data[i]; }

        T operator () (int i, int j, int k) const { return data[i][j][k]; }

        Field_& operator = (double s) {
            for (int i = 0; i < ni; i++)
                for (int j = 0; j < nj; j++)
                    for (int k = 0; k < nk; k++)
                        data[i][j][k] = s;

            return (*this);
        }

        void operator /= (const Field_ &other) {
            for (int i = 0; i < ni; i++)
                for (int j = 0; j < nj; j++)
                    for (int k = 0; k < nk; k++) {
                        if (other.data[i][j][k] != 0) data[i][j][k] /= other(i, j, k);
                        else data[i][j][k] = 0;
                    }
        }

        Field_& operator += (const Field_ &other) {
            for (int i = 0; i < ni; i++)
                for (int j = 0; j < nj; j++)
                    for (int k = 0; k < nk; k++)
                        data[i][j][k] += other(i, j, k);

            return (*this);
        }

        Field_& operator *= (double s) {
            for (int i = 0; i < ni; i++)
                for (int j = 0; j < nj; j++)
                    for (int k = 0; k < nk; k++)
                        data[i][j][k] *= s;

            return (*this);
        }

        friend Field_<T> operator * (double s, const Field_<T>&f) {
            Field_<T> r(f);
            return std::move(r *= s);
        }

    	template<typename S>
    	friend std::ostream& operator << (std::ostream &out, Field_<S> &f);

        void scatter(double3 lc, T value) {
            double3 dl;
            for (int i = 0; i < 3; i++) dl[i] = lc[i] - (int) lc[i];

            int i = (int) lc[0];
            int j = (int) lc[1];
            int k = (int) lc[2];

            data[i][j][k] += value*(1-dl[0])*(1-dl[1])*(1-dl[2]);
            data[i+1][j][k] += value*(dl[0])*(1-dl[1])*(1-dl[2]);
            data[i][j+1][k] += value*(1-dl[0])*(dl[1])*(1-dl[2]);
            data[i+1][j+1][k] += value*(dl[0])*(dl[1])*(1-dl[2]);
            data[i][j][k+1] += value*(1-dl[0])*(1-dl[1])*(dl[2]);
            data[i+1][j][k+1] += value*(dl[0])*(1-dl[1])*(dl[2]);
            data[i][j+1][k+1] += value*(1-dl[0])*(dl[1])*(dl[2]);
            data[i+1][j+1][k+1] += value*(dl[0])*(dl[1])*(dl[2]);
        }

        T gather(double3 lc) {
            int i = (int) lc[0];
            int j = (int) lc[1];
            int k = (int) lc[2];

            double di = lc[0] - i;
            double dj = lc[1] - j;
            double dk = lc[2] - k;

            return (
                data[i][j][k]*(1-di)*(1-dj)*(1-dk) + 
                data[i+1][j][k]*(di)*(1-dj)*(1-dk) + 
                data[i][j+1][k]*(1-di)*(dj)*(1-dk) + 
                data[i+1][j+1][k]*(di)*(dj)*(1-dk) + 
                data[i][j][k+1]*(1-di)*(1-dj)*(dk) + 
                data[i+1][j][k+1]*(di)*(1-dj)*(dk) + 
                data[i][j+1][k+1]*(1-di)*(dj)*(dk) + 
                data[i+1][j+1][k+1]*(di)*(dj)*(dk) 
            );
        }

        const int ni, nj, nk;

    protected:
        T ***data;
};

template<typename T>
std::ostream& operator<<(std::ostream &out, Field_<T> &f)
{
	for (int k=0;k<f.nk;k++,out<<"\n")
		for (int j=0;j<f.nj;j++)
			for (int i=0;i<f.ni;i++) out<<f.data[i][j][k]<<" ";
	return out;
}

using Field = Field_<double>;
using Field3 = Field_<double3>;

#endif
