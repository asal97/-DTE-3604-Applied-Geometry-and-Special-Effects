#ifndef BSPLINE_H
#define BSPLINE_H

#include <parametrics/gmpcurve.h>

#include "vector"

namespace custom {
using namespace GMlib;
template <typename T>
class Bspline : public PCurve<T,3>{
    GM_SCENEOBJECT(Bspline)
public:
     Bspline();
     Bspline(const DVector<Vector<T, 3>>& c);
     Bspline(const DVector<Vector<T,3>>& p, int n);


     virtual ~Bspline();

     void generateKnotVector(int n,T start,T end);
     int findIndex(T t) const;

    T W(int d, int i, T t) const;
    Vector<T,3> B(T t, int i) const;

    bool  isClosed() const override;
protected:
  // Virtual functions from PCurve, which have to be implemented locally
  void            eval(T t, int d, bool l) const override;
  T               getStartP() const override;
  T               getEndP()   const override;
  // control points
  DVector<Vector<T, 3>> _c;
  // knot vector
  std::vector<T>  _t;

  int      _k=3; //order
  int      _n;
  int      _d=2; // degree

private:


};



template <typename T>
Bspline<T>::Bspline() {}

template <typename T>
Bspline<T>::~Bspline(){}


template <typename T>
Bspline<T>::Bspline(const DVector<Vector<T, 3>>& c)
    : GMlib::PCurve<T,3>(0, 0, 0), _c(c)
   {
    generateKnotVector(_c.getDim(),T(0),T(1));
   }

template <typename T>
Bspline<T>::Bspline(const DVector<Vector<T,3>>& p, int n)
   : PCurve<T,3>(0, 0, 0) ,_n(n)
{

    auto x = std::vector<T>{T(0)};

    for (int i = 1; i < p.getDim(); i++)
       x.push_back((p[i] - p[i-1]).getLength() + x.back());



     generateKnotVector(n,x[0],x.back());


//    // Initialize A and filled with zeros
    auto A = DMatrix<T>(p.getDim(), n, T(0));

    for (int i = 0; i < x.size(); i++) {
      int index = findIndex(x[i]);
      auto b = B(x[i], index);
      A[i][index-2] = b[0];
      A[i][index-1] = b[1];
      A[i][index]   = b[2];
    }
    auto AT = A;
    AT.transpose();

    auto _B = AT * A;

    _B.invert();


    auto y = AT * p;


    _c = _B*y;

}

template <typename T>
Vector<T, 3> Bspline<T>::B(T t, int i) const {
    T w_1_i = W(1,i,t);
    T w_2_i_1 = W(2, i - 1, t);
    T w_2_i = W(2, i, t);
    // formula on page 95 of the book
    T b1 = (1 - w_1_i) * (1 - w_2_i_1);
    T b2 = w_1_i * (1 - w_2_i) + w_2_i_1 * (1 -w_1_i);
    T b3 = w_1_i * w_2_i;

    return Vector<T,3>{b1, b2, b3};
}
template <typename T>
T Bspline<T>::W(int d, int i, T t) const
{
    return (t - _t[i]) / (_t[i + d] - _t[i]);
}

template <typename T>
void Bspline<T>::eval(T t, int d, bool /*l*/) const
{

    this->_p.setDim( d + 1 );

    int i = findIndex(t);

    //get basis
    auto a = (1 - W(1, i, t)) * (1 - W(2, i - 1, t));
    auto b = (1 - W(1, i, t)) * W(2, i - 1, t) + (W(1, i, t) * (1 - W(2, i, t)));
    auto c =  W(1, i, t) * W(2, i, t);


    this->_p[0] = a * _c[i - 2] + b * _c[i - 1] + c * _c[i];
}

template <typename T>
inline void Bspline<T>::generateKnotVector(int n,T start, T end) {

    _t.clear();
    _t.push_back(start);

    _t.push_back(start);

    n = n - 2;


    T delta = (end - start) /n;

    for (int i = 0; i <= n; i++) {
    _t.push_back(start + i * delta);
    }


     _t.push_back(end);
     _t.push_back(end);

}
template <typename T>
int Bspline<T>::findIndex(T t) const
 {
    if (t >= this->getEndP()){
        return _t.size() - 4;
    }
    else{
        return std::distance(_t.begin(), std::upper_bound(_t.begin(), _t.end(), t)) -1;
    }
 }



template <typename T>
bool Bspline<T>::isClosed() const {return false;}

template <typename T>
T Bspline<T>::getStartP() const
{
    return _t[_d];
}

template <typename T>
T Bspline<T>::getEndP() const
{
    return _t[_t.size() - 3];
}




}

#endif // BSPLINE_H
