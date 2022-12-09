//#include "modelcurve.h"


//    template <typename T>
//    ModelCurve<T>::ModelCurve(){ }

//    template <typename T>
//    ModelCurve<T>::ModelCurve(T radius, T c) : PCurve<T,3>(20, 0, 0), _r(radius), _c(c)
//    {
//    }

//    template <typename T>
//    ModelCurve<T>::~ModelCurve()
//    {
//    }

//    template <typename T>
//    bool ModelCurve<T>::isClosed() const
//    {
//        return true;
//    }

//    template <typename T>
//    void ModelCurve<T>::eval(T t, int d, bool /*l*/) const
//    {
//        this->_p.setDim( d + 1 );
//        const T x = _r * cos(t) * cos(_c * t);
//        const T y = _r * sin(t) * cos(_c * t);
//        this->_p[0][0] = x;
//        this->_p[0][1] = y;
//        this->_p[0][2] = T(0);
//    }

//    template <typename T>
//    T ModelCurve<T>::getStartP() const
//    {
//        return T(0);
//    }

//    template <typename T>
//    T ModelCurve<T>::getEndP() const
//    {
//        return T(5* M_2PI);
//    }

