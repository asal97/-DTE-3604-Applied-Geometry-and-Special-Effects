#ifndef MODELCURVE_H
#define MODELCURVE_H

#include "../../gmlib/modules/parametrics/gmpcurve.h"
#include "vector"

using namespace GMlib;
template <typename T>
class ModelCurve : public PCurve<T,3>
{
    GM_SCENEOBJECT(ModelCurve)
public:
    ModelCurve();
     ModelCurve(T a);
    virtual ~ModelCurve();

    bool  isClosed() const override;
    protected:
      // Virtual functions from PCurve, which have to be implemented locally
      void            eval(T t, int d, bool l) const override;
      T               getStartP() const override;
      T               getEndP()   const override;
      // control points

      T _a;

};

    template <typename T>
    ModelCurve<T>::ModelCurve(){ }

    template <typename T>
    ModelCurve<T>::ModelCurve(T a) : PCurve<T,3>(0, 0, 0), _a(a)
    {
    }

    template <typename T>
    ModelCurve<T>::~ModelCurve()
    {
    }

    template <typename T>
    bool ModelCurve<T>::isClosed() const
    {
        return true;
    }

    template <typename T>
    void ModelCurve<T>::eval(T t, int d, bool /*l*/) const
    {
        this->_p.setDim( d + 1 );
        const T x = _a * cos(t) * cos(2 * t);
        const T y = _a * sin(t) * (2+cos(2 * t));
        this->_p[0][0] = x;
        this->_p[0][1] = y;
        this->_p[0][2] = T(0);
    }

    template <typename T>
    T ModelCurve<T>::getStartP() const
    {
        return T(0);
    }

    template <typename T>
    T ModelCurve<T>::getEndP() const
    {
        return T(M_2PI);
    }

#endif // MODELCURVE_H
