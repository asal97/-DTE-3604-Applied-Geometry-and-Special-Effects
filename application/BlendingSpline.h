#ifndef BLENDINGSPLINE_H
#define BLENDINGSPLINE_H

#include <parametrics/curves/gmpsubcurve.h>
#include "../../gmlib/modules/parametrics/gmpcurve.h"


using namespace GMlib;

template <typename T>
class BlendingSpline : public PCurve<T, 3> {

  GM_SCENEOBJECT(BlendingSpline)

  public:

  BlendingSpline(PCurve<T, 3>* curve, int n);
  virtual ~BlendingSpline(){};


protected:
    void eval(T t, int d, bool left) const override;
    T getStartP() const override;
    T getEndP() const override;

    bool isClosed() const override;


    T W(int d, int i, T t) const;
    Vector<T, 2> B(T t, int i) const;
    T blend(T w) const;
    int findIndex(T t) const;
    void generateKnotVector(int n,T start, T end);
    void makeLocalCurve(PCurve<T,3>* curve);

    void localSimulate(double dt) override;


private:

    std::vector<PSubCurve<T> *> _c;
    std::vector<T> _t;

    T _t_start;
    T _t_end;
    bool _closed;
    double _angle=0.0;
    double _rot = 0.01;
    int flag = 5;
    int flip = -1;
    int count = 0;
    T rythm = 13;
    float   _animSpeed = 0.02;
    const int range_from  = 0;
    const int range_to    = 255;

};

template <typename T>
BlendingSpline<T>::BlendingSpline(PCurve<T, 3> *curve, int n)
    : PCurve<T, 3>(0, 0, 0), _t_start{curve->getParStart()},
      _t_end{curve->getParEnd()}, _closed{curve->isClosed()} {

    generateKnotVector(n, _t_start, _t_end);
    makeLocalCurve(curve);

}
template <typename T>
void BlendingSpline<T>::makeLocalCurve(PCurve<T,3>* curve){


    int n = _t.size() - 2;
    _c.resize(n+1);

    for(int i=0;i<n;i++){
        _c[i] = new PSubCurve<T>(curve,_t[i],_t[i+2],_t[i+1]);
        _c[i]->toggleDefaultVisualizer();
//        _c[i]->sample(5, 5);
        _c[i]->setCollapsed(true);
        _c[i]->setParent(this);
//        this->insert(_c[i]);


    }

    _c[n] = _c[0];



}

template <typename T>
void BlendingSpline<T>::eval(T t, int d, bool /*l*/) const
{

    this->_p.setDim( d + 1 );

    int i = findIndex(t);

    //get basis
    auto b = B(t,i);


    this->_p = b[0] *( _c[i - 1]->evaluateParent(t,d)) + b[1] * (_c[i]->evaluateParent(t,d));
}
template <typename T>
void BlendingSpline<T>::generateKnotVector(int n, T start, T end) {

    std::cout<<"n"<<n<<std::endl;
    _t.push_back(start);

    T delta = (end - start) /n;

    for (int i = 0; i <= n; i++) {
    _t.push_back(start + i * delta);
    }


    _t[0] = start - delta;

    std::cout<<_t<<std::endl;

}
template <typename T>
int BlendingSpline<T>::findIndex(T t) const
 {
    if (t == this->getEndP()) {
        return _t.size() - 2;
      }
    else {
        return std::distance(_t.begin(),
                             std::upper_bound(
                                 _t.begin(),
                                 _t.end(),
                                 t)) - 1;
      }
 }


template <typename T>
T BlendingSpline<T>::W(int d, int i, T t) const
{
    return (t - _t[i]) / (_t[i + d] - _t[i]);
}

template <typename T>
Vector<T, 2> BlendingSpline<T>::B(T t, int i) const {

    T w_1_i = W(1,i,t);
    //symmetric TB-functions order 2 B2(t) = t - (1/(2pi))* sin(2pi * t)
    T B = w_1_i - T(1.0 / (M_2PI)) * sin(M_2PI * w_1_i);
    Vector<T,2> Bf = {1-B,B};
    return Bf;
}

template <typename T>
void BlendingSpline<T>::localSimulate(double dt){

        float ang = 0.05;

       for (int i = 0; i < _c.size() - 1; i+= 5) {
         auto dir = _c[i]->getGlobalPos() - Vector<T,3>(0, 0, 0);
         _c[i]->rotateParent(ang, dir.normalize());
       }

       count ++;
       if (count%45==0)
                 flip = -1*flip;

       if (count%100==0){
           std::random_device                  rand_dev;
           std::mt19937                        generator(rand_dev());
           std::uniform_int_distribution<int>  distr(range_from, range_to);
           this->setColor(Color(distr(generator), distr(generator) ,distr(generator)));
       }

       double r = 0.000001 + (1.0 - (1.0/ M_PI)*(2*cos(_angle) - 2*cos(2*_angle))) * dt;

       for (int i = 0; i < _c.size() - 1; i+= 1) {
         auto dir = Vector<T,3>(flip*cos(rythm * M_PI ),flip* sin(rythm * M_PI), 0);
         _c[i]->translate( _animSpeed * dir);
         flag -= 1;

         _c[i]->scale(Point<T,3>(1, 1, 1) * (1 + 0.01 * sin(_angle)));
       }
       _angle += r;
       if (_angle >= M_2PI){_angle -= M_2PI;}
        this->resample();
        this->setEditDone();
}

template <typename T> T BlendingSpline<T>::getStartP() const {
  return _t_start;
}

template <typename T> T BlendingSpline<T>::getEndP() const { return _t_end; }

template <typename T> bool BlendingSpline<T>::isClosed() const {
  return _closed;
}

#endif // BLENDINGSPLINE_H
