#ifndef BLENDINGSPLINESURFACE_H
#define BLENDINGSPLINESURFACE_H

using namespace GMlib;


#include <parametrics/curves/gmpsubcurve.h>
#include "../../gmlib/modules/parametrics/gmpcurve.h"
#include "PSimpleSubSurf.h"

template <typename T>
class BlendingSplineSurface : public PSurf<T, 3> {


    GM_SCENEOBJECT(BlendingSplineSurface)

    public:

    BlendingSplineSurface(PSurf<T, 3> *surface, int n_u,int n_v);
    virtual ~BlendingSplineSurface(){};




  protected:
    void eval(T u,T v, int du,int dv, bool /*l*/,bool /*l*/) const override;

    T getStartPU() const override;
    T getEndPU() const override;

    T getStartPV() const override;
    T getEndPV() const override;

    bool isClosedU() const override;
    bool isClosedV() const override;

    T W(int d, int i, T t,const std::vector<T>& knotVec) const;
    Vector<T, 2> B(const std::vector<T> & knotVec,int i, T t) const;

    T dB(const std::vector<T> & knotVec,int i, T t) const;
    T dW(int d, int i,const std::vector<T>& knotVec) const;

    int findIndex(T t,bool u) const;
    void generateKnotVector(int n, std::vector<T> &t,T start, T end,bool closed);
    void makeLocalSurface(PSurf<T, 3>* surface);



    void localSimulate(double dt) override;


  private:
    DMatrix<PSurf<T,3> *> _s;
    std::vector<T> _tu;
    std::vector<T> _tv;

    T u_start;
    T u_end;

    T v_start;
    T v_end;

    bool _closedU;
    bool _closedV;
  };

template <typename T>
BlendingSplineSurface<T>::BlendingSplineSurface(PSurf<T, 3> *surface, int n_u,int n_v)
    :  u_start(surface->getParStartU()),
      u_end(surface->getParEndU()), v_start(surface->getParStartV()),
      v_end(surface->getParEndV()),
      _closedU(surface->isClosedU()), _closedV(surface->isClosedV()) {



  generateKnotVector(n_u, _tu ,u_start, u_end,_closedU);
  generateKnotVector(n_v, _tv ,v_start, v_end,_closedV);

  makeLocalSurface(surface);

}

template <typename T>
void BlendingSplineSurface<T>::makeLocalSurface(PSurf<T, 3>* surface){

    int nu = _tu.size() - 2 ;
    int nv = _tv.size() - 2  ;
    _s.setDim((_closedU ? (nu+1) : nu), (_closedV ? (nv+1) : nv));


    for (int i = 0; i < nu; i++) {
        for (int j = 0; j < nv; j++){
            _s[i][j] = new PSimpleSubSurf<T>(surface,_tu[i],_tu[i+2],_tu[i+1],_tv[j],_tv[j+2],_tv[j+1]);
            this->insert(_s[i][j]);
            _s[i][j]->toggleDefaultVisualizer();
            _s[i][j]->sample(5, 5,1,1);
            _s[i][j]->setCollapsed(true);
        }

        if (_closedV)
            _s[i][nv] = _s[i][0];
    }
    if (_closedU)
        for (int j  = 0; j < _s.getDim2(); j++)
            _s[nu][j] = _s[0][j];




}

template <typename T>
void BlendingSplineSurface<T>::generateKnotVector(int n, std::vector<T> &t,T start, T end,bool closed) {


        t.push_back(start);

        n = n - (closed ? 0 : 1);


        T delta = (end - start) /n;

        for (int i = 0; i <= n; i++) {
        t.push_back(start + i * delta);
        }

        if(closed){
          t[0] = start - delta;
        }
        if(!closed){
            t.push_back(end);
        }

    }

template <typename T> T BlendingSplineSurface<T>::getStartPU() const
{
  return u_start;
}

template <typename T> T BlendingSplineSurface<T>::getEndPU() const
{ return u_end; }

template <typename T> T BlendingSplineSurface<T>::getStartPV() const
{
  return v_start;
}

template <typename T> T BlendingSplineSurface<T>::getEndPV() const
{ return v_end; }


template <typename T>
void BlendingSplineSurface<T>::eval(T u,T v, int du,int dv, bool l1,  bool l2) const
{
    this->_p.setDim( du + 1 ,dv + 1);

    int i_u = findIndex(u,true);
    int i_v = findIndex(v,false);


    auto B_u = B(_tu,i_u,u);

    auto B_v = B(_tv,i_v,v);

    auto Bpart1 = B_u[0] * B_v[0];
    auto Bpart2 = B_u[1] * B_v[0];
    auto Bpart3 = B_u[0] * B_v[1];
    auto Bpart4 = B_u[1] * B_v[1];


    auto p1 = Bpart1 * ((_s(i_u-1)(i_v-1))->evaluateParent(u,v, du, dv))(0)(0);
    auto p2 = Bpart2 * (_s(i_u)(i_v-1)->evaluateParent(u,v, du, dv))(0)(0);
    auto p3 =Bpart3* (_s(i_u-1)(i_v)->evaluateParent(u,v, du, dv))(0)(0);
    auto p4 = Bpart4 * (_s(i_u)(i_v)->evaluateParent(u,v, du, dv))(0)(0);

    this->_p[0][0] = p1 + p2 + p3 + p4;


    auto B_up = dB(_tu,i_u,u);
    auto p1PrimeU = -1* B_up * B_v[0] * ((_s(i_u-1)(i_v-1)->evaluateParent(u,v, du, dv))(0)(0)) + Bpart1* ((_s(i_u-1)(i_v-1)->evaluateParent(u,v, du, dv))(0)(1));
    auto p2PrimeU =  B_up * B_v[0] * ((_s(i_u)(i_v-1)->evaluateParent(u,v, du, dv))(0)(0)) +Bpart2 * ((_s(i_u)(i_v-1)->evaluateParent(u,v, du, dv))(0)(1));
    auto p3PrimeU = -1 *  B_up * B_v[1] * ((_s(i_u-1)(i_v)->evaluateParent(u,v, du, dv))(0)(0)) + Bpart3 *((_s(i_u-1)(i_v)->evaluateParent(u,v, du, dv))(0)(1));
    auto p4PrimeU =  B_up * B_v[1] * ((_s(i_u)(i_v)->evaluateParent(u,v, du, dv))(0)(0)) + Bpart4 *((_s(i_u)(i_v)->evaluateParent(u,v, du, dv))(0)(1));

    this->_p[0][1] = p1PrimeU + p2PrimeU + p3PrimeU + p4PrimeU;


    auto B_vp = dB(_tv,i_v,v);
    auto p1PrimeV = B_u[0] * -1 * B_vp * ((_s(i_u-1)(i_v-1)->evaluateParent(u,v, du, dv))(0)(0)) + Bpart1 * ((_s(i_u-1)(i_v-1)->evaluateParent(u,v, du, dv))(1)(0));
    auto p2PrimeV =   B_u[1] *  -1 * B_vp * ((_s(i_u)(i_v-1)->evaluateParent(u,v, du, dv))(0)(0)) +Bpart2 * ((_s(i_u)(i_v-1)->evaluateParent(u,v, du, dv))(1)(0));
    auto p3PrimeV =  B_u[0]* B_vp * ((_s(i_u-1)(i_v)->evaluateParent(u,v, du, dv))(0)(0)) +Bpart3 *((_s(i_u-1)(i_v)->evaluateParent(u,v, du, dv))(1)(0));
    auto p4PrimeV =   B_u[1]* B_vp * ((_s(i_u)(i_v)->evaluateParent(u,v, du, dv))(0)(0)) + Bpart4 * ((_s(i_u)(i_v)->evaluateParent(u,v, du, dv))(1)(0));

    this->_p[1][0] = p1PrimeV + p2PrimeV + p3PrimeV + p4PrimeV;

}

template <typename T>
int BlendingSplineSurface<T>::findIndex(T t,bool u) const
 {
    auto _t = ((u)? _tu : _tv);

    int _n = (u ?(_closedU ? _tu.size()-1: _tu.size()-2) : (_closedV ? _tv.size()-1: _tv.size()-2));

    if(t == _t[_n])
        return _n - 1;

    for(int i = 1; i < _n; i++){
        if(t < _t[i + 1])
            return i;
    }


 }


template <typename T>
T BlendingSplineSurface<T>::W(int d, int i, T t,const std::vector<T>& knotVec) const
{
    return (t - knotVec[i]) / (knotVec[i + d] - knotVec[i]);
}
template <typename T>
Vector<T, 2> BlendingSplineSurface<T>::B(const std::vector<T> & knotVec,int i, T t) const
{
    T w_1_i = W(1,i,t,knotVec);
    T B = w_1_i - T(1.0 / (M_2PI)) * sin(M_2PI * w_1_i);
    Vector<T,2> Bf = {1-B,B};
    return Bf;
}
template <typename T>
T BlendingSplineSurface<T>::dB(const std::vector<T> & knotVec,int i, T t) const
{
    T w_1_i = dW(1,i,knotVec);
    T Bprime = 1 - cos(M_2PI * w_1_i);
    return Bprime;
}
template <typename T>
T BlendingSplineSurface<T>::dW(int d, int i,const std::vector<T>& knotVec) const
{
    return 1/ (knotVec[i + d] - knotVec[i]);
}

template <typename T>
bool BlendingSplineSurface<T>::isClosedU() const {
  return _closedU;
}

template <typename T>
bool BlendingSplineSurface<T>::isClosedV() const {
  return _closedV;
}




template <typename T>
void BlendingSplineSurface<T>::localSimulate(double dt){


        this->sample(this->_visu[0][0][0],this->_visu[0][0][1],1,1);
        this->setEditDone();
}
#endif // BLENDINGSPLINESURFACE_H
