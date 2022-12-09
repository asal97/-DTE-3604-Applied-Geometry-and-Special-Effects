#ifndef LANERIESENFELD_H
#define LANERIESENFELD_H

#include <vector>
#include <parametrics/gmpcurve.h>

using namespace GMlib;
template <typename T>
class LaneRiesenfeld : public PCurve<T,3>
{
   GM_SCENEOBJECT(LaneRiesenfeld)
public:
    LaneRiesenfeld(const DVector<Vector<T, 3>> &p, int d);
    ~LaneRiesenfeld();

    bool isClosed() const override ;

protected:
  void eval(T t, int d, bool l) const override ;
  T getStartP() const override ;
  T getEndP() const override ;
  void resample(std::vector<DVector<Vector<T,3>>>& p, Sphere<T,3>& s, const std::vector<T>& t, int d) const;
  int  doublePart(std::vector<DVector<Vector<T, 3>>>& p, int& n) const;
  void  smoothPartClosed(std::vector<DVector<Vector<T, 3>>>& p, int& n, int d) const;
  void          localSimulate(double dt) override;
private:
  DVector<Vector<T, 3>> _p;
  int _d;
  int _k;
};



template <typename T>
inline
LaneRiesenfeld<T>::LaneRiesenfeld(const DVector<Vector<T, 3>> &p, int d)
: PCurve<T,3>(0, 0, 0), _p(p), _d(d),_k(d+1){
    _dm = GM_DERIVATION_DD;
}

template <typename T>
LaneRiesenfeld<T>::~LaneRiesenfeld() {}

template <typename T>
inline
void LaneRiesenfeld<T>::resample(std::vector<DVector<Vector<T,3>>>& p, Sphere<T,3>& s, const std::vector<T>& t, int d) const {
    int n=_p.getDim(); // Number of intervals
    int k = this->_visu.no_sample;
     int m = n * (2 << k - 1) + 1;//The final number of points
    p.resize(m);
      s.reset();

    for(int i = 0; i < p.size(); i++)
        p[i].setDim(1);

    for(int i = 0; i < n; i++)
        p[i][0] = _p[i]; //inserting the initial points

    p[n][0] = _p[0]; //closing the curve

    for(int i = 0; i < k; i++){ //for each level of refinement
        n = doublePart(p, n);
        smoothPartClosed(p, n, d);

    }

    computeSurroundingSphere(p, s);
      if (d > _der_implemented || (d > 0 && this->_dm == GMlib::GM_DERIVATION_DD))
        GMlib::DD::compute1D(p, t, isClosed(), d, _der_implemented);
}




template <typename T>
inline
int  LaneRiesenfeld<T>::doublePart(std::vector<DVector<Vector<T, 3>>>& p, int& n) const
{
    for(int i = n; i > 0; i--) {
       p[2*i][0] = p[i][0];
       p[2*i-1][0] = 0.5*(p[i][0] + p[i-1][0]);
         }
    return 2 * n;

}
template <typename T>
inline
void  LaneRiesenfeld<T>::smoothPartClosed(std::vector<DVector<Vector<T, 3>>>& p, int& n, int d) const
{
    for (int j = 1; j < d; j++) {
      for (int i = 0; i < n; i++) {
        p[i][0] = T(0.5) * (p[i][0] + p[i + 1][0]);
      }
      p[n][0] = p[0][0];
    }
}

template <typename T>
void LaneRiesenfeld<T>::eval(T t, int d, bool l) const
{

}

template <typename T>
bool LaneRiesenfeld<T>::isClosed() const
{
    return true;
}

template <typename T>
T LaneRiesenfeld<T>::getStartP() const
{
    return 0;
}

template <typename T>
T LaneRiesenfeld<T>::getEndP() const
{
    return 1;
}
template <typename T>
void LaneRiesenfeld<T>::localSimulate( double dt ) {
    //    this->resample();
}


#endif // LANERIESENFELD_H
