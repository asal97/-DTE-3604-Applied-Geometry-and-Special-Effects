
#include <iostream>


#include "scenario.h"
#include "testtorus.h"
#include "parametrics/curves/gmpcircle.h"

#include "parametrics/surfaces/gmpcylinder.h"
#include "parametrics/surfaces/gmptorus.h"
#include "parametrics/surfaces/gmpplane.h"
#include "application/Bspline.h"
#include "application/laneriesenfeld.h"
#include "application/modelcurve.h"
#include "application/BlendingSpline.h"
#include "application/BlendingSplineSurface.h"
// hidmanager
#include "hidmanager/defaulthidmanager.h"

// gmlib
#include <scene/light/gmpointlight.h>
#include <scene/sceneobjects/gmpathtrack.h>
#include <scene/sceneobjects/gmpathtrackarrows.h>


// qt
#include <QQuickItem>

template <typename T>
inline std::ostream &operator<<(std::ostream &out, const std::vector<T> &v) {
  out << v.size() << std::endl;
  for (uint i = 0; i < v.size(); i++)
    out << " " << v[i];
  out << std::endl;
  return out;
}

void Scenario::initializeScenario() {

  // Insert a light
  GMlib::Point<GLfloat, 3> init_light_pos(2.0, 4.0, 10);
  GMlib::PointLight *light =
      new GMlib::PointLight(GMlib::GMcolor::white(), GMlib::GMcolor::white(),
                            GMlib::GMcolor::white(), init_light_pos);
  light->setAttenuation(0.8f, 0.002f, 0.0008f);
  this->scene()->insertLight(light, false);

  // Insert Sun
  this->scene()->insertSun();

  // Default camera parameters
  int init_viewport_size = 600;
  GMlib::Point<float, 3> init_cam_pos(0.0f, 0.0f, 0.0f);
  GMlib::Vector<float, 3> init_cam_dir(0.0f, 1.0f, 0.0f);
  GMlib::Vector<float, 3> init_cam_up(1.0f, 0.0f, 0.0f);

  // Projection cam
  auto proj_rcpair = createRCPair("Projection");
  proj_rcpair.camera->set(init_cam_pos, init_cam_dir, init_cam_up);
  proj_rcpair.camera->setCuttingPlanes(1.0f, 8000.0f);
  proj_rcpair.camera->rotateGlobal(GMlib::Angle(-45),
                                   GMlib::Vector<float, 3>(1.0f, 0.0f, 0.0f));
  proj_rcpair.camera->translateGlobal(
      GMlib::Vector<float, 3>(0.0f, -20.0f, 20.0f));
  scene()->insertCamera(proj_rcpair.camera.get());
  proj_rcpair.renderer->reshape(
      GMlib::Vector<int, 2>(init_viewport_size, init_viewport_size));

  /***************************************************************************
   *                                                                         *
   * Standar example, including path track and path track arrows             *
   *                                                                         *
   ***************************************************************************/

  GMlib::Material mm(GMlib::GMmaterial::polishedBronze());
  mm.set(45.0);

  /// Bspline
//  GMlib::DVector<GMlib::Vector<float, 3>> c(10);
//  c[0] = {0, 3, -5};
//  c[1] = {1, 0, 0};
//  c[2] = {3, 5, 0};
//  c[3] = {5, 1, 0};
//  c[4] = {3, -3, 5};
//  c[5] = {1, -5, 0};
//  c[6] = {3, 5, -5};
//  c[7] = {10, 9, 0};
//  c[8] = {7, -3, 0};
//  c[9] = {13, -5, 5};


//  auto* b_spline_curve_control = new custom::Bspline<float>(c);
//  b_spline_curve_control->toggleDefaultVisualizer();
//  b_spline_curve_control->sample(200, 0);
//  this->scene()->insert(b_spline_curve_control);
//  b_spline_curve_control->setLineWidth(4);

//  auto spline3 = new custom::Bspline<float>(c, 5);
//  spline3->setColor(GMlib::GMcolor::yellow());
//  spline3->toggleDefaultVisualizer();
//  spline3->sample(200, 0);
//  this->scene()->insert(spline3);
//  spline3->setLineWidth(4);
//  spline3->translate({0, 20, 0});




   /// LaneRiesenfeld
   ///
//  GMlib::DVector<GMlib::Vector<float,3>> P(8);
//  P[0] = {0, 0, 0}; P[1] = {1, 1, 0};
//  P[2] = {2, 0, 0}; P[3] = {3, 2, 0};
//  P[4] = {4, 1, 0}; P[5] = {5, 1, 0};
//  P[6] = {6, 2, 0}; P[7] = {7, 0 ,0};

//  auto* subdivisionCurve = new LaneRiesenfeld<float>(P, 2);
//  subdivisionCurve->toggleDefaultVisualizer();
//  subdivisionCurve->sample(20,2);
//  this->scene()->insert(subdivisionCurve);


  /// BlendingSpline

  /// ModelCurve

   auto CORNOID = new ModelCurve<float>(10);
//   CORNOID->toggleDefaultVisualizer();
//   CORNOID->sample(200);
//   this->scene()->insert(CORNOID);




   auto* Blending = new BlendingSpline<float>(CORNOID,10);
   Blending->toggleDefaultVisualizer();
   Blending->sample(200, 0);
   this->scene()->insert(Blending);
//   Blending->translate({0, 40, 0});
   Blending->setLineWidth(4);










//  auto plane = new GMlib::PPlane<float>(GMlib::Point<float, 3>(0.0f, 0.0f,
//        0.0f),
//                                   GMlib::Vector<float, 3>(6.0f, 0.0f,
//                                   0.0f), GMlib::Vector<float, 3>(0.0f,
//                                   0.0f, 4.0f));

//    auto Torus = new PTorus<float>(3,1,1);
//  auto cylinder = new PCylinder<float>(5,5,5);
//    auto plane = new PPlane<float>(3,3,3);
//    auto* Blending = new BlendingSplineSurface<float>(Torus,5,5);
//    Blending->toggleDefaultVisualizer();
//    Blending->sample(20, 20,1,1);
//    Blending->setMaterial(GMlib::GMmaterial::emerald());
//    this->scene()->insert(Blending);



//auto bspline = new Bspline<float>(c);
//bspline->toggleDefaultVisualizer();
//bspline->sample(50, 0);
//this->scene()->insert(bspline);

//auto fem = new FEMObject();
}

void Scenario::cleanupScenario() {}

void Scenario::callDefferedGL() {

  GMlib::Array<const GMlib::SceneObject *> e_obj;
  this->scene()->getEditedObjects(e_obj);

  for (int i = 0; i < e_obj.getSize(); i++){
    if (e_obj(i)->isVisible()){
      e_obj[i]->replot();}
  }

//  fem->replot();
}
