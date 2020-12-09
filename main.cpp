#include "cubic_style_precomputation.h"
#include "cubic_style_single_iteration.h"
#include <cubic_style_data.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/snap_points.h>
#include <Eigen/Core>
#include <iostream>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

// The code is ref to main file of deformation in class
// https://github.com/alecjacobson/geometry-processing-deformation 
// Undoable
struct State
{
    // Rest and transformed control points
    Eigen::MatrixXd CV;
    bool placing_handles = true;
} s;

int main(int argc, char *argv[])
{
  Eigen::MatrixXd V,U;
  Eigen::MatrixXi F;

  // Load input meshes
  igl::readOBJ((argc>1?argv[1]:"../data/spot.obj"),V,F);
  U = V;
  Eigen::Vector3d m = V.colwise().minCoeff();
  // Eigen::Vector3d M = V.colwise().maxCoeff();

  igl::opengl::glfw::Viewer viewer;
  std::cout<<R"(
  [space]  Toggle whether placing control points or deforming
  >,.      increase lambda by 0.02 
  <,,      decrease lambda by 0.02
  )";
  
  cubic_style_data data;
  data.lambda = 2e-1;

  // https://libigl.github.io/tutorial/
  // Attach a menu plugin
  igl::opengl::glfw::imgui::ImGuiMenu menu;
  viewer.plugins.push_back(&menu);

  // Customize the menu
  menu.callback_draw_viewer_window = [](){};
  menu.callback_draw_custom_window = [&]()
  {
    // Define next window position + size
    ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(200, 160), ImGuiCond_FirstUseEver);
    ImGui::Begin(
        "Output", nullptr,
        ImGuiWindowFlags_NoSavedSettings
    );

    // Expose the same variable directly ...
    static std::string str = argc>1?argv[1]:"../data/spot.obj";
    ImGui::InputText("Name", str);

    ImGui::PushItemWidth(-80);
    ImGui::DragScalar("Lambda", ImGuiDataType_Double, &data.lambda, 0.1, 0, 0, "%.4f");
    ImGui::PopItemWidth();
    ImGui::End();
  };

  const auto & update = [&]()
  {
    // predefined colors
    const Eigen::RowVector3d yellow(1.0,0.9,0.2);
    const Eigen::RowVector3d blue(0.2,0.3,0.8);
    const Eigen::RowVector3d grey(0.8,0.8,0.8);
    if(s.placing_handles)
    {
      viewer.data().set_vertices(V);
      viewer.data().set_colors(grey);
    }else
    {
      // SOLVE FOR DEFORMATION
      cubic_style_single_iteration(data, U);
      viewer.data().set_vertices(U);
      viewer.data().set_colors(yellow);
    }
    viewer.data().compute_normals();
  };

  viewer.callback_key_pressed = 
    [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)
  {
    switch(key)
    {
      case '>':
      case '.':
      {
        // increase lambda
        data.lambda += 0.02;
        break;
      }
      case '<':
      case ',':
      {
        // decrease lambda
        data.lambda -= 0.02;
        break;
      }
      case ' ':
        // push_undo();
        s.placing_handles ^= 1;
        // set constrained F(0,0)
        s.CV = Eigen::MatrixXd();
        s.CV.resize(1,3);
        s.CV.row(0) = V.row(F(0,0));
        if(!s.placing_handles)
        {
          // Switching to deformation mode
          igl::snap_points(s.CV,V,data.b);
          data.bc.setZero(data.b.size(),3);
          for (int i=0; i<data.b.size(); i++) {
            data.bc.row(i) = V.row(data.b(i));
          }
          // PRECOMPUTATION FOR DEFORMATION
          U = V;
          cubic_style_precomputation(V,F,data);
        }
        break;
      default:
        return false;
    }
    update();
    return true;
  };
  
  viewer.callback_pre_draw = 
    [&](igl::opengl::glfw::Viewer &)->bool
  {
    if(viewer.core().is_animating && !s.placing_handles)
    {
      update();
    }
    return false;
  };
  viewer.data().set_mesh(V,F);
  viewer.data().show_lines = false;
  viewer.core().is_animating = true;
  viewer.data().face_based = true;
  update();
  viewer.launch();
  return EXIT_SUCCESS;
}
