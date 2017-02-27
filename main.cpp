#include "biharmonic_precompute.h"
#include "biharmonic_solve.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/read_triangle_mesh.h>
#include <igl/viewer/Viewer.h>
#include <igl/project.h>
#include <igl/unproject.h>
#include <igl/snap_points.h>
#include <igl/slice.h>
#include <igl/harmonic.h>
#include <igl/unproject_onto_mesh.h>
#include <Eigen/Core>
#include <string>
#include <iostream>
#include <stack>

struct State
{
  Eigen::MatrixXd CV;
  Eigen::MatrixXd CU;
  bool placing_handles = true;
} s;

int main(int argc, char *argv[])
{
  // Undo Management
  std::stack<State> undo_stack,redo_stack;
  const auto push_undo = [&](State & _s=s)
  {
    undo_stack.push(_s);
    // clear
    redo_stack = std::stack<State>();
  };
  const auto undo = [&]()
  {
    if(!undo_stack.empty())
    {
      redo_stack.push(s);
      s = undo_stack.top();
      undo_stack.pop();
    }
  };
  const auto redo = [&]()
  {
    if(!redo_stack.empty())
    {
      undo_stack.push(s);
      s = redo_stack.top();
      redo_stack.pop();
    }
  };

  Eigen::MatrixXd V,U;
  Eigen::MatrixXi F;
  long sel = -1;
  Eigen::RowVector3f last;
  igl::min_quad_with_fixed_data<double> biharmonic_data;

  // Load input meshes
  igl::read_triangle_mesh(
    (argc>1?argv[1]:"../shared/data/bunny.off"),V,F);
  U = V;
  // Load data into MatrixXd rather than VectorXd for simpler `smooth` API
  // Just use y-coordinates as data to be smoothed
  // Create a libigl Viewer object to toggle between point cloud and mesh
  igl::viewer::Viewer viewer;
  std::cout<<R"(
[space]  Toggle whether placing handles or deforming
)";

  bool plot_parameterization = false;
  const auto & update = [&]()
  {
    using namespace igl;
    if(s.placing_handles)
    {
      viewer.data.set_vertices(V);
      viewer.data.set_colors(
        Eigen::RowVector3d(GOLD_DIFFUSE[0], GOLD_DIFFUSE[1], GOLD_DIFFUSE[2]));
      viewer.data.set_points(s.CV,Eigen::RowVector3d(1.0,0.5,0.1));
    }else
    {

      Eigen::MatrixXd D;
      biharmonic_solve(biharmonic_data,s.CU-s.CV,D);
      U = V+D;
      viewer.data.set_vertices(U);
      viewer.data.set_colors(
        Eigen::RowVector3d(SILVER_DIFFUSE[0], SILVER_DIFFUSE[1], SILVER_DIFFUSE[2]));
      viewer.data.set_points(s.CU,Eigen::RowVector3d(1.0,0.5,0.1));
    }
    viewer.data.compute_normals();
  };
  viewer.callback_mouse_down = 
    [&](igl::viewer::Viewer&, int, int)->bool
  {
    last = Eigen::RowVector3f(
      viewer.current_mouse_x,
      viewer.core.viewport(3) - viewer.current_mouse_y,
      0);
    if(s.placing_handles)
    {
      int fid;
      Eigen::Vector3f bc;
      // Cast a ray in the view direction starting from the mouse position
      if(igl::unproject_onto_mesh(
        last.head(2),
        viewer.core.view * viewer.core.model,
        viewer.core.proj, 
        viewer.core.viewport, 
        V, F, 
        fid, bc))
      {
        push_undo();
        s.CV.conservativeResize(s.CV.rows()+1,3);
        s.CV.row(s.CV.rows()-1) = 
          V.row(F(fid,0))*bc(0)+ V.row(F(fid,1))*bc(1)+ V.row(F(fid,2))*bc(2);
        update();
        return true;
      }
    }else
    {
      // Get closest control point
      Eigen::MatrixXf CP;
      igl::project(
        Eigen::MatrixXf(s.CU.cast<float>()),
        viewer.core.view * viewer.core.model, 
        viewer.core.proj, viewer.core.viewport, CP);
      Eigen::VectorXf D = (CP.rowwise()-last).rowwise().norm();
      sel = (D.minCoeff(&sel) < 30)?sel:-1;
      if(sel != -1)
      {
        last(2) = CP(sel,2);
        push_undo();
        update();
        return true;
      }
    }
    return false;
  };

  viewer.callback_mouse_move = [&](igl::viewer::Viewer &, int,int)->bool
  {
    if(sel!=-1)
    {
      Eigen::RowVector3f drag(
        viewer.current_mouse_x,
        viewer.core.viewport(3) - viewer.current_mouse_y,
        last(2));
      Eigen::RowVector3f drag_scene,last_scene;
      igl::unproject(drag, viewer.core.view * viewer.core.model, viewer.core.proj,viewer.core.viewport,drag_scene);
      igl::unproject(last, viewer.core.view * viewer.core.model, viewer.core.proj,viewer.core.viewport,last_scene);
      s.CU.row(sel) += (drag_scene-last_scene).cast<double>();
      last = drag;

      update();
      return true;
    }
    return false;
  };
  viewer.callback_mouse_up = [&](igl::viewer::Viewer&, int, int)->bool
  {
    sel = -1;
    return false;
  };
  viewer.callback_key_pressed = 
    [&](igl::viewer::Viewer &, unsigned int key, int mod)
  {
    switch(key)
    {
      case 'R':
      case 'r':
      {
        push_undo();
        s.CU = s.CV;
        break;
      }
      case ' ':
        push_undo();
        s.placing_handles ^= 1;
        if(!s.placing_handles)
        {
          // Switching to deformation mode
          s.CU = s.CV;
          // PRECOMPUTATION
          Eigen::VectorXi b;
          igl::snap_points(s.CV,V,b);
          biharmonic_precompute(V,F,b,biharmonic_data);
        }
        break;
      default:
        return false;
    }
    update();
    return true;
  };

  viewer.callback_key_down = 
    [&](igl::viewer::Viewer &, unsigned char key, int mod)->bool
  {
    switch(key)
    {
      case 'Z':
        if(mod & GLFW_MOD_SUPER)
        {
          if(mod & GLFW_MOD_SHIFT)
          {
            redo();
          }else
          {
            undo();
          }
        }
        update();
        break;
    }
    return false;
  };
  viewer.data.set_mesh(V,F);
  viewer.core.show_lines = false;
  viewer.core.is_animating = true;
  update();
  viewer.launch();

  return EXIT_SUCCESS;
}
