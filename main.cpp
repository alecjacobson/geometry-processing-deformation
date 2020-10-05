#include "arap_precompute.h"
#include "arap_single_iteration.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/project.h>
#include <igl/unproject.h>
#include <igl/snap_points.h>
#include <igl/unproject_onto_mesh.h>
#include <Eigen/Core>
#include <iostream>
#include <string>
#include <stack>
#include <chrono>


using namespace Eigen;
using namespace std;


std::ofstream outputFile;


Vector3d count_time;



using namespace Eigen;

// Undoable
struct State
{
  // Rest and transformed control points
  Eigen::MatrixXd CV, CU;
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


  Eigen::MatrixXd R_last;
  Eigen::MatrixXd V,U;
  Eigen::MatrixXi F;
  long sel = -1;
  Eigen::RowVector3f last_mouse;
  igl::min_quad_with_fixed_data<double> arap_data;
  Eigen::SparseMatrix<double> arap_K;

  // test IO
  outputFile.open("./axis_angle_decimated_knight.txt", std::ios_base::app);


  // Load input meshes
  igl::read_triangle_mesh(
    (argc>1?argv[1]:"../data/decimated-knight.off"),V,F);

  // init R_last
  R_last.resize(3*V.rows(), 3);
  for (int i = 0; i < V.rows(); i++) {
    R_last.block(3 * i, 0, 3, 3) = MatrixXd::Identity(3, 3);
  }

  int num_of_group = V.rows()/8+1;
  int num_of_ele = num_of_group*72;


  // initialize rotation matrices
  void *R = NULL;
  posix_memalign(&R, 32, num_of_ele * sizeof(float));
  float *Rf = (float *)R;
  for (int i = 0; i < num_of_group; i++) {
    for (int j = 0; j < 8; j++) {
        Rf[i*72+0*8+j] = 1.0; // R(0,0)
        Rf[i*72+1*8+j] = 0.0; // R(0,1)
        Rf[i*72+2*8+j] = 0.0; // R(0,2)
        Rf[i*72+3*8+j] = 0.0; // R(1,0)
        Rf[i*72+4*8+j] = 1.0; // R(1,1)
        Rf[i*72+5*8+j] = 0.0; // R(1,2)
        Rf[i*72+6*8+j] = 0.0; // R(2,0)
        Rf[i*72+7*8+j] = 0.0; // R(2,1)
        Rf[i*72+8*8+j] = 1.0; // R(2,2)
    }
  }


  void *M = NULL;
  posix_memalign(&M, 32, num_of_ele * sizeof(float));
  float *Mf = (float *)M;


  U = V;
  igl::opengl::glfw::Viewer viewer;
  std::cout<<R"(
[click]  To place new control point
[drag]   To move control point
[space]  Toggle whether placing control points or deforming
U,u      Update deformation (i.e., run another iteration of solver)
R,r      Reset control points 
⌘ Z      Undo
⌘ ⇧ Z    Redo
)";
  enum Method
  {
    ARAP = 0,
  } method = ARAP;

  const auto & update = [&]()
  {
    // predefined colors
    const Eigen::RowVector3d orange(1.0,0.7,0.2);
    const Eigen::RowVector3d yellow(1.0,0.9,0.2);
    const Eigen::RowVector3d blue(0.2,0.3,0.8);
    const Eigen::RowVector3d green(0.2,0.6,0.3);
    if(s.placing_handles)
    {
      viewer.data().set_vertices(V);
      viewer.data().set_colors(blue);
      viewer.data().set_points(s.CV,orange);
    }else
    {
      // SOLVE FOR DEFORMATION
      switch(method)
      {
        default:
        case ARAP:
        {
          count_time = arap_single_iteration(arap_data,arap_K,s.CU,outputFile,R_last,U,Rf,Mf,num_of_group);
          break;
        }
      }
      viewer.data().set_vertices(U);
      viewer.data().set_colors(blue);
      viewer.data().set_points(s.CU,orange);
    }
    viewer.data().compute_normals();
  };
  viewer.callback_mouse_down = 
    [&](igl::opengl::glfw::Viewer&, int, int)->bool
  {
    last_mouse = Eigen::RowVector3f(
      viewer.current_mouse_x,viewer.core().viewport(3)-viewer.current_mouse_y,0);
    if(s.placing_handles)
    {
      // Find closest point on mesh to mouse position
      int fid;
      Eigen::Vector3f bary;
      if(igl::unproject_onto_mesh(
        last_mouse.head(2),
        viewer.core().view,
        viewer.core().proj, 
        viewer.core().viewport, 
        V, F, 
        fid, bary))
      {
        long c;
        bary.maxCoeff(&c);
        Eigen::RowVector3d new_c = V.row(F(fid,c));
        if(s.CV.size()==0 || (s.CV.rowwise()-new_c).rowwise().norm().minCoeff() > 0)
        {
          push_undo();
          s.CV.conservativeResize(s.CV.rows()+1,3);
          // Snap to closest vertex on hit face
          s.CV.row(s.CV.rows()-1) = new_c;
          update();
          return true;
        }
      }
    }else
    {
      // Move closest control point
      Eigen::MatrixXf CP;
      igl::project(
        Eigen::MatrixXf(s.CU.cast<float>()),
        viewer.core().view,
        viewer.core().proj, viewer.core().viewport, CP);
      Eigen::VectorXf D = (CP.rowwise()-last_mouse).rowwise().norm();
      sel = (D.minCoeff(&sel) < 30)?sel:-1;
      if(sel != -1)
      {
        last_mouse(2) = CP(sel,2);
        push_undo();
        update();
        return true;
      }
    }
    return false;
  };

  viewer.callback_mouse_move = [&](igl::opengl::glfw::Viewer &, int,int)->bool
  {
    if(sel!=-1)
    {
      Eigen::RowVector3f drag_mouse(
        viewer.current_mouse_x,
        viewer.core().viewport(3) - viewer.current_mouse_y,
        last_mouse(2));
      Eigen::RowVector3f drag_scene,last_scene;
      igl::unproject(
        drag_mouse,
        viewer.core().view,
        viewer.core().proj,
        viewer.core().viewport,
        drag_scene);
      igl::unproject(
        last_mouse,
        viewer.core().view,
        viewer.core().proj,
        viewer.core().viewport,
        last_scene);
      s.CU.row(sel) += (drag_scene-last_scene).cast<double>();
      last_mouse = drag_mouse;
      update();
      return true;
    }
    return false;
  };
  viewer.callback_mouse_up = [&](igl::opengl::glfw::Viewer&, int, int)->bool
  {
    sel = -1;
    return false;
  };
  viewer.callback_key_pressed = 
    [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)
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
      case 'U':
      case 'u':
      {
        // Just trigger an update
        break;
      }
      case ' ':
        push_undo();
        s.placing_handles ^= 1;
        if(!s.placing_handles && s.CV.rows()>0)
        {
          // Switching to deformation mode
          s.CU = s.CV;
          Eigen::VectorXi b;
          igl::snap_points(s.CV,V,b);
          // PRECOMPUTATION FOR DEFORMATION
          arap_precompute(V,F,b,arap_data,arap_K);
        }
        break;
      default:
        return false;
    }
    update();
    return true;
  };

  // Special callback for handling undo
  viewer.callback_key_down = 
    [&](igl::opengl::glfw::Viewer &, unsigned char key, int mod)->bool
  {
    if(key == 'Z' && (mod & GLFW_MOD_SUPER))
    {
      (mod & GLFW_MOD_SHIFT) ? redo() : undo();
      update();
      return true;
    }
    return false;
  };
  viewer.callback_pre_draw = 
    [&](igl::opengl::glfw::Viewer &)->bool
  {
    if(viewer.core().is_animating && !s.placing_handles && method == ARAP)
    {
      count_time = arap_single_iteration(arap_data,arap_K,s.CU,outputFile,R_last,U,Rf,Mf,num_of_group);
      update();
    }
    return false;
  };


  viewer.core().background_color = Eigen::Vector4f(1,1,1,0);

  viewer.data().set_mesh(V,F);
  viewer.data().show_lines = false;
  viewer.core().is_animating = true;
  viewer.data().face_based = true;
  update();
  viewer.launch();



  outputFile.close();

  return EXIT_SUCCESS;


}
