/* ========================================================================= *
 *                                                                           *
 *                               OpenMesh                                    *
 *           Copyright (c) 2001-2015, RWTH-Aachen University                 *
 *           Department of Computer Graphics and Multimedia                  *
 *                          All rights reserved.                             *
 *                            www.openmesh.org                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * This file is part of OpenMesh.                                            *
 *---------------------------------------------------------------------------*
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  *
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 *                                                                           *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 *                                                                           *
 * 3. Neither the name of the copyright holder nor the names of its          *
 *    contributors may be used to endorse or promote products derived from   *
 *    this software without specific prior written permission.               *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       *
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED *
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           *
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER *
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 *                                                                           *
 * ========================================================================= */

//=============================================================================
//
//  CLASS ModQuadricT
//
//=============================================================================

#ifndef OSG_MODQUADRIC_HH
#define OSG_MODQUADRIC_HH


//== INCLUDES =================================================================

#include <float.h>
#include <OpenMesh/Tools/Decimater/ModBaseT.hh>
#include <OpenMesh/Core/Utils/Property.hh>
#include <OpenMesh/Core/Utils/vector_cast.hh>
#include <OpenMesh/Core/Geometry/QuadricT.hh>


//== NAMESPACE ================================================================

namespace OpenMesh  {
namespace Decimater {


//== CLASS DEFINITION =========================================================


/** \brief Mesh decimation module computing collapse priority based on error quadrics.
 *
 *  This module can be used as a binary and non-binary module.
 */
template <class MeshT>
class ModQuadricT : public ModBaseT<MeshT>
{
public:

  // Defines the types Self, Handle, Base, Mesh, and CollapseInfo
  // and the memberfunction name()
  DECIMATING_MODULE( ModQuadricT, MeshT, Quadric );

public:

  /** Constructor
   *  \internal
   */
  ModQuadricT( MeshT &_mesh )
    : Base(_mesh, false)
  {
    unset_max_err();
    Base::mesh().add_property( quadrics_ );

    //Kevin Coopman Added these
    SetNewVertex(0);
    quadric_max_error_ = 0;
    quadric_min_error_ = 10000;
    quadric_sum_error_ = 0;
    quadric_count      = 0;
  }


  /// Destructor
  virtual ~ModQuadricT()
  {
    Base::mesh().remove_property(quadrics_);
  }


public: // inherited

  /// Initalize the module and prepare the mesh for decimation.
  virtual void initialize(void);

  // -------------------------------------------------------------------
  double det(double a11, double a12, double a13,
	     double a21, double a22, double a23,
	     double a31, double a32, double a33)
  {
    double det =  a11*a22*a33 + a13*a21*a32 + a12*a23*a31
      - a13*a22*a31 - a11*a23*a32- a12*a21*a33;
    return det;
  }

  // ---------------------------------------------------------------------
  
  /** Compute collapse priority based on error quadrics.
   *
   *  \see ModBaseT::collapse_priority() for return values
   *  \see set_max_err()
   */
  virtual float collapse_priority(const CollapseInfo& _ci)
  {
    using namespace OpenMesh;

    typedef Geometry::QuadricT<double> Q;

    //OK we add Q1 and Q2 together
    Q q = Base::mesh().property(quadrics_, _ci.v0);
    q += Base::mesh().property(quadrics_, _ci.v1);

    // evaluate quadric Q at (3D or 4D) vector v: v*Q*v
    double err = q(_ci.p1);

    if(create_new_vertex_)
      {
	// calculate 
	double a0 = q.a();
	double b1 = q.b();
	double c2 = q.c();
	double d3 = q.d();
	double e4 = q.e();
	double f5 = q.f();
	double g6 = q.g();
	double h7 = q.h();
	double i8 = q.i();
	
	double determinate = det(a0, b1, c2, b1, e4, f5, c2, f5, h7);
	if ( determinate != 0)
	  {
	    double x = -1/determinate*(det(b1, c2, d3, e4, f5, g6, f5, h7, i8));
	    double y =  1/determinate*(det(a0, c2, d3, b1, f5, g6, c2, h7, i8));
	    double z = -1/determinate*(det(a0, b1, d3, b1, e4, g6, c2, f5, i8));
	    
	    typename Mesh::Point new_point;
	    new_point[0] = x;
	    new_point[1] = y;
	    new_point[2] = z;
	    err = q(new_point);	    
	  }
      }
    
    return float( (err < max_err_) ? err : float( Base::ILLEGAL_COLLAPSE ) );
  }


  // Kevin added this
 /// Pre-process halfedge collapse (accumulate quadrics)
  virtual void preprocess_collapse(const CollapseInfo& _ci) override
  {
    using namespace OpenMesh;
    typedef Geometry::QuadricT<double> Q;
    double error = 0;
    
    //OK we add Q1 and Q2 together, the edge we are going to Collapse 
    Q q = Base::mesh().property(quadrics_, _ci.v0);
    q  += Base::mesh().property(quadrics_, _ci.v1);
    error = q(_ci.p1);
	
    if(create_new_vertex_)
      {    	
	// calculate 
	double a0 = q.a();
	double b1 = q.b();
	double c2 = q.c();
	double d3 = q.d();
	double e4 = q.e();
	double f5 = q.f();
	double g6 = q.g();
	double h7 = q.h();
	double i8 = q.i();
	
	double determinate = det(a0, b1, c2, b1, e4, f5, c2, f5, h7);
	if ( determinate != 0)
	  {
	    double x = -1/determinate*(det(b1, c2, d3, e4, f5, g6, f5, h7 , i8));	// vx = A41/det(q_delta)
	    double y =  1/determinate*(det(a0, c2, d3, b1, f5, g6, c2, h7 , i8));	// vy = A42/det(q_delta)
	    double z = -1/determinate*(det(a0, b1, d3, b1, e4, g6, c2, f5,  i8));	// vz = A43/det(q_delta)
	    
	    typename Mesh::Point new_point;
	    new_point[0] = x;
	    new_point[1] = y;
	    new_point[2] = z;
	    
	    //lets get cost/error so we can keep track of max and min error
	    error  = q(new_point);

	    // get the vertex to remain
	    typename Mesh::VertexHandle v1(Base::mesh().to_vertex_handle(_ci.v0v1));
	    // set the point in the mesh
	    Base::mesh().set_point(v1,new_point);	    
	  }    
      }

    //let keep some stats on error here
    if(error > quadric_max_error_)
      {
	quadric_max_error_ = error;
	//std::cout << "max error occured = " << error << std::endl;
      }
    if(error < quadric_min_error_) quadric_min_error_ = error;
    
    quadric_sum_error_ += error;
    quadric_count++;
        
  }
  
  /// Post-process halfedge collapse (accumulate quadrics)
  virtual void postprocess_collapse(const CollapseInfo& _ci)
  {
    Base::mesh().property(quadrics_, _ci.v1) +=
      Base::mesh().property(quadrics_, _ci.v0);
  } 

  /// set the percentage of maximum quadric error
  void set_error_tolerance_factor(double _factor);



public: // specific methods

  /** Set maximum quadric error constraint and enable binary mode.
   *  \param _err    Maximum error allowed
   *  \param _binary Let the module work in non-binary mode in spite of the
   *                 enabled constraint.
   *  \see unset_max_err()
   */
  void set_max_err(double _err, bool _binary=true)
  {
    max_err_ = _err;
    Base::set_binary(_binary);
  }

  /// Unset maximum quadric error constraint and restore non-binary mode.
  /// \see set_max_err()
  void unset_max_err(void)
  {
    max_err_ = DBL_MAX;
    Base::set_binary(false);
  }

  /// Return value of max. allowed error.
  double max_err() const { return max_err_; }

  double quadric_avg_error() { return quadric_sum_error_/(double)quadric_count; }
  double quadric_max_error() { return quadric_max_error_; }
  double quadric_min_error() { return quadric_min_error_; }

  // turn off and on new vertex method
  void SetNewVertex( bool true_false)
  {
    create_new_vertex_ = true_false;
    //if(true_false)
    //  std::cout << "New Vertex is ON \n";
    //else
    //  std::cout << "New Vertex is OFF \n";
  }
  
private:

  // maximum quadric error
  double max_err_;
  
  // Kevin added these
  bool   create_new_vertex_;
  double quadric_max_error_;
  double quadric_min_error_;
  double quadric_sum_error_;
  int    quadric_count;
  
  // this vertex property stores a quadric for each vertex
  VPropHandleT< Geometry::QuadricT<double> >  quadrics_;
};

//=============================================================================
} // END_NS_DECIMATER
} // END_NS_OPENMESH
//=============================================================================
#if defined(OM_INCLUDE_TEMPLATES) && !defined(OPENMESH_DECIMATER_MODQUADRIC_CC)
#define OSG_MODQUADRIC_TEMPLATES
#include "ModQuadricT_impl.hh"
#endif
//=============================================================================
#endif // OSG_MODQUADRIC_HH defined
//=============================================================================

