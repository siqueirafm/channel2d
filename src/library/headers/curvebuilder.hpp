#pragma once

/** 
 * \file curvebuilder.hpp
 *
 * \brief  Definition of  a class  for threading  a b-spline  curve of
 * degree 3  through a planar channel  defined by a pair  of polygonal
 * chains.
 *
 * \author
 * Marcelo Ferreira Siqueira \n
 * Universidade Federal do Rio Grande do Norte, \n
 * Departamento de Matem&aacute;tica, \n
 * mfsiqueira at mat (dot) ufrn (dot) br
 *
 * \version 1.0
 * \date May 2016
 *
 * \attention This program is distributed WITHOUT ANY WARRANTY, and it
 *            may be freely redistributed under the condition that the
 *            copyright notices  are not removed,  and no compensation
 *            is received. Private, research, and institutional use is
 *            free. Distribution of this  code as part of a commercial
 *            system  is permissible ONLY  BY DIRECT  ARRANGEMENT WITH
 *            THE AUTHOR.
 */

#include "exceptionobject.hpp"           // ExceptionObject
#include "tabulatedfunction.hpp"         // TabulatedFunction
#include "coefficient.hpp"               // Coefficient
#include "bound.hpp"                     // Bound

extern "C" {
#include "glpk.h"                        // For all GLPK functions
}

#include <vector>                        // std::vector
#include <string>                        // std::string
#include <sstream>                       // std::stringstream
#include <cassert>                       // assert
#include <cstdlib>                       // size_t


/** 
 * \defgroup ChannelNameSpace Namespace channel.
 * @{
 */

/**
 * \namespace channel
 *
 * \brief   The  namespace   channel  contains   the  definition   and
 * implementation of a  set of classes for threading  a cubic b-spline
 * curve  into  a given  planar  channel  delimited by  two  polygonal
 * chains.
 *
 */

namespace channel {


  /**
   * \class CurveBuilder
   *
   * \brief This class provides methods for threading a cubic b-spline
   * curve through a  planar channel delimited by a  pair of polygonal
   * chains.
   *
   * \attention This class is based on a particular case (i.e., planar
   *             and cubic curves) of the  method described by Myles &
   *             Peters in
   *
   *             A. Myles and J. Peters,
   *             Threading splines through 3d channels
   *             Computer-Aided Design, 37(2), 139-148, 2005.
   *
   */
  
  class CurveBuilder {
  private:
    // ---------------------------------------------------------------
    //
    // Private data members
    //
    // ---------------------------------------------------------------
  
    size_t _np ; ///< The number of b-spline segments per channel segment.
    size_t _nc ; ///< The number of segments of the channel.
    bool _closed ; ///< A flag to indicate whether the channel is closed.
    std::vector< double > _lxcoords ; ///< X coordinates of the lower polygonal chain of the channel.
    std::vector< double > _lycoords ; ///< Y coordinates of the lower polygonal chain of the channel.  
    std::vector< double > _uxcoords ; ///< X coordinates of the upper polygonal chain of the channel.
    std::vector< double > _uycoords ; ///< Y coordinates of the upper polygonal chain of the channel.
    TabulatedFunction* _tf ;  ///< A pointer to the lower and upper \f$a\f$ functions.
    std::vector< std::vector< Coefficient > > _coefficients ;  ///< Coefficients of the constraints of the linear program.
    std::vector< Bound > _bounds ;  ///< Type of the constraints and their bounds.
    std::vector< double > _ctrlpts ; ///< X and Y coordinates of the control points of the resulting b-spline.
    std::vector< double > _secdiff ; ///< Lower and upper bounds on the second difference values.
    double _ofvalue ;  ///< Optimal value (i.e., minimum) of the objective function.
    
  public: 
    // ---------------------------------------------------------------
    //
    // Public methods
    //
    // ---------------------------------------------------------------

    
    /**
     * \fn CurveBuilder( size_t np , size_t nc , bool closed , double* lx , double* ly , double* ux , double* uy ) 
     *
     * \brief Creates an instance of this class.
     *
     * \param np The number of b-spline segments.
     * \param nc The number of c-segments of the channel.
     * \param closed A flag to indicate whether the channel is closed. 
     * \param lx A  pointer to an array with the  x-coordinates of the
     * lower envelope of the channel.
     * \param ly A  pointer to an array with the  y-coordinates of the
     * lower envelope of the channel.
     * \param ux A  pointer to an array with the  x-coordinates of the
     * upper envelope of the channel.
     * \param uy A  pointer to an array with the  y-coordinates of the
     * upper envelope of the channel.
     *
     */
    CurveBuilder(
                 size_t np ,
                 size_t nc ,
                 bool closed ,
                 double* lx ,
                 double* ly ,
                 double* ux ,
                 double* uy
                ) ;
    
    
    /**
     * \fn CurveBuilder( const CurveBuilder& b )
     *
     * \brief Clones an instance of this class.
     *
     * \param b A reference to another instance of this class.
     *
     */
    CurveBuilder( const CurveBuilder& b ) ;


    /**
     * \fn bool build( int& error )
     *
     * \brief Solves the channel problem by solving a linear program.
     *
     * \param error Code returned by the LP solver whenever a solution
     * could not be  found. If a solution is found,  this parameter is
     * ignored.
     *
     * \return The logic  value true if the LP solver  is able to find
     * an  optimal solution  for the  channel problem;  otherwise, the
     * logic value false is returned.
     *
     */
    bool build( int& error ) ;


    /**
     * \fn size_t get_degree() const
     *
     * \brief Returns the degree of the bspline curve.
     *
     * \return The degree of the bspline curve.
     *
     */
    size_t get_degree() const
    {
      return 3 ;
    }


    /**
     * \fn size_t get_number_of_segments() const
     *
     * \brief Returns the number of b-spline segments.
     *
     * \return The number of b-spline segments.
     *
     */
    size_t get_number_of_segments() const
    {
      return _np ;
    }
    
    
    /**
     * \fn size_t get_number_of_csegments() const
     *
     * \brief Returns the number of c-segments of the channel.
     *
     * \return The number of c-segments of the channel.
     *
     */
    size_t get_number_of_csegments() const
    {
      return _nc ;
    }
    
    
    /**
     * \fn bool is_curve_closed() const
     *
     * \brief Returns  the logic value  true if the b-spline  curve is
     * closed, and the logic value false otherwise.
     *
     * \return The logic  value true if the b-spline  curve is closed,
     * and the logic value false otherwise.
     *
     */
    bool is_curve_closed() const
    {
      return _closed ;
    }


    /**
     * \fn size_t get_number_of_control_points() const
     *
     * \brief Returns the number of control points of the b-spline.
     *
     * \return The number of control points of the b-spline.
     *
     */
    size_t get_number_of_control_points() const
    {
      return _np + get_degree() ;
    }


    /**
     * \fn size_t get_number_of_constraints() const
     *
     * \brief Returns the number of constraints of the instance of the
     * linear program  corresponding to the channel  problem solved by
     * this class.
     *
     * \return The number of constraints of the instance of the linear
     * program  corresponding to  the channel  problem solved  by this
     * class.
     *
     */
    size_t get_number_of_constraints() const
    {
      if ( _coefficients.empty() ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "No constraint has been created so far" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }
	    
      return _coefficients.size() ;
    }

 
    /**
     * \fn double get_control_value( const size_t i , const size_t v ) const 
     *
     * \brief  Returns the  \f$v\f$-th  coordinate  of the  \f$i\f$-th
     * control point of the b-spline curve threaded into the channel.
     *
     * \param  i The  index of  the  \f$i\f$-th control  point of  the
     * b-spline curve.

     * \param v The \f$v\f$-th  Cartesian coordinate of the \f$i\f$-th
     * control point of the b-spline curve.
     *
     * \return  The \f$v\f$-th  coordinate of  the \f$i\f$-th  control
     * point of the b-spline curve threaded into the channel.
     *
     */
    double get_control_value(
                             const size_t i ,
                             const size_t v
                            )
      const 
    {        
      if ( i >= ( _np + 3 ) ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "Index of the control point is out of range" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }
        
      if ( v >= 2 ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "Index of the Cartesian coordinate is out of range" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }
        
      if ( _ctrlpts.empty() ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "Control points have not been computed" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }

      size_t index = ( 2 * i ) + v ;
      
      return _ctrlpts[ index ] ;
    }


    /**
     * \fn size_t get_number_of_coefficients_in_the_ith_constraint( const size_t i ) const
     *
     * \brief  Returns the  number of  coefficients of  the \f$i\f$-th
     * constraint of the instance of the linear program.
     *
     * \param i The index of a constraint. 
     *
     * \return The number of coefficients of the \f$i\f$-th constraint
     * of the instance of the linear program.
     *
     */
    size_t get_number_of_coefficients_in_the_ith_constraint( const size_t i ) const
    {
      if ( _coefficients.empty() ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "No constraint has been created so far" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }

      if ( i >= _coefficients.size() ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "Constraint index is out of range" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }

      return _coefficients[ i ].size() ;
    }
    
      
    /**
     * \fn size_t get_coefficient_identifier( const size_t i , const size_t j ) const
     *
     * \brief Returns the index of  the column that corresponds to the
     * \f$j\f$-th  coefficient of  the  \f$i\f$-th  constraint in  the
     * matrix associated with the linear program (LP) instance.
     *
     * \param i The index of a constraint.
     * \param j The \f$j\f$-th (nonzero) coefficient of the \f$i\f$-th
     * constraint.
     *
     * \return  The  index  of  the column  that  corresponds  to  the
     * \f$j\f$-th  coefficient of  the  \f$i\f$-th  constraint in  the
     * matrix associated with the linear program (LP) instance.
     *
     */
    size_t get_coefficient_identifier( const size_t i , const size_t j ) const
    {
      if ( _coefficients.empty() ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "No constraint has been created so far" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }

      if ( i >= _coefficients.size() ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "Constraint index is out of range" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }

      if ( j >= _coefficients[ i ].size() ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "Coefficient index is out of range" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }
	   
      return _coefficients[ i ][ j ].get_col() ;
    }


    /**
     * \fn size_t get_coefficient_value( const size_t i , const size_t j ) const
     *
     * \brief  Returns  the \f$(  i  ,  j  )\f$  entry of  the  matrix
     * associated with the instance of the linear program.
     *
     * \param i The index of a constraint.
     * \param j The \f$j\f$-th (nonzero) coefficient of the \f$i\f$-th
     * constraint.
     *
     * \return The \f$( i , j )\f$ entry of the matrix associated with
     * the instance of the linear program.
     *
     */
    double get_coefficient_value( const size_t i , const size_t j ) const
    {
      if ( _coefficients.empty() ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "No constraint has been created so far" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }

      if ( i >= _coefficients.size() ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "Constraint index is out of range" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }

      if ( j >= _coefficients[ i ].size() ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "Coefficient index is out of range" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }
	   
      return _coefficients[ i ][ j ].get_value() ;
    }


    /**
     * \fn double get_bound_of_ith_constraint( const size_t i ) const
     *
     * \brief Returns  the real  value on the  right-hand side  of the
     * equality   or  inequality   corresponding  to   the  \f$i\f$-th
     * constraint.
     *
     * \param i The index of a constraint. 
     *
     * \return The real  value on the right-hand side  of the equality
     * or inequality corresponding to the \f$i\f$-th constraint.
     *
     */
    double get_bound_of_ith_constraint( const size_t i ) const
    {
      if ( _coefficients.empty() ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "No constraint has been created so far" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }

      if ( i >= _coefficients.size() ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "Constraint index is out of range" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }

#ifdef DEBUGMODE
      assert( _bounds.size() == _coefficients.size() ) ;
      assert( _bounds.size() > std::vector< std::vector< Bound > >::size_type( i ) ) ;
#endif

      return _bounds[ i ].get_value() ;
    }


    /**
     * \fn bool is_equality( const size i ) const
     *
     * \brief Returns  the logic value  true if  the type of  the i-th
     * constraint  is equality;  otherwise,  returns  the logic  value
     * false.
     *
     * \param i The index of a constraint. 
     *
     * \return The logic value true if the type of the i-th constraint
     * is equality; otherwise, the logic value false is returned.
     *
     */
    bool is_equality( const size_t i ) const
    {
      if ( _coefficients.empty() ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "No constraint has been created so far" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }

      if ( i >= _coefficients.size() ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "Constraint index is out of range" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }

#ifdef DEBUGMODE
      assert( _bounds.size() == _coefficients.size() ) ;
      assert( _bounds.size() > i ) ;
#endif

      return _bounds[ i ].get_type() == Bound::EQT ;
    }

    /**
     * \fn bool is_greater_than_or_equal_to( const size_t i ) const
     *
     * \brief Returns the  logic value true if the  i-th constraint is
     * an inequality of the type  greater than or equal to; otherwise,
     * returns the logic value false.
     *
     * \param i The index of a constraint. 
     *
     * \return  The logic  value true  if  the i-th  constraint is  an
     * inequality of the type greater than or equal to; otherwise, the
     * logic value false is returned.
     *
     */
    bool is_greater_than_or_equal_to( const size_t i ) const
    {
      if ( _coefficients.empty() ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "No constraint has been created so far" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }

      if ( i >= _coefficients.size() ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "Constraint index is out of range" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }

#ifdef DEBUGMODE
      assert( _bounds.size() == _coefficients.size() ) ;
      assert( _bounds.size() > i ) ;
#endif

      return _bounds[ i ].get_type() == Bound::GTE ;
    }


    /**
     * \fn bool is_less_than_or_equal_to( const size_t i ) const
     *
     * \brief Returns the  logic value true if the  i-th constraint is
     * an inequality  of the  type less than  or equal  to; otherwise,
     * returns the logic value false.
     *
     * \param i The index of a constraint. 
     *
     * \return  The logic  value true  if  the i-th  constraint is  an
     * inequality of  the type less  than or equal to;  otherwise, the
     * logic value false is returned.
     *
     */
    bool is_less_than_or_equal_to( const size_t i ) const
    {
      if ( _coefficients.empty() ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "No constraint has been created so far" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }

      if ( i >= _coefficients.size() ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "Constraint index is out of range" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }

#ifdef DEBUGMODE
      assert( _bounds.size() == _coefficients.size() ) ;
      assert( _bounds.size() > i ) ;
#endif

      return _bounds[ i ].get_type() == Bound::LTE ;
    }

    
    /**
     * \fn double get_lower_bound_on_second_difference_value( const size_t p , const size_t i , const size_t v ) const 
     *
     * \brief Returns the lower bound (found  by the LP solver) on the
     * \f$v\f$-th coordinate  of the  \f$i\f$-th second  difference of
     * the  \f$i\f$-th curve  segment of  the b-spline  curve threaded
     * into the channel.
     *
     * \param p The index of a curve segment of the b-spline.
     * \param i The  index of the \f$i\f$-th second  difference of the
     * \f$p\f$-th curve segment of the b-spline.
     * \param v The \f$v\f$-th  Cartesian coordinate of the \f$i\f$-th
     * control point of the \f$p\f$-th curve segment of the b-spline.
     *
     * \return  The  lower bound  (found  by  the  LP solver)  on  the
     * \f$v\f$-th  coordinate  of  the  \f$i\f$-th  second  difference
     * vector of  the \f$p\f$-th curve  segment of the  b-spline curve
     * threaded into the channel.
     *
     */
    double get_lower_bound_on_second_difference_value(
                                                      const size_t p ,
                                                      const size_t i ,
                                                      const size_t v
                                                     )
      const 
    {
      if ( p >= _np ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "Index of the curve segment is out of range" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }

      if ( ( i < 1 ) || ( i > 3 ) ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "Index of the second difference vector is out of range" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }

      if ( v >= 2 ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "Index of the Cartesian coordinate is out of range" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }

      if ( _secdiff.empty() ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "Second differences have not been computed" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }

      size_t index = ( 4 * 2 * p ) + ( 4 * ( i - 1 ) ) + v ;
      
      return _secdiff[ index ] ;
    }


    /**
     * \fn double get_upper_bound_on_second_difference_value( const size_t p , const size_t i , const size_t v ) const 
     *
     * \brief Returns the upper bound (found  by the LP solver) on the
     * \f$v\f$-th  coordinate  of  the  \f$i\f$-th  second  difference
     * vector of  the \f$p\f$-th curve  segment of the  b-spline curve
     * threaded into the channel.
     *
     * \param p The index of the curve segment of the b-spline.
     * \param i The  index of the \f$i\f$-th second  difference of the
     * \f$p\f$-th segment of the b-spline.
     * \param v The \f$v\f$-th  Cartesian coordinate of the \f$i\f$-th
     * control point of the \f$p\f$-th curve segment of the b-spline.
     *
     * \return  The  upper bound  (found  by  the  LP solver)  on  the
     * \f$v\f$-th  coordinate  of  the  \f$i\f$-th  second  difference
     * vector of  the \f$p\f$-th curve  segment of the  b-spline curve
     * threaded into the channel.
     *
     */
    double get_upper_bound_on_second_difference_value(
                                                      const size_t p ,
                                                      const size_t i ,
                                                      const size_t v
                                                     )
      const
    {
      if ( p >= _np ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "Index of the curve is out of range" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }

      if ( ( i < 1 ) || ( i > 3 ) ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "Index of the second difference vector is out of range" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }

      if ( v >= 2 ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "Index of the Cartesian coordinate is out of range" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }
      
      if ( _secdiff.empty() ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "Second differences have not been computed" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }

      size_t index = ( 4 * 2 * p ) + ( 4 * ( i - 1 ) ) + 2 + v ;
      
      return _secdiff[ index ] ;
    }


    /**
     * \fn double minimum_value() const
     *
     * \brief  Returns the  optimal (minimum)  value of  the objective
     * function of the instance of the channel problem as found by the
     * LP solver.
     *
     * \return The  optimal (minimum) value of  the objective function
     * of  the instance  of the  channel problem  as found  by the  LP
     * solver.
     *
     */
    double minimum_value() const
    {
      return _ofvalue ;
    }
        
        
    /**
     * \fn std::string get_solver_error_message( int error )
     *
     * \brief Returns the error message  of the GLPK solver associated
     * with a given error code.
     *
     * \param error Error code returned by the LP solver.
     *
     * \return The error message of  the GLPK solver associated with a
     * given error code.
     *
     */
     std::string get_solver_error_message( int error )
     {
       std::string message ;
       switch ( error ) {
         case GLP_EBADB :
           message = "Unable to start the search because the number of basic variables is not the same as the number of rows in the problem object." ;
           break ;
         case GLP_ESING :
           message = "Unable to start the search because the basis matrix corresponding to the initial basis is singular within the working precision." ;
           break ;
         case GLP_ECOND :
           message = "Unable to start the search because the basis matrix corresponding to the initial basis is ill-conditioned." ;
           break ;
         case GLP_EBOUND :
           message = "Unable to start the search because some double-bounded variables have incorrect bounds." ;
           break ;
         case GLP_EFAIL :
           message = "The search was prematurely terminated due to the solver failure." ;
           break ;
         case GLP_EOBJLL :
           message = "The search was prematurely terminated because the objective function being maximized has reached its lower limit and continues decreasing." ;
           break ;
         case GLP_EOBJUL :
           message = "The search was prematurely terminated because the objective function being minimized has reached its upper limit and continues increasing." ;
           break ;
         case GLP_EITLIM :
           message = "The search was prematurely terminated because the simplex iteration limit has been exceeded." ;
           break ;
         case GLP_ETMLIM :
           message = "The search was prematurely terminated because the time limit has been exceeded." ;
           break ;
         case GLP_ENOPFS :
           message = "The LP problem instance has no primal feasible solution." ;
           break ;
         case GLP_ENODFS :
           message = "The LP problem instance has no dual feasible solution." ;
           break ;
         default :
           message = "Unknown reason." ;
           break ;
       }
       
       return message ;
     }

  
  private:
  
    // ---------------------------------------------------------------
    //
    // Private methods
    //
    // ---------------------------------------------------------------

    /**
     * \fn void compute_min_max_constraints( size_t& eqline )
     *
     * \brief Computes the equations defining the min-max constraints.
     *
     * \param eqline A reference to the counter of equations.
     *
     */
    void compute_min_max_constraints( size_t& eqline ) ;


    /**
     * \fn void compute_correspondence_constraints( size_t& eqline )
     *
     * \brief Computes  the equations defining the  constraints on the
     * location of the  endpoints of the b-spline  curve threaded into
     * the channel.
     *
     * \param eqline A reference to the counter of equations. 
     *
     */
    void compute_correspondence_constraints( size_t& eqline ) ;


    /**
     * \fn void compute_sleeve_corners_in_channel_constraints( size_t& eqline )
     *
     * \brief  Computes the  equations defining  the constraints  that
     * ensure  that the  breakpoints  of the  sleeves  are inside  the
     * channel.
     *
     * \param eqline A reference to the counter of equations. 
     *
     */
    void compute_sleeve_corners_in_channel_constraints( size_t& eqline ) ;


    /**
     * \fn void compute_channel_corners_outside_sleeve_constraints( size_t& eqline )
     *
     * \brief  Computes the  equations defining  the constraints  that
     * ensure  that the  corners of  the  channel are  located on  the
     * boundary or outside the sleeve.
     *
     * \param eqline A reference to the counter of equations. 
     *
     */
    void compute_channel_corners_outside_sleeve_constraints( size_t& eqline ) ;


    /**
     * \fn void compute_sleeve_inside_csegment_constraints( size_t& eqline )
     *
     * \brief  Computes the  equations defining  the constraints  that
     * ensure the bspline segments  associated with a c-segment remain
     * inside it.
     *
     * \param eqline A reference to the counter of equations. 
     *
     */
    void compute_sleeve_inside_csegment_constraints( size_t& eqline ) ;

    
    /**
     * \fn void compute_normal_to_lower_envelope( const size_t s , double& nx , double& ny ) const
     *
     * \brief  Computes  an  outward  normal to  the  \f$s\f$-th  line
     * segment of the lower envelope of the channel.
     *
     * \param s Index of a line segment of the lower channel envelope.
     * \param nx A reference to  the first Cartesian coordinate of the
     * normal.
     * \param ny A reference to the second Cartesian coordinate of the
     * normal.
     *
     */
    void compute_normal_to_lower_envelope(
                                          const size_t s ,
                                          double& nx ,
                                          double& ny
                                         )
      const ;


    /**
     * \fn void compute_normal_to_upper_envelope( const size_t s , double& nx , double& ny ) const
     *
     * \brief  Computes  an  outward  normal to  the  \f$s\f$-th  line
     * segment of the upper envelope of the channel.
     *
     * \param s Index of a line segment of the upper channel envelope.
     * \param nx A reference to  the first Cartesian coordinate of the
     * normal.
     * \param ny A reference to the second Cartesian coordinate of the
     * normal.
     *
     */
    void compute_normal_to_upper_envelope(
                                          const size_t s ,
                                          double& nx ,
                                          double& ny
                                         )
      const ;
        
        
    /**
     * \fn void compute_normal_to_csection( const size_t s , double& nx , double& ny ) const
     *
     * \brief Computes  a normal  to the  \f$s\f$-th c-section  of the
     * channel.
     *
     * \param s Index of a c-section of the channel.
     * \param nx  A reference to  the first Cartesian coordinate  of the
     * normal.
     * \param ny A  reference to the Second Cartesian  coordinate of the
     * normal.
     *
     */
    void compute_normal_to_csection(
                                    const size_t s ,
                                    double& nx ,
                                    double& ny
                                   )
      const ;


    /**
     * \fn size_t compute_control_value_column_index( const size_t p , const size_t i , const size_t v ) const
     *
     * \brief Computes the  index of the linear  program matrix column
     * corresponding to  the x-  or y-coordinate  of the  i-th control
     * point of the  p-th segment of the b-spline to  be threaded into
     * the channel.
     *
     * \param p Index of the b-spline segment.
     * \param i Index of a control point of the p-th b-spline segment.
     * \param v Index of the x- or y-coordinate of the control point.
     *
     * \return  The   index  of  the  linear   program  matrix  column
     * corresponding to  the x-  or y-coordinate  of the  i-th control
     * point of the  p-th segment of the b-spline to  be threaded into
     * the channel.
     *
     */
    size_t compute_control_value_column_index(
                                              const size_t p ,
                                              const size_t i ,
                                              const size_t v
                                             )
      const ;


    /**
     * \fn void insert_coefficient( const size_t eqline , const size_t index , const double value )
     *
     * \brief Assigns  a value to the  coefficient of an unknown  of a
     * given constraint  of the linear  program (LP).  The  unknown is
     * identified by its corresponding  column index in the associated
     * matrix of the LP.
     *
     * \param eqline Matrix row corresponding to the constraint.
     * \param index  Matrix column index corresponding to the unknown.
     * \param value  Value to be assigned to the unknown coefficient.
     *
     */
    void insert_coefficient(
                            const size_t eqline ,
                            const size_t index ,
                            const double value
                           )
    {
      _coefficients[ eqline ].push_back(
                                        Coefficient(
                                                    eqline ,
                                                    index ,
                                                    value
                                                   )
                                       ) ;
      
      return ;
    }


    /**
     * \fn void insert_bound( const size_t eqline , const Bound::CONSTRAINTYPE type , const double value )
     *
     * \brief  Assigns  a real  value  to  the  right-hand side  of  a
     * constraint  (equality  or inequality)  of  an  instance of  the
     * linear program associated with the channel problem.
     *
     * \param eqline Matrix row corresponding to the constraint.
     * \param type   Type of the bound (==, >= or <=).
     * \param value Bound value (right-hand side of the constraint)
     *
     */
    void insert_bound(
                      const size_t eqline ,
                      const Bound::CONSTRAINTYPE type ,
                      const double value
                     )
    {
      _bounds[ eqline ] = Bound(
                                type ,
                                value ,
                                eqline
                               ) ;
      
      return ;
    }

    
    /**
     * \fn size_t compute_second_difference_column_index( const size_t p , const size_t i , const size_t l , const size_t v ) const
     *
     * \brief Computes the  index of the linear  program matrix column
     * corresponding to  the x- or  y-coordinate of the l-th  bound of
     * the i-th second difference of  the p-th segment of the b-spline
     * to be threaded into the channel.
     *
     * \param p Index of the b-spline segment.
     * \param i  Index of the  second difference of the  p-th b-spline
     * segment.
     * \param l Index of the l-th  bound of the second difference (0 -
     * lower bound; 1 - upper bound).
     * \param  v  Index  of  the  x- or  y-coordinate  of  the  second
     * difference bound.
     *
     * \return  The   index  of  the  linear   program  matrix  column
     * corresponding to  the x- or  y-coordinate of the l-th  bound of
     * the i-th second difference of  the p-th segment of the b-spline
     * to be threaded into the channel.
     *
     */
    size_t compute_second_difference_column_index(
                                                  const size_t p ,
                                                  const size_t i ,
                                                  const size_t l ,
                                                  const size_t v
                                                 )
      const ;


    /**
     * \fn size_t compute_index_of_endpoint_barycentric_coordinate( const size_t i ) const
     *
     * \brief Computes the  index of the linear  program matrix column
     * corresponding  to  the   barycentric  coordinate  defining  the
     * \f$i\f$-th endpoint of the b-spline.
     *
     * \param i Index of the \f$i\f$-th barycentric coordinate.
     *
     * \return  The   index  of  the  linear   program  matrix  column
     * corresponding  to  the   barycentric  coordinate  defining  the
     * \f$i\f$-th  endpoint  of  the  b-spline.
     *
     */
    size_t compute_index_of_endpoint_barycentric_coordinate( const size_t i ) const ;


    /**
     * \fn size_t compute_index_of_corner_barycentric_coordinate( const size_t i ) const
     *
     * \brief Computes the  index of the linear  program matrix column
     * corresponding to  the barycentric coordinate associated  with a
     * channel corner.
     *
     * \param i Index of a channel corner.
     *
     * \return  The   index  of  the  linear   program  matrix  column
     * corresponding to  the barycentric coordinate associated  with a
     * channel corner.
     *
     */
    size_t compute_index_of_corner_barycentric_coordinate( const size_t i ) const ;
    

    /**
     * \fn void insert_min_max_constraints( const size_t eqline , const size_t lo , const size_t up , const size_t b0 , const size_t b1 , const size_t b2 )
     *
     * \brief Inserts  the coefficients of the  equations defining the
     * three min-max  constraints into the matrix  associated with the
     * linear  program  (LP), and  sets  the  right-hand side  of  the
     * constraints as well.
     *
     * \param eqline A reference to the counter of equations.
     * \param  lo  Column  index  of  the lower  bound  for  a  second
     * difference.
     * \param  up  Column  index  of  the upper  bound  for  a  second
     * difference.
     * \param b0 Column index of  the first control value defining the
     * second difference.
     * \param b1 Column index of the second control value defining the
     * second difference.
     * \param b2 Column index of  the third control value defining the
     * second difference.
     *
     */
    void insert_min_max_constraints(
                                    const size_t eqline ,
                                    const size_t lo ,
                                    const size_t up ,
                                    const size_t b0 ,
                                    const size_t b1 ,
                                    const size_t b2
                                   ) ;


    /**
     * \fn void insert_extreme_point_correspondence_constraint( const size_t eqline , const std::vector< size_t >& col , const std::vector< double >& val , const double rhs )
     *
     * \brief  Inserts  into  the   linear  program  (LP)  matrix  the
     * coefficients of the unknowns and the right-hand side value of a
     * constraint  corresponding to  the location  of the  starting or
     * final point of the b-spline curve.
     *
     * \param eqline A reference to the counter of equations. 
     * \param  col  An  array  with   the  LP  matrix  column  indices
     * corresponding to the unknowns of the correspondence constraint.
     * \param  val  An array  with  the  values corresponding  to  the
     * unknowns of the correspondence constraint.
     * \param rhs The right-hand side value of the constraint.
     *
     */
    void insert_extreme_point_correspondence_constraint(
                                                        const size_t eqline ,
                                                        const std::vector< size_t >& col ,
                                                        const std::vector< double >& val ,
                                                        const double rhs 
                                                       ) ;


    /**
     * \fn void insert_periodic_correspondence_constraints( const size_t eqline , const std::vector< size_t >& strx , const std::vector< size_t >& stry , const std::vector< size_t >& endx , const std::vector< size_t >& endy )
     *
     * \brief  Inserts  into  the   linear  program  (LP)  matrix  the
     * coefficients of the unknowns and  the right-hand side values of
     * the constraints that ensure that the first three control points
     * are the same as the last three control points (in this order).
     *
     * \param eqline A reference to the counter of equations. 
     * \param strx An  array with the column indices of  the LP matrix
     * corresponding to  the first Cartesian coordinates  of the first
     * three control points.
     * \param stry An  array with the column indices of  the LP matrix
     * corresponding to the second  Cartesian coordinates of the first
     * three control points.
     * \param endx An  array with the column indices of  the LP matrix
     * corresponding to  the first  Cartesian coordinates of  the last
     * three control points.
     * \param endy An  array with the column indices of  the LP matrix
     * corresponding to  the second Cartesian coordinates  of the last
     * three control points.
     *
     */
    void insert_periodic_correspondence_constraints(
                                                    const size_t eqline ,
                                                    const std::vector< size_t >& strx ,
                                                    const std::vector< size_t >& stry ,
                                                    const std::vector< size_t >& endx ,
                                                    const std::vector< size_t >& endy
                                                   ) ;


    /**
     * \fn void insert_nonlinear_terms_of_epiece_point_lower_bound( const size_t eqline , const double s , const size_t c , const std::vector< std::vector< std::vector< size_t > > >& sd , const std::vector< std::vector< double > >& nl , const std::vector< std::vector< double > >& nu )
     *
     * \brief Inserts the coefficients  of the second difference terms
     * of the  equation defining lower  bounds for the  e-piece points
     * into  the matrix  associated  with an  instance  of the  Linear
     * Program (LP).  The terms belong  to the constraint  that forces
     * the  e-piece points  to be  inside a  certain c-section  of the
     * channel.
     *
     * \param eqline A counter for the number of constraints.
     * \param s A parameter value identifying a point on the e-piece.
     * \param c An index identifying a c-segment of the channel.
     * \param sd Array with the LP matrix column indices corresponding
     * to the lower  and upper bounds on  second differences occurring
     * in the  equation defining the  e-piece points belonging  to the
     * c-segment.
     * \param nl Array of Cartesian  coordinates of outward normals to
     * the  supporting lines  of the  lower envelope  segments of  the
     * channel.
     * \param nu Array of Cartesian  coordinates of outward normals to
     * the  supporting lines  of the  upper envelope  segments of  the
     * channel.
     * 
     */
    void insert_nonlinear_terms_of_epiece_point_lower_bound(
                                                            const size_t eqline ,
                                                            const double s ,
                                                            const size_t c ,
                                                            const std::vector< std::vector< std::vector< size_t > > >& sd ,
                                                            const std::vector< std::vector< double > >& nl ,
                                                            const std::vector< std::vector< double > >& nu
                                                           ) ;


    /**
     * \fn void insert_nonlinear_terms_of_epiece_point_lower_bound( const size_t eqline , const double s , const size_t c , const std::vector< std::vector< std::vector< size_t > > >& sd , const std::vector< std::vector< double > >& ncsec )
     *
     * \brief Inserts the coefficients  of the second difference terms
     * of the  equation defining lower  bounds for the  e-piece points
     * into  the matrix  associated  with an  instance  of the  Linear
     * Program (LP).  The  terms belong to the  constraint that forces
     * one e-piece point to be on the  right or left side of a channel
     * c-section.
     *
     * \param eqline A counter for the number of constraints.
     * \param s A parameter value identifying a point on the e-piece.
     * \param c An index identifying a c-segment of the channel.
     * \param sd Array with the LP matrix column indices corresponding
     * to the lower  and upper bounds on  second differences occurring
     * in the  equation defining the  e-piece points belonging  to the
     * c-segment.
     * \param  ncsec   Array  of  Cartesian  coordinates   of  normals
     * (pointing  to  the  left)  to   the  supporting  lines  of  the
     * c-sections of the channel.
     * 
     */
    void insert_nonlinear_terms_of_epiece_point_lower_bound(
                                                            const size_t eqline ,
                                                            const double s ,
                                                            const size_t c ,
                                                            const std::vector< std::vector< std::vector< size_t > > >& sd ,
                                                            const std::vector< std::vector< double > >& ncsec
                                                           ) ;

    
    /**
     * \fn void insert_nonlinear_terms_of_epiece_point_upper_bound( const size_t eqline , const double s , const size_t c , const std::vector< std::vector< std::vector< size_t > > >& sd , const std::vector< std::vector< double > >& nl , const std::vector< std::vector< double > >& nu  )
     *
     * \brief  Inserts  into the  matrix  associated  with the  Linear
     * Program (LP) the coefficients of  the lower and upper bounds of
     * the  second difference  terms  of the  equation defining  upper
     * bounds  for  the e-piece  points.   These  terms occur  in  the
     * constraint that keep  the sleeve inside a  certain c-section of
     * the channel.
     *
     * \param eqline A counter for the number of constraints.
     * \param s A parameter value identifying a point on the e-piece.
     * \param c An index identifying a c-segment of the channel.
     * \param sd Array with the LP matrix column indices corresponding
     * to the lower  and upper bounds on  second differences occurring
     * in the equation defining the points on the e-piece matched with
     * the c-segment.
     * \param nl Array of Cartesian  coordinates of outward normals to
     * the  supporting lines  of the  lower envelope  segments of  the
     * channel.
     * \param nu Array of Cartesian  coordinates of outward normals to
     * the  supporting lines  of the  upper envelope  segments of  the
     * channel.
     * 
     */
    void insert_nonlinear_terms_of_epiece_point_upper_bound(
                                                            const size_t eqline ,
                                                            const double s ,
                                                            const size_t c ,
                                                            const std::vector< std::vector< std::vector< size_t > > >& sd ,
                                                            const std::vector< std::vector< double > >& nl ,
                                                            const std::vector< std::vector< double > >& nu
                                                           ) ;


    /**
     * \fn void insert_nonlinear_terms_of_epiece_point_upper_bound( const size_t eqline , const double s , const size_t c , const std::vector< std::vector< std::vector< size_t > > >& sd , const std::vector< std::vector< double > >& ncsec )
     *
     * \brief Inserts the coefficients  of the second difference terms
     * of the  equation defining upper  bounds for the  e-piece points
     * into  the matrix  associated  with an  instance  of the  Linear
     * Program (LP).  The  terms belong to the  constraint that forces
     * one e-piece point to be on the  right or left side of a channel
     * c-section.
     *
     * \param eqline A counter for the number of constraints.
     * \param s A parameter value identifying a point on the e-piece.
     * \param c An index identifying a c-segment of the channel.
     * \param sd Array with the LP matrix column indices corresponding
     * to the lower  and upper bounds on  second differences occurring
     * in the  equation defining the  e-piece points belonging  to the
     * c-segment.
     * \param  ncsec   Array  of  Cartesian  coordinates   of  normals
     * (pointing  to  the  left)  to   the  supporting  lines  of  the
     * c-sections of the channel.
     * 
     */
    void insert_nonlinear_terms_of_epiece_point_upper_bound(
                                                            const size_t eqline ,
                                                            const double s ,
                                                            const size_t c ,
                                                            const std::vector< std::vector< std::vector< size_t > > >& sd ,
                                                            const std::vector< std::vector< double > >& ncsec
                                                           ) ;
        
        
    /**
     * \fn void insert_linear_terms_of_epiece_point_bounds( const size_t eqline , const double s , const double t , const size_t , const size_t c , const std::vector< std::vector< size_t > >& cp , const std::vector< std::vector< double > >& nl , const std::vector< std::vector< double > >& nu )
     *
     * \brief  Inserts the  coefficients of  the linear  terms of  the
     * equation defining lower and upper bounds for the e-piece points
     * into the matrix associated with the Linear Program (LP).  These
     * terms occur in the constraint that enforces an e-piece point to
     * stay inside channel.
     *
     * \param eqline A counter for the number of constraints.
     * \param s A parameter value identifying a point on the e-piece.
     * \param t A parameter value  identifying the b-spline point that
     * corresponds to the point on the e-piece at parameter \e s.
     * \param p Index of the  b-spline segment containing the b-spline
     * point at parameter \e t.
     * \param c An  index identifying the c-segment  the e-piece point
     * belongs to.
     * \param cp Array with the LP matrix column indices corresponding
     * to the control  points of the b-spline  defining the \f$p\f$-th
     * piece of the curve.
     * \param nl Array of Cartesian  coordinates of outward normals to
     * the  supporting lines  of the  lower envelope  segments of  the
     * channel.
     * \param nu Array of Cartesian  coordinates of outward normals to
     * the  supporting lines  of the  upper envelope  segments of  the
     * channel.
     *
     */
    void insert_linear_terms_of_epiece_point_bounds(
                                                    const size_t eqline ,
                                                    const double s ,
                                                    const double t ,
                                                    const size_t p ,
                                                    const size_t c ,
                                                    const std::vector< std::vector< size_t > >& cp ,
                                                    const std::vector< std::vector< double > >& nl ,
                                                    const std::vector< std::vector< double > >& nu
                                                   ) ;


    /**
     * \fn void insert_linear_terms_of_epiece_point_bounds( const size_t eqline , const double s , const double t , const size_t , const size_t c , const std::vector< std::vector< size_t > >& cp , const std::vector< std::vector< double > >& ncsec )
     *
     * \brief  Inserts the  coefficients of  the linear  terms of  the
     * equation defining lower and upper bounds for the e-piece points
     * into the matrix associated with the Linear Program (LP).  These
     * terms occur in the constraint that enforces an e-piece point to
     * stay  either  on  the  right  or on  left  side  of  a  channel
     * c-section.
     *
     * \param eqline A counter for the number of constraints.
     * \param s A parameter value identifying a point on the e-piece.
     * \param t A parameter value  identifying the b-spline point that
     * corresponds to the point on the e-piece at parameter \e s.
     * \param p Index of the  b-spline segment containing the b-spline
     * point at parameter \e t.
     * \param c An  index identifying the c-segment  the e-piece point
     * belongs to.
     * \param cp Array with the LP matrix column indices corresponding
     * to the control  points of the b-spline  defining the \f$p\f$-th
     * piece of the curve.
     * \param  ncsec   Array  of  Cartesian  coordinates   of  normals
     * (pointing  to  the  left)  to   the  supporting  lines  of  the
     * c-sections of the channel.
     *
     */
    void insert_linear_terms_of_epiece_point_bounds(
                                                    const size_t eqline ,
                                                    const double s ,
                                                    const double t ,
                                                    const size_t p ,
                                                    const size_t c ,
                                                    const std::vector< std::vector< size_t > >& cp ,
                                                    const std::vector< std::vector< double > >& ncsec
                                                   ) ;


    /**
     * \fn void insert_rhs_of_sleeve_corners_in_channel_constraints( const size_t eqline , const size_t c , const std::vector< std::vector< double > >& nl , const std::vector< std::vector< double > >& nu )
     *
     * \brief  Inserts  into the  matrix  associated  with the  Linear
     * Program (LP) the right-hand side values of the constraints that
     * enforce  a sleeve  point  to  stay inside  a  c-segment of  the
     * channel. The  type of each constraint  (equality or inequality:
     * ==, >= or <=) is also set here.
     *
     * \param eqline A counter for the number of constraints.
     * \param c An  index identifying the c-segment  the e-piece point
     * belongs to.
     * \param nl Array of Cartesian  coordinates of outward normals to
     * the  supporting lines  of the  lower envelope  segments of  the
     * channel.
     * \param nu Array of Cartesian  coordinates of outward normals to
     * the  supporting lines  of the  upper envelope  segments of  the
     * channel.
     * 
     */
    void insert_rhs_of_sleeve_corners_in_channel_constraints(
                                                             const size_t eqline ,
                                                             const size_t c ,
                                                             const std::vector< std::vector< double > >& nl ,
                                                             const std::vector< std::vector< double > >& nu
                                                            ) ;


    /**
     * \fn void insert_rhs_of_sleeve_inside_csegment_constraints( const size_t eqline , const size_t c , const std::vector< std::vector< double > >& ncsec )
     *
     * \brief  Inserts  into the  matrix  associated  with the  Linear
     * Program (LP) the right-hand side values of the constraints that
     * enforce one e-piece  breakpoint to stay on the right  side of a
     * c-section  of the  channel, and  another e-piece  breakpoint to
     * stay on the left side of the same c-section.
     *
     * \param eqline A counter for the number of constraints.
     * \param c An index identifying  the c-segment the e-piece points
     * belongs to.
     * \param  ncsec   Array  of  Cartesian  coordinates   of  normals
     * (pointing  to  the  left)  to   the  supporting  lines  of  the
     * c-sections of the channel.
     * 
     */
    void insert_rhs_of_sleeve_inside_csegment_constraints(
                                                          const size_t eqline ,
                                                          const size_t c ,
                                                          const std::vector< std::vector< double > >& ncsec
                                                         ) ;


    /**
     * \fn void evaluate_bounding_polynomial( const size_t j , const double t , double& lower , double& upper )
     *
     * \brief Obtains a  lower bound and an upper bound  for the value
     * of  a  precomputed, bounding polynomial  at  a given  parameter
     * value.
     *
     * \param j An index for the precomputed, bounding polynomial.
     * \param t A parameter value.
     * \param lower A reference to the lower bound.
     * \param upper A reference to the upper bound.
     *
     */
    void evaluate_bounding_polynomial(
                                      const size_t j ,
                                      const double t ,
                                      double& lower ,
                                      double& upper
                                     ) ;


    /**
     * \fn void insert_csegment_constraint( const size_t eqline , const double lower , const double upper , const size_t sdlo , const size_t sdup , const double normal )
     *
     * \brief Inserts the  coefficients of the lower  and upper bounds
     * of  a  constraint  second   difference  term  into  the  matrix
     * associated with  an instance of  the linear program  (LP).  The
     * term  belongs to  the equation  defining the  upper (or  lower)
     * bound of  a point of  an e-piece.  The constraint  ensures that
     * the point lies  on a specific side of  the oriented suppporting
     * line of one of the four line segments delimiting a c-segment of
     * the channel.
     *
     * \param eqline A counter for the number of constraints. 
     * \param lower  Coefficient of the second  difference lower bound
     * term.
     * \param upper  Coefficient of the second  difference upper bound
     * term.
     * \param sdlo The index of  the LP matrix column corresponding to
     * the second difference lower bound term.
     * \param sdup The index of  the LP matrix column corresponding to
     * the second difference upper bound term.
     * \param normal A normal to a supporting, oriented line of one of
     * the four line  segments delimiting a specific  c-segment of the
     * channel.
     *
     */
    void insert_csegment_constraint(
                                    const size_t eqline ,
                                    const double lower ,
                                    const double upper ,
                                    const size_t sdlo ,
                                    const size_t sdup ,
                                    const double normal
                                   ) ;


    /**
     * \fn void insert_csegment_constraint( const size_t eqline , const double c0 , const double c1 , const double c2 , const double c3 , const size_t b0 , const size_t b1 , const size_t b2 , const size_t b3 , const double normal )
     *
     * \brief  Inserts the  coefficients of  the linear  terms of  the
     * upper and  lower bounds of  an e-piece point equation  into the
     * matrix associated with an instance  of the linear program (LP).
     * The  constraint ensures  that  the point  of  the e-piece  lies
     * inside a c-segment of the channel.
     *
     * \param eqline A counter for the number of constraints.
     * \param  c0  Coefficient  of  the first  control  point  of  the
     * b-spline segment  containing the curve point  associated to the
     * e-piece point.
     * \param  c1  Coefficient of  the  second  control point  of  the
     * b-spline segment  containing the curve point  associated to the
     * e-piece point.
     * \param  c2  Coefficient  of  the third  control  point  of  the
     * b-spline segment  containing the curve point  associated to the
     * e-piece point.
     * \param  c3  Coefficient of  the  fourth  control point  of  the
     * b-spline segment  containing the curve point  associated to the
     * e-piece point.
     * \param b0  Index of the  LP matrix column corresponding  to the
     * first  control point  of  the b-spline  segment containing  the
     * curve point associated to the e-piece point.
     * \param b1  Index of the  LP matrix column corresponding  to the
     * second  control point  of the  b-spline segment  containing the
     * curve point associated to the e-piece point.
     * \param b2  Index of the  LP matrix column corresponding  to the
     * third  control point  of  the b-spline  segment containing  the
     * curve point associated to the e-piece point.
     * \param b3  Index of the  LP matrix column corresponding  to the
     * fourth  control point  of the  b-spline segment  containing the
     * curve point associated to the e-piece point.
     * \param normal A normal to a supporting, oriented line of one of
     * the four line  segments delimiting a specific  c-segment of the
     * channel.
     *
     */
    void insert_csegment_constraint(
                                    const size_t eqline ,
                                    const double c0 ,
                                    const double c1 ,
                                    const double c2 ,
                                    const double c3 ,
                                    const size_t b0 ,
                                    const size_t b1 ,
                                    const size_t b2 ,
                                    const size_t b3 ,
                                    const double normal
                                   ) ;


    /**
     * \fn void insert_csegment_constraint( const size_t eqline , const double c0 , const double c1 , const double c2 , const double c3 , const double c4 , const size_t b0 , const size_t b1 , const size_t b2 , const size_t b3 , const size_t b4 , const double normal )
     *
     * \brief  Inserts the  coefficients of  the linear  terms of  the
     * upper and  lower bounds of  an e-piece point equation  into the
     * matrix associated with an instance  of the linear program (LP).
     * This point belongs to an  e-piece segment whose endpoints bound
     * b-spline curve  points in  two distinct, but  consecutive curve
     * segments.  The  constraint ensures that the  e-piece point lies
     * inside a c-segment of the channel.
     *
     * \param eqline A counter for the number of constraints.
     * \param  c0  Coefficient  of  the first  control  point  of  the
     * b-spline segment containing the curve point associated with the
     * right endpoint of the e-piece segment.
     * \param  c1  Coefficient of  the  second  control point  of  the
     * b-spline segment containing the curve point associated with the
     * right endpoint of the e-piece segment.
     * \param  c2  Coefficient  of  the third  control  point  of  the
     * b-spline segment containing the curve point associated with the
     * right endpoint of the e-piece segment.
     * \param  c3  Coefficient of  the  fourth  control point  of  the
     * b-spline segment containing the curve point associated with the
     * right endpoint of the e-piece segment.
     * \param  c4  Coefficient of  the  fourth  control point  of  the
     * b-spline segment containing the curve point associated with the
     * left endpoint of the e-piece segment.
     * \param b0  Index of the  LP matrix column corresponding  to the
     * first  control point  of  the b-spline  segment containing  the
     * curve point associated  with the right endpoint  of the e-piece
     * segment.
     * \param b1  Index of the  LP matrix column corresponding  to the
     * second  control point  of the  b-spline segment  containing the
     * curve point associated  with the right endpoint  of the e-piece
     * segment.
     * \param b2  Index of the  LP matrix column corresponding  to the
     * third  control point  of  the b-spline  segment containing  the
     * curve point associated  with the right endpoint  of the e-piece
     * segment.
     * \param b3  Index of the  LP matrix column corresponding  to the
     * fourth  control point  of the  b-spline segment  containing the
     * curve point associated  with the right endpoint  of the e-piece
     * segment.
     * \param b4  Index of the  LP matrix column corresponding  to the
     * fourth  control point  of the  b-spline segment  containing the
     * curve point  associated with the  left endpoint of  the e-piece
     * segment.
     * \param normal A normal to a supporting, oriented line of one of
     * the four line  segments delimiting a specific  c-segment of the
     * channel.
     *
     */
    void insert_csegment_constraint(
                                    const size_t eqline ,
                                    const double c0 ,
                                    const double c1 ,
                                    const double c2 ,
                                    const double c3 ,
                                    const double c4 ,
                                    const size_t b0 ,
                                    const size_t b1 ,
                                    const size_t b2 ,
                                    const size_t b3 ,
                                    const size_t b4 ,
                                    const double normal
                                   ) ;


    /**
     * \fn int solve_lp( const size_t rows , const size_t cols )
     *
     * \brief Solves  the linear program corresponding  to the channel
     * problem.
     *
     * \param rows The number of constraints of the linear program.
     * \param cols The number of unknowns of the linear program.
     *
     * \return  The code  returned by  the LP  solver to  indicate the
     * status  of  the  computation  of the  solution  of  the  linear
     * program.
     */
    int solve_lp(
                 const size_t rows ,
                 const size_t cols
                )  ;


    /**
     * \fn void set_up_lp_constraints( glp_prob* lp ) const
     *
     * \brief Assemble the matrix of  constraints of the linear program,
     * and define  the type (equality  or inequality) and bounds  on the
     * constraints.
     *
     * \param lp A pointer to the instance of the LP program.
     *
     */
    void set_up_lp_constraints( glp_prob* lp ) const ;


    /**
     * \fn void set_up_structural_variables( glp_prob* lp ) const
     *
     * \brief  Define  lower and/or  upper  bounds  on the  structural
     * variables of  the linear  program corresponding to  the channel
     * problem.
     *
     * \param lp A pointer to the instance of the LP program.
     *
     */
    void set_up_structural_variables( glp_prob* lp ) const ;


    /**
     * \fn void set_up_objective_function( glp_prob* lp ) const
     *
     * \brief  Define the  objective  function of  the linear  program
     * corresponding to  the channel problem, which  is a minimization
     * problem.
     *
     * \param lp A pointer to the instance of the LP program.
     *
     */
    void set_up_objective_function( glp_prob* lp ) const ;


    /**
     * \fn void get_lp_solver_result_information( glp_prob* lp )
     *
     * \brief Obtain the optimal values found by the LP solver for the
     * structural values  of the  linear programming  corresponding to
     * the channel problem.
     *
     * \param lp A pointer to the instance of the LP program.
     *
     */
    void get_lp_solver_result_information( glp_prob* lp ) ;
    

  } ;
  

}

/** @} */ //end of group class.
