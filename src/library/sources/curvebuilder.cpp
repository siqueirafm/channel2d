/** 
 * \file curvebuilder.cpp
 *
 * \brief Implementation of a class  for threading a b-spline curve of
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

#include "curvebuilder.hpp"       // CurveBuilder

#include "exceptionobject.hpp"    // ExceptionObject
#include "a3.hpp"                 // a3

extern "C" {
#include "glpk.h"                 // For all GLPK functions
}

#include <cmath>                  // sqrt, fabs, ceil, floor
#include <cassert>                // assert
#include <sstream>                // std::stringstream
#include <iostream>               // std::cerr, std::endl
#include <vector>                 // std::vector
#include <cfloat>                 // DBL_MAX
#include <cstdlib>                // exit, EXIT_FAILURE


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
   * \fn CurveBuilder::CurveBuilder( size_t np , size_t nc , bool closed , double* lx , double* ly , double* ux , double* uy )
   *
   * \brief Creates an instance of this class.
   *
   * \param np The number of b-spline segments.
   * \param nc The number of c-segments of the channel.
   * \param closed A flag to indicate whether the channel is closed.
   * \param lx  A pointer to  an array  with the x-coordinates  of the
   * lower envelope of the channel.
   * \param ly  A pointer to  an array  with the y-coordinates  of the
   * lower envelope of the channel.
   * \param ux  A pointer to  an array  with the x-coordinates  of the
   * upper envelope of the channel.
   * \param uy  A pointer to  an array  with the y-coordinates  of the
   * upper envelope of the channel.
   *
   */
  CurveBuilder::CurveBuilder(
                             size_t np ,
                             size_t nc ,
                             bool closed ,
                             double* lx ,
                             double* ly ,
                             double* ux ,
                             double* uy 
                            )
  :
    _np( np ) ,
    _nc( nc ) ,
    _closed ( closed )
  {
    if ( _closed ) {
      if ( _np < 4 ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "The number of curve segments must be at least 4 for a closed curve" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }
      if ( _nc < 3 ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "The number of segments of a closed channel must be at least 3" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }
    }
    else {
      if ( _np < 1 ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "The number of curve segments must be at least 1" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }
      if ( _nc < 1 ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "The number of segments of an open channel must be at least 1" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }
    }
    
    size_t nn = ( _closed ) ? _nc : ( _nc + 1 ) ;

    _lxcoords.resize( nn ) ;
    _lycoords.resize( nn ) ;
    _uxcoords.resize( nn ) ;
    _uycoords.resize( nn ) ;

    for ( unsigned i = 0 ; i < nn ; i++ ) {
      _lxcoords[ i ] = lx[ i ] ;
      _lycoords[ i ] = ly[ i ] ;
      _uxcoords[ i ] = ux[ i ] ;
      _uycoords[ i ] = uy[ i ] ;
    }

    _tf = new a3() ;

    _ofvalue = DBL_MAX ;

    return ;
  }


  /**
   * \fn CurveBuilder::CurveBuilder( const CurveBuilder& b )
   *
   * \brief Clones an instance of this class.
   *
   * \param b A reference to another instance of this class.
   *
   */
  CurveBuilder::CurveBuilder( const CurveBuilder& b )
  :
    _np( b._np ) ,
    _nc( b._nc ) ,
    _closed( b._closed ) ,
    _lxcoords( b._lxcoords ) ,
    _lycoords( b._lycoords ) ,
    _uxcoords( b._uxcoords ) ,
    _uycoords( b._uycoords ) ,
    _tf( b._tf ) ,
    _coefficients( b._coefficients ) ,
    _bounds( b._bounds ) ,
    _ctrlpts( b._ctrlpts ) ,
    _secdiff( b._secdiff ) ,
    _ofvalue( b._ofvalue )
  {
  }

  /**
   * \fn bool CurveBuilder::build( int& error )
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
  bool
  CurveBuilder::build( int& error )
  {
    // Compute the number  of linear constraints (i.e.,  the number of
    // rows of the matrix) of the linear program whose solution yields
    // the curve.
    size_t rows = (
                     ( 6 * ( _np + 1 ) )         // min-max
                   + ( ( _closed ) ? 8 : 4 )     // correspondence
                   + ( 2 * ( _nc - 1 ) )         // channel corners
                   + ( ( 3 * 4 * _np ) - 4 )     // sleeve corners
		               + ( 4 * ( _nc - 1 ) )         // sleeve in csegments
                  ) ;

    // Compute  the  unknowns (i.e.,  the  number  of columns  of  the
    // matrix)  of  the  linear  program  whose  solution  yields  the
    // b-spline curve.
    size_t cols = ( 6 * _np ) + 10 + ( ( _closed ) ? 1 : 2 ) + ( _nc - 1 ) ;

    //
    // Allocate memory for the array of coefficients and bounds.
    //
    _coefficients.resize( rows ) ;
    _bounds.resize( rows ) ;

    //
    // Initialize the equation counter.
    //
    size_t eqline = 0 ;

    //
    // Compute the min-max constraints.
    //
    compute_min_max_constraints( eqline ) ;

    //
    // Compute the correspondence constraints.
    //
    compute_correspondence_constraints( eqline ) ;

    //
    // Compute channel corners outside sleeve constraints.
    //
    compute_channel_corners_outside_sleeve_constraints( eqline ) ;

    //
    // Compute the sleeve corners in channel constraints.
    //
    compute_sleeve_corners_in_channel_constraints( eqline ) ;

    //
    // Compute the sleeve inside csegment constraints.
    //
    compute_sleeve_inside_csegment_constraints( eqline ) ;
    
    //
    // Solve the LP and get the solution.
    //
    error = solve_lp( rows , cols ) ;
    
    return ( error == 0 ) ;
  }

  
  // -----------------------------------------------------------------
  //
  // Private methods
  //
  // -----------------------------------------------------------------


  /**
   * \fn void CurveBuilder::compute_min_max_constraints( size_t& eqline )
   *
   * \brief Computes the equations defining the min-max constraints.
   *
   * \param eqline A reference to the counter of equations. 
   *
   */
  void
  CurveBuilder::compute_min_max_constraints( size_t& eqline )
  {
    //
    // Obtain the min-max constraints for each second difference.
    //
    
    for ( size_t j = 1 ; j < 3 ; j++ ) {
      for ( size_t v = 0 ; v < 2 ; v++ ) {
        // Get  the column  indices of  the lower  bound and  of the
        // upper bound  of the  v-th coordinate  of the  j-th second
        // difference.
        size_t jl = compute_second_difference_column_index(
                                                           0 ,
                                                           j ,
                                                           0 ,
                                                           v
                                                          ) ;
          
        size_t ju = compute_second_difference_column_index(
                                                           0 ,
                                                           j ,
                                                           1 ,
                                                           v
                                                          ) ;
          
        // Get the column indices of  the v-th coordinates that define
        // the j-th second difference of the p-th curve segment.
        size_t c1 = compute_control_value_column_index(
                                                       0 ,
                                                       j - 1 ,
                                                       v
                                                      ) ;
        size_t c2 = compute_control_value_column_index(
                                                       0,
                                                       j ,
                                                       v
                                                      ) ;
        size_t c3 = compute_control_value_column_index(
                                                       0 ,
                                                       j + 1 ,
                                                       v
                                                      ) ;
          
        //
        // Set the nonzero coefficients of the next three equations.
        //
        insert_min_max_constraints(
                                   eqline ,
                                   jl ,
                                   ju ,
                                   c1 ,
                                   c2 ,
                                   c3
                                  ) ;

        //
        // Increment equation counter
        //
        eqline += 3 ;
      }
    }

    for ( size_t p = 1 ; p < _np ; p++ ) {
      for ( size_t v = 0 ; v < 2 ; v++ ) {
        // Get the column indices of the  lower bound and of the upper
        // bound of the v-th coordinate of the 2nd second difference.
        size_t jl = compute_second_difference_column_index(
                                                           p ,
                                                           2 ,
                                                           0 ,
                                                           v
                                                          ) ;
          
        size_t ju = compute_second_difference_column_index(
                                                           p ,
                                                           2 ,
                                                           1 ,
                                                           v
                                                          ) ;
	  
        // Get the column indices of  the v-th coordinates that define
        // the i-th second difference of the p-th curve segment.
        size_t c1 = compute_control_value_column_index(
                                                       p ,
                                                       1 ,
                                                       v
                                                      ) ;
        size_t c2 = compute_control_value_column_index(
                                                       p ,
                                                       2 ,
                                                       v
                                                      ) ;
        size_t c3 = compute_control_value_column_index(
                                                       p ,
                                                       3 ,
                                                       v
                                                      ) ;

        //
        // Set the nonzero coefficients of the next three equations.
        //
        insert_min_max_constraints(
                                   eqline ,
                                   jl ,
                                   ju ,
                                   c1 ,
                                   c2 ,
                                   c3
                                  ) ;
        
        //
        // Increment equation counter
        //
        eqline += 3 ;
      }
    }

    return ;
  }


  /**
   * \fn void CurveBuilder::compute_correspondence_constraints( size_t& eqline )
   *
   * \brief  Computes the  equations defining  the constraints  on the
   * location of the endpoints of the b-spline curve threaded into the
   * channel.
   *
   * \param eqline A reference to the counter of equations. 
   *
   */
  void
  CurveBuilder::compute_correspondence_constraints( size_t& eqline )
  {
    // Get the column  index of the x- and y-coordinates  of the first
    // three control  points of the  b-spline to be threaded  into the
    // channel.

    std::vector< size_t > strx( 4 ) ;

    strx[ 0 ] = compute_control_value_column_index( 0 , 0 , 0 ) ;
    strx[ 1 ] = compute_control_value_column_index( 0 , 1 , 0 ) ;
    strx[ 2 ] = compute_control_value_column_index( 0 , 2 , 0 ) ;

    // Get  the column  index of  the barycentric  coordinate defining
    // the first  endpoint of the  b-spline with respect to  the first
    // channel points.
    
    strx[ 3 ] = compute_index_of_endpoint_barycentric_coordinate( 0 ) ;

    //
    // Get the coefficients of the unknowns of the constraint.
    //
    std::vector< double > vals( 4 ) ;

    vals[ 0 ] = double( 1 ) / double( 6 ) ;
    vals[ 1 ] = double( 2 ) / double( 3 ) ;
    vals[ 2 ] = vals[ 0 ] ;
    vals[ 3 ] = _lxcoords[ 0 ] - _uxcoords[ 0 ]  ;

    //
    // Constraint corresponding to first Cartesian coordinate.
    //
    insert_extreme_point_correspondence_constraint(
                                                   eqline ,
                                                   strx ,
                                                   vals ,
                                                   _lxcoords[ 0 ]
                                                  ) ;

    ++eqline ;  // increment the equation counter.

    //
    // Constraint corresponding to second Cartesian coordinate.
    //

    std::vector< size_t > stry( 4 ) ;

    stry[ 0 ] = compute_control_value_column_index( 0 , 0 , 1 ) ;
    stry[ 1 ] = compute_control_value_column_index( 0 , 1 , 1 ) ;
    stry[ 2 ] = compute_control_value_column_index( 0 , 2 , 1 ) ;
    stry[ 3 ] = strx[ 3 ] ;

    vals[ 3 ] = _lycoords[ 0 ] - _uycoords[ 0 ]  ;

    insert_extreme_point_correspondence_constraint(
                                                   eqline ,
                                                   stry ,
                                                   vals ,
                                                   _lycoords[ 0 ]
                                                  ) ;

    ++eqline ;  // increment the equation counter.
  
    // ---------------------------------------------------------------
    //
    // If the curve is closed, then the last three control points must
    // match the  first three control  points. Otherwise, we  must fix
    // the position of the final of  the curve, which differs from the
    // starting one.
    //
    // ---------------------------------------------------------------

    // Get the  column index of the  x- and y-coordinates of  the last
    // three control  points of the  b-spline to be threaded  into the
    // channel.

    std::vector< size_t > endx( 4 ) ;

    endx[ 0 ] = compute_control_value_column_index( _np - 1 , 1 , 0 ) ;
    endx[ 1 ] = compute_control_value_column_index( _np - 1 , 2 , 0 ) ;
    endx[ 2 ] = compute_control_value_column_index( _np - 1 , 3 , 0 ) ;

    std::vector< size_t > endy( 4 ) ;

    endy[ 0 ] = compute_control_value_column_index( _np - 1 , 1 , 1 ) ;
    endy[ 1 ] = compute_control_value_column_index( _np - 1 , 2 , 1 ) ;
    endy[ 2 ] = compute_control_value_column_index( _np - 1 , 3 , 1 ) ;

    if ( _closed ) {
      //
      // Compute the  equations that  match the  first three  and last
      // three control points of the b-spline: last is equal to third,
      // ...
      //

      insert_periodic_correspondence_constraints(
                                                 eqline ,
                                                 strx ,
                                                 stry ,
                                                 endx ,
                                                 endy 
                                                ) ;

      eqline += 6 ;   // increment the equation counter.
    }
    else {
      //
      // Compute the equations determining the b-spline final point.
      //
      
      // Get the  column index  of the  barycentric coordinate  of the
      // final point of the b-spline  with respect to the final points
      // of the channel.
    
      endx[ 3 ] = compute_index_of_endpoint_barycentric_coordinate( 1 ) ;
      vals[ 3 ] = _lxcoords[ _nc ] - _uxcoords[ _nc ] ;

      //
      // Constraint corresponding to first Cartesian coordinate.
      //
      insert_extreme_point_correspondence_constraint(
                                                     eqline ,
                                                     endx ,
                                                     vals ,
                                                     _lxcoords[ _nc ]
                                                    ) ;

      endy[ 3 ] = endx[ 3 ] ;    
      vals[ 3 ] = _lycoords[ _nc ] - _uycoords[ _nc ] ;

      ++eqline ;  // increment the equation counter.
      
      insert_extreme_point_correspondence_constraint(
                                                     eqline ,
                                                     endy ,
                                                     vals ,
                                                     _lycoords[ _nc ]
                                                    ) ;

      ++eqline ;  // increment the equation counter.
    }
      
    return ;
  }


  /**
   * \fn void CurveBuilder::compute_sleeve_corners_in_channel_constraints( size_t& eqline )
   *
   * \brief  Computes  the  equations defining  the  constraints  that
   * ensure  that  the  breakpoints  of the  sleeves  are  inside  the
   * channel.
   *
   * \param eqline A reference to the counter of equations. 
   *
   */
  void
  CurveBuilder::compute_sleeve_corners_in_channel_constraints( size_t& eqline )
  {
    // 
    // Compute outward normals to the line segments of the channel.
    //
    std::vector< std::vector< double > > nl( _nc , std::vector< double >( 2 , 0 ) ) ;
    std::vector< std::vector< double > > nu( _nc , std::vector< double >( 2 , 0 ) ) ;
    
    for ( size_t c = 0 ; c < _nc ; c++ ) {
      compute_normal_to_lower_envelope( c , nl[ c ][ 0 ] , nl[ c ][ 1 ] ) ;
      compute_normal_to_upper_envelope( c , nu[ c ][ 0 ] , nu[ c ][ 1 ] ) ;
    }

    // Each segment of the b-spline must  be enclosed by a sleeve with
    // four breakpoints, two of which are shared with the previous and
    // next segment  (if any).  Each  breakpoint is constrained  to be
    // bounded by a pair of parallel segments (lower and upper) of the
    // channel.

    const size_t lo = ( 3 *   3 ) + 1 ;
    const size_t up = ( 3 * _np ) + 8 ;
    const double NcOverNp = _nc / double( _np ) ;

    for ( size_t u = lo ; u <= up ; u++ ) {
      //
      // Find the index of the channel segment corresponding to \e u.
      //
      double t = u / double( 3 ) ;
      double s = t - floor( t ) ;
      size_t p = ( size_t ) floor( t - 3 ) ;
      size_t c = ( size_t ) ( ( t - 3 ) * NcOverNp ) ;

      // Compute  the  column indices  of  the  linear program  matrix
      // corresponding to  the four  control points defining  the p-th
      // segment of the b-spline curve.

      std::vector< std::vector< size_t > > cp( 4 , std::vector< size_t >( 2 , 0 ) ) ;

      for ( size_t i = 0 ; i < 4 ; i++ ) {
        for  ( size_t j = 0 ; j < 2 ; j++ ) {
          cp[ i ][ j ] = compute_control_value_column_index( p , i , j ) ;
        }
      }

      // Compute  the  column indices  of  the  linear program  matrix
      // of the values of the second difference bounds associated with
      // the p-th segment of the b-spline curve.

      std::vector< std::vector< std::vector< size_t > > > sd(
                                                             2 ,
                                                             std::vector< std::vector< size_t > >
                                                             (
                                                               2 ,
                                                               std::vector< size_t >( 2 , 0 )
                                                             )
                                                            ) ;
        
      for ( size_t j = 1 ; j < 3 ; j++ ) {
        for ( size_t l = 0 ; l < 2 ; l++ ) {
          for ( size_t v = 0 ; v < 2 ; v++ ) {
            sd[ j - 1 ][ l ][ v ] = compute_second_difference_column_index( p , j , l , v ) ;
          }
        }
      }

      // ------------------------------------------------------------
      //
      // Process nonlinear terms of Constraint (3a).
      //
      // ------------------------------------------------------------

      //
      // Nonlinear terms of \f$\stackrel{e}{\sim}^p\f$
      //

      insert_nonlinear_terms_of_epiece_point_lower_bound(
                                                         eqline ,
                                                         s ,
                                                         c ,
                                                         sd ,
                                                         nl ,
                                                         nu
                                                        ) ;

      //
      // Nonlinear terms of \f$\stackrel{\sim}{e}^p\f$.
      //

      insert_nonlinear_terms_of_epiece_point_upper_bound(
                                                         eqline + 2 ,
                                                         s ,
                                                         c ,
                                                         sd ,
                                                         nl ,
                                                         nu
                                                        ) ;


      // ------------------------------------------------------------
      //
      // Process linear terms of Constraint (3a).
      //
      // ------------------------------------------------------------

      //
      // Linear terms of \f$\stackrel{e}{\sim}^p\f$
      //

      insert_linear_terms_of_epiece_point_bounds(
                                                 eqline ,
                                                 s ,
                                                 t ,
                                                 p ,
                                                 c ,
                                                 cp ,
                                                 nl ,
                                                 nu
                                                ) ;

      //
      // Linear terms of \f$\stackrel{\sim}{e}^p\f$.
      //

      insert_linear_terms_of_epiece_point_bounds(
                                                 eqline + 2 ,
                                                 s ,
                                                 t ,
                                                 p ,
                                                 c ,
                                                 cp ,
                                                 nl ,
                                                 nu
                                                ) ;

      // -------------------------------------------------------------
      //
      // Compute right-hand side of the constraints.
      //
      // -------------------------------------------------------------

      insert_rhs_of_sleeve_corners_in_channel_constraints(
                                                          eqline ,
                                                          c ,
                                                          nl ,
                                                          nu
                                                         ) ;

      //
      // Increment equation counter.
      //
      eqline += 4 ;
    }

    return ;
  }


  /**
   * \fn void CurveBuilder::compute_channel_corners_outside_sleeve_constraints( size_t& eqline )
   *
   * \brief  Computes  the  equations defining  the  constraints  that
   * ensure  that  the corners  of  the  channel  are located  on  the
   * boundary or outside the sleeve.
   *
   * \param eqline A reference to the counter of equations. 
   *
   */
  void
  CurveBuilder::compute_channel_corners_outside_sleeve_constraints( size_t& eqline )
  {
    //
    // This restriction applies to channels with at least 3 c-sections
    // only.
    //
    if ( _nc < 2 ) {
      return ;
    }

    // For  each  inner  corner  of   the  given  channel,  compute  a
    // constraint that ensures that the  channel corner is outside the
    // sleeve.
    const double NpOverNc = _np / double( _nc ) ;
    const double onesixth = 1 / double( 6 ) ;
  
    for ( size_t c = 1 ; c < _nc ; c++ ) {
      //
      // Find the parameter \e t corresponding to the \e c corner.
      //
      double t = ( c * NpOverNc ) + 3 ;

      //
      // Find the curve segment \e  p containing point at parameter \e
      // t.
      size_t p = ( size_t ) floor( t - 3 ) ;

      // Compute  the  column indices  of  the  linear program  matrix
      // corresponding to  the four  control points defining  the p-th
      // segment of the b-spline curve.

      std::vector< std::vector< size_t > > cp( 4 , std::vector< size_t >( 2 , 0 ) ) ;

      for ( size_t i = 0 ; i < 4 ; i++ ) {
        for  ( size_t j = 0 ; j < 2 ; j++ ) {
          cp[ i ][ j ] = compute_control_value_column_index( p , i , j ) ;
        }
      }

      //
      // Compute the coefficients of the control points.
      //
      double s = t - floor( t ) ;
      std::vector< double > coeffs( 4 ) ;

      coeffs[ 0 ] = onesixth * ( 1 + s * ( -3 + s * ( 3 - s ) ) ) ;
      coeffs[ 1 ] = onesixth * ( 4 + s * s * ( - 6 + 3 * s ) ) ;
      coeffs[ 2 ] = onesixth * ( 1 + s * ( 3 + s * ( 3 - 3 * s ) ) ) ;
      coeffs[ 3 ] = onesixth * ( s * s * s ) ;
      
      // Get  the   LP  matrix  column  index   corresponding  to  the
      // barycentric coordinate associated with the c-th corner of the
      // channel.
      size_t k = compute_index_of_corner_barycentric_coordinate( c ) ;

      //
      // Insert the constraints into the LP program.
      //
      for ( size_t i = 0 ; i < 4 ; i++ ) {
        insert_coefficient( eqline     , cp[ i ][ 0 ] , coeffs[ i ] ) ;
        insert_coefficient( eqline + 1 , cp[ i ][ 1 ] , coeffs[ i ] ) ;
      }
	
      insert_coefficient( eqline     , k , ( _lxcoords[ c ] - _uxcoords[ c ] ) ) ;
      insert_coefficient( eqline + 1 , k , ( _lycoords[ c ] - _uycoords[ c ] ) ) ;
	
      insert_bound( eqline     , Bound::EQT , _lxcoords[ c ] ) ;
      insert_bound( eqline + 1 , Bound::EQT , _lycoords[ c ] ) ;
	
      //
      // Increment equation counter.
      //
      eqline += 2 ;
    }

    return ;
  }


  /**
   * \fn void CurveBuilder::compute_sleeve_inside_csegment_constraints( size_t& eqline )
   *
   * \brief  Computes  the  equations defining  the  constraints  that
   * ensure the  bspline segments  associated with a  c-segment remain
   * inside it.
   *
   * \param eqline A reference to the counter of equations. 
   *
   */
  void
  CurveBuilder::compute_sleeve_inside_csegment_constraints( size_t& eqline )
  {
    //
    // This restriction applies to channels with at least 3 c-sections
    // only.
    //
    if ( _nc < 2 ) {
      return ;
    }

    //
    // Compute normals to the c-sections of the channel.
    //
    
    const size_t NumberOfCSections = ( _closed ) ? _nc : _nc + 1 ;
    
    std::vector< std::vector< double > > ncsec(
                                               NumberOfCSections ,
                                               std::vector< double >( 2 , 0 )
                                              ) ;
    
    for ( size_t c = 0 ; c < NumberOfCSections ; c++ ) {
      compute_normal_to_csection(
                                 c ,
                                 ncsec[ c ][ 0 ] ,
                                 ncsec[ c ][ 1 ]
                                ) ;
    }

    // For  each  inner  corner  of the  given  channel,  compute  two
    // constraints   which  ensure   that   the  e-piece   breakpoints
    // immediately  on the  right  (resp. left)  of the  corresponding
    // c-section remain  in the right  (resp.  left) c-segment  of the
    // channel.

    const double NpOverNc = _np / double( _nc ) ;
    const double onethird = 1 / double( 3 ) ;
    const double twothird = 2 * onethird ; 
  
    for ( size_t c = 1 ; c < _nc ; c++ ) {
      //
      // Find the parameter \e t corresponding to the \e c corner.
      //
      double t = ( c * NpOverNc ) + 3 ;

      //
      // Find the curve segment \e  p containing point at parameter \e
      // t.
      size_t p = ( size_t ) floor( t - 3 ) ;

      // Compute the  indices of  the curve segments  corresponding to
      // the e-piece breakpoints immediately to  the right and left of
      // point \c p(t).

      double s = t - floor( t ) ;

      size_t p1 , p2 ;
      double s1 , s2 ;
      
      if ( s == 0 ) {
        p1 = p - 1 ;
        p2 = p ;
        s1 = twothird ;
        s2 = onethird ;
      }
      else if ( s < onethird ) {
        p1 = p ;
        p2 = p ;
        s1 = 0 ;
        s2 = onethird ;
      }
      else if ( s < twothird ) {
        p1 = p ;
        p2 = p ;
        s1 = onethird ;
        s2 = twothird ;
      }
      else {
        p1 = p ;
        p2 = p ;
        s1 = twothird ;
        s2 = 1 ;
      }

      double t1 = p1 + 3 + s1 ;
      double t2 = p2 + 3 + s2 ;

      // Compute the column indices of  the LP matrix corresponding to
      // the  second  difference  bounds  associated  with  the  p1-th
      // segment.

      std::vector< std::vector< std::vector< size_t > > > sd1(
                                                              2 ,
                                                              std::vector< std::vector< size_t > >
                                                              (
                                                               2 ,
                                                               std::vector< size_t >( 2 , 0 )
                                                              )
                                                             ) ;
      
      for ( size_t j = 1 ; j < 3 ; j++ ) {
        for ( size_t l = 0 ; l < 2 ; l++ ) {
          for ( size_t v = 0 ; v < 2 ; v++ ) {
            sd1[ j - 1 ][ l ][ v ] = compute_second_difference_column_index(
                                                                            p1 ,
                                                                            j ,
                                                                            l ,
                                                                            v
                                                                           ) ;
          }
        }
      }

      // Compute the column indices of  the LP matrix corresponding to
      // the  second  difference  bounds  associated  with  the  p2-th
      // segment.

      std::vector< std::vector< std::vector< size_t > > > sd2(
                                                              2 ,
                                                              std::vector< std::vector< size_t > >
                                                              (
                                                               2 ,
                                                               std::vector< size_t >( 2 , 0 )
                                                              )
                                                             ) ;
      
      for ( size_t j = 1 ; j < 3 ; j++ ) {
        for ( size_t l = 0 ; l < 2 ; l++ ) {
          for ( size_t v = 0 ; v < 2 ; v++ ) {
            sd2[ j - 1 ][ l ][ v ] = compute_second_difference_column_index(
                                                                            p2 ,
                                                                            j ,
                                                                            l ,
                                                                            v
                                                                           ) ;
          }
        }
      }

      // Compute the column indices of  the LP matrix corresponding to
      // the four  control points  defining the  p1-th segment  of the
      // curve.

      std::vector< std::vector< size_t > > cp1( 4 , std::vector< size_t >( 2 , 0 ) ) ;
      
      for ( size_t i = 0 ; i < 4 ; i++ ) {
        for  ( size_t j = 0 ; j < 2 ; j++ ) {
          cp1[ i ][ j ] = compute_control_value_column_index( p1 , i , j ) ;
        }
      }

      // Compute the column indices of  the LP matrix corresponding to
      // the four  control points  defining the  p2-th segment  of the
      // curve.

      std::vector< std::vector< size_t > > cp2( 4 , std::vector< size_t >( 2 , 0 ) ) ;
      
      for ( size_t i = 0 ; i < 4 ; i++ ) {
        for  ( size_t j = 0 ; j < 2 ; j++ ) {
          cp2[ i ][ j ] = compute_control_value_column_index( p2 , i , j ) ;
        }
      }

      // -------------------------------------------------------------
      //
      // Process nonlinear terms of Constraint (3c).
      //
      // -------------------------------------------------------------
      
      //
      // Nonlinear terms of \f$\stackrel{e}{\sim}^p( s_1 )\f$
      //      
      insert_nonlinear_terms_of_epiece_point_lower_bound(
                                                         eqline ,
                                                         s1 ,
                                                         c ,
                                                         sd1 ,
                                                         ncsec
                                                        ) ;

      //
      // Nonlinear terms of \f$\stackrel{\sim}{e}^p( s_1 )\f$.
      //
      insert_nonlinear_terms_of_epiece_point_upper_bound(
                                                         eqline + 1 ,
                                                         s1 ,
                                                         c ,
                                                         sd1 ,
                                                         ncsec
                                                        ) ;

      //
      // Nonlinear terms of \f$\stackrel{e}{\sim}^p( s_2 )\f$
      //      
      insert_nonlinear_terms_of_epiece_point_lower_bound(
                                                         eqline + 2 ,
                                                         s2 ,
                                                         c ,
                                                         sd2 ,
                                                         ncsec
                                                        ) ;

      //
      // Nonlinear terms of \f$\stackrel{\sim}{e}^p( s_2 )\f$.
      //
      insert_nonlinear_terms_of_epiece_point_upper_bound(
                                                         eqline + 3 ,
                                                         s2 ,
                                                         c ,
                                                         sd2 ,
                                                         ncsec
                                                        ) ;

      // -------------------------------------------------------------
      //
      // Process linear terms of Constraint (3c).
      //
      // -------------------------------------------------------------

      //
      // Linear terms of \f$\stackrel{e}{\sim}^p( s_1 )\f$
      //     
      insert_linear_terms_of_epiece_point_bounds(
                                                 eqline ,
                                                 s1 ,
                                                 t1 ,
                                                 p1 ,
                                                 c ,
                                                 cp1 ,
                                                 ncsec
                                                ) ;

      //
      // Linear terms of \f$\stackrel{\sim}{e}^p( s_2 )\f$.
      //
      insert_linear_terms_of_epiece_point_bounds(
                                                 eqline + 2 ,
                                                 s2 ,
                                                 t2 ,
                                                 p2 ,
                                                 c ,
                                                 cp2 ,
                                                 ncsec
                                                ) ;

      // -------------------------------------------------------------
      //
      // Compute right-hand side of the constraints.
      //
      // -------------------------------------------------------------

      insert_rhs_of_sleeve_inside_csegment_constraints(
                                                       eqline ,
                                                       c ,
                                                       ncsec
                                                      ) ;	
      //
      // Increment equation counter.
      //
      eqline += 4 ;
    }

    return ;
  }
  
  
  /**
   * \fn void CurveBuilder::compute_normal_to_lower_envelope( const size_t s , double& nx , double& ny ) const
   *
   * \brief Computes an outward normal  to the \f$s\f$-th line segment
   * of the lower envelope of the channel.
   *
   * \param s Index of a line segment of the lower channel envelope.
   * \param nx  A reference to  the first Cartesian coordinate  of the
   * normal.
   * \param ny A  reference to the second Cartesian  coordinate of the
   * normal.
   *
   */
  void
  CurveBuilder::compute_normal_to_lower_envelope(
                                                 const size_t s ,
                                                 double& nx ,
                                                 double& ny
                                                )
    const
  {
#ifdef DEBUGMODE
    assert( s < _nc ) ;
#endif

    size_t t = s + 1 ;
    
    if ( _closed ) {
       t %= _nc ;
    }
    
    nx = _lycoords[ s ] - _lycoords[ t ] ;
    ny = _lxcoords[ t ] - _lxcoords[ s ] ;

    return ;
  }


  /**
   * \fn void CurveBuilder::compute_normal_to_upper_envelope( const size_t s , double& nx , double& ny ) const
   *
   * \brief Computes an outward normal  to the \f$s\f$-th line segment
   * of the upper envelope of the channel.
   *
   * \param s Index of a line segment of the upper channel envelope.
   * \param nx  A reference to  the first Cartesian coordinate  of the
   * normal.
   * \param ny A  reference to the second Cartesian  coordinate of the
   * normal.
   *
   */
  void
  CurveBuilder::compute_normal_to_upper_envelope(
                                                 const size_t s ,
                                                 double& nx ,
                                                 double& ny
                                                )
    const
  {
#ifdef DEBUGMODE
    assert( s < _nc ) ;
#endif
    
    size_t t = s + 1 ;
    
    if ( _closed ) {
      t %= _nc ;
    }
    
    nx = _uycoords[ t ] - _uycoords[ s ] ;
    ny = _uxcoords[ s ] - _uxcoords[ t ] ;

    return ;
  }

  
  /**
   * \fn void CurveBuilder::compute_normal_to_csection( const size_t s , double& nx , double& ny ) const
   *
   * \brief  Computes a  normal  to the  \f$s\f$-th  c-section of  the
   * channel.
   *
   * \param s Index of a c-section of the channel.
   * \param nx  A reference to  the first Cartesian coordinate  of the
   * normal.
   * \param ny A  reference to the Second Cartesian  coordinate of the
   * normal.
   *
   */
  void
  CurveBuilder::compute_normal_to_csection(
                                           const size_t s ,
                                           double& nx ,
                                           double& ny
                                          )
  const
  {
#ifdef DEBUGMODE
    assert( s <= _nc ) ;
#endif
  
    size_t t = ( _closed ) ? ( s % _nc ) : s ;
    
    nx = _lycoords[ t ] - _uycoords[ t ] ;
    ny = _uxcoords[ t ] - _lxcoords[ t ] ;
    
    return ;
  }
  

  /**
   * \fn size_t CurveBuilder::compute_control_value_column_index( const size_t p , const size_t i , const size_t v ) const
   *
   * \brief Computes  the index  of the  linear program  matrix column
   * corresponding to the x- or y-coordinate of the i-th control point
   * of  the p-th  segment of  the b-spline  to be  threaded into  the
   * channel.
   *
   * \param p Index of the b-spline segment.
   * \param i Index of a control point of the p-th b-spline segment.
   * \param v Index of the x- or y-coordinate of the control point.
   *
   * \return   The  index   of  the   linear  program   matrix  column
   * corresponding to the x- or y-coordinate of the i-th control point
   * of  the p-th  segment of  the b-spline  to be  threaded into  the
   * channel.
   *
   */
  size_t
  CurveBuilder::compute_control_value_column_index(
                                                   const size_t p ,
                                                   const size_t i ,
                                                   const size_t v 
                                                  )
    const
  {
#ifdef DEBUGMODE
    assert( p < _np ) ;
    assert( i <=  3 ) ;
    assert( v <=  1 ) ;
#endif

    return 2 * ( p + i ) + v ;
  }


  /**
   * \fn size_t CurveBuilder::compute_second_difference_column_index( const size_t p , const size_t i , const size_t l , const size_t v ) const
   *
   * \brief Computes  the index  of the  linear program  matrix column
   * corresponding to the  x- or y-coordinate of  the \f$l\f$-th bound
   * of the \f$i\f$-th second difference  of the \f$p\f$-th segment of
   * the b-spline to be threaded into the channel.
   *
   * \param p Index of the b-spline segment.
   * \param  i  Index  of  the second  difference  of  the  \f$p\f$-th
   * b-spline segment.
   * \param l Index  of the \f$l\f$-th bound of  the second difference
   * (0 - lower bound; 1 - upper bound).
   * \param  v Index  of  the \f$x\f$-  or  \f$y\f$-coordinate of  the
   * second difference bound.
   *
   * \return   The  index   of  the   linear  program   matrix  column
   * corresponding  to  the  \f$x\f$-  or  \f$y\f$-coordinate  of  the
   * \f$l\f$-th  bound  of the  \f$i\f$-th  second  difference of  the
   * \f$p\f$-th  segment  of the  b-spline  to  be threaded  into  the
   * channel.
   *
   */
  size_t
  CurveBuilder::compute_second_difference_column_index(
                                                       const size_t p ,
                                                       const size_t i ,
                                                       const size_t l ,
                                                       const size_t v 
                                                      )
    const
  {
#ifdef DEBUGMODE
    assert( p < _np ) ;
    assert( i >= 1 ) ;
    assert( i <= 2 ) ;
    assert( l <= 1 ) ;
    assert( v <= 1 ) ;
#endif
  
    size_t offset = ( 2 * _np ) + 6 ; 

    return offset + ( 4 * ( p + i - 1 ) ) + ( 2 * l ) + v ;
  }


  /**
   * \fn size_t CurveBuilder::compute_index_of_endpoint_barycentric_coordinate( const size_t i ) const
   *
   * \brief Computes  the index  of the  linear program  matrix column
   * corresponding   to  the   barycentric  coordinate   defining  the
   * \f$i\f$-th endpoint of the b-spline.
   *
   * \param i Index of the \f$i\f$-th barycentric coordinate.
   *
   * \return  The   index  of  the  linear   program  matrix  column
   * corresponding  to  the   barycentric  coordinate  defining  the
   * \f$i\f$-th  endpoint  of  the  b-spline.
   *
   */
  size_t
  CurveBuilder::compute_index_of_endpoint_barycentric_coordinate( const size_t i ) const
  {
    #ifdef DEBUGMODE
    if ( _closed ) {
      assert( i == 0 ) ;
    }
    else {
      assert( i <= 1 ) ;
    }
#endif
  
    size_t offset =  ( 6 * _np ) + 10 ;

    return offset + i ;
  }


  /**
   * \fn size_t CurveBuilder::compute_index_of_corner_barycentric_coordinate( const size_t i ) const
   *
   * \brief Computes  the index  of the  linear program  matrix column
   * corresponding  to the  barycentric coordinate  associated with  a
   * channel corner.
   *
   * \param i Index of a channel corner.
   *
   * \return   The  index   of  the   linear  program   matrix  column
   * corresponding  to the  barycentric coordinate  associated with  a
   * channel corner.
   *
   */
  size_t
  CurveBuilder::compute_index_of_corner_barycentric_coordinate( const size_t i ) const
  {
#ifdef DEBUGMODE
    assert( i >   0 ) ;
    assert( i < _nc ) ;
#endif
  
    size_t offset =  ( 6 * _np ) + 10 + ( ( _closed ) ? 0 : 1 ) ;

    return offset + i ;
    
  }


  /**
   * \fn void CurveBuilder::insert_min_max_constraints( const size_t eqline , const size_t lo , const size_t up , const size_t b0 , const size_t b1 , const size_t b2 )
   *
   * \brief  Inserts the  coefficients of  the equations  defining the
   * three  min-max constraints  into the  matrix associated  with the
   * linear  program  (LP),  and  sets  the  right-hand  side  of  the
   * constraints as well.
   *
   * \param eqline A reference to the counter of equations.
   * \param  lo  Column  index  of   the  lower  bound  for  a  second
   * difference.
   * \param  up  Column  index  of   the  upper  bound  for  a  second
   * difference.
   * \param b0  Column index of  the first control value  defining the
   * second difference.
   * \param b1 Column  index of the second control  value defining the
   * second difference.
   * \param b2  Column index of  the third control value  defining the
   * second difference.
   *
   */
  void
  CurveBuilder::insert_min_max_constraints(
                                           const size_t eqline ,
                                           const size_t lo ,
                                           const size_t up ,
                                           const size_t b0 ,
                                           const size_t b1 ,
                                           const size_t b2
                                          )
  {
    // First  min-max  constraint:  the  upper  bound  of  the  second
    // difference must  be greater than or  equal to the value  of the
    // second difference:
          
    const double onesixth = double( 1 ) / double( 6 ) ;

    insert_coefficient( eqline , up ,             1 ) ;
    insert_coefficient( eqline , b0 , -1 * onesixth ) ;
    insert_coefficient( eqline , b1 ,  2 * onesixth ) ;
    insert_coefficient( eqline , b2 , -1 * onesixth ) ;
    
    insert_bound( eqline , Bound::GTE , 0 ) ;
    
    // Second  min-max  constraint:  the  upper bound  on  the  second
    // difference must be greater than or equal to zero (i.e., must be
    // non-negative).

    insert_coefficient( eqline + 1 , up , 1 ) ;    

    insert_bound( eqline + 1 , Bound::GTE , 0 ) ;
    
    // Third min-max constraint: the sum of the upper and lower bounds
    // of the second difference must be  equal to the value the second
    // difference.
    
    insert_coefficient( eqline + 2 , up ,             1 ) ;
    insert_coefficient( eqline + 2 , lo ,             1 ) ;
    insert_coefficient( eqline + 2 , b0 , -1 * onesixth ) ;
    insert_coefficient( eqline + 2 , b1 ,  2 * onesixth ) ;
    insert_coefficient( eqline + 2 , b2 , -1 * onesixth ) ;
    
    insert_bound( eqline + 2 , Bound::EQT , 0 ) ;
 
    return ;
  }


  /**
   * \fn void CurveBuilder::insert_extreme_point_correspondence_constraint( const size_t eqline , const std::vector< size_t >& col , const std::vector< double >& val , const double rhs  )
   *
   * \brief  Inserts   into  the   linear  program  (LP)   matrix  the
   * coefficients of the  unknowns and the right-hand side  value of a
   * constraint corresponding to the location of the starting or final
   * point of the b-spline curve.
   *
   * \param eqline A reference to the counter of equations. 
   * \param  col   An  array  with   the  LP  matrix   column  indices
   * corresponding to the unknowns of the correspondence constraint.
   * \param val An array with the values corresponding to the unknowns
   * of the correspondence constraint.
   * \param rhs The right-hand side value of the constraint.
   *
   */
  void
  CurveBuilder::insert_extreme_point_correspondence_constraint(
                                                               const size_t eqline ,
                                                               const std::vector< size_t >& col ,
                                                               const std::vector< double >& val ,
                                                               const double rhs
                                                              )
  {
    for ( size_t i = 0 ; i < 4 ; i++ ) {
      insert_coefficient( eqline , col[ i ] , val[ i ] ) ;
    }
    
    insert_bound( eqline , Bound::EQT , rhs ) ;

    return ;
  }


  /**
   * \fn void CurveBuilder::insert_periodic_correspondence_constraints( const size_t eqline , const std::vector< size_t >& strx , const std::vector< size_t >& stry , const std::vector< size_t >& endx , const std::vector< size_t >& endy )
   *
   * \brief  Inserts   into  the   linear  program  (LP)   matrix  the
   * coefficients of  the unknowns and  the right-hand side  values of
   * the constraints that  ensure that the first  three control points
   * are the same as the last three control points (in this order).
   *
   * \param eqline A reference to the counter of equations. 
   * \param strx  An array with  the column  indices of the  LP matrix
   * corresponding  to the  first Cartesian  coordinates of  the first
   * three control points.
   * \param stry  An array with  the column  indices of the  LP matrix
   * corresponding  to the second Cartesian  coordinates of  the first
   * three control points.
   * \param endx  An array with  the column  indices of the  LP matrix
   * corresponding  to the  first  Cartesian coordinates  of the  last
   * three control points.
   * \param endy  An array with  the column  indices of the  LP matrix
   * corresponding  to the  second Cartesian  coordinates of  the last
   * three control points.
   *
   */
  void
  CurveBuilder::insert_periodic_correspondence_constraints(
                                                           const size_t eqline ,
                                                           const std::vector< size_t >& strx ,
                                                           const std::vector< size_t >& stry ,
                                                           const std::vector< size_t >& endx ,
                                                           const std::vector< size_t >& endy
                                                          )
  {
    for ( size_t j = 0 ; j < 3 ; j++ ) {
      insert_coefficient( eqline + 2 * j , strx[ j ] ,  1 ) ;
      insert_coefficient( eqline + 2 * j , endx[ j ] , -1 ) ;

      insert_bound( eqline + 2 * j , Bound::EQT , 0 ) ;

      insert_coefficient( eqline + 2 * j + 1 , stry[ j ] ,  1 ) ;
      insert_coefficient( eqline + 2 * j + 1 , endy[ j ] , -1 ) ;

      insert_bound( eqline + 2 * j + 1 , Bound::EQT , 0 ) ;
    }

    return ;
  }

  
  /**
   * \fn void CurveBuilder::insert_nonlinear_terms_of_epiece_point_lower_bound( const size_t eqline , const double s , const size_t c , const std::vector< std::vector< std::vector< size_t > > >& sd , const std::vector< std::vector< double > >& nl , const std::vector< std::vector< double > >& nu )
   *
   * \brief Inserts the coefficients of the second difference terms of
   * the equation  defining lower bounds  for the e-piece  points into
   * the  matrix associated  with an  instance of  the Linear  Program
   * (LP).  The terms belong to the constraint that forces the e-piece
   * points to be inside a certain c-section of the channel.
   *
   * \param eqline A counter for the number of constraints.
   * \param s A parameter value identifying a point on the e-piece.
   * \param c An index identifying a c-segment of the channel.
   * \param sd Array  with the LP matrix  column indices corresponding
   * to the lower and upper  bounds on second differences occurring in
   * the  equation  defining  the  e-piece  points  belonging  to  the
   * c-segment.
   * \param nl  Array of Cartesian  coordinates of outward  normals to
   * the  supporting  lines of  the  lower  envelope segments  of  the
   * channel.
   * \param nu  Array of Cartesian  coordinates of outward  normals to
   * the  supporting  lines of  the  upper  envelope segments  of  the
   * channel.
   * 
   */
  void
  CurveBuilder::insert_nonlinear_terms_of_epiece_point_lower_bound(
                                                                   const size_t eqline ,
                                                                   const double s ,
                                                                   const size_t c ,
                                                                   const std::vector< std::vector< std::vector< size_t > > >& sd ,
                                                                   const std::vector< std::vector< double > >& nl ,
                                                                   const std::vector< std::vector< double > >& nu
                                                                  )
  {
    // Insert into the matrix associated  with the linear program (LP)
    // the  coefficients  of the  second  differences  of the  e-piece
    // breakpoint lower bound \f$\stackrel{e}{\sim}^p\f$ in constraint
    // (3a).

    //
    // The computation is performed for each second difference j.
    //
 
    for ( size_t j = 1 ; j < 3 ; j++ ) {
      //
      // Get lower and upper bounds for the special polynomial.
      //
      double dl ;
      double du ;
      evaluate_bounding_polynomial(
                                   j  ,
                                   s  ,
                                   du ,   // switch lower and upper bounds.
                                   dl     // switch lower and upper bounds.
                                  ) ;

      //
      // The coefficients are the same for each Cartesian coordinate.
      //

      for ( size_t v = 0 ; v < 2 ; v++ ) {
        // Point \f$\stackrel{e}{\sim}^p( s )  \f$ of the e-piece must
        // be above  the lower envelope  of the c-th c-segment  of the
        // channel.
        insert_csegment_constraint(
                                   eqline ,
                                   dl ,
                                   du ,
                                   sd[ j - 1 ][ 0 ][ v ] ,
                                   sd[ j - 1 ][ 1 ][ v ] ,
                                   nl[ c ][ v ] 
                                  ) ;

        // Point \f$\stackrel{e}{\sim}^p( s )  \f$ of the e-piece must
        // be below  the upper envelope  of the c-th c-segment  of the
        // channel.
        insert_csegment_constraint(
                                   eqline + 1 ,
                                   dl ,
                                   du ,
                                   sd[ j - 1 ][ 0 ][ v ] ,
                                   sd[ j - 1 ][ 1 ][ v ] ,
                                   nu[ c ][ v ] 
                                  ) ;
      }
    }

    return ;
  }


  /**
   * \fn void CurveBuilder::insert_nonlinear_terms_of_epiece_point_lower_bound( const size_t eqline , const double s , const size_t c , const std::vector< std::vector< std::vector< size_t > > >& sd , const std::vector< std::vector< double > >& ncsec )
   *
   * \brief Inserts the coefficients of the second difference terms of
   * the equation  defining lower bounds  for the e-piece  points into
   * the  matrix associated  with an  instance of  the Linear  Program
   * (LP).  The terms belong to the constraint that forces one e-piece
   * point to be on the right or left side of a channel c-section.
   *
   * \param eqline A counter for the number of constraints.
   * \param s A parameter value identifying a point on the e-piece.
   * \param c An index identifying a c-segment of the channel.
   * \param sd Array  with the LP matrix  column indices corresponding
   * to the lower and upper  bounds on second differences occurring in
   * the  equation  defining  the  e-piece  points  belonging  to  the
   * c-segment.
   * \param ncsec Array of  Cartesian coordinates of normals (pointing
   * to the  left) to the  supporting lines  of the c-sections  of the
   * channel.
   * 
   */
  void
  CurveBuilder::insert_nonlinear_terms_of_epiece_point_lower_bound(
								   const size_t eqline ,
								   const double s ,
								   const size_t c ,
								   const std::vector< std::vector< std::vector< size_t > > >& sd ,
								   const std::vector< std::vector< double > >& ncsec
								  )
  {
    // Insert into the matrix associated  with the linear program (LP)
    // the  coefficients  of the  second  differences  of the  e-piece
    // breakpoint lower bound \f$\stackrel{e}{\sim}^p\f$ in constraint
    // (3c).

    //
    // The computation is performed for each second difference j.
    //
 
    for ( size_t j = 1 ; j < 3 ; j++ ) {
      //
      // Get lower and upper bounds for the special polynomial.
      //
      double dl ;
      double du ;
      evaluate_bounding_polynomial(
                                   j  ,
                                   s  ,
                                   du ,   // switch lower and upper bounds.
                                   dl     // switch lower and upper bounds.
                                  ) ;

      //
      // The coefficients are the same for each Cartesian coordinate.
      //

      for ( size_t v = 0 ; v < 2 ; v++ ) {
        // Point \f$\stackrel{e}{\sim}^p( s )  \f$ of the e-piece must
        // be either  on the  right side  or on the  left side  of the
        // channel c-section.
        insert_csegment_constraint(
                                   eqline ,
                                   dl ,
                                   du ,
                                   sd[ j - 1 ][ 0 ][ v ] ,
                                   sd[ j - 1 ][ 1 ][ v ] ,
                                   ncsec[ c ][ v ] 
                                  ) ;
      }
    }

    return ;
  }


  /**
   * \fn void CurveBuilder::insert_nonlinear_terms_of_epiece_point_upper_bound( const size_t eqline , const double s , const size_t c , const std::vector< std::vector< std::vector< size_t > > >& sd , const std::vector< std::vector< double > >& nl , const std::vector< std::vector< double > >& nu )
   *
   * \brief Inserts the coefficients of the second difference terms of
   * the equation  defining lower bounds  for the e-piece  points into
   * the  matrix associated  with an  instance of  the Linear  Program
   * (LP).  The terms belong to the constraint that forces the e-piece
   * points to be inside a certain c-section of the channel.
   *
   * \param eqline A counter for the number of constraints.
   * \param s A parameter value identifying a point on the e-piece.
   * \param c An index identifying a c-segment of the channel.
   * \param sd Array  with the LP matrix  column indices corresponding
   * to the lower and upper  bounds on second differences occurring in
   * the  equation  defining  the  e-piece  points  belonging  to  the
   * c-segment.
   * \param nl  Array of Cartesian  coordinates of outward  normals to
   * the  supporting  lines of  the  lower  envelope segments  of  the
   * channel.
   * \param nu  Array of Cartesian  coordinates of outward  normals to
   * the  supporting  lines of  the  upper  envelope segments  of  the
   * channel.
   * 
   */
  void
  CurveBuilder::insert_nonlinear_terms_of_epiece_point_upper_bound(
                                                                   const size_t eqline ,
                                                                   const double s ,
                                                                   const size_t c ,
                                                                   const std::vector< std::vector< std::vector< size_t > > >& sd ,
                                                                   const std::vector< std::vector< double > >& nl ,
                                                                   const std::vector< std::vector< double > >& nu
                                                                  )
  {
    // Insert into the matrix associated  with the linear program (LP)
    // the  coefficients  of the  second  differences  of the  e-piece
    // breakpoint lower bound \f$\stackrel{\sim}{e}^p\f$ in constraint
    // (3a).

    //
    // The computation is performed for each second difference j.
    //
 
    for ( size_t j = 1 ; j < 3 ; j++ ) {
      //
      // Get lower and upper bounds for the special polynomial.
      //
      double dl ;
      double du ;
      evaluate_bounding_polynomial(
                                   j  ,
                                   s  ,
                                   dl ,    // DON't switch lower and upper bounds.
                                   du      // DON't switch lower and upper bounds.
                                  ) ;

      //
      // The coefficients are the same for each Cartesian coordinate.
      //

      for ( size_t v = 0 ; v < 2 ; v++ ) {
        // Point \f$\stackrel{e}{\sim}^p( s )  \f$ of the e-piece must
        // be above  the lower envelope  of the c-th c-segment  of the
        // channel.
        insert_csegment_constraint(
                                   eqline ,
                                   dl ,
                                   du ,
                                   sd[ j - 1 ][ 0 ][ v ] ,
                                   sd[ j - 1 ][ 1 ][ v ] ,
                                   nl[ c ][ v ]
                                  ) ;
        
        // Point \f$\stackrel{e}{\sim}^p( s )  \f$ of the e-piece must
        // be below  the upper envelope  of the c-th c-segment  of the
        // channel.
        insert_csegment_constraint(
                                   eqline + 1 ,
                                   dl ,
                                   du ,
                                   sd[ j - 1 ][ 0 ][ v ] ,
                                   sd[ j - 1 ][ 1 ][ v ] ,
                                   nu[ c ][ v ]
                                  ) ;
      }
    }

    return ;
  }


  /**
   * \fn void CurveBuilder::insert_nonlinear_terms_of_epiece_point_upper_bound( const size_t eqline , const double s , const size_t c , const std::vector< std::vector< std::vector< size_t > > >& sd , const std::vector< std::vector< double > >& ncsec )
   *
   * \brief Inserts the coefficients of the second difference terms of
   * the equation  defining upper bounds  for the e-piece  points into
   * the  matrix associated  with an  instance of  the Linear  Program
   * (LP).  The terms belong to the constraint that forces one e-piece
   * point to be on the right or left side of a channel c-section.
   *
   * \param eqline A counter for the number of constraints.
   * \param s A parameter value identifying a point on the e-piece.
   * \param c An index identifying a c-segment of the channel.
   * \param sd Array  with the LP matrix  column indices corresponding
   * to the lower and upper  bounds on second differences occurring in
   * the  equation  defining  the  e-piece  points  belonging  to  the
   * c-segment.
   * \param ncsec Array of  Cartesian coordinates of normals (pointing
   * to the  left) to the  supporting lines  of the c-sections  of the
   * channel.
   * 
   */
  void
  CurveBuilder::insert_nonlinear_terms_of_epiece_point_upper_bound(
								   const size_t eqline ,
								   const double s ,
								   const size_t c ,
								   const std::vector< std::vector< std::vector< size_t > > >& sd ,
								   const std::vector< std::vector< double > >& ncsec
								  )
  {
    // Insert into the matrix associated  with the linear program (LP)
    // the  coefficients  of the  second  differences  of the  e-piece
    // breakpoint lower bound \f$\stackrel{\sim}{e}^p\f$ in constraint
    // (3c).

    //
    // The computation is performed for each second difference j.
    //
 
    for ( size_t j = 1 ; j < 3 ; j++ ) {
      //
      // Get lower and upper bounds for the special polynomial.
      //
      double dl ;
      double du ;
      evaluate_bounding_polynomial(
                                   j  ,
                                   s  ,
                                   dl ,   // DON't switch lower and upper bounds.
                                   du     // DON't switch lower and upper bounds.
                                  ) ;

      //
      // The coefficients are the same for each Cartesian coordinate.
      //

      for ( size_t v = 0 ; v < 2 ; v++ ) {
        // Point \f$\stackrel{\sim}{e}^p( s )  \f$ of the e-piece must
        // be either  on the  right side  or on the  left side  of the
        // channel c-section.
        insert_csegment_constraint(
                                   eqline ,
                                   dl ,
                                   du ,
                                   sd[ j - 1 ][ 0 ][ v ] ,
                                   sd[ j - 1 ][ 1 ][ v ] ,
                                   ncsec[ c ][ v ] 
                                  ) ;
      }
    }

    return ;
  }


  /**
   * \fn void CurveBuilder::insert_linear_terms_of_epiece_point_bounds( const size_t eqline , const double s , const double t , const size_t p , const size_t c , const std::vector< std::vector< size_t > >& cp , const std::vector< std::vector< double > >& nl , const std::vector< std::vector< double > >& nu )
   *
   * \brief  Inserts  the coefficients  of  the  linear terms  of  the
   * equation defining lower  and upper bounds for  the e-piece points
   * into the matrix  associated with the Linear  Program (LP).  These
   * terms occur in  the constraint that enforces an  e-piece point to
   * stay inside channel.
   *
   * \param eqline A counter for the number of constraints.
   * \param s A parameter value identifying a point on the e-piece.
   * \param t  A parameter value  identifying the b-spline  point that
   * corresponds to the point on the e-piece at parameter \e s.
   * \param p  Index of the  b-spline segment containing  the b-spline
   * point at parameter \e t.
   * \param c  An index  identifying the  c-segment the  e-piece point
   * belongs to.
   * \param cp Array  with the LP matrix  column indices corresponding
   * to the  control points  of the  b-spline defining  the \f$p\f$-th
   * piece of the curve.
   * \param nl  Array of Cartesian  coordinates of outward  normals to
   * the  supporting  lines of  the  lower  envelope segments  of  the
   * channel.
   * \param nu  Array of Cartesian  coordinates of outward  normals to
   * the  supporting  lines of  the  upper  envelope segments  of  the
   * channel.
   * 
   */
  void
  CurveBuilder::insert_linear_terms_of_epiece_point_bounds(
                                                           const size_t eqline ,
                                                           const double s ,
                                                           const double t ,
                                                           const size_t p ,
                                                           const size_t c ,
                                                           const std::vector< std::vector< size_t > >& cp ,
                                                           const std::vector< std::vector< double > >& nl ,
                                                           const std::vector< std::vector< double > >& nu
                                                          )
  {
    //
    // The coefficients are the same for each Cartesian coordinate.
    //
    const double onesixth = double( 1 ) / double( 6 ) ;
    
    // The upper and  lower bounds on the e-piece points  must be on
    // or  above the  lower envelope  of the  c-th c-segment  of the
    // channel.
    
    const double c0 = onesixth * ( 1 - s ) ;
    const double c1 = ( ( -2 + 3 * s ) * onesixth ) + ( p + 4 - t ) ;
    const double c2 = ( (  1 - 3 * s ) * onesixth ) + ( t - p - 3 ) ;
    const double c3 = onesixth * s ;
    
    for ( size_t v = 0 ; v < 2 ; v++ ) {
      //
      // Compute constraints for the v-th Cartesian coordinate.
      //
      insert_csegment_constraint(
                                 eqline ,
                                 c0 ,
                                 c1 ,
                                 c2 ,
                                 c3 ,
                                 cp[ 0 ][ v ] ,
                                 cp[ 1 ][ v ] ,
                                 cp[ 2 ][ v ] ,
                                 cp[ 3 ][ v ] ,
                                 nl[ c ][ v ]
                                ) ;
      
      // The upper and  lower bounds on the e-piece points  must be on
      // or  below the  upper envelope  of the  c-th c-segment  of the
      // channel.
      insert_csegment_constraint(
                                 eqline + 1 ,
                                 c0 ,
                                 c1 ,
                                 c2 ,
                                 c3 ,
                                 cp[ 0 ][ v ] ,
                                 cp[ 1 ][ v ] ,
                                 cp[ 2 ][ v ] ,
                                 cp[ 3 ][ v ] ,
                                 nu[ c ][ v ]
                                ) ;
    }

    return ;
  }


  /**
   * \fn void CurveBuilder::insert_linear_terms_of_epiece_point_bounds( const size_t eqline , const double s , const double t , const size_t , const size_t c , const std::vector< std::vector< size_t > >& cp , const std::vector< std::vector< double > >& ncsec )
   *
   * \brief  Inserts  the coefficients  of  the  linear terms  of  the
   * equation defining lower  and upper bounds for  the e-piece points
   * into the matrix  associated with the Linear  Program (LP).  These
   * terms occur in  the constraint that enforces an  e-piece point to
   * stay either on the right or on left side of a channel c-section.
   *
   * \param eqline A counter for the number of constraints.
   * \param s A parameter value identifying a point on the e-piece.
   * \param t  A parameter value  identifying the b-spline  point that
   * corresponds to the point on the e-piece at parameter \e s.
   * \param p  Index of the  b-spline segment containing  the b-spline
   * point at parameter \e t.
   * \param c  An index  identifying the  c-segment the  e-piece point
   * belongs to.
   * \param cp Array  with the LP matrix  column indices corresponding
   * to the  control points  of the  b-spline defining  the \f$p\f$-th
   * piece of the curve.
   * \param ncsec Array of  Cartesian coordinates of normals (pointing
   * to the  left) to the  supporting lines  of the c-sections  of the
   * channel.
   *
   */
  void
  CurveBuilder::insert_linear_terms_of_epiece_point_bounds(
                                                           const size_t eqline ,
                                                           const double s ,
                                                           const double t ,
                                                           const size_t p ,
                                                           const size_t c ,
                                                           const std::vector< std::vector< size_t > >& cp ,
                                                           const std::vector< std::vector< double > >& ncsec
                                                          )
  {
    //
    // The coefficients are the same for each Cartesian coordinate.
    //
    const double onesixth = double( 1 ) / double( 6 ) ;
    
    // The upper and lower bounds on  the e-piece point must be either
    // on the left or on the right side of a c-section of the channel.
    
    const double c0 = onesixth * ( 1 - s ) ;
    const double c1 = ( ( -2 + 3 * s ) * onesixth ) + ( p + 4 - t ) ;
    const double c2 = ( (  1 - 3 * s ) * onesixth ) + ( t - p - 3 ) ;
    const double c3 = onesixth * s ;
    
    for ( size_t v = 0 ; v < 2 ; v++ ) {
      //
      // Compute constraints for the v-th Cartesian coordinate.
      //

      //
      // Lower bound --> Equation eqline
      //
      insert_csegment_constraint(
                                 eqline ,
                                 c0 ,
                                 c1 ,
                                 c2 ,
                                 c3 ,
                                 cp[ 0 ][ v ] ,
                                 cp[ 1 ][ v ] ,
                                 cp[ 2 ][ v ] ,
                                 cp[ 3 ][ v ] ,
                                 ncsec[ c ][ v ]
                                ) ;

      //
      // Upper bound --> Equation eqline + 1
      //
      insert_csegment_constraint(
                                 eqline + 1 ,
                                 c0 ,
                                 c1 ,
                                 c2 ,
                                 c3 ,
                                 cp[ 0 ][ v ] ,
                                 cp[ 1 ][ v ] ,
                                 cp[ 2 ][ v ] ,
                                 cp[ 3 ][ v ] ,
                                 ncsec[ c ][ v ]
                                ) ;
    }

    return ;
  }


  /**
   * \fn void CurveBuilder::insert_rhs_of_sleeve_corners_in_channel_constraints( const size_t eqline , const size_t c , const std::vector< std::vector< double > >& nl , const std::vector< std::vector< double > >& nu )
   *
   * \brief Inserts into the matrix associated with the Linear Program
   * (LP) the right-hand side values of the constraints that enforce a
   * sleeve point to stay inside a  c-segment of the channel. The type
   * of each constraint (equality or inequality: ==, >= or <=) is also
   * set here.
   *
   * \param eqline A counter for the number of constraints.
   * \param c  An index  identifying the  c-segment the  e-piece point
   * belongs to.
   * \param nl  Array of Cartesian  coordinates of outward  normals to
   * the  supporting  lines of  the  lower  envelope segments  of  the
   * channel.
   * \param nu  Array of Cartesian  coordinates of outward  normals to
   * the  supporting  lines of  the  upper  envelope segments  of  the
   * channel.
   * 
   */
  void
  CurveBuilder::insert_rhs_of_sleeve_corners_in_channel_constraints(
                                                                    const size_t eqline ,
                                                                    const size_t c ,
                                                                    const std::vector< std::vector< double > >& nl ,
                                                                    const std::vector< std::vector< double > >& nu
                                                                   )
  {
    insert_bound( eqline     , Bound::LTE , _lxcoords[ c ] * nl[ c ][ 0 ] + _lycoords[ c ] * nl[ c ][ 1 ] ) ;
    insert_bound( eqline + 1 , Bound::LTE , _uxcoords[ c ] * nu[ c ][ 0 ] + _uycoords[ c ] * nu[ c ][ 1 ] ) ;
    
    insert_bound( eqline + 2 , Bound::LTE , _lxcoords[ c ] * nl[ c ][ 0 ] + _lycoords[ c ] * nl[ c ][ 1 ] ) ;
    insert_bound( eqline + 3 , Bound::LTE , _uxcoords[ c ] * nu[ c ][ 0 ] + _uycoords[ c ] * nu[ c ][ 1 ] ) ;
    
    return ;
  }

  
  /**
   * \fn void CurveBuilder::insert_rhs_of_sleeve_inside_csegment_constraints( const size_t eqline , const size_t c , const std::vector< std::vector< double > >& ncsec )
   *
   * \brief Inserts into the matrix associated with the Linear Program
   * (LP) the right-hand  side values of the  constraints that enforce
   * one e-piece breakpoint  to stay on the right side  of a c-section
   * of the  channel, and  another e-piece breakpoint  to stay  on the
   * left side of the same c-section.
   *
   * \param eqline A counter for the number of constraints.
   * \param c  An index identifying  the c-segment the  e-piece points
   * belongs to.
   * \param ncsec Array of  Cartesian coordinates of normals (pointing
   * to the  left) to the  supporting lines  of the c-sections  of the
   * channel.
   * 
   */
  void
  CurveBuilder::insert_rhs_of_sleeve_inside_csegment_constraints(
								 const size_t eqline ,
								 const size_t c ,
								 const std::vector< std::vector< double > >& ncsec
								)
  {
    insert_bound( eqline     , Bound::LTE , _lxcoords[ c ] * ncsec[ c ][ 0 ] + _lycoords[ c ] * ncsec[ c ][ 1 ] ) ;
    insert_bound( eqline + 1 , Bound::LTE , _uxcoords[ c ] * ncsec[ c ][ 0 ] + _uycoords[ c ] * ncsec[ c ][ 1 ] ) ;
    
    insert_bound( eqline + 2 , Bound::GTE , _lxcoords[ c ] * ncsec[ c ][ 0 ] + _lycoords[ c ] * ncsec[ c ][ 1 ] ) ;
    insert_bound( eqline + 3 , Bound::GTE , _uxcoords[ c ] * ncsec[ c ][ 0 ] + _uycoords[ c ] * ncsec[ c ][ 1 ] ) ;
  }


  /**
   * \fn void CurveBuilder::evaluate_bounding_polynomial( const size_t j , const double t , double& lower , double& upper )
   *
   * \brief Obtains a lower bound and  an upper bound for the value of
   * a precomputed, bounding polynomial at a given parameter value.
   *
   * \param j An index for the precomputed, bounding polynomial.
   * \param t A parameter value.
   * \param lower A reference to the lower bound.
   * \param upper A reference to the upper bound.
   *
   */
  void
  CurveBuilder::evaluate_bounding_polynomial(
                                             const size_t j ,
                                             const double t ,
                                             double& lower ,
                                             double& upper
                                            )
  {
    try {
      lower = _tf->alower( j , t ) ;
      upper = _tf->aupper( j , t ) ;
    }
    catch ( const ExceptionObject& xpt ) {
      treat_exception( xpt ) ;
      exit( EXIT_FAILURE ) ;
    }   

    return ;
  }

  
  /**
   * \fn void CurveBuilder::insert_csegment_constraint( const size_t eqline , const double lower , const double upper , const size_t sdlo , const size_t sdup , const double normal )
   *
   * \brief Inserts the coefficients of  the lower and upper bounds of
   * a constraint  second difference  term into the  matrix associated
   * with an instance of the linear program (LP).  The term belongs to
   * the equation defining the upper (or lower) bound of a point of an
   * e-piece.   The  constraint  ensures  that the  point  lies  on  a
   * specific side of the oriented suppporting line of one of the four
   * line segments delimiting a c-segment of the channel.
   *
   * \param eqline A counter for the number of constraints. 
   * \param  lower Coefficient  of the  second difference  lower bound
   * term.
   * \param  upper Coefficient  of the  second difference  upper bound
   * term.
   * \param sdlo The  index of the LP matrix  column corresponding to
   * the second difference lower bound term.
   * \param sdup The  index of the LP matrix  column corresponding to
   * the second difference upper bound term.
   * \param normal A  normal to a supporting, oriented line  of one of
   * the four  line segments  delimiting a  specific c-segment  of the
   * channel.
   *
   */
  void
  CurveBuilder::insert_csegment_constraint(
                                           const size_t eqline ,
                                           const double lower ,
                                           const double upper ,
                                           const size_t sdlo ,
                                           const size_t sdup ,
                                           const double normal
                                          )
  {
    double temp ;

    temp = lower * normal ;
    if ( temp != 0 ) {
      insert_coefficient( eqline , sdlo , temp ) ;
    }

    temp = upper * normal ;
    if ( temp != 0 ) {
      insert_coefficient( eqline , sdup , temp ) ;
    }

    return ;
  }


  /**
   * \fn void CurveBuilder::insert_csegment_constraint( const size_t eqline , const double c0 , const double c1 , const double c2 , const double c3 , const size_t b0 , const size_t b1 , const size_t b2 , const size_t b3 , const double normal )
   *
   * \brief Inserts the coefficients of  the linear terms of the upper
   * and lower  bounds of  an e-piece point  equation into  the matrix
   * associated  with an  instance  of the  linear  program (LP).  The
   * constraint ensures  that the point  of the e-piece lies  inside a
   * c-segment of the channel.
   *
   * \param eqline A counter for the number of constraints.
   * \param c0 Coefficient of the  first control point of the b-spline
   * segment  containing the  curve  point associated  to the  e-piece
   * point.
   * \param c1 Coefficient of the second control point of the b-spline
   * segment  containing the  curve  point associated  to the  e-piece
   * point.
   * \param c2 Coefficient of the  third control point of the b-spline
   * segment  containing the  curve  point associated  to the  e-piece
   * point.
   * \param c3 Coefficient of the fourth control point of the b-spline
   * segment  containing the  curve  point associated  to the  e-piece
   * point.
   * \param  b0 Index  of the  LP matrix  column corresponding  to the
   * first control point of the  b-spline segment containing the curve
   * point associated to the e-piece point.
   * \param  b1 Index  of the  LP matrix  column corresponding  to the
   * second control point of the b-spline segment containing the curve
   * point associated to the e-piece point.
   * \param  b2 Index  of the  LP matrix  column corresponding  to the
   * third control point of the  b-spline segment containing the curve
   * point associated to the e-piece point.
   * \param  b3 Index  of the  LP matrix  column corresponding  to the
   * fourth control point of the b-spline segment containing the curve
   * point associated to the e-piece point.
   * \param normal A  normal to a supporting, oriented line  of one of
   * the four  line segments  delimiting a  specific c-segment  of the
   * channel.
   *
   */
  void
  CurveBuilder::insert_csegment_constraint(
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
                                          )
  {
    double temp = c0 * normal ;
    if ( temp != 0 ) {
      insert_coefficient( eqline , b0 , temp ) ;
    }

    temp = c1 * normal ;
    if ( temp != 0 ) {
      insert_coefficient( eqline , b1 , temp ) ;
    }

    temp = c2 * normal ;
    if ( temp != 0 ) {
      insert_coefficient( eqline , b2 , temp ) ;
    }

    temp = c3 * normal ;
    if ( temp != 0 ) {
      insert_coefficient( eqline , b3 , temp ) ;
    }

    return ;
  }

  
  /**
   * \fn void CurveBuilder::insert_csegment_constraint( const size_t eqline , const double c0 , const double c1 , const double c2 , const double c3 , const double c4 , const size_t b0 , const size_t b1 , const size_t b2 , const size_t b3 , const size_t b4 , const double normal )
   *
   * \brief Inserts the coefficients of  the linear terms of the upper
   * and lower  bounds of  an e-piece point  equation into  the matrix
   * associated with  an instance  of the  linear program  (LP).  This
   * point  belongs  to  an  e-piece  segment  whose  endpoints  bound
   * b-spline  curve points  in  two distinct,  but consecutive  curve
   * segments.   The constraint  ensures that  the e-piece  point lies
   * inside a c-segment of the channel.
   *
   * \param eqline A counter for the number of constraints.

   * \param c0 Coefficient of the  first control point of the b-spline
   * segment  containing the  curve  point associated  with the  right
   * endpoint of the e-piece segment.
   * \param c1 Coefficient of the second control point of the b-spline
   * segment  containing the  curve  point associated  with the  right
   * endpoint of the e-piece segment.
   * \param c2 Coefficient of the  third control point of the b-spline
   * segment  containing the  curve  point associated  with the  right
   * endpoint of the e-piece segment.
   * \param c3 Coefficient of the fourth control point of the b-spline
   * segment  containing the  curve  point associated  with the  right
   * endpoint of the e-piece segment.
   * \param c4 Coefficient of the fourth control point of the b-spline
   * segment  containing  the curve  point  associated  with the  left
   * endpoint of the e-piece segment.
   * \param  b0 Index  of the  LP matrix  column corresponding  to the
   * first control point of the  b-spline segment containing the curve
   * point associated with the right endpoint of the e-piece segment.
   * \param  b1 Index  of the  LP matrix  column corresponding  to the
   * second control point of the b-spline segment containing the curve
   * point associated with the right endpoint of the e-piece segment.
   * \param  b2 Index  of the  LP matrix  column corresponding  to the
   * third control point of the  b-spline segment containing the curve
   * point associated with the right endpoint of the e-piece segment.
   * \param  b3 Index  of the  LP matrix  column corresponding  to the
   * fourth control point of the b-spline segment containing the curve
   * point associated with the right endpoint of the e-piece segment.
   * \param  b4 Index  of the  LP matrix  column corresponding  to the
   * fourth control point of the b-spline segment containing the curve
   * point associated with the left endpoint of the e-piece segment.
   * \param normal A  normal to a supporting, oriented line  of one of
   * the four  line segments  delimiting a  specific c-segment  of the
   * channel.
   *
   */
  void
  CurveBuilder::insert_csegment_constraint(
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
                                          )
  {
    insert_csegment_constraint(
                               eqline ,
                               c0 ,
                               c1 ,
                               c2 ,
                               c3 ,
                               b0 ,
                               b1 ,
                               b2 ,
                               b3 ,
                               normal
                              ) ;
    
    double temp = c4 * normal ;
    if ( temp != 0 ) {
      insert_coefficient( eqline , b4 , temp ) ;
    }
    
    return ;
  }


  /**
   * \fn int CurveBuilder::solve_lp( const size_t rows , const size_t cols )
   *
   * \brief  Solves the  linear program  corresponding to  the channel
   * problem.
   *
   * \param rows The number of constraints of the linear program.
   * \param cols The number of unknowns of the linear program.
   *
   * \return The code returned by the LP solver to indicate the status
   * of the computation of the solution of the linear program.
   */
  int
  CurveBuilder::solve_lp(
                         const size_t rows ,
                         const size_t cols
                        )
  {
    /*
     * Create the LP problem.
     */
    glp_prob* lp = glp_create_prob() ;
    
    /*
     * Set up the number of constraints and structural variables.
     */
    glp_add_rows( lp , int( rows ) ) ;
    glp_add_cols( lp , int( cols ) ) ;
    
    /*
     * Set the problem as a minimization one.
     */
    glp_set_obj_dir( lp , GLP_MIN ) ;
    
    /*
     * Set up the constraints of the problem.
     */
    set_up_lp_constraints( lp ) ;
    
    /*
     * Define bounds on the structural variables of the problem.
     */
    set_up_structural_variables( lp ) ;
    
    /*
     * Define objective function.
     */
    set_up_objective_function( lp ) ;
    
    /*
     * Set parameters of the solver.
     */
    glp_smcp param ;
    glp_init_smcp( &param ) ;
    
    param.msg_lev  = GLP_MSG_OFF ;
    param.presolve = GLP_ON ;
    
    /*
     * Call the solver.
     */
    
    int res = glp_simplex( lp , &param ) ;
    
    if ( res == 0 ) {
      /*
       * Get the solver result information.
       */
      get_lp_solver_result_information( lp ) ;
    }
    
    /*
     * Release memory held by the solver.
     */
    glp_delete_prob( lp ) ;
    
    return res ;
  }
  
  
  /**
   * \fn void CurveBuilder::set_up_lp_constraints( glp_prob* lp ) const
   *
   * \brief Assemble the matrix of  constraints of the linear program,
   * and define  the type (equality  or inequality) and bounds  on the
   * constraints.
   *
   * \param lp A pointer to the instance of the LP program.
   *
   */
  void
  CurveBuilder::set_up_lp_constraints( glp_prob* lp ) const
  {
    /*
     * Set up the bounds on the constraints of the problem.
     */
    
    for ( size_t j = 0 ; j < _bounds.size() ; j++ ) {
#ifdef DEBUGMODE
      assert( j == _bounds[ j ].get_row() ) ;
#endif
      
      int i = int( _bounds[ j ].get_row() + 1 ) ;
      
      std::stringstream ss ( std::stringstream::in | std::stringstream::out ) ;
      ss << "c" << i ;
      glp_set_row_name( lp , i , ss.str().c_str() ) ;
      
      double val = _bounds[ j ].get_value() ;
      if ( _bounds[ j ].get_type() == Bound::LTE ) {
        glp_set_row_bnds( lp , i , GLP_UP ,   0 , val ) ;
      }
      else if ( _bounds[ j ].get_type() == Bound::GTE ) {
        glp_set_row_bnds( lp , i , GLP_LO , val ,   0 ) ;
      }
      else {
        glp_set_row_bnds( lp , i , GLP_FX , val , val ) ;
      }
    }
    
    
    /*
     * Obtain the coefficients of the constraints of the problem.
     */
    
    std::vector< int    > ia ; ia.push_back( 0 ) ;  // GLPK starts indexing array \e ia at 1
    std::vector< int    > ja ; ja.push_back( 0 ) ;  // GLPK starts indexing array \e ja at 1
    std::vector< double > ar ; ar.push_back( 0 ) ;  // GLPK starts indexing array \e ar at 1
    
    int h = 0 ;
    for ( size_t j = 0 ; j < _coefficients.size() ; j++ ) {
      for ( size_t k = 0 ; k < _coefficients[ j ].size() ; k++ ) {
#ifdef DEBUGMODE
        assert( _coefficients[ j ][ k ].get_row() == j ) ;
#endif
        ia.push_back( int( _coefficients[ j ][ k ].get_row() + 1 ) ) ;
        ja.push_back( int( _coefficients[ j ][ k ].get_col() + 1 ) ) ;
        ar.push_back( _coefficients[ j ][ k ].get_value() ) ;
        ++h ;
      }
    }
    
    glp_load_matrix(
                    lp ,
                    h ,
                    &ia[ 0 ] ,
                    &ja[ 0 ] ,
                    &ar[ 0 ]
                   ) ;
    
    return ;
  }
  
  
  /**
   * \fn void CurveBuilder::set_up_structural_variables( glp_prob* lp ) const
   *
   * \brief  Define  lower  and/or  upper  bounds  on  the  structural
   * variables  of the  linear  program corresponding  to the  channel
   * problem.
   *
   * \param lp A pointer to the instance of the LP program.
   *
   */
  void
  CurveBuilder::set_up_structural_variables( glp_prob* lp ) const
  {
    //
    // Set up bounds for the first two second differences.
    //
    for ( size_t i = 1 ; i <= 2 ; i++ ) {
      for ( size_t l = 0 ; l < 2 ; l++ ) {
        for ( size_t v = 0 ; v < 2 ; v++ ) {
          size_t c = compute_second_difference_column_index(
                                                            0 ,
                                                            i ,
                                                            l ,
                                                            v
                                                           ) ;
          if ( l == 0 ) {
            std::stringstream ss ( std::stringstream::in | std::stringstream::out ) ;
            if ( v == 0 ) {
              ss << "mx" << i ;
            }
            else {
              ss << "my" << i ;
            }
            glp_set_col_name( lp , int( c ) + 1 , ss.str().c_str() ) ;
            glp_set_col_bnds( lp , int( c ) + 1 , GLP_UP , 0 , 0 ) ;
          }
          else {
            std::stringstream ss ( std::stringstream::in | std::stringstream::out ) ;
            if ( v == 0 ) {
              ss << "px" << i ;
            }
            else {
              ss << "py" << i ;
            }
            glp_set_col_name( lp , int( c ) + 1 , ss.str().c_str() ) ;
            glp_set_col_bnds( lp , int( c ) + 1 , GLP_LO , 0 , 0 ) ;
          }
        }
      }
    }
    
    //
    // Set up bounds for the remaining second differences.
    //
    for ( size_t p = 1 ; p < _np ; p++ ) {
      for ( size_t l = 0 ; l < 2 ; l++ ) {
        for ( size_t v = 0 ; v < 2 ; v++ ) {
          size_t c = compute_second_difference_column_index(
                                                            p ,
                                                            2 ,
                                                            l ,
                                                            v
                                                           ) ;
          if ( l == 0 ) {
            std::stringstream ss ( std::stringstream::in | std::stringstream::out ) ;
            if ( v == 0 ) {
              ss << "mx" << p + 2 ;
            }
            else {
              ss << "my" << p + 2 ;
            }
            glp_set_col_name( lp , int( c ) + 1 , ss.str().c_str() ) ;
            glp_set_col_bnds( lp , int( c ) + 1 , GLP_UP , 0 , 0 ) ;
          }
          else {
            std::stringstream ss ( std::stringstream::in | std::stringstream::out ) ;
            if ( v == 0 ) {
              ss << "px" << p + 2 ;
            }
            else {
              ss << "py" << p + 2 ;
            }
            glp_set_col_name( lp , int( c ) + 1 , ss.str().c_str() ) ;
            glp_set_col_bnds( lp , int( c ) + 1 , GLP_LO , 0 , 0 ) ;
          }
        }
      }
    }
    
    //
    // Set up bounds for the first four control points.
    //
    for ( size_t i = 0 ; i < 4 ; i++ ) {
      for ( size_t v = 0 ; v < 2 ; v++ ) {
        size_t c = compute_control_value_column_index(
                                                      0 ,
                                                      i ,
                                                      v
                                                     ) ;
          
        std::stringstream ss ( std::stringstream::in | std::stringstream::out ) ;
        if ( v == 0 ) {
          ss << "x" << i + 1 ;
        }
        else {
          ss << "y" << i + 1 ;
        }
        glp_set_col_name( lp , int( c ) + 1 , ss.str().c_str() ) ;
        glp_set_col_bnds( lp , int( c ) + 1 , GLP_FR , 0 , 0 ) ;
      }
    }

    for ( size_t p = 1 ; p < _np ; p++ ) {
      for ( size_t v = 0 ; v < 2 ; v++ ) {
        size_t c = compute_control_value_column_index(
                                                      p ,
                                                      3 ,
                                                      v
                                                     ) ;
          
        std::stringstream ss ( std::stringstream::in | std::stringstream::out ) ;
        if ( v == 0 ) {
          ss << "x" << p + 4 ;
        }
        else {
          ss << "y" << p + 4 ;
        }
        glp_set_col_name( lp , int( c ) + 1 , ss.str().c_str() ) ;
        glp_set_col_bnds( lp , int( c ) + 1 , GLP_FR , 0 , 0 ) ;
      }
    }

    size_t s = compute_index_of_endpoint_barycentric_coordinate( 0 ) ;
    std::stringstream ss ( std::stringstream::in | std::stringstream::out ) ;
    ss << "st" ;
    glp_set_col_name( lp , int( s ) + 1 , ss.str().c_str() ) ;
    glp_set_col_bnds( lp , int( s ) + 1 , GLP_DB , 0.40 , 0.60 ) ;
    
    if ( !_closed ) {
      size_t e = compute_index_of_endpoint_barycentric_coordinate( 1 ) ;
      std::stringstream ss2 ( std::stringstream::in | std::stringstream::out ) ;
      ss2 << "en" ;
      glp_set_col_name( lp , int( e ) + 1 , ss2.str().c_str() ) ;
      glp_set_col_bnds( lp , int( e ) + 1 , GLP_DB , 0.40 , 0.60 ) ;
    }

    for ( size_t i = 1 ; i < _nc ; i++ ) {
      size_t corner_coord = compute_index_of_corner_barycentric_coordinate( i ) ;
      std::stringstream ss2 ( std::stringstream::in | std::stringstream::out ) ;
      ss2 << "co" << i ;
      glp_set_col_name( lp , int( corner_coord ) + 1 , ss2.str().c_str() ) ;
      glp_set_col_bnds( lp , int( corner_coord ) + 1 , GLP_DB , 0.40 , 0.60 ) ;
    }

    return ;
  }
  
  /**
   * \fn void CurveBuilder::set_up_objective_function( glp_prob* lp ) const
   *
   * \brief  Define  the  objective  function of  the  linear  program
   * corresponding  to the  channel problem,  which is  a minimization
   * problem.
   *
   * \param lp A pointer to the instance of the LP program.
   *
   */
  void
  CurveBuilder::set_up_objective_function( glp_prob* lp ) const
  {
    //
    // Add the first two second difference bounds to the function.
    //
    for ( size_t i = 1 ; i < 3 ; i++ ) {
      for ( size_t l = 0 ; l < 2 ; l++ ) {
        for ( size_t v = 0 ; v < 2 ; v++ ) {
          size_t c = compute_second_difference_column_index(
                                                            0 ,
                                                            i ,
                                                            l ,
                                                            v
                                                           ) ;
            
          if ( l == 0 ) {
            glp_set_obj_coef( lp , int( c ) + 1 , -1 ) ;
          }
          else {
            glp_set_obj_coef( lp , int( c ) + 1 ,  1 ) ;
          }
        }
      }
    }
    
    //
    // Add the remaining second difference bounds to the function.
    //
    for ( size_t p = 1 ; p < _np ; p++ ) {
      for ( size_t l = 0 ; l < 2 ; l++ ) {
        for ( size_t v = 0 ; v < 2 ; v++ ) {
          size_t c = compute_second_difference_column_index(
                                                            p ,
                                                            2 ,
                                                            l ,
                                                            v
                                                           ) ;
            
          if ( l == 0 ) {
            glp_set_obj_coef( lp , int( c ) + 1 , -1 ) ;
          }
          else {
            glp_set_obj_coef( lp , int( c ) + 1 ,  1 ) ;
          }
        }
      }
    }
    
    return ;
  }
  
  
  /**
   * \fn void CurveBuilder::get_lp_solver_result_information( glp_prob* lp )
   *
   * \brief Obtain the  optimal values found by the LP  solver for the
   * structural values of the  linear programming corresponding to the
   * channel problem.
   *
   * \param lp A pointer to the instance of the LP program.
   *
   */
  void
  CurveBuilder::get_lp_solver_result_information( glp_prob* lp )
  {
    //
    // Obtain the control points of the spline curve.
    //
    for ( size_t i = 0 ; i < 4 ; i++ ) {
      for ( size_t v = 0 ; v < 2 ; v++ ) {
        size_t c = compute_control_value_column_index(
                                                      0 ,
                                                      i ,
                                                      v
                                                     ) ;
        _ctrlpts.push_back(
                           glp_get_col_prim(
                                            lp ,
                                            int( c ) + 1
                                           )
                          ) ;
      }
    }
    
    for ( size_t p = 1 ; p < _np ; p++ ) {
      for ( size_t v = 0 ; v < 2 ; v++ ) {
        size_t c = compute_control_value_column_index(
                                                      p ,
                                                      3 ,
                                                      v
                                                     ) ;
        _ctrlpts.push_back(
                           glp_get_col_prim(
                                            lp ,
                                            int( c ) + 1
                                           )
                          ) ;
      }
    }
    
    //
    // Obtain the lower and upper bounds of the second differences.
    //
    for ( size_t i = 1 ; i < 3 ; i++ ) {
      for ( size_t l = 0 ; l < 2 ; l++ ) {
        for ( size_t v = 0 ; v < 2 ; v++ ) {
          size_t c = compute_second_difference_column_index(
                                                            0 ,
                                                            i ,
                                                            l ,
                                                            v
                                                           ) ;
            
          _secdiff.push_back(
                             glp_get_col_prim(
                                              lp ,
                                              int( c ) + 1
                                             )
                            ) ;
        }
      }
    }
    
    for ( size_t p = 1 ; p < _np ; p++ ) {
      for ( size_t l = 0 ; l < 2 ; l++ ) {
        for ( size_t v = 0 ; v < 2 ; v++ ) {
          size_t c = compute_second_difference_column_index(
                                                            p ,
                                                            2 ,
                                                            l ,
                                                            v
                                                           ) ;
            
          _secdiff.push_back(
                             glp_get_col_prim(
                                              lp ,
                                              int( c ) + 1
                                             )
                            ) ;
        }
      }
    }
        
    //
    // Obtain the minimum value of the objective function.
    //
    _ofvalue = glp_get_obj_val( lp ) ;
    
    return ;
  }


}

/** @} */ //end of group class.
