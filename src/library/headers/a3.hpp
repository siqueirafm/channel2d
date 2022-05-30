#pragma once

/** 
 * \file a3.hpp
 *
 * \brief  Definition of  a  class for  representing piecewise  linear
 * enclosures of  certain cubic polynomial functions  in B&eacute;zier
 * form.
 *
 * \author
 * Marcelo Ferreira Siqueira \n
 * Universidade Federal do Rio Grande do Norte, \n
 * Departamento de Matem&aacute;tica, \n
 * mfsiqueira at mat (dot) ufrn (dot) br
 *
 * \version 1.0
 * \date March 2016
 *
 * \attention This program is distributed WITHOUT ANY WARRANTY, and it
 *            may be freely redistributed under the condition that the
 *            copyright notices  are not removed,  and no compensation
 *            is received. Private, research, and institutional use is
 *            free. Distribution of this  code as part of a commercial
 *            system  is permissible ONLY  BY DIRECT  ARRANGEMENT WITH
 *            THE AUTHOR.
 */

#include "exceptionobject.hpp"               // ExceptionObject
#include "tabulatedfunction.hpp"             // TabulatedFunction

#include <cmath>                             // sqrt
#include <cassert>                           // assert
#include <sstream>                           // std::stringstream
#include <cstdlib>                           // size_t


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
   * \class a3
   *
   * \brief   This  class   represents  two-sided,   piecewise  linear
   * enclosures  for   two  polynomial   functions  of  degree   3  in
   * B&eacute;zier form.
   *
   * \attention This class is based on the work described in
   *
   *     
   *             J. Peters and X. Wu.
   *             On  the  optimality   of  piecewise  linear  max-norm
   *             enclosures  based on  slefes. In  Proceedings of  the
   *             2002 St Malo conference on Curves and Surfaces, 2003.
   *
   */
  
  class a3 : public TabulatedFunction {
  protected:

    // ---------------------------------------------------------------
    //
    // Protected data
    //
    // ---------------------------------------------------------------

    double _l0 ;  ///< 1st control value of the lower enclosure of the polynomial \f$a_1\f$. 
    double _l1 ;  ///< 2nd control value of the lower enclosure of the polynomial \f$a_1\f$. 
    double _l2 ;  ///< 3rd control value of the lower enclosure of the polynomial \f$a_1\f$.
    double _l3 ;  ///< 4th control value of the lower enclosure of the polynomial \f$a_1\f$. 

  public:
  
    // ---------------------------------------------------------------
    //
    // Public methods
    //
    // ---------------------------------------------------------------
    
    /**
     * \fn a3( )
     *
     * \brief Creates an instance of this class.
     *
     */
    a3( ) 
    {
      const double B1 = sqrt( 57 ) ;
      const double B2 = sqrt( -10 + 2 * B1 ) ;
      const double B3 = ( double( 38 ) / 9 ) * B1 ;
      const double B4 = ( double( 1 ) / 27 ) ;

      _l0 = B4 * ( 30 - B3 ) ;
      _l1 = B4 * ( 20.0 - B3 ) ;
      _l2 = B4 * ( 25.0 + ( ( ( B1 - 9 ) * B2 ) / 2 ) - B3 ) ;

      _l3 = B4 * ( ( double( 261 ) / 8 )
          + ( ( ( B1 - 9 ) * B2 ) / 4 )
          + ( ( ( 3 * B2 - B1 ) / 8 ) * ( sqrt( 11 - 12 * B2 - 2 * B1 * ( 1 - B2 ) ) ) ) - B3 ) ;

      return ;
    }
    
    /**
     * \fn virtual double alower( const size_t i , const double u ) const
     *
     * \brief Evaluates the piecewise linear function corresponding to
     * the lower enclosure  of the \f$i\f$-th tabulated  function at a
     * point in \f$[0,1]\f$.
     *
     * \param i The index of the i-th polynomial function.
     * \param u A value in the interval \f$[0,1]\f$.
     *
     * \return   The   value   of  the   piecewise   linear   function
     * corresponding  to   the  lower  enclosure  of   the  \f$i\f$-th
     * tabulated function at a point in \f$[0,1]\f$.
     *
     */
    double alower(
                  const size_t i ,
                  const double u
                 )
      const override
    {
      if ( ( i != 1 ) && ( i != 2 ) ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "Index of the polynomial function is out of range" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }
        
      if ( ( u < 0 ) || ( u > 1 ) ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "Parameter value must belong to the interval [0,1]" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }

      return ( i == 1 ) ? a1lower( u ) : a1lower( 1 - u ) ; 
    }


    /**
     * \fn virtual double aupper( const size_t i , const double u ) const throw( ExceptionObject )
     *
     * \brief Evaluates the piecewise linear function corresponding to
     * the upper enclosure  of the \f$i\f$-th tabulated  function at a
     * point in \f$[0,1]\f$.
     *
     * \param i The index of the \f$i\f$-th polynomial function.
     * \param u A value in the interval \f$[0,1]\f$.
     *
     * \return   The   value   of  the   piecewise   linear   function
     * corresponding  to   the  upper  enclosure  of   the  \f$i\f$-th
     * tabulated function at a point in \f$[0,1]\f$.
     *
     */
    double aupper(
                  const size_t i ,
                  const double u
                 )
      const override
    {
      if ( ( i != 1 ) && ( i != 2 ) ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "Index of the polynomial function is out of range" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }

      if ( ( u < 0 ) || ( u > 1 ) ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "Parameter value must belong to the interval [0,1]" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }

      return ( i == 1 ) ? a1upper( u ) : a1upper( 1 - u ) ;
    }
    

    /**
     * \fn virtual double a( const size_t i , const double u ) const
     *
     * \brief Computes the value of the \f$i\f$-th polynomial function
     * \f$a\f$ at  a given  point of the  interval \f$[0,1]\f$  of the
     * real line.
     *
     * \param i Index of the i-th polynomial function.
     * \param u A parameter point in the interval \f$[0,1]\f$.
     *
     * \return The value of the \f$i\f$-th polynomial function \f$a\f$
     * at a  given point  \f$u\f$ of the  interval \f$[0,1]\f$  of the
     * real line.
     *
     */
    double a(
             const size_t i ,
             const double u
            )
      const override
    {
      if ( ( i != 1 ) && ( i != 2 ) ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "Index of the polynomial function is out of range" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }

      if ( ( u < 0 ) || ( u > 1 ) ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "Parameter value must belong to the interval [0,1]" ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }

      return ( i == 1 ) ? a1( u ) : a1( 1 - u ) ;
    }


    /**
     * \fn virtual unsigned degree() const
     *
     * \brief Returns the degree of tabulated functions.
     *
     * \return The degree of the tabulated functions.
     *
     */
    unsigned degree() const override
    {
      return 3 ;
    }


  protected:
    
    // ---------------------------------------------------------------
    //
    // Protected methods
    //
    // ---------------------------------------------------------------

    /**
     * \fn double a1lower( const double u ) const
     *
     * \brief  Compute the  image of  a  given point  of the  interval
     * \f$[0,1]\f$ under the lower enclosure of function \f$a_1\f$.
     *
     * \param u A point in the interval \f$[0,1]\f$.
     *
     * \return The image of a  given point of the interval \f$[0,1]\f$
     * under the lower enclosure of function \f$a_1\f$.
     *
     */
    double a1lower( const double u ) const
    {
      const double onethird = double( 1 ) / 3 ;

      double res = _l0 * h( u )
	               + _l1 * h( u -     onethird )
	               + _l2 * h( u - 2 * onethird )
	               + _l3 * h( u - 1 ) ;

      return res ;
    }


    /**
     * \fn double a1upper( const double u ) const
     *
     * \brief  Compute the  image of  a  given point  of the  interval
     * \f$[0,1]\f$ under the upper enclosure of function \f$a_1\f$.
     *
     * \param u A point in the interval \f$[0,1]\f$.
     *
     * \return The image of a  given point of the interval \f$[0,1]\f$
     * under the upper enclosure of function \f$a_1\f$.
     *
     */
    double a1upper( const double u ) const
    {
      const double onethird = double( 1 ) / 3 ;

      double res = -10 * h( u - onethird ) - 8 * h( u - 2 * onethird ) ;

      res *= ( double(1) / 27 ) ;
    
      return res ;
    }


    /**
     * \fn double a1( const double u ) const
     *
     * \brief  Computes the  value  of the  cubic polynomial  function
     * \f$a_1\f$ at a  given point of the interval  \f$[0,1]\f$ of the
     * real line.
     *
     * \param u A parameter point in the interval \f$[0,1]\f$.
     *
     * \return The value of the cubic polynomial function \f$a_1\f$ at
     * a given point of the interval \f$[0,1]\f$ of the real line.
     *
     */
    double a1( const double u ) const
    {
#ifdef DEBUGMODE
      assert( u >= 0 ) ;
      assert( u <= 1 ) ;
#endif

      return -u * ( 2 - u * ( 3 - u ) ) ;
    }


    /**
     * \fn double h( const double u ) const
     *
     * \brief Computes the value of a piecewise linear hat function at
     * a given point of the real line.
     *
     * \param u A parameter point of the real line.
     *
     * \return The value of a piecewise linear hat function at a given
     * point of the real line.
     *
     */
    double h( const double u ) const
    {
      const double onethird = 1.0 / 3.0 ;
      
      if ( u <= -onethird ) {
        return 0 ;
      }
      else if ( u <= 0 ) {
        return 3 * u + 1 ;
      }
      else if ( u <= onethird ) {
        return 1 - 3 * u ;
      }
      
      return 0 ;
    }


  } ;

}

/** @} */ //end of group class.
