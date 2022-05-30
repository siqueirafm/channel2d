#pragma once

/** 
 * \file tabulatedfunction.hpp
 *
 * \brief Definition  of an abstract class  for representing piecewise
 * linear  enclosures of  certain  polynomial  functions of  arbitrary
 * degree.
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
   * \class TabulatedFunction
   *
   * \brief   This  class   represents  two-sided,   piecewise  linear
   * enclosures of a set of \f$(d-1)\f$ polynomial functions of degree
   * \f$d\f$  in  B&eacute;zier form.   The  enclosures  must be  made
   * available  by implementating  a  pure virtual  method in  derived
   * classes.
   *
   * \attention This class is based on several papers surveyed in
   *
   *             J.   Peters.
   *             Efficient one-sided linearization of spline geometry.
   *             Proceeding  of the  10th International  Conference on
   *             Mathematics of Surfaces,  Leeds, UK, September 15-17,
   *             2003,  p.   297-319.   (Lecture  Notes   in  Computer
   *             Science,   volume  2768,   Eds.   M.J.    Wilson  and
   *             R.R. Martin).
   *
   */
  
  class TabulatedFunction {
  public:
  
    // ---------------------------------------------------------------
    //
    // Public methods
    //
    // ---------------------------------------------------------------
    
    /**
     * \fn TabulatedFunction( )
     *
     * \brief Creates an instance of this class.
     *
     */
    TabulatedFunction( ) 
    {
    }
     
    
    /**
     * \fn virtual ~TabulatedFunction( )
     *
     * \brief Releases the memory held by an instance of this class.
     *
     */
    virtual ~TabulatedFunction( )
    {
    }

    
    /**
     * \fn virtual double alower( const size_t i , const double u ) const = 0
     *
     * \brief Evaluates the piecewise linear function corresponding to
     * the lower enclosure  of the \f$i\f$-th tabulated  function at a
     * point in \f$[0,1]\f$.
     *
     * \param i The index of the \f$i\f$-th polynomial function.
     * \param u A value in the interval \f$[0,1]\f$.
     *
     * \return   The   value   of  the   piecewise   linear   function
     * corresponding  to   the  lower  enclosure  of   the  \f$i\f$-th
     * tabulated function at a point in \f$[0,1]\f$.
     *
     */
    virtual double alower( const size_t i , const double u ) const = 0 ;


    /**
     * \fn virtual double aupper( const size_t i , const double u ) const = 0
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
    virtual double aupper( const size_t i , const double u ) const = 0 ;


    /**
     * \fn virtual double a( const size_t i , const double u ) const = 0
     *
     * \brief Computes the value of the \f$i\f$-th polynomial function
     * \f$a\f$ at  a given  point of the  interval \f$[0,1]\f$  of the
     * real line.
     *
     * \param i The index of the \f$i\f$-th polynomial function.
     * \param u A parameter point in the interval \f$[0,1]\f$.
     *
     * \return The value of the \f$i\f$-th polynomial function \f$a\f$
     * at a  given point  \f$u\f$ of the  interval \f$[0,1]\f$  of the
     * real line.
     *
     */
    virtual double a( const size_t i , const double u ) const = 0 ;


    /**
     * \fn virtual unsigned degree() const = 0
     *
     * \brief Returns the degree of the tabulated functions.
     *
     * \return The degree of the tabulated functions.
     *
     */
    virtual unsigned degree() const = 0 ;
       
  } ;

}

/** @} */ //end of group class.
