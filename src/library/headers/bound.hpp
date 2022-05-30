#pragma once

/**
 * \file bound.hpp
 *
 * \brief Definition of a class for  representing the type of a linear
 * constraint (i.e., equality or  inequality) and its right-hand side:
 * a real number.
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

#include <cstdlib>         // size_t


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
   * \class Bound
   *
   * \brief  This class  represents the  type of  a constraint  (i.e.,
   * equality or inequality)  and the value of its  right-hand side: a
   * real number.
   *
   */
  
  class Bound {
  public:

    // -----------------------------------------------------------------
    //
    // Type definitions
    //
    // -----------------------------------------------------------------

    /**
     * \typedef CONSTRAINTYPE
     *
     * \brief Defines a type for the type of a constraint.
     *
     */
    typedef enum { EQT , LTE , GTE } CONSTRAINTYPE ;

    
  protected:

    // ---------------------------------------------------------------
    //
    // Private data
    //
    // ---------------------------------------------------------------

    CONSTRAINTYPE _ctype ; ///< The type of the constraint associated with this bound.
    double _value ;        ///< The bound value.
    size_t _row ;          ///< The identifier of the constraint associated with this bound.
    

  public:
  
    // ---------------------------------------------------------------
    //
    // Public methods
    //
    // ---------------------------------------------------------------
    
    /**
     * \fn Bound()
     *
     * \brief Creates an instance of this class.
     *
     */
    Bound()
    :
      _ctype( EQT ) ,
      _value( 0 ) ,
      _row( 0 )
    {
    }

    
    /**
     * \fn Bound( const CONSTRAINTYPE type , const double value , const size_t row ) 
     *
     * \brief Creates an instance of this class.
     *
     * \param type  The type  of the  constraint associated  with this
     * bound.
     * \param value The value of the bound.
     * \param  row The  identifier of  the constraint  associated with
     * this bound.
     *
     */
    Bound( const CONSTRAINTYPE type , const double value , const size_t row )
    :
      _ctype( type ) ,
      _value( value ) ,
      _row( row )
    {
    }
    
    /**
     * \fn CONSTRAINTYPE get_type() const 
     *
     * \brief Returns the type of  the constraint associated with this
     * bound.
     *
     *
     * \return The type of the constraint associated with this bound.
     *
     */
    CONSTRAINTYPE get_type() const
    {
      return _ctype ;
    }


    /**
     * \fn double get_value() const 
     *
     * \brief Returns the value of this bound.
     *
     * \return The value of this bound.
     *
     */
    double get_value() const
    {
      return _value ;
    }


    /**
     * \fn size_t get_row() const 
     *
     * \brief Returns the identifier of the constraint associated with
     * this bound. This identifier corresponds  to the number of a row
     * in the coefficient  matrix associated with of  a linear program
     * instance.
     *
     *
     * \return The  identifier of the constraint  associated with this
     * bound.
     *
     */
    size_t get_row() const
    {
      return _row ;
    }

  } ;

}

/** @} */ //end of group class.
