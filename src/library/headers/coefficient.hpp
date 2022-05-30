#pragma once

/** 
 * \file coefficient.hpp
 *
 * \brief Definition of a class for representing a nonzero coefficient
 * of an unknown of a constraint  (inequality or equality) of a linear
 * program instance.
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
   * \class Coefficient
   *
   * \brief This class represents a  nonzero coefficient of an unknown
   * of  a constraint  (inequality or  equality) of  a linear  program
   * instance.
   *
   */
  
  class Coefficient {
  protected:

    // ---------------------------------------------------------------
    //
    // Private data
    //
    // ---------------------------------------------------------------

    size_t _row ;    ///< The identifier of the constraint this coefficient belongs to.
    size_t _col ;    ///< The identifier of the unknown multiplied by this coefficient.
    double _value ;  ///< The coefficient value.


  public:
  
    // ---------------------------------------------------------------
    //
    // Public methods
    //
    // ---------------------------------------------------------------
    
    /**
     * \fn Coefficient()
     *
     * \brief Creates an instance of this class.
     *
     */
    Coefficient()
    :
      _row( 0 ) ,
      _col( 0 ) ,
      _value( 0 )
    {
    }

    
    /**
     * \fn Coefficient( const size_t row , const size_t col , const double value )
     *
     * \brief Creates an instance of this class.
     *
     * \param row The identifier of the  constraint  this  coefficient
     * belongs to.
     * \param col  The identifier  of the  unknown multiplied  by this
     * coefficient.
     * \param value The value of the coefficient.
     *
     */
    Coefficient( const size_t row , const size_t col , double value )
    :
      _row( row ) ,
      _col( col ) ,
      _value( value )
    {
    }

    
    /**
     * \fn size_t get_row() const 
     *
     * \brief Returns the identifier of the constraint the coefficient
     * is associated  with. This identifier corresponds  to the number
     * of  a row  in the  constraint  coefficient matrix  of a  linear
     * program.
     *
     *
     * \return The  identifier of  the constraint this  coefficient is
     * associated with.
     *
     */
    size_t get_row() const
    {
      return _row ;
    }


    /**
     * \fn size_t get_col() const 
     *
     * \brief Returns the identifier of the unknown multiplied by this
     * coefficient in a constraint of a linear program instance.  This
     * identifier  corresponds  to  the  number of  a  column  in  the
     * coefficient matrix of the linear program instance.
     *
     *
     * \return  The  identifier  of  the unknown  multiplied  by  this
     * coefficient in a constraint of a linear program instance.
     *
     */
    size_t get_col() const
    {
      return _col ;
    }


    /**
     * \fn double get_value() const 
     *
     * \brief Returns the value of this coefficient.
     *
     * \return The value of this coefficient.
     *
     */
    double get_value() const
    {
      return _value ;
    }
    
  } ;

}

/** @} */ //end of group class.
