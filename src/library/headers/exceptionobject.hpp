#pragma once

/** 
 * \file exceptionobject.hpp
 *
 * \brief Definition of a class for handling exceptions.
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

#include <string>                   // std::string
#include <stdexcept>                // std::exception
#include <iostream>                 // std::cerr, std::endl


/**
 * \def treat_exception( e )
 *
 * \brief  Prints out  the description  of  the error  that caused  an
 * exception as well as the file containing the instruction that threw
 * the exception and the line of the instruction in the file.
 *
 * \param e An exception.
 *
 */
#define treat_exception( e ) \
  std::cerr << std::endl \
            << "Exception: " << e.get_description() << std::endl \
            << "File: "      << e.get_file()        << std::endl \
            << "Line: "      << e.get_line()        << std::endl \
            << std::endl ;


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
   * \class ExceptionObject
   *
   * \brief This class extends class  \e exception of STL and provides
   * us with a customized way of handling exceptions and showing error
   * messages.
   *
   */
  
  class ExceptionObject : public std::exception {
  protected:

    // ---------------------------------------------------------------
    //
    // Protected data
    //
    // ---------------------------------------------------------------

    std::string  _location ;     ///< Location of the error in the line that caused the exception.
    std::string  _description ;  ///< Description of the error.
    std::string  _file ;         ///< File where the error occured.
    unsigned _line ;             ///< Line of the file where the error occurred.

  public:
  
    // ---------------------------------------------------------------
    //
    // Public methods
    //
    // ---------------------------------------------------------------
    
    
    /**
     * \fn ExceptionObject()
     *
     * \brief Creates an instance of this class.
     *
     */
    ExceptionObject()
    :
      _location( "Unknown" ) ,
      _description( "Unknown" ) ,
      _file( "Unknown" ) ,
      _line( 0 )
    {
    }

    
    /**
     * \fn ExceptionObject( const char* file , unsigned ln )
     *
     * \brief Creates an instance of this class.
     *
     * \param  file A  pointer  to  the name  of  the  file where  the
     * exception occurred.
     * \param ln  Number of the  line containing the  instruction that
     * caused the exception.
     *
     */
    ExceptionObject( const char* file , unsigned ln )
    :
      _location( "Unknown" ) ,
      _description( "Unknown" ) ,
      _file( file ) ,
      _line( ln )
    {
    }

    
    /**
     * \fn ExceptionObject( const char* file , unsigned int ln , const char* desc )
     *
     * \brief Creates an instance of this class.
     *
     * \param  file A  pointer  to  the name  of  the  file where  the
     * exception occurred.
     * \param ln  Number of the  line containing the  instruction that
     * caused the exception.
     * \param desc A pointer to a description of the error that caused
     * the exception.
     *
     */
    ExceptionObject( const char* file , unsigned int ln , const char* desc )
    :
      _location( "Unknown" ) ,
      _description( desc ) ,
      _file( file ) ,
      _line( ln )
    {
    }

    
    /**
     * \fn ExceptionObject( const char* file , unsigned ln , const char* desc , const char* loc  )
     *
     * \brief Creates an instance of this class.
     *
     * \param  file A  pointer  to  the name  of  the  file where  the
     * exception occurred.
     * \param ln  Number of the  line containing the  instruction that
     * caused the exception.
     * \param desc A pointer to a description of the error that caused
     * the exception.
     * \param loc  A pointer to  the location of the  exception inside
     * the line where it occurred.
     *
     */
    ExceptionObject( const char* file , unsigned ln , const char* desc , const char* loc )
    :
      _location( loc ) ,
      _description( desc ) ,
      _file( file ) ,
      _line( ln )
    {
    }


    /**
     * \fn ExceptionObject( const ExceptionObject& xpt ) : exception()
     *
     * \brief Clones an instance of this class.
     *
     * \param xpt A reference to another instance of this class.
     *
     */
    ExceptionObject( const ExceptionObject& xpt ) : exception()
    {
      _location = xpt._location ;
      _description = xpt._description ;
      _file = xpt._file ;
      _line = xpt._line ;
    }


    /**
     * \fn virtual ~ExceptionObject() throw()
     *
     * \brief Releases the memory held by an instance of this class.
     *
     */
    virtual ~ExceptionObject() throw()
    {
    }


    /**
     * \fn ExceptionObject& operator=( const ExceptionObject& xpt )
     *
     * \brief Overloads the assignment operator. 
     *
     */
    ExceptionObject& operator=( const ExceptionObject& xpt )
    {
      _location = xpt._location ;
      _description = xpt._description ;
      _file = xpt._file ;
      _line = xpt._line ;
      
      return *this ;
    }


    /**
     * \fn virtual const char* get_name_of_class() const 
     *
     * \brief Returns the name of this class.
     *
     * \return The name of this class.
     *
     */
    virtual const char* get_name_of_class() const 
    {
      return "ExceptionObject" ;
    }


    /**
     * \fn virtual void set_location( const std::string& s )   
     *
     * \brief Assigns a location to this exception.
     *
     * \param s A string containing the location.
     *
     */
    virtual void set_location( const std::string& s )
    {
      _location = s ;
    }


    /**
     * \fn virtual void set_location( const char* s )   
     *
     * \brief Assigns a location to this exception.
     *
     * \param s A pointer to a string containing the location.
     *
     */
    virtual void set_location( const char* s )
    {
      _location = s ;
    }


    /**
     * \fn virtual void set_description( const std::string& s )   
     *
     * \brief Assigns a description to this exception.
     *
     * \param s A string containing the description.
     *
     */
    virtual void set_description( const std::string& s )
    {
      _description = s ;
    }


    /**
     * \fn virtual void set_description( const char* s )   
     *
     * \brief Assigns a description to this exception.
     *
     * \param s A pointer to a string containing the description.
     *
     */
    virtual void set_description( const char* s )
    {
      _description = s ;
    }


    /**
     * \fn virtual const char* get_location() const 
     *
     * \brief Returns the location where this exception occurs.
     *
     * \return The location where this exception occurs.
     *
     */
    virtual const char* get_location() const 
    {
      return _location.c_str() ;
    }

    
    /**
     * \fn virtual const char* get_description() const 
     *
     * \brief  Returns a  description of  the error  that caused  this
     * exception.
     *
     * \return A description of the error that caused this exception.
     *
     */
    virtual const char* get_description() const 
    {
      return _description.c_str() ;
    }


    /**
     * \fn virtual const char* get_file() const 
     *
     * \brief Returns  the name of  the file containing the  line that
     * caused the exception.
     *
     * \return The  name of the  file containing the line  that caused
     * the exception.
     *
     */
    virtual const char* get_file() const 
    {
      return _file.c_str() ;
    }


    /**
     * \fn virtual unsigned get_line() const
     *
     * \brief Returns the line that caused this exception.
     *
     * \return The line that caused this exception.
     *
     */
    virtual unsigned get_line() const 
    {
      return _line ;
    }


    /**
     * \fn virtual const char* what() const throw()
     *
     * \brief  Returns a  description of  the error  that caused  this
     * exception.
     *
     * \return A description of the error that caused this exception.
     *
     */
    virtual const char* what() const throw()
    {
      return _description.c_str() ;
    }
 
  } ;

}

/** @} */ //end of group class.
