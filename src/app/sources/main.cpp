/** 
 * \file main.cpp
 *
 * \brief A simple program for testing the channel-2d library.
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

#include <exceptionobject.hpp>      // channel::ExceptionObject
#include <curvebuilder.hpp>         // channel::CurveBuilder

#include <iostream>                 // std::cout, std::endl, std::cerr
#include <fstream>                  // std::ifstream, std::ofstream
#include <sstream>                  // std::sstream
#include <string>                   // std::string
#include <cstdlib>                  // exit, EXIT_SUCCESS, EXIT_FAILURE, size_t
#include <iomanip>                  // std::setprecision
#include <cassert>                  // assert
#include <ctime>                    // time, clock, CLOCKS_PER_SEC, clock_t
#include <cmath>                    // fabs

using std::cin  ;
using std::cout ;
using std::cerr ;
using std::endl ;
using std::string ;
using channel::CurveBuilder ;
using channel::ExceptionObject ;



/**
 * \fn void read_input( const string& fn , size_t& np , size_t& nc , bool& closed , double*& lx , double*& ly , double*& ux , double*& uy )
 *
 * \brief Read in a file describing a polygonal channel.
 *
 * \param fn The name of a file describing a polygonal channel.
 * \param np A reference to the number of b-spline segments.
 * \param nc A reference to the number of c-segments of the channel.
 * \param closed A reference to a flag to indicate whether the channel
 * is closed.
 * \param  lx  A  reference  to  a   pointer  to  an  array  with  the
 * x-coordinates of the lower polygonal chain of the channel.
 * \param  ly  A  reference  to  a   pointer  to  an  array  with  the
 * y-coordinates of the lower polygonal chain of the channel.
 * \param  ux  A  reference  to  a   pointer  to  an  array  with  the
 * x-coordinates of the upper polygonal chain of the channel.
 * \param  uy  A  reference  to  a   pointer  to  an  array  with  the
 * y-coordinates of the upper polygonal chain of the channel.
 *
 */
void read_input(
                const string& fn ,
                size_t& np ,
                size_t& nc ,
                bool& closed ,
                double*& lx ,
                double*& ly ,
                double*& ux ,
                double*& uy
               );


/**
 * \fn void write_solution( const string& fn , const CurveBuilder& b )
 *
 * \brief Write the control points of  the b-spline curve to an output
 * file.
 *
 * \param fn The name of the output file.
 * \param b An instance of the b-spline curve builder.
 *
 */
void write_solution(
                    const string& fn ,
                    const CurveBuilder& b
                   ) ;


/**
 * \fn void write_lp( const string& fn , const CurveBuilder& b )
 *
 * \brief Write the  instance of the linear program  problem solved by
 * this program in  CPLEX format. The output file can  be given to the
 * \e gpsolve function of the GNU GLPK or to debug the assembly of the
 * constraints.
 *
 * \param fn The name of the output file.
 * \param b An instance of the spline curve builder.
 *
 */
void write_lp(
              const string& fn ,
              const CurveBuilder& b
             ) ;


/**
 * \fn int main( int argc , char* argv[]  )
 *
 * \brief A simple program for testing the bc2d library.
 *
 * \param argc The number of command-line arguments.
 * \param argv An array with the command-line arguments.
 *
 * \return An integer number.
 */
int main( int argc , char* argv[] ) {
  //
  // Check command-line arguments.
  //
  
  if ( ( argc != 3 ) && ( argc != 4 ) ) {
    cerr << "Usage: "
         << endl
         << "\t\t channel2d arg1 arg2 [ arg3 ]"
         << endl
         << "\t\t arg1: name of the file describing the polygonal channel"
         << endl
         << "\t\t arg2: name of the output file describing the computed cubic b-spline curve"
         << endl
         << "\t\t arg3: name of an output file to store a CPLEX format definition of the LP solved by this program (OPTIONAL)"
         << endl
         << endl ;
    cerr.flush() ;
    
    return EXIT_FAILURE ;
  }
  
  //
  // Read in the input file.
  //
  
  clock_t start, end ;
  
  cerr << endl
       << "Reading file describing a polygonal channel..."
       << endl ;
  cerr.flush() ;
  
  string fn1( argv[ 1 ] ) ;

  size_t np ;
  size_t nc ;
  bool closed ;
  double* lx ;
  double* ly ;
  double* ux ;
  double* uy ;

  start = clock() ;
  try {
    read_input( fn1 , np , nc , closed , lx , ly , ux , uy ) ;
  }
  catch ( const ExceptionObject& xpt ) {
    treat_exception( xpt ) ;
    exit( EXIT_FAILURE ) ;
  }  
  end = clock() ;
  
  cerr << ( (double) ( end - start ) ) / CLOCKS_PER_SEC
       << " seconds."
       << endl
       << endl ;
  cerr.flush() ;
  
  //
  // Compute a cubic b-spline curve that passes through the channel.
  //
  
  cerr << "Computing a cubic b-spline curve that passes through the channel... "
       << endl ;
  cerr.flush() ;
  
  start = clock() ;
  CurveBuilder* builder = 0 ;
  try {
    builder = new CurveBuilder(
                               np ,
                               nc ,
                               closed ,
                               &lx[ 0 ] ,
                               &ly[ 0 ] ,
                               &ux[ 0 ] ,
                               &uy[ 0 ] 
                              ) ;
  }
  catch ( const ExceptionObject& xpt ) {
    treat_exception( xpt ) ;
    exit( EXIT_FAILURE ) ;
  }
  
  int error ;
  bool res = builder->build( error ) ;
  end = clock() ;
  
  cerr << ( (double) ( end - start ) ) / CLOCKS_PER_SEC
       << " seconds."
       << endl
       << endl ;
  cerr.flush() ;

  if ( res ) {
    //
    // Write the control points of the b-spline to a file.
    //
    cerr << "Writing out the control points of the b-spline to a file..."
	       << endl ;
    cerr.flush() ;

    start = clock() ;
    string fn2( argv[ 2 ] ) ;
    write_solution(
                   fn2 ,
                   *builder
                  ) ;
    end = clock() ;
  
    cerr << ( (double) ( end - start ) ) / CLOCKS_PER_SEC
         << " seconds."
         << endl
         << endl ;
    cerr.flush() ;
  }
  else {
    //
    // Print the error message returned by the LP solver.
    //
    cerr << endl
         << "ATTENTION: "
         << endl
         << builder->get_solver_error_message( error )
         << endl
         << endl ;
  }

  //
  // Generate a description of the linear program in CPLEX format.
  //
  if ( argc == 4 ) {
    cerr << "Writing out a description of the linear program in CPLEX format..."
	       << endl ;
    cerr.flush() ;
  
    start = clock() ;  
    string fn3( argv[ 3 ] ) ;
    write_lp(
             fn3 ,
             *builder
            ) ;
    end = clock() ;
  
    cerr << ( (double) ( end - start ) ) / CLOCKS_PER_SEC
         << " seconds."
         << endl
         << endl ;
    cerr.flush() ;
  }
  
  //
  // Release memory
  //
  
  cerr << "Releasing memory..."
       << endl ;
  cerr.flush() ;
  
  start = clock() ;
  if ( lx != 0 ) delete[ ] lx ;
  if ( ly != 0 ) delete[ ] ly ;
  if ( ux != 0 ) delete[ ] ux ;
  if ( uy != 0 ) delete[ ] uy ;
  if ( builder != 0 ) delete builder ;
  end = clock() ;
  
  cerr << ( (double) ( end - start ) ) / CLOCKS_PER_SEC
  << " seconds."
  << endl
  << endl ;
  cerr.flush() ;
  
  //
  // Done.
  //
  
  cerr << "Finished."
       << endl
       << endl
       << endl ;
  cerr.flush() ;
  
  return EXIT_SUCCESS ;
}


/**
 * \fn void read_input( const string& fn , size_t& np , size_t& nc , bool& closed , double*& lx , double*& ly , double*& ux , double*& uy )
 *
 * \brief Read in a file describing a polygonal channel.
 *
 * \param fn The name of a file describing a polygonal channel.
 * \param  np A  reference  to  the number  of  b-spline segments.
 * \param nc A reference to the number of c-segments of the channel.
 * \param closed A reference to a flag to indicate whether the channel
 * is closed.
 * \param  lx  A  reference  to  a   pointer  to  an  array  with  the
 * x-coordinates of the lower polygonal chain of the channel.
 * \param  ly  A  reference  to  a   pointer  to  an  array  with  the
 * y-coordinates of the lower polygonal chain of the channel.
 * \param  ux  A  reference  to  a   pointer  to  an  array  with  the
 * x-coordinates of the upper polygonal chain of the channel.
 * \param  uy  A  reference  to  a   pointer  to  an  array  with  the
 * y-coordinates of the upper polygonal chain of the channel.
 *
 */
void read_input(
                const string& fn ,
                size_t& np ,
                size_t& nc ,
                bool& closed ,
                double*& lx ,
                double*& ly ,
                double*& ux ,
                double*& uy
               )
{
  //
  // Open the input file
  //
  std::ifstream in( fn.c_str() ) ;
  
  if ( in.is_open() ) {
    //
    // Read in the number of segments of the b-spline.
    //
    in >> np ;
    
    //
    // Read in the number of c-segments of the channel.
    //
    in >> nc ;
    
    //
    // Read in the flag indicating whether the channel is closed.
    //
    unsigned flag ;
    in >> flag ;
    
    if ( ( flag != 0 ) && ( flag != 1 ) ) {
      std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
      ss << "Flag value indicating whether the channel is closed or open is invalid" ;
      in.close() ;
      throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
    }
      
    closed = ( flag == 1 ) ;
      
    if ( closed ) {
      if ( np < 4 ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "The number of curve segments must be at least 4 for a closed curve" ;
        in.close() ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }
      if ( nc < 3 ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "The number of segments of a closed channel must be at least 3" ;
        in.close() ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }
    }
    else {
      if ( np < 1 ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "The number of curve segments must be at least 1" ;
        in.close() ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }
      if ( nc < 1 ) {
        std::stringstream ss( std::stringstream::in | std::stringstream::out ) ;
        ss << "The number of segments of an open channel must be at least 1" ;
        in.close() ;
        throw ExceptionObject( __FILE__ , __LINE__ , ss.str().c_str() ) ;
      }
    }
      
    //
    // Read in the channel vertex coordinates.
    //
    const size_t nn = ( closed ) ? nc : ( nc + 1 ) ;
      
    lx = new double[ nn ] ;
    ly = new double[ nn ] ;
      
    for ( size_t i = 0 ; i < nn ; i++ ) {
      //
      // Read in the X and Y coordinates of the i-th vertex.
      //
      in >> lx[ i ] ;
      in >> ly[ i ] ;
    }
      
    ux = new double[ nn ] ;
    uy = new double[ nn ] ;
      
    for ( size_t i = 0 ; i < nn ; i++ ) {
      //
      // Read in the X and Y coordinates of the i-th vertex.
      //
      in >> ux[ i ] ;
      in >> uy[ i ] ;
    }
      
    //
    // Close file
    //
      
    in.close() ;
  }
  
  return ;
}


/**
 * \fn void write_solution( const string& fn , const CurveBuilder& b )
 *
 * \brief Write the control points of  the b-spline curve to an output
 * file.
 *
 * \param fn The name of the output file.
 * \param b An instance of the b-spline curve builder.
 *
 */
void write_solution(
                    const string& fn ,
                    const CurveBuilder& b
                   )
{
  using std::endl ;
  
  std::ofstream ou( fn.c_str() ) ;
  
  if ( ou.is_open() ) {
    //
    // Set the precision of the floating-point numbers.
    //
    
    ou << std::setprecision( 6 ) << std::fixed ;
    
    //
    // Write out the number of control points and the degree of the b-spline.
    //
    
    size_t ncp = b.get_number_of_control_points() ;
    
    ou << ncp
       << '\t'
       << 3
       << endl ;
    
    for ( size_t i = 0 ; i < ncp ; i++ ) {
      double x ;
      double y ;
      try {
        x = b.get_control_value( i , 0 ) ;
        y = b.get_control_value( i , 1 ) ;
      }
      catch ( const ExceptionObject& xpt ) {
        treat_exception( xpt ) ;
        ou.close() ;
        exit( EXIT_FAILURE ) ;
      }
      ou << x << '\t' << y << endl ;
    }
    
    //
    // Close file
    //
    
    ou.close() ;
  }

  return ;
}


/**
 * \fn void write_lp( const string& fn , const CurveBuilder& b )
 *
 * \brief Write the  instance of the linear program  problem solved by
 * this program in  CPLEX format. The output file can  be given to the
 * \e gpsolve function of the GNU GLPK or to debug the assembly of the
 * constraints.
 *
 * \param fn The name of the output file.
 * \param b An instance of the spline curve builder.
 *
 */
void write_lp(
              const string& fn ,
              const CurveBuilder& b
             )
{
  //
  // Create the output file
  //

  const size_t np = b.get_number_of_segments() ;
  const size_t dg = b.get_degree() ;

  const size_t NumberOfControlPoints = np + dg ;
  const size_t NumberOfSecondDifferences = ( np - 1 ) * ( dg - 2 ) + ( dg - 1 ) ;

  std::ofstream ou( fn.c_str() ) ;
  
  if ( ou.is_open() ) {
    //
    // Set the precision of the floating-point numbers.
    //
    
    ou << std::setprecision( 6 ) << std::fixed ;

    //
    // Write the objective function
    //
    
    ou << "Minimize"
       << std::endl ;
    ou << '\t'
       << "obj: " ;
    
    size_t j = 1 ;
    for ( size_t i = 0 ; i < NumberOfSecondDifferences ; i++ ) {
      ou << "-mx" << j << " " ;
      ou << "-my" << j << " " ;
      ou << "+px" << j << " " ;
      ou << "+py" << j << " " ;
      ++j ;
    }

    ou << std::endl ;
    
    //
    // Write the constraints
    //

    ou << "Subject To" << std::endl ;
    
    try {
      
      for( size_t i = 0 ; i < b.get_number_of_constraints() ; i++ ) {
        //
        // Write out the number of the constraint.
        //
        ou << '\t' << "c" << i + 1 << ": " ;
        
        //
        // Get the coefficients of the i-th constraint.
        //
        for( j = 0 ; j < b.get_number_of_coefficients_in_the_ith_constraint( i ) ; ++j ) {
          //
          // Get the column index of the coefficient.
          //
          size_t col = b.get_coefficient_identifier( i , j ) ;
          
          //
          // Get the value of the coefficient.
          //
          double value = b.get_coefficient_value( i , j ) ;
          
          //
          // Compute the  index of the  curve piece associated  with the
          // coefficient, and  find the type of  structural variable the
          // coefficient is.
          //
          if ( col < ( 2 * NumberOfControlPoints ) ) {
            if ( ( col % 2 ) == 0 ) {
              if ( value >= 0 ) {
                ou << "+" << value << "x" << ( col >> 1 ) + 1 << " " ;
              }
              else {
                ou << value << "x" << ( col >> 1 ) + 1 << " " ;
              }
            }
            else {
              if ( value >= 0 ) {
                ou << "+" << value << "y" << ( col >> 1 ) + 1 << " " ;
              }
              else {
                ou << value << "y" << ( col >> 1 ) + 1 << " " ;
              }
            }
          }
          else if ( col < ( ( 2 * NumberOfControlPoints ) + ( 4 * NumberOfSecondDifferences ) ) ) {
            col -= ( 2 * NumberOfControlPoints ) ;
            size_t ro = col % 4 ;
            if ( ro == 0 ) {
              if ( value >= 0 ) {
                ou << "+" << value << "mx" << ( col / 4 ) + 1 << " " ;
              }
              else {
                ou << value << "mx" << ( col / 4 ) + 1 << " " ;
              }
            }
            else if ( ro == 1 ) {
              if ( value >= 0 ) {
                ou << "+" << value << "my" << ( col / 4 ) + 1 << " " ;
              }
              else {
                ou << value << "my" << ( col / 4 ) + 1 << " " ;
              }
            }
            else if ( ro == 2 ) {
              if ( value >= 0 ) {
                ou << "+" << value << "px" << ( col / 4 ) + 1 << " " ;
              }
              else {
                ou << value << "px" << ( col / 4 ) + 1 << " " ;
              }
            }
            else if ( ro == 3 ) {
              if ( value >= 0 ) {
                ou << "+" << value << "py" << ( col / 4 ) + 1 << " " ;
              }
              else {
                ou << value << "py" << ( col / 4 ) + 1 << " " ;
              }
            }
          }
          else {
            col -= ( ( 2 * NumberOfControlPoints ) + ( 4 * NumberOfSecondDifferences ) ) ;
            
            if ( col == 0 ) {
              if ( value >= 0 ) {
                ou << "+" << value << "as" << " " ;
              }
              else if ( value < 0 ) {
                ou << value << "as" << " " ;
              }
            }
            else if ( ( col == 1 ) && !b.is_curve_closed() ) {
              if ( value >= 0 ) {
                ou << "+" << value << "ae" << " " ;
              }
              else if ( value < 0 ) {
                ou << value << "ae" << " " ;
              }
            }
            else if ( !b.is_curve_closed() ) {
              if ( value >= 0 ) {
                ou << "+" << value << "co" << ( col - 1 ) << " " ;
              }
              else if ( value < 0 ) {
                ou << value << "co" << ( col - 1 ) << " " ;
              }
            }
            else {
              if ( value >= 0 ) {
                ou << "+" << value << "co" << col << " " ;
              }
              else if ( value < 0 ) {
                ou << value << "co" << col << " " ;
              }
            }
          }
        }
        
        if ( b.is_equality( i ) ) {
          ou << " = " ;
        }
        else if ( b.is_less_than_or_equal_to( i ) ) {
          ou << " <= " ;
        }
        else {
#ifdef DEBUGMODE
          assert( b.is_greater_than_or_equal_to( i ) ) ;
#endif
          ou << " >= " ;
        }
        
        ou << b.get_bound_of_ith_constraint( i ) << std::endl ;
      }
      
    }
    catch ( const ExceptionObject& xpt ) {
      treat_exception( xpt ) ;
      exit( EXIT_FAILURE ) ;
    }
  
    //
    // Write the bounds
    //
    
    ou << "Bounds" << std::endl ;
    
    for ( unsigned k = 0 ; k < NumberOfControlPoints ; k++ ) {
      ou << '\t' << "x" << k + 1 << " free" << std::endl ;
      ou << '\t' << "y" << k + 1 << " free" << std::endl ;
    }
    
    for ( unsigned k = 0 ; k < NumberOfSecondDifferences ; k++ ) {
      ou << '\t' << "-inf <= mx" << k + 1 <<    " <= 0" << std::endl ;
      ou << '\t' << "-inf <= my" << k + 1 <<    " <= 0" << std::endl ;
      ou << '\t' <<    "0 <= px" << k + 1 << " <= +inf" << std::endl ;
      ou << '\t' <<    "0 <= py" << k + 1 << " <= +inf" << std::endl ;
    }

    ou << '\t' << 0.40 << " <= as <= " << 0.60 << std::endl ;
    
    if ( !b.is_curve_closed() ) {
      ou << '\t' << 0.40 << " <= ae <= " << 0.60 << std::endl ;
    }

    const size_t NumberOfCSegments = b.get_number_of_csegments() ;

    for ( unsigned k = 1 ; k < NumberOfCSegments ; k++ ) {
      ou << '\t'
      << 0.40
      << " <= co"
      << k
      << " <= "
      << 0.60
      << std::endl ;
    }

    ou << "End" << std::endl ;
    
    //
    // Close file
    //
    
    ou.close() ;
  }
  
  return ;
}
