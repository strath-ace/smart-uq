# include "Non-Intrusive/utils.h"

//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the smaller of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
  int value;

  if ( i1 < i2 ) 
  {
    value = i1;
  }
  else 
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
  int value;

  if ( i2 < i1 ) 
  {
    value = i1;
  }
  else 
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

string i4_to_string ( int i4, string format )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_STRING converts an I4 to a C++ string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I4, an integer.
//
//    Input, string FORMAT, the format string.
//
//    Output, string I4_TO_STRING, the string.
//
{
  char i4_char[80];
  string i4_string;

  sprintf ( i4_char, format.c_str ( ), i4 );

  i4_string = string ( i4_char );

  return i4_string;
}
//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  const double value = 2.220446049250313E-016;

  return value;
}

//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
//****************************************************************************80

void print_double_vector( vector<double> v, string filnam, string header )

//****************************************************************************80
//
//  Purpose:
//
//    PRINT_DOUBLE_VECTOR print a vector if filname is a non-empty string
//
//  Discussion
//
//    It works only for double. If filname is a non-empty string, it opens the file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//   14 October 2015
//
//  Author:
//
//    Chiara Tardioli
//
//  Parameters:
//
//    Input, vector<int> V, vector of doubles
//
//    Input, string FILNAM, name of the output file to open (empty for the standard output)
//
{
  vector<double>::iterator it;

  if( filnam.length() )
    {
      ofstream out;
      const char *fil = filnam.c_str();
      out.open(fil);

      out << "    " ;
      for(it = v.begin(); it != v.end(); ++it)
	out << *it << " ";
      out << endl;
      out.close();
    }
  else
    {
      cout << endl << header << endl;
      for(it = v.begin(); it != v.end(); ++it)
	cout << *it << " ";
      cout << endl;
    }
}
//****************************************************************************80

void print_double_matrix(vector< vector<double> > vv, string filnam, string header)

//****************************************************************************80
//
//  Purpose:
//
//    PRINT_DOUBLE_MATRIX print a matrix if filname is a non-empty string
//
//  Discussion
//
//    It works only for double. If filname is a non-empty string, it opens the file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//   14 October 2015
//
//  Author:
//
//    Chiara Tardioli
//
//  Parameters:
//
//    Input, vector< vector<int> > VV, matrix of doubles
//
//    Input, string FILNAM, name of the output file to open (empty for the standard output)
//
{
  vector<double>::iterator it;

  if( filnam.length() )
    {
      ofstream out;
      const char *fil = filnam.c_str();
      out.open(fil);

      for( unsigned int i=0; i<vv.size(); ++i)
	{
	  out << "    " ;
	  for(it = vv[i].begin(); it != vv[i].end(); ++it)
	    out << *it << " ";
	  out << endl;
	}
      out.close();
    }
  else
    {
      cout << endl << header << endl;
      for( unsigned int i=0; i<vv.size(); ++i)
	{
	  for(it = vv[i].begin(); it != vv[i].end(); ++it)
	    cout << *it << " ";
	  cout << endl;
	}
    }
}
