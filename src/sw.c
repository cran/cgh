// Smith-Waterman algorithm
// to detect copy number changes
// in microarray CGH data

// T.S.Price 2004-2008

#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Random.h>




#define MAX( x, y ) ( ( (x) > (y) ) ? (x) : (y) )




//////////////////////////////////////////////////////////////////////



// module for Iarray data structure

struct island {
  double score;   // island score
  int start;      // start of island
  int length;     // length of island
};
typedef struct island Island;

struct iarray {
  Island* island;   // ptr to array of islands
  int length;       // length (capacity) of array
  int top;          // top of array (next position to be filled)
};
typedef struct iarray *Iarray;



extern inline Iarray iarray_init( int length );
extern inline void iarray_pushback( Iarray i, Island j );
extern inline void iarray_print( Iarray i );
extern inline int iarray_top( Iarray i );
extern inline Island iarray_element( Iarray i, int element );
extern inline Island island_new( double score, int start, int length );
extern inline double island_score( Island j );
extern inline int island_start( Island j );
extern inline int island_length( Island j );



Iarray
iarray_init( int length )
{
  Iarray i;
  if ( length <= 0 ) {
    length = 1;
  }
  i = ( Iarray ) S_alloc( 1, sizeof( struct iarray ) );
  i->island = ( Island * ) S_alloc( length, sizeof( Island ) );
  i->length = length;
  i->top = 0;
  return i;
}

void
iarray_pushback( Iarray i, Island j )
{
  // check for overflow and reallocate memory if necessary
  if ( i->top >= i->length ) {
    i->island = ( Island * ) S_realloc( ( char * ) i->island, i->length * 2, i->length, sizeof( Island ) );
    i->length *= 2;
  }
  i->island[ i->top++ ] = j;
}

void
iarray_print( Iarray i )
{
  int j;
  Island k;
  Rprintf( "  Score  Start  Length\n" );
  for ( j = 0; j < i->top; j++ ) {
    k = i->island[ j ];
    Rprintf( "%7.4f %6d %7d\n",
      island_score( k ),
      island_start( k ),
      island_length( k ) );
  }
}

int
iarray_top( Iarray i )
{
  return i->top;
}

Island
iarray_element( Iarray i, int element )
{
  return i->island[ element ];
}

Island
island_new( double score, int start, int length )
{
  Island j = { score, start, length };
  return j;
}

double
island_score( Island j )
{
  return j.score;
}

int
island_start( Island j )
{
  return j.start;
}

int
island_length( Island j )
{
  return j.length;
}


//////////////////////////////////////////////////////////////////////



// module for Vector data structure


struct vector {
  double *element;
  int length;
};
typedef struct vector *Vector;


extern inline double vector_element( Vector v, int index );
extern inline void vector_setElement( Vector v, int index, double element );
extern inline int vector_length( Vector v );
extern inline void vector_setLength( Vector v, int len );
extern inline void vector_copy( Vector x, Vector y );
extern inline void vector_REALcopy( double *r, Vector v );
extern inline void vector_copyArray( Vector v, double *array );
extern inline void vector_copyLongIntArray( Vector v, long int *array );
extern inline Vector vector_init( int len );
extern inline void vector_print( Vector v );
extern inline void vector_printInt( Vector v );

double
vector_element( Vector v, int index )
{
  return v->element[ index ];
}

void
vector_setElement( Vector v, int index, double element )
{
  v->element[ index ] = element;
}

// get length of vector
int
vector_length( Vector v )
{
  return v->length;
}

// set length of vector and reallocate memory
void
vector_setLength( Vector v, int len )
{
  v->element = ( double * ) S_realloc(
    ( char * ) v, len, v->length, sizeof( double ) );
  v->length = len;
}

// let vector x = vector y
void
vector_copy( Vector x, Vector y )
{
  int len = vector_length( y );
  vector_setLength( x, len );
  while ( len-- ) {
    vector_setElement( x, len, vector_element( y, len ) );
  }
}

// coerce Vector to REAL
void
vector_REALcopy( double *r, Vector v )
{
  int len = vector_length( v );
  while ( len-- ) {
    r[ len ] = vector_element( v, len );
  }
}

// coerce double * array to Vector
void
vector_copyArray( Vector v, double *a )
{
  int len = vector_length( v );
  while ( len-- ) {
    vector_setElement( v, len , a[ len ] );
  }
}

// coerce long int * array to Vector
void
vector_copyLongIntArray( Vector v, long int *a )
{
  int len = vector_length( v );
  while ( len-- ) {
    vector_setElement( v, len , ( double ) a[ len ] );
  }
}

Vector
vector_init( int len )
{
  Vector v;
  if ( len <= 0 ) {
    len = 1;
  }
  v = ( Vector ) S_alloc( 1, sizeof( struct vector ) );
  v->element = ( double * ) S_alloc( len, sizeof( double ) );
  v->length = len;
  return v;
}

// print vector to 4dp
void
vector_print( Vector v )
{
  int i;
  int len = vector_length( v );
  for ( i = 0; i < len; i++ ) {
    if ( i == ( i / 10 ) * 10 ) {
      Rprintf( "%4d> ", i );
    }
    Rprintf( "%8.4f", vector_element( v, i ) );
    if ( i == ( i / 10 ) * 10 + 9 ) {
      Rprintf( "\n" );
    }
  }
  Rprintf( "\n" );
}

// print vector to 0dp
void
vector_printInt( Vector v )
{
  int i;
  int len = vector_length( v );
  for ( i = 0; i < len; i++ ) {
    if ( i == ( i / 10 ) * 10 ) {
      Rprintf( "%4d> ", i );
    }
    Rprintf( "%8.0f", vector_element( v, i ) );
    if ( i == ( i / 10 ) * 10 + 9 ) {
      Rprintf( "\n" );
    }
  }
  Rprintf( "\n" );
}



//////////////////////////////////////////////////////////////////////



// Smith-Waterman algorithm

extern inline void sw_scores( Vector s, Vector x );
extern inline Island sw_top_island( Vector x );
extern inline double sw_top_island_score( Vector x );
extern inline void sw_loop( Iarray results, Vector x, int max_i );


// calculate partial sums s from scores x
void
sw_scores( Vector s, Vector x )
{
  int i;
  int len = vector_length( x );
  vector_setElement( s, 0, MAX( vector_element( x, 0 ), 0.0 ) );
  for ( i = 1; i < len; i++ ) {
    vector_setElement( s, i,
      MAX( vector_element( s, i - 1 ) + vector_element( x, i ), 0.0 ) );
  }
}

// identify top-scoring island in input vector
Island
sw_top_island( Vector x )
{
  int z = -1;             // location of last zero
  int b = -1;             // location of last zero before maximum
  int e = -1;             // location of maximum
  double m = 0.0;         // maximum
  int len = vector_length( x );
  int i;
  double temp;
  Vector s = vector_init( len );

  // calculate partial sums
  sw_scores( s, x );

  // find highest scoring island
  for ( i = 0; i < len; i++ ) {
    temp = vector_element( s, i );
    if ( temp == 0.0 ) {
      z = i;
    }
    if ( temp > m ) {
      b = z;
      e = i;
      m = temp;
      }
    }

  // save details of island: score = m == 0.0 if no island
  return island_new( m, b + 1, e - b );
}

// identify top island score in input vector
double
sw_top_island_score( Vector x )
{
  double m = 0.0;         // maximum
  int len = vector_length( x );
  int i;
  double temp;
  Vector s = vector_init( len );

  // calculate partial sums
  sw_scores( s, x );

  // find highest scoring island
  for ( i = 0; i < len; i++ ) {
    temp = vector_element( s, i );
    if ( temp > m ) {
      m = temp;
      }
    }

  // return score = m == 0.0 if no island
  return m;
}

void
sw_loop( Iarray results, Vector x, int max_i )
{
  if ( max_i <= 0 ) {
    max_i = INT_MAX;
  }
  int i, loop;
  Island j;

  // loop to find at maximum `max_i' islands
  for ( loop = 0; loop++ < max_i; ) {

    // find top-scoring island
    j = sw_top_island( x );

    // if there is an island:
    if ( island_score( j ) > 0.0 ) {

      // save the details of island in results
      iarray_pushback( results, j );

      // remove the island from x
      // by replacing it with zeros
      for ( i = island_length( j ); i; ) {
        vector_setElement( x, --i + island_start( j ), 0.0 );
      }
    }

    // repeat until no more islands
    else break;
  }
}



//////////////////////////////////////////////////////////////////////



// R interface for Smith-Waterman algorithm

extern SEXP sw( SEXP xR, SEXP max_iR, SEXP traceR );

SEXP
sw( SEXP xR, SEXP max_iR, SEXP traceR )
{
  //Rprintf( "Start of C function sw\n" );

  // coerce input variables to correct data types
  // xR = floating point vector of scores
  if ( !isReal( xR ) ) {
    xR = coerceVector( xR, INTSXP );
    Rprintf( "`xR' coerced to real vector\n" );
  }
  // max_i = max number of islands to extract
  // if max_iR = NULL or max_iR < 0 all islands are extracted
  int mi;
  if ( max_iR == R_NilValue ) {
    mi = 0;
  }
  else {
    if ( !isInteger( max_iR ) ) {
      max_iR = coerceVector( max_iR, INTSXP );
    }
    mi = INTEGER( max_iR )[ 0 ];
    if ( ( ISNAN( mi ) ) || ( !R_FINITE( mi ) ) )
      mi = 0;
    }

  // declare and initialise variables
  int i, top;
  int len = LENGTH( xR );
  int trace = INTEGER( coerceVector( traceR, INTSXP ) )[ 0 ];
  Island j;
  Iarray sw_results = iarray_init( 1 );
  Vector x = vector_init( len );
  Vector s = vector_init( len );
  double *r, *scores, temp;
  int *starts, *lengths;
  SEXP sR, swR, swR_names, i_scores, i_starts, i_lengths;
  PROTECT( swR = allocVector( VECSXP, 5 ) );
  PROTECT( swR_names = allocVector( STRSXP, 5 ) ); //VECSXP before
  SET_STRING_ELT( swR_names, 0, mkChar( "x" ) );
  SET_STRING_ELT( swR_names, 1, mkChar( "s" ) );
  SET_STRING_ELT( swR_names, 2, mkChar( "score" ) );
  SET_STRING_ELT( swR_names, 3, mkChar( "start" ) );
  SET_STRING_ELT( swR_names, 4, mkChar( "length" ) );
  SET_NAMES( swR, swR_names );
  UNPROTECT( 1 );

  // copy xR into Vector x
  // checking for missing or infinite values
  for ( i = 0; i < len; i++ ) {
    temp = REAL( xR )[ i ];
    if ( ( ISNAN( temp ) ) || ( !R_FINITE( temp ) ) ) {
      REprintf( "NaN / infinite value in input vector\n" );
    }
    vector_setElement( x, i, temp );
  }

  // calculate running totals and copy into sR
  sw_scores( s, x );
  PROTECT( sR = allocVector( REALSXP, len ) );
  r = REAL( sR );
  vector_REALcopy( r, s );

  // print input and running total vectors
  if ( trace ) {
    Rprintf( "Input vector:\n" );
    vector_print( x );
    Rprintf( "Running totals:\n" );
    vector_print( s );
  }

  // copy input vector and running totals into output
  SET_VECTOR_ELT( swR, 0, xR );
  SET_VECTOR_ELT( swR, 1, sR );
  UNPROTECT( 1 );
  SET_VECTOR_ELT( swR, 2, R_NilValue );
  SET_VECTOR_ELT( swR, 3, R_NilValue );
  SET_VECTOR_ELT( swR, 4, R_NilValue );

  // run Smith-Waterman algorithm and print results
  sw_loop( sw_results, x, mi );
  /* Rprintf( "\nSmith-Waterman algorithm results:\n" );    */
  /* iarray_print( sw_results );                            */

  // convert results to R format and return
  top = iarray_top( sw_results );
  if ( top ) {
    PROTECT( i_scores = allocVector( REALSXP, top ) );
    PROTECT( i_starts = allocVector( INTSXP, top ) );
    PROTECT( i_lengths = allocVector( INTSXP, top ) );
    scores = REAL( i_scores );
    starts = INTEGER( i_starts );
    lengths = INTEGER( i_lengths );
    for ( i = 0; i < top; i++ ) {
      j = iarray_element( sw_results, i );
      scores[ i ] = island_score( j );
      starts[ i ] = island_start( j ) + 1;               // R vectors are numbered beginning at 1 not 0
      lengths[ i ] = island_length( j );
    }
    SET_VECTOR_ELT( swR, 2, i_scores );
    SET_VECTOR_ELT( swR, 3, i_starts );
    SET_VECTOR_ELT( swR, 4, i_lengths );
    UNPROTECT( 3 );
  }
  UNPROTECT( 1 );
  return swR;
}



//////////////////////////////////////////////////////////////////////



// R interface for permutation test
// optimized for speed and memory usage in version 1.0-5

extern SEXP sw_permTest( SEXP xR, SEXP max_iR, SEXP nItR, SEXP seedR,
  SEXP traceR, SEXP envR );


SEXP
sw_permTest( SEXP xR, SEXP max_iR, SEXP nItR, SEXP seedR, SEXP traceR,
  SEXP envR )
{
  // declare variables
  int i, m, n, ni, nullseed = 0;
  int len = LENGTH( xR );
  int trace = INTEGER( coerceVector( traceR, INTSXP ) )[ 0 ];
  long int loop;
  long int nIter = INTEGER( coerceVector( nItR, INTSXP ) )[ 0 ];
  long int *p;
  double temp, *x_perm, *score, *i_scores, *pVal;
  Vector q;
  double j;
  SEXP swR, pValR, tempR, runifR, s, t;
  if ( TYPEOF( seedR ) == NILSXP ) {
    nullseed = 1;
  }
  if ( trace ) {
    Rprintf( "Running permutation test\n\n" );
  }

  // initialise variables
  PROTECT( swR = allocVector( VECSXP, 5 ) );
  swR = sw( xR, max_iR, traceR );
  x_perm = ( double * ) S_alloc( len, sizeof( double ) );
  score = ( double * ) S_alloc( len, sizeof( double ) );
  for ( i = 0; i < len; i++ ) {
    x_perm[ i ] = REAL( xR )[ i ];
  }
  ni = LENGTH( VECTOR_ELT( swR, 2 ) );                // number of islands extracted
  p = ( long int * ) S_alloc( ni, sizeof( double ) );
  q = vector_init( ni );
  i_scores = ( double * ) S_alloc( ni, sizeof( double ) );
  pVal = ( double * ) S_alloc( ni, sizeof( double ) );
  for ( i = 0; i < ni; i++ ) {
    p[ i ] = 0.0;
    i_scores[ i ] = REAL( VECTOR_ELT( swR, 2 ) )[ i ];
  }
  UNPROTECT( 1 ); // swR

  if ( trace ) {
    Rprintf( "`nIter' == %d\n", nIter );
  }
  if ( nIter <= 0 ) {
    PROTECT( pValR = allocVector( NILSXP, 1 ) );
    UNPROTECT( 1 ); // pValR
    return pValR;
  }

  // initialise random number generator
  if ( nullseed ) {
    GetRNGstate();
  }
  else {
    SEXP kindR;
    PROTECT( t = s = allocList( 3 ) );
    SET_TYPEOF( s, LANGSXP );
    SETCAR( t, install( "set.seed" ) );
    t = CDR( t );
    SETCAR( t, seedR );
    SET_TAG( t, install( "seed" ) );
    t = CDR( t );
    PROTECT( kindR = allocVector( STRSXP, 1 ) );
    SET_STRING_ELT( kindR, 0, mkChar( "Mersenne-Twister" ) );
    SETCAR( t, kindR );
    SET_TAG( t, install( "kind" ) );
    tempR = eval( s, envR );
    UNPROTECT( 2 ); // t, s, kindR
  }

  // loop over `nIter' permutations of xR
  PROTECT( t = s = allocList( 3 ) );
  PROTECT( runifR = allocVector( REALSXP, len ) );
  for ( loop = 0; loop++ < nIter; ) {

    if ( trace ) {
      Rprintf( "%9d ", loop );
    }
    if ( nullseed ) {
      // generate permutation of xR
      // by swapping each element of `x_perm' in turn
      // with another randomly chosen element
      for ( m = 0; m < len; ) {
        n = ( int ) ( unif_rand() * ( double ) len );
        temp = x_perm[ m ];
        x_perm[ m++ ] = x_perm[ n ];
        x_perm[ n ] = temp;
      }
    }
    else {
      // generate vector of random numbers
      SET_TYPEOF( s, LANGSXP );
      SETCAR( t, install( "runif" ) );
      t = CDR( t );
      SETCAR( t, allocVector( INTSXP, 1 ) );
      INTEGER( CAR( t ) )[ 0 ] = len;
      SET_TAG( t, install( "n" ) );
      t = CDR( t );
      SETCAR( t, allocVector( INTSXP, 1 ) );
      INTEGER( CAR( t ) )[ 0 ] = len;
      SET_TAG( t, install( "max" ) );
      runifR = eval( s, envR );

      // generate permutation of xR
      // by swapping each element of `x_perm' in turn
      // with another randomly chosen element
      for ( m = 0; m < len; ) {
        n = ( int ) REAL( runifR )[ m ];
        temp = x_perm[ m ];
        x_perm[ m++ ] = x_perm[ n ];
        x_perm[ n ] = temp;
      }
    }

    // do Smith-Waterman algorithm on permuted data
    // calculate partial sums
    score[ 0 ] = MAX( x_perm[ 0 ], 0.0 );
    for ( i = 1; i < len; i++ ) {
      score[ i ] = MAX( score[ i - 1 ] + x_perm[ i ], 0.0 );
    }

    // find highest scoring island
    j = 0.0; // maximum
    for ( i = 0; i < len; i++ ) {
      temp = score[ i ];
      if ( temp > j ) {
        j = temp;
      }
    }
    // j == 0.0 if no island

    if ( trace ) {
      Rprintf( "%9f", j );
    }

    // compare score for top island in permuted data
    // with scores from the extracted islands
    for ( i = 0; i < ni; i++ ) {
      if ( j > i_scores[ i ] ) {
        for ( m = i; m < ni; ) {
          ++p[ m++ ];
        }
        break;
      }
    }
    if ( trace ) {
      vector_copyLongIntArray( q, p );
      vector_printInt( q );
    }
  }
  UNPROTECT( 2 ); // t, s, runifR

  if ( nullseed ) {
    PutRNGstate();
  }

  // divide ( p + 1 ) by ( nIter + 1 ) to get permutation p-values
  // convert to R numeric vector `pValR' and return
  PROTECT( pValR = allocVector( REALSXP, ni ) );
  pVal = REAL( pValR );
  for ( i = 0; i < ni; i++ ) {
    pVal[ i ] = ( ( double ) p[ i ] + 1.0 ) / ( ( double ) nIter + 1.0 );
  }
  UNPROTECT( 1 ); // pValR
  return pValR;
}



//////////////////////////////////////////////////////////////////////



extern void R_init_sw( DllInfo *info );
extern void R_unload_sw( DllInfo *info );

// Registration information for DLL
static R_CallMethodDef callMethods[] = {
  { "sw", ( DL_FUNC ) &sw, 3,
    /* { REALSXP, INTSXP, LGLSXP } */ },
  { "sw_permTest", ( DL_FUNC ) &sw_permTest, 6,
    /* { REALSXP, INTSXP, INTSXP, INTSXP, LGLSXP, ENVSXP } */ },
  { NULL, NULL, 0 }
};

void
R_init_sw( DllInfo *info )
{
  //R_useDynamicSymbols( dll, FALSE );                   // don't know what this command does
  R_registerRoutines( info, NULL, callMethods, NULL, NULL );
}

void
R_unload_sw( DllInfo *info )
{
  /* Release resources. */
}
