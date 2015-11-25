/*
  File:            RSW.h

  Copyright 2015 Jeff Kinne, Yongheng Bai, Brandon Donham.
  Permission to use for academic, non-profit purposes is granted, in which case this
  copyright notice should be maintained and original authors acknowledged.

  Author:          Jeff Kinne, jkinne@cs.indstate.edu

  Version history: Not tracked in this file, see splitPairs.cpp and compare.cpp

  Contents:        Definitions of classes used for keeping track of data read
                   in for the RSR splitPairs.cpp and compare.cpp programs.  RSW.h
                   is shared by those two files.
*/

#ifndef RSW_H_
#define RSW_H_

/*
  Used to store an RSW record when it's read off of disk into memory.
  Only the fields that are needed for determining matched pairs and supporting
  reads are kept.
*/
class RSW {
 public:
  const char *id; // pointer to id of strand
  char side;      // from left or right side of splitting a read
  int length;     // number of base pairs of this split
  char direction; // which strand/direction the match of this split is on
  const char *chromosome; // pointer to which chromosome
  long int position;      // position on chromosome where split matches
  int count; // not currently used
  int hash;  // not currently used
  const char *sequence; // sequence
};



/*
  Used to store a record from refFlat file that lists genes and the
  chromosomes/positions the genes are at.
*/
class RSW_Known {
 public:
  const char *id1; // pointer to name of gene
  const char *id2; // unused
  const char *chromosome; // chromosome gene is on
  char direction;         // which direction
  long int position1;     // start of gene
  long int position2;     // end of gene
};

/*
  Used to store a record from refFlat introns/extrons boundaries
  file.  Same fields as RSW_Known
*/
class RSW_Boundaries {
 public:
  const char *id1;
  const char *id2;
  const char *chromosome;
  char direction;
  long int length;
  long int position1;
  long int position2;
};



/*
  Used to store a information about a possible splice junction.
*/
class RSW_splice {
 public:
  const char * id;       // id of a read this junction is based on
  const char * geneName; //
  int geneUnknown;       // is this junction within a gene or outside of known genes
  const char * chromosome;
  char direction;
  long int positionSmaller; // one end of the junction
  long int positionLarger;  // other end of the junction
  long int minSmallSupport; // smallest position of the split of a supporting read
  long int maxLargeSupport; // largest position of the split of a supporting read
  unordered_map <const char *, RSW_splice *> supported_reads; // list of supporting reads
  bool novel;               // does this show up in refFlatBoundary file as already known

  const char *sequenceSmaller; // bracketed sequence from smaller position
  const char *sequenceLarger;  // bracketed sequence from larger position
};

/*
  Used in compare.cpp to read in from the output of splitPairs.cpp - information about
  a possible splice junction.
*/
class RSW_result {
 public:
  const char * geneName;
  const char * chromosome;
  long int supportCount; // how many reads supported this junction
  long int spliceLength; // number of base pairs spliced out
  long int minSmallSupport; // from RSW_splice record
  long int maxLargeSupport; // from RSW_splice record
  bool novel; // does this show up in refFlatBoundary file as already known
};



/*
  Function: get_line, read a line from a file, removing newline character

  Parameters: f - opened file handle, s - char array to read into, maxChars - size of s

  Return: 1 if line received, 0 if EOF encountered, -1 if error

 return 1 if line received, return 0 if EOF was encountered, -1 if error
*/
int get_line(FILE *f, char s[], int maxChars) {
  int i = 0;
  int ch;
  do {
    ch = fgetc(f);
    if (ch == EOF) 
      break;
    if (i < maxChars) {
      s[i] = ch; i++;
    }
  } while (ch != '\n');
  if (i > 0) {
    if (s[i-1] == '\n')
      s[i-1] = '\0';
    else 
      s[i] = '\0';
  }
  else 
    s[0] = '\0';

  if (ch == EOF) return 0;
  else if (i >= maxChars) return -1; // error, not enough room in s for line
  else return 1;
}

// These should be long enough so that no line in any of the input files is
// longer than this amount ...  It is a bit trusting of us to assume this, but
// if there is a problem an error message will be printed.
#define MAX_LINE 10000
#define MAX_STR_LEN 10000


#endif
