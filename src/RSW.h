#ifndef RSW_H_
#define RSW_H_

class RSW {
 public:
  const char *id;
  char side;
  int length;
  char direction;
  const char *chromosome;
  long int position;
  int count;
  int hash;
};



class RSW_Known {
 public:
  const char *id1;
  const char *id2;
  const char *chromosome;
  char direction;
  long int position1;
  long int position2;
};

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



class RSW_splice {
 public:
  const char * id;
  const char * geneName;
  int geneUnknown;
  const char * chromosome;
  char direction;
  long int positionSmaller;
  long int positionLarger;
  long int minSmallSupport;
  long int maxLargeSupport;
  unordered_map <const char *, RSW_splice *> supported_reads;
  //vector <RSW_splice *> supported_reads;
  bool reported; // has this result already been reported as supporting another read
  bool novel; // does this show up in refFlatBoundary file as already known
};

class RSW_result {
 public:
  const char * geneName;
  const char * chromosome;
  long int supportCount;
  long int spliceLength;
  long int minSmallSupport;
  long int maxLargeSupport;
  bool novel; // does this show up in refFlatBoundary file as already known
};



// return 1 if line received, return 0 if EOF was encountered, -1 if error
// remove newline character from string.
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

#define MAX_LINE 10000
#define MAX_STR_LEN 10000


#endif
