// old file, apparently no longer used...

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHR "chr"
#define SUFFIX ".txt"
#define MAX_CHROMOSOMES (308 + 1)
#define MAX_FILES  ( MAX_CHROMOSOMES + 5 )

#define chrX "chrX"
#define chrY "chrY"
#define chrM "chrM"
#define UC "chr_uc"
#define UNKNOWN "chr_unknown"

FILE *chrFiles[MAX_FILES];
int openFiles = 0;
char *prefix;

FILE *pick_file(char *chr) {
    int id;
    char *fn_chr = chr;
    if (!chr) id = MAX_FILES -1;
    else {
        id = atoi((chr+3));
        if (0 == id) {
            if (strcmp(chr, chrX) == 0)  id = MAX_CHROMOSOMES +1;
            else if (strcmp(chr, chrY) == 0) id = MAX_CHROMOSOMES +2;
            else if (strcmp(chr, chrM) == 0) id = MAX_CHROMOSOMES +3;
            else if (strncmp(chr, "uc_", 3) == 0) {
                id = MAX_CHROMOSOMES +4;
                fn_chr = UC;
            } else {
                //id = 0
                fn_chr = UNKNOWN;
            }
       }
    }
    //fprintf(stderr,"pick_file:: guessing file %i for %s\n",id,chr+3);
    if (!chrFiles[id]) {
        char *fn;
        int sz = sizeof(char) * (strlen(fn_chr)+strlen(prefix)+strlen(SUFFIX)+1+1); // bug fix 6 jan 2016, added +1 to include space for terminating NULL character
        fn = (char *)malloc(sizeof(char) * sz);
        memset(fn, 0, sz);
        sprintf(fn,"%s.%s%s",prefix,fn_chr,SUFFIX);
        if (!(chrFiles[id] = fopen(fn,"w"))) {
            fprintf(stderr,"Panic! Cannot open chrom file %s for writing!\n",fn);
            perror("");
            fprintf(stderr,"Ignoring.");
            return 0;
        }
    }
    return chrFiles[id];
}

char *get_chrom(char *line) {
    char *tok;
    tok = strtok(line,"\t");
    tok = strtok(NULL,"\t");
    tok = strtok(NULL,"\t");
    return tok;
}

void parse(FILE *in) {
    char *line = 0;
    char *tmp, *chrom;
    size_t n;
    FILE *outfile;

    tmp = 0;
    while(!(feof(in)) && (getline(&line,&n,in) != -1)) {
        if (tmp) { free(tmp); tmp=0; }
        tmp = strdup(line);
        chrom = get_chrom(tmp);
        outfile = pick_file(chrom);
        if (!outfile) continue;
        fprintf(outfile,"%s",line);
    }
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <reads file>\n", argv[0]);
        return 1;
    }
    FILE *fr;

    memset(chrFiles, 0, MAX_CHROMOSOMES * sizeof(FILE *));

    prefix = strdup(argv[1]);
    if (!(fr = fopen(argv[1], "r"))) {
        fprintf(stderr, "Could not open %s for reading", argv[1]);
        perror("");
        return 1;
    }
    
    parse(fr);
    fcloseall();
    return 0;
}
