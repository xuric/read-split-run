#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define CHR "Chr"
#define SUFFIX ".txt"
#define MAX_CHROMOSOMES 308


FILE *chrFiles[MAX_CHROMOSOMES + 2];
int openFiles = 0;
char *prefix;
const char *const file310 = "Chr_unknown";

FILE *pick_file(char *chr) {
    int id;
    if (!chr) id = MAX_CHROMOSOMES+1;
    else id = atoi((chr+3));
    if (0 == id) id = MAX_CHROMOSOMES+1;
    //fprintf(stderr,"pick_file:: guessing file %i for %s\n",id,chr+3);
    if (!chrFiles[id]) {
        char *fn;
        if (MAX_CHROMOSOMES+1 == id) {
            fn = (char *)malloc(strlen(file310)+strlen(prefix)+strlen(SUFFIX)+3);
            sprintf(fn,"%s.%s%s",prefix,file310,SUFFIX);
        }
        else {
            fn = (char *)malloc(strlen(chr)+strlen(prefix)+strlen(SUFFIX)+3);
            sprintf(fn,"%s.%s%s",prefix,chr,SUFFIX);
        }
        if (!(chrFiles[id] = fopen(fn,"w"))) {
 	    fprintf(stderr,"Panic! Cannot open chrom file %s for writing!\n", fn);
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
    while(tok && strncmp(tok,CHR,strlen(CHR)) != 0)
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
    clock_t after, before = clock();

    memset(chrFiles, 0, MAX_CHROMOSOMES * sizeof(FILE *));

    prefix = strdup(argv[1]);
    if (!(fr = fopen(argv[1], "r"))) {
        fprintf(stderr, "Could not open %s for reading", argv[1]);
        perror("");
        return 1;
    }
    
    parse(fr);
    fcloseall();
    after = clock();
    printf("Elapsed time: %lf seconds\n", (after-before)*(1.0/CLOCKS_PER_SEC));
    return 0;
}
