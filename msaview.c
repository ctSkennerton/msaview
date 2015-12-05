#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <getopt.h>
#include <math.h>
#include <stdarg.h>

#include "easel/include/esl_config.h"
#include "easel/include/easel.h"
#include "easel/include/esl_msa.h"
#include "easel/include/esl_msafile.h"
#include "easel/include/esl_alphabet.h"

#include "termbox/include/termbox.h"

// These define bit flags for amino acid residues usful for OR-ing together in consensus or coloring rules

#define RES_G 0x2       // Glycine         Gly                                 
#define RES_A 0x4       // Alanine         Ala                                 
#define RES_L 0x8       // Leucine         Leu                                 
#define RES_M 0x10      // Methionine      Met                             
#define RES_F 0x20      // Phenylalanine   Phe                         
#define RES_W 0x40      // Tryptophan      Trp                             
#define RES_K 0x80      // Lysine          Lys                                 
#define RES_Q 0x100     // Glutamine       Gln                             
#define RES_E 0x200     // Glutamic Acid   Glu                         
#define RES_S 0x400     // Serine          Ser                                 
#define RES_P 0x800     // Proline         Pro
#define RES_V 0x1000    // Valine          Val
#define RES_I 0x2000    // Isoleucine      Ile
#define RES_C 0x4000    // Cysteine        Cys
#define RES_Y 0x8000    // Tyrosine        Tyr
#define RES_H 0x10000   // Histidine       His
#define RES_R 0x20000   // Arginine        Arg
#define RES_N 0x40000   // Asparagine      Asn
#define RES_D 0x80000   // Aspartic Acid   Asp
#define RES_T 0x100000  // Threonine       Thr


// lookup table which stores the integer representation of the residue
const int res_lookup_table[] = {
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
/*     A      B  C      D      E      F      G      H      I      J  K      L      M      N      O*/
    0, RES_A, 0, RES_C, RES_D, RES_E, RES_F, RES_G, RES_H, RES_I, 0, RES_K, RES_L, RES_M, RES_N, 0,
/*  P      Q      R      S      T      U  V      W      X  Y      Z*/
    RES_P, RES_Q, RES_R, RES_S, RES_T, 0, RES_V, RES_W, 0, RES_Y, 0, 0,  0,  0,  0,  0,
/*     a      b  c      d      e      f      g      h      i      j  k      l      m      n      o*/
    0, RES_A, 0, RES_C, RES_D, RES_E, RES_F, RES_G, RES_H, RES_I, 0, RES_K, RES_L, RES_M, RES_N, 0,
/*  p      q      r      s      t      u  v      w      x  y      z*/
    RES_P, RES_Q, RES_R, RES_S, RES_T, 0, RES_V, RES_W, 0, RES_Y, 0, 0,  0,  0,  0,  0,
};

typedef struct {
    char  name;         // Single character name of the rule
    float perc;         // The percent id of residues needed in the column of the alignment for this rule to take effect
    int   residue_list; // The residues that count toward this rule. An OR-ed int based on the defines above.
}ConsensusRules_t;

// ClustalX defines 29 consensus rules so lets just hard code this as an array. 
ConsensusRules_t consensusRules[29];

void init_consensus_rules()
{
    /* These are the default rules as they appear in clustalX parameter file
     * They get translated below as an array of structs
    % = 60% w:l:v:i:m:a:f:c:y:h:p
    # = 80% w:l:v:i:m:a:f:c:y:h:p
    - = 50% e:d
    + = 60% k:r
    g = 50% g
    n = 50% n
    q = 50% q:e
    p = 50% p
    t = 50% t:s
    A = 85% a
    C = 85% c
    D = 85% d
    E = 85% e
    F = 85% f
    G = 85% g
    H = 85% h
    I = 85% i
    K = 85% k
    L = 85% l
    M = 85% m
    N = 85% n
    P = 85% p
    Q = 85% q
    R = 85% r
    S = 85% s
    T = 85% t
    V = 85% v
    W = 85% w
    Y = 85% y
    */
    consensusRules[0].name = '%';
    consensusRules[0].perc = 0.60;
    consensusRules[0].residue_list = RES_W | RES_L | RES_V | RES_I | RES_M | RES_A | RES_F | RES_C | RES_Y | RES_H | RES_P;
    consensusRules[1].name = '#';
    consensusRules[1].perc = .80;
    consensusRules[1].residue_list = RES_W | RES_L | RES_V | RES_I | RES_M | RES_A | RES_F | RES_C | RES_Y | RES_H | RES_P;
    consensusRules[2].name = '-';
    consensusRules[2].perc = .50;
    consensusRules[2].residue_list = RES_E | RES_D;
    consensusRules[3].name = '+';
    consensusRules[3].perc = .60;
    consensusRules[3].residue_list = RES_K | RES_R;
    consensusRules[4].name = 'g';
    consensusRules[4].perc = .50;
    consensusRules[4].residue_list = RES_G; 
    consensusRules[5].name = 'n';
    consensusRules[5].perc = .50;
    consensusRules[5].residue_list = RES_N; 
    consensusRules[6].name = 'q';
    consensusRules[6].perc = .50;
    consensusRules[6].residue_list = RES_Q | RES_E;
    consensusRules[7].name = 'p';
    consensusRules[7].perc = .50;
    consensusRules[7].residue_list = RES_P; 
    consensusRules[8].name = 't';
    consensusRules[8].perc = .50;
    consensusRules[8].residue_list = RES_T | RES_S;
    consensusRules[9].name = 'A';
    consensusRules[9].perc = .85;
    consensusRules[9].residue_list = RES_A; 
    consensusRules[10].name = 'C';
    consensusRules[10].perc = .85;
    consensusRules[10].residue_list = RES_C; 
    consensusRules[11].name = 'D';
    consensusRules[11].perc = .85;
    consensusRules[11].residue_list = RES_D; 
    consensusRules[12].name = 'E';
    consensusRules[12].perc = .85;
    consensusRules[12].residue_list = RES_E; 
    consensusRules[13].name = 'F';
    consensusRules[13].perc = .85;
    consensusRules[13].residue_list = RES_F; 
    consensusRules[14].name = 'G';
    consensusRules[14].perc = .85;
    consensusRules[14].residue_list = RES_G; 
    consensusRules[15].name = 'H';
    consensusRules[15].perc = .85;
    consensusRules[15].residue_list = RES_H; 
    consensusRules[16].name = 'I';
    consensusRules[16].perc = .85;
    consensusRules[16].residue_list = RES_I; 
    consensusRules[17].name = 'K';
    consensusRules[17].perc = .85;
    consensusRules[17].residue_list = RES_K; 
    consensusRules[18].name = 'L';
    consensusRules[18].perc = .85;
    consensusRules[18].residue_list = RES_L; 
    consensusRules[19].name = 'M';
    consensusRules[19].perc = .85;
    consensusRules[19].residue_list = RES_M; 
    consensusRules[20].name = 'N';
    consensusRules[20].perc = .85;
    consensusRules[20].residue_list = RES_N; 
    consensusRules[21].name = 'P';
    consensusRules[21].perc = .85;
    consensusRules[21].residue_list = RES_P; 
    consensusRules[22].name = 'Q';
    consensusRules[22].perc = .85;
    consensusRules[22].residue_list = RES_Q; 
    consensusRules[23].name = 'R';
    consensusRules[23].perc = .85;
    consensusRules[23].residue_list = RES_R; 
    consensusRules[24].name = 'S';
    consensusRules[24].perc = .85;
    consensusRules[24].residue_list = RES_S; 
    consensusRules[25].name = 'T';
    consensusRules[25].perc = .85;
    consensusRules[25].residue_list = RES_T; 
    consensusRules[26].name = 'V';
    consensusRules[26].perc = .85;
    consensusRules[26].residue_list = RES_V; 
    consensusRules[27].name = 'W';
    consensusRules[27].perc = .85;
    consensusRules[27].residue_list = RES_W; 
    consensusRules[28].name = 'Y';
    consensusRules[28].perc = .85;
    consensusRules[28].residue_list = RES_Y; 
}

typedef enum {
    red = 0,
    blue,
    green,
    cyan,
    pink,
    magenta,
    yellow,
    orange
} clustalx_colors_t;

typedef struct {
    char               residue;
    clustalx_colors_t  color;
    char *             rules_list;
} ColorRules_t;

ColorRules_t colorRules[21];

void init_color_rules()
{
    /*
       g = ORANGE
       p = YELLOW
       t = GREEN if t:S:T:%:#
       s = GREEN if t:S:T:#
       n = GREEN if n:N:D
       q = GREEN if q:Q:E:+:K:R
       w = BLUE if %:#:A:C:F:H:I:L:M:V:W:Y:P:p
       l = BLUE if %:#:A:C:F:H:I:L:M:V:W:Y:P:p
       v = BLUE if %:#:A:C:F:H:I:L:M:V:W:Y:P:p
       i = BLUE if %:#:A:C:F:H:I:L:M:V:W:Y:P:p
       m = BLUE if %:#:A:C:F:H:I:L:M:V:W:Y:P:p
       a = BLUE if %:#:A:C:F:H:I:L:M:V:W:Y:P:p:T:S:s:G
       f = BLUE if %:#:A:C:F:H:I:L:M:V:W:Y:P:p
       c = BLUE if %:#:A:F:H:I:L:M:V:W:Y:S:P:p
       c = PINK if C
       h = CYAN if %:#:A:C:F:H:I:L:M:V:W:Y:P:p
       y = CYAN if %:#:A:C:F:H:I:L:M:V:W:Y:P:p
       e = MAGENTA if -:D:E:q:Q
       d = MAGENTA if -:D:E:n:N
       k = RED if +:K:R:Q
       r = RED if +:K:R:Q
       */
    colorRules[0].residue = 'G';
    colorRules[0].color = orange;
    colorRules[0].rules_list = NULL;
    colorRules[1].residue = 'P';
    colorRules[1].color = yellow;
    colorRules[1].rules_list = NULL;
    colorRules[2].residue = 'T';
    colorRules[2].color = green;
    colorRules[2].rules_list = "tST%#";
    colorRules[3].residue = 'S';
    colorRules[3].color = green;
    colorRules[3].rules_list = "tST#";
    colorRules[4].residue = 'N';
    colorRules[4].color = green;
    colorRules[4].rules_list = "nND";
    colorRules[5].residue = 'Q';
    colorRules[5].color = green;
    colorRules[5].rules_list = "qQE+KR";
    colorRules[6].residue = 'W';
    colorRules[6].color = blue;
    colorRules[6].rules_list = "%#ACFHILMVWYPp";
    colorRules[7].residue = 'L';
    colorRules[7].color = blue;
    colorRules[7].rules_list = "%#ACFHILMVWYPp";
    colorRules[8].residue = 'V';
    colorRules[8].color = blue;
    colorRules[8].rules_list = "%#ACFHILMVWYPp";
    colorRules[9].residue = 'I';
    colorRules[9].color = blue;
    colorRules[9].rules_list = "%#ACFHILMVWYPp";
    colorRules[10].residue = 'M';
    colorRules[10].color = blue;
    colorRules[10].rules_list = "%#ACFHILMVWYPp";
    colorRules[11].residue = 'A';
    colorRules[11].color = blue;
    colorRules[11].rules_list = "%#ACFHILMVWYPpTSsG";
    colorRules[12].residue = 'F';
    colorRules[12].color = blue;
    colorRules[12].rules_list = "%#ACFHILMVWYPp";
    colorRules[13].residue = 'C';
    colorRules[13].color = blue;
    colorRules[13].rules_list = "%#AFHILMVWYSPp";
    colorRules[14].residue = 'C';
    colorRules[14].color = pink;
    colorRules[14].rules_list = "C";
    colorRules[15].residue = 'H';
    colorRules[15].color = cyan;
    colorRules[15].rules_list = "%#ACFHILMVWYPp";
    colorRules[16].residue = 'Y';
    colorRules[16].color = cyan;
    colorRules[16].rules_list = "%#ACFHILMVWYPp";
    colorRules[17].residue = 'E';
    colorRules[17].color = magenta;
    colorRules[17].rules_list = "-DEqQ";
    colorRules[18].residue = 'D';
    colorRules[18].color = magenta;
    colorRules[18].rules_list = "-DEnN";
    colorRules[19].residue = 'K';
    colorRules[19].color = red;
    colorRules[19].rules_list = "+KRQ";
    colorRules[20].residue = 'R';
    colorRules[20].color = red;
    colorRules[20].rules_list = "+KRQ";

}

void usage()
{
fprintf(stderr, "msaview [-f <format>] <msafile>\n\
  Input format choices:   \n\
                           a2m        \n\
                           afa        \n\
                           clustal    \n\
                           clustallike\n\
                           pfam       \n\
                           phylip     \n\
                           phylips    \n\
                           psiblast   \n\
                           selex      \n\
                           stockholm  \n\
\nThe defult is to guess the format\n");
exit(1);
}

// returns the length of characters required to print a number
int numLen(int n)
{
    if(n == 0) return 1;
    else return floor(log10(abs(n))) + 1;
}


typedef struct {
    uint16_t fg;
    uint16_t bg;
} ColorPair_t;

ColorPair_t custom_colors[8];

// maps characters to their corresponding colors
char color_table[] = {
    -1,   -1, -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,   -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,   -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,   -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    /*  A      B     C        D     E     F      G        H    I      J      K      L     M     N        O*/
    -1, green, blue, magenta, cyan, cyan, green, magenta, red, green, green, red, green, green, orange, yellow,
    /*  P    Q      R    S      T        U        V        W       X       Y    Z*/
    magenta, orange, red, orange, orange, magenta, green, green, orange, green, blue, -1,  -1,  -1,  -1,  -1,
    -1, red, blue, green, cyan, pink, magenta, yellow, orange, red, blue, green, cyan, pink, magenta, yellow,
    orange, red, blue, green, cyan, pink, magenta, yellow, orange, red, blue, -1, -1, -1, -1, -1
};
void init_clustalx_colors()
{
    //red
    custom_colors[0].bg = 196;
    custom_colors[0].fg = 232;
    //blue
    custom_colors[1].bg = 33;
    custom_colors[1].fg = 232;
    //green
    custom_colors[2].bg = 29;
    custom_colors[2].fg = 232;
    // cyan
    custom_colors[3].bg = 80;
    custom_colors[3].fg = 232;
    //pink
    custom_colors[4].bg = 213;
    custom_colors[4].fg = 232;
    //magenta
    custom_colors[5].bg = 165;
    custom_colors[5].fg = 232;
    //yellow
    custom_colors[6].bg = 220;
    custom_colors[6].fg = 232;
    //orange
    custom_colors[7].bg = 208;
    custom_colors[7].fg = 232;
/*
    if(!has_colors())
    {
        fprintf(stderr, "Cannot change colors\n");
    }

    // map enumerated names onto the terminal colors
    // https://lh3.googleusercontent.com/-JBl1Qa6UoBo/USXe5Wzw5uI/AAAAAAAAEeI/f0tyZjXBiyw/s800/2013-02-21--15%253A03%253A58.png
    // http://misc.flogisoft.com/_media/bash/colors_format/256-colors.sh.png
*/
}


void print_tb(const char *str, int x, int y, uint16_t fg, uint16_t bg)
{
    while (*str) {
        uint32_t uni;
        str += tb_utf8_char_to_unicode(&uni, str);
        tb_change_cell(x, y, uni, fg, bg);
        x++;
    }
}

void printf_tb(int x, int y, uint16_t fg, uint16_t bg, const char *fmt, ...)
{
    char buf[4096];
    va_list vl;
    va_start(vl, fmt);
    vsnprintf(buf, sizeof(buf), fmt, vl);
    va_end(vl);
    print_tb(buf, x, y, fg, bg);
}

void write_position(int rows, int cols, int sidebar, int start_col){
    int i, j;
    for(i = 0; i < cols; i++)
    { tb_change_cell(i, 0, ' ', 0, 255); }         // write a black bar across the top
    for(i = sidebar, j = start_col + 1; i < cols; ++i, ++j)
    {
        if(j % 10 == 0)
        {
            int n = numLen(j);
            //only print the number if it doesn't overflow the end of the screen
            if(i + n + 1 <= cols)
            {
                printf_tb(i, 0, 0, 255, "|%d", j);
            }
        }
    }
}

// compare array in qsort returning the index
int *cmp_array;
int index_cmp(const void *a, const void *b){
    int ia = *(int *)a;
    int ib = *(int *)b;
    return cmp_array[ia] > cmp_array[ib] ? -1 : cmp_array[ia] < cmp_array[ib];
}


void determine_consensus_character(ESL_MSA * msa)
{
    int col, seq_idx;
    for(col = 0; col <= msa->alen - 1; ++col)
    {
        // lookup table which stores the count of characters per column
        int consensus_lookup_table[] = {
            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /*     A  B  C  D  E  F  G  H  I  J  K  L  M  N  O*/
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        /*  P  Q  R  S  T  U  V  W  X  Y  Z*/
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0,  0,
        /*     a  b  c  d  e  f  g  h  i  j  k  l  m  n  o*/
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        /*  p  q  r  s  t  u  v  w  x  y  z*/
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        };
        for(seq_idx = 0; seq_idx <= msa->nseq - 1; ++seq_idx)
        {
            char res = msa->aseq[seq_idx][col];
            consensus_lookup_table[res]++;
        }

        // now see if any of the residues meet the requirments of the consensus rules
        int index[128];
        int i;
        for(i=0;i<128;i++){
            index[i] = i;
        }
        // copy the local table on to the globally defined pointer
        cmp_array = consensus_lookup_table;
        // qsort will sort the index based on the values in the globally defined array
        // as I've defined a custom comparison function
        qsort(index, 128, sizeof(*index), index_cmp);

        
        // see if the top residue is one that we are interested in
        // gap characters and anything that isn't an amino acid will be 0
        int consensus_res = res_lookup_table[index[0]];
        //printf("%d\t%d\t%d\t%.2f\t%c\t%#08x\n",col, consensus_lookup_table[index[0]], index[0], (float) consensus_lookup_table[index[0]] / (float) msa->nseq, (char) index[0], consensus_res);
        if(consensus_res)
        {
            int i;
            float perc = (float) consensus_lookup_table[index[0]] / (float) msa->nseq;
            //printf("\ttesting against rules:\n");
            for(i = 0; i < 29; ++i)
            {
                //printf("\t%c\t%.2f\t%#08x", consensusRules[i].name, consensusRules[i].perc, (consensus_res & consensusRules[i].residue_list));
                if(consensusRules[i].perc <= perc && (consensus_res & consensusRules[i].residue_list))
                {
                    // this rule is a winner assign this column to this rule
                    //printf("\tmatched rule %c at position %d with character %c\n", consensusRules[i].name, col, index[0]);
                    msa->rf[col] = consensusRules[i].name;
                    //break;
                }
                //printf("\n");
            }
        }
        //printf("\n");
    }
}

int main(int argc, char * argv[])
{
    //FILE * logfile = fopen("log", "w");
    ESLX_MSAFILE *afp;
    ESL_MSA      *msa;
    int opterr = 0;
    int c;
    int esl_format = eslMSAFILE_UNKNOWN; 
    while ((c = getopt (argc, argv, "hf:")) != -1)
    {
        switch (c)
        {
            case 'f':
                esl_format = esl_sqio_EncodeFormat(optarg);
                break;
            case 'h':
                usage();
                break;
            case '?':
                if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf (stderr,
                            "Unknown option character `\\x%x'.\n",
                            optopt);
                usage();
            default:
                usage();
        }
    }
    if(optind >= argc)
    {
        // no input file
        fprintf(stderr, "No input file provided\n");
        usage();
    }
    int status = eslx_msafile_Open(NULL, argv[optind], NULL, esl_format, NULL, &afp);
    if (status != eslOK) eslx_msafile_OpenFailure(afp, status);

    status = eslx_msafile_Read( afp, &msa);
    if (status != eslOK) eslx_msafile_ReadFailure(afp, status);
    eslx_msafile_Close(afp);

    init_consensus_rules();
    init_color_rules();
    // prepare the RF field for consensus information
    msa->rf = (char*) malloc(msa->alen+1);
    int y;
    for(y = 0; y < msa->alen; ++y)
    {
        msa->rf[y] = '\0';
    }
    /*int alphabet;
    esl_msa_GuessAlphabet(msa, &alphabet);
    if(msa->abc == NULL)
    {
        msa->abc = esl_alphabet_Create(alphabet);
    }

    if(msa->rf == NULL)
    {
        // calculate the consensus from the alignment, if there isn't one already
        esl_msa_ReasonableRF(msa, 0.3, 1, msa->rf); 
    }
    int z;
    for(z = 0; z < msa->alen; ++z)
    {
        printf("%c", msa->rf[z]);
    }
    printf("\n");
    exit(1);*/
    determine_consensus_character(msa);
    /*for(y = 0; y < msa->alen; ++y)
    {
        printf("%c", (msa->rf[y] == '\0') ? '.' : msa->rf[y]);
    }
    printf("\n");
    exit(1);*/
    // it's time to get to work actually displaying things  
    int phys_row, phys_col;

    /*
     * The ncurses library works by setting up one or more "windows"
     * in which to display your desired text. To pull this off, it
     * needs to be able to take over control of the terminal, thereby
     * overriding the normal i/o system (stdio or iostream)
     *
     * That being the case, we'll need to do some work initializing
     * and configuring ncurses before we can get down to the business
     * of displaying stuff.
     *
     * initscr() start_rows ncurses and takes over control of the terminal.
     * once it's been called, the window stdscr becomes available, and
     * corresponds to the size of the terminal when ncurses start_rowed
     */

    tb_init();
    tb_select_input_mode(TB_INPUT_ESC);
    tb_select_output_mode(TB_OUTPUT_256);
    init_clustalx_colors();
    phys_row = tb_height();
    phys_col = tb_width();

    // this struct holds keyboard events for us to process
    struct tb_event ev;

    /*
     * time to do the actual paging
     *
     * we'll do this by keeping track of the first line we should print
     * and then printing from there until we run out of screen or file,
     * whichever happens first
     */
    int i, j;
    unsigned int start_row = 0;		// defines first line in the window
    unsigned int start_col = 0;		// defines first column in the window
    unsigned int sidebar   = 15;     // the length of the sidebar
    unsigned int max_sidebar = phys_col * 0.2 ; // the sidebar should be at most 1/5 of the screen 
    for(j = 0; j < msa->nseq; ++j)
    {
        unsigned int sqname_len = strlen(msa->sqname[j]);
        if(sqname_len > sidebar && sqname_len <= max_sidebar)
        {
            sidebar = sqname_len;
        }
    }

    /*
     * do loops are effectively upsidedown while loops.
     * they always run at least once, which is handy
     */
    do {
        // we'll only be here if a keypress occurred, so process that first
        switch (ev.type) {
            case TB_EVENT_KEY:
            {
                switch(ev.key) {
                    case 'q':
                    case 'Q':
                    case TB_KEY_CTRL_X:
                        goto CLEANUP;
                        break;
                        // KEY_UP is defined by ncurses, and corresponds to the up arrow
                    case TB_KEY_ARROW_UP:			
                        if(start_row > 0) start_row--;	// move the start line up if not at the top
                        break;
                    case TB_KEY_ARROW_DOWN:			// move down if not at end of file
                        if(start_row < (msa->nseq - (phys_row - 1))) start_row++;
                        break;
                    case TB_KEY_ARROW_LEFT:          // move back one column if not at the begining
                        if(start_col > 0) start_col--;
                        break;
                    case TB_KEY_ARROW_RIGHT:         // move forward one column
                        if(start_col < (msa->alen - (phys_col - 1) + sidebar)) start_col++;
                        break;
                }
            }
        }

        /*
         * clear the internal buffer ready for more rendering 
         */
        tb_clear();

        //we'll loop from the starting row until we run out of screen or file

        for(i = 0; (i< phys_row - 1) && (i + start_row < msa->nseq); i++) {
            /*
             * attron and attroff turn attributes on and off, respectively.
             * pass the attribute you want to fiddle with in, and it'll get set/cleared
             */
            //attron(A_REVERSE);		// set the text mode to reverse (bg on fg)

            /*
             * ncurses provides three output functions in the style of printf
             * mvprintw takes the row to print in, the column to start at,
             * a printf-style format string, and a list of params
             */
            printf_tb(0, i+1, 0, 255, "%-*s", sidebar, msa->sqname[i + start_row]); 
            //mvprintw(i+1, 0, "%-*s", sidebar, msa->sqname[i + start_row]); 	// print line number
            //attroff(A_REVERSE);		// back to normal

            //mvprintw(i+1, sidebar, "%s", msa->aseq[i+start_row] + start_col); // print line contents
            
            int j;
            for(j = 0; j < phys_col - sidebar; j++)
            {
                char c = toupper(msa->aseq[i+start_row][start_col + j]);
                int k;
                clustalx_colors_t d = -1;
                for(k = 0; k < 21; ++k)
                {
                   if(c == colorRules[k].residue)
                   {
                       if(colorRules[k].rules_list == NULL)
                       {
                           d = colorRules[k].color;
                           //fprintf(logfile, "1: %d\t%c\t%c\t%d\t%d\t%c\n", start_col + j, c, colorRules[k].residue, d, k, msa->rf[start_col + j]);
                           break;
                       } 
                       if(msa->rf[start_col + j] == NULL)
                       {
                           //fprintf(logfile, "3: %d\t%c\t%c\t%d\t%d\t%c\n", start_col + j, c, colorRules[k].residue, d, k, msa->rf[start_col + j]);
                           break;
                       }
                       
                       if(strchr(colorRules[k].rules_list, msa->rf[start_col + j]))
                       {
                           d = colorRules[k].color;
                           //fprintf(logfile, "2: %d\t%c\t%c\t%d\t%d\t%c\n", start_col + j, c, colorRules[k].residue, d, k, msa->rf[start_col + j]);
                           break;
                       }
                   }
                }
                if(d != -1)
                {
                    tb_change_cell(j+sidebar, i+1, c, custom_colors[d].fg, custom_colors[d].bg );
                }
                else
                {
                    tb_change_cell(j+sidebar, i+1, c, 7, 232 );
                }

            }
            
        }

        // draw the status bars
        write_position(phys_row, phys_col, sidebar, start_col);


        tb_present();
    }
    while(tb_poll_event(&ev));

CLEANUP:
    tb_shutdown();
    esl_msa_Destroy(msa);
    //fclose(logfile);

    return 0;
}

//    printf("%s\n", msa->aseq[0]);
//    return 0;
//}
