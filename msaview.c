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

#include "termbox/include/termbox.h"

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
/*
void write_position(int rows, int cols, int sidebar, int start_col){
    attron(A_REVERSE);			// set the text mode to reverse (bg on fg)
    move(0,0);          		// move to the upper-left corner
    int i, j;
    for(i = 0; i < cols; i++)
    { addch(' '); }			// write a black bar across the top
    for(i = sidebar, j = start_col + 1; i < cols; ++i, ++j)
    {
        if(j % 10 == 0)
        {
            //only print the number if it doesn't overflow the end of the screen
            if(i + numLen(j) + 1 <= cols)
            {
                mvprintw(0,i, "|%d", j);
            }
        }
    }
    attroff(A_REVERSE);			// switch back to white on black
}
*/
// the default terminal colors go from 1 to 16. I don't want to mess with them
// so I'm going to start my colors at 17.
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
    -1, red, blue, green, cyan, pink, magenta, yellow, orange, red, blue, green, cyan, pink, magenta, yellow,
    orange, red, blue, green, cyan, pink, magenta, yellow, orange, red, blue, -1,  -1,  -1,  -1,  -1,
    -1, red, blue, green, cyan, pink, magenta, yellow, orange, red, blue, green, cyan, pink, magenta, yellow,
    orange, red, blue, green, cyan, pink, magenta, yellow, orange, red, blue, -1, -1, -1, -1, -1
};
void init_clustalx_colors()
{
    
    custom_colors[0].fg = 196;
    custom_colors[0].bg = 232;
    custom_colors[1].fg = 21;
    custom_colors[1].bg = 232;
    custom_colors[2].fg = 22;
    custom_colors[2].bg = 232;
    custom_colors[3].fg = 80;
    custom_colors[3].bg = 232;
    custom_colors[4].fg = 213;
    custom_colors[4].bg = 232;
    custom_colors[5].fg = 165;
    custom_colors[5].bg = 232;
    custom_colors[6].fg = 220;
    custom_colors[6].bg = 232;
    custom_colors[7].fg = 172;
    custom_colors[7].bg = 232;
/*
    if(!has_colors())
    {
        fprintf(stderr, "Cannot change colors\n");
    }

    // map enumerated names onto the terminal colors
    // https://lh3.googleusercontent.com/-JBl1Qa6UoBo/USXe5Wzw5uI/AAAAAAAAEeI/f0tyZjXBiyw/s800/2013-02-21--15%253A03%253A58.png
    // http://misc.flogisoft.com/_media/bash/colors_format/256-colors.sh.png
    init_pair(red, -1, 196);
    init_pair(blue, -1, 21);
    init_pair(green, -1, 22);
    init_pair(cyan, 242, 80);
    init_pair(pink, 242, 213);
    init_pair(magenta, -1, 165);
    init_pair(yellow, 242, 220);
    init_pair(orange, -1, 172);
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

int main(int argc, char * argv[])
{
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
    tb_select_input_mode(TB_INPUT_ESC | TB_INPUT_MOUSE);
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
                char c = msa->aseq[i+start_row][start_col + j];
                char d = color_table[c];
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

    return 0;
}

//    printf("%s\n", msa->aseq[0]);
//    return 0;
//}
