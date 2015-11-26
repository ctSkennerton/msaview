#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <getopt.h>
#include <ncurses.h>
#include <math.h>

#include "easel/include/esl_config.h"
#include "easel/include/easel.h"
#include "easel/include/esl_msa.h"
#include "easel/include/esl_msafile.h"
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

    initscr();

    /*
     * we don't know the size of stdscr as yet, but that would probably
     * be good to know. getmaxyx takes the window to size and variables
     * in which to place its dimensions
     */
    getmaxyx(stdscr, phys_row, phys_col);

    /*
     * normally, calls to get input will wait until they see a line break
     * (exactly as stdio or iostream would work). switching on raw
     * mode means we get data as soon as it's available.
     */
    raw();

    /*
     * switching on keypad mode means that we'll receive things like
     * arrow and function keypresses, in addition to normal text data
     */
    keypad(stdscr, TRUE);

    /*
     * by default, ncurses will "echo" keypresses - that is, print the
     * character on screen. we don't want that, so we'll shut echo off
     */
    noecho();

    /*
     * time to do the actual paging
     *
     * we'll do this by keeping track of the first line we should print
     * and then printing from there until we run out of screen or file,
     * whichever happens first
     */
    int i, ch = 0;
    unsigned int start_row = 0;		// defines first line in the window
    unsigned int start_col = 0;		// defines first column in the window
    unsigned int sidebar   = 15;     // the length of the sidebar
    unsigned int max_sidebar = phys_col * 0.2 ; // the sidebar should be at most 1/5 of the screen 
    int j;
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
        switch (ch) {
            // KEY_UP is defined by ncurses, and corresponds to the up arrow
            case KEY_UP:			
                if(start_row > 0) start_row--;	// move the start line up if not at the top
                break;
            case KEY_DOWN:			// move down if not at end of file
                if(start_row < (msa->nseq - (phys_row - 1))) start_row++;
                break;
            case KEY_LEFT:          // move back one column if not at the begining
                if(start_col > 0) start_col--;
                break;
            case KEY_RIGHT:         // move forward one column
                if(start_col < (msa->alen - (phys_col - 1) + sidebar)) start_col++;
                break;
        }

        /*
         * for the most part, ncurses behaves the way you'd expect a graphics
         * library to work. the window retains contents until cleared, so we
         * should clear it before we start laying down new data
         */
        clear();

        //we'll loop from the starting row until we run out of screen or file

        for(i = 0; (i< phys_row - 1) && (i + start_row < msa->nseq); i++) {
            /*
             * attron and attroff turn attributes on and off, respectively.
             * pass the attribute you want to fiddle with in, and it'll get set/cleared
             */
            attron(A_REVERSE);		// set the text mode to reverse (bg on fg)

            /*
             * ncurses provides three output functions in the style of printf
             * mvprintw takes the row to print in, the column to start at,
             * a printf-style format string, and a list of params
             */
            mvprintw(i+1, 0, "%-*s", sidebar, msa->sqname[i + start_row]); 	// print line number
            attroff(A_REVERSE);		// back to normal

            mvprintw(i+1, sidebar, "%s", msa->aseq[i+start_row] + start_col); // print line contents
        }

        // draw the status bars
        write_position(phys_row, phys_col, sidebar, start_col);
        //write_status(phys_row, phys_col, argv[1], start_row + 1, start_row + i);


        /*
         * you can think of ncurses as though it double-buffers. changes
         * are made to a buffer behind the scenes and then written, all
         * at once, to screen when refresh() is called
         */
        refresh();
    }
    while((ch = getch()) != 'q');		// continue until user hits q

    /*
     * as with anything you have to explicitely set up, you have to
     * tell ncurses to tear itself down before your program exits.
     * exitwin() does this.
     */
    endwin();

    return 0;
}

//    printf("%s\n", msa->aseq[0]);
//    return 0;
//}
