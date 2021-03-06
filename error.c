/* ------- file: -------------------------- Error.c -----------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Dec 23 14:47:16 2008 --

       --------------------------                      ----------RH-- */

/* --- Print warning, or print error string and exit when error level
       is above preset treshold --                     -------------- */

#include <stdlib.h>
#include <errno.h>

#include "rh.h"
#include "error.h"
#include "inputs.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

extern CommandLine commandline;

/* ------- begin -------------------------- Error.c ----------------- */

void Error(enum errorlevel level, const char *routineName,
           const char *messageStr) {
  char errorStr[MAX_MESSAGE_LENGTH];
  enum errorlevel defaultLevel = ERROR_LEVEL_1;

  switch (level) {
  case MESSAGE:
    if (!commandline.quiet)
      fprintf(commandline.logfile, "%s", (messageStr) ? messageStr : "");
      fflush(commandline.logfile);
    return;
  case WARNING:
    fprintf(commandline.logfile, "\n-WARNING in routine %s\n %s\n", routineName,
            (messageStr) ? messageStr : " (Undocumented)\n");
    fflush(commandline.logfile);
    return;
  default:
    if (level < defaultLevel) {
      fprintf(commandline.logfile, "\a\n-ERROR in routine %s\n %s \n %s\n",
              routineName, (messageStr) ? messageStr : " (Undocumented)\n",
              "Trying to continue.....");
      fflush(commandline.logfile);
      return;
    } else {
      sprintf(errorStr, "\a\n\n-FATAL_ERROR in routine %s\n %s \n %s\n",
              routineName, (messageStr) ? messageStr : " (Undocumented)\n",
              "Exiting.....");

      fprintf(commandline.logfile, "%s", errorStr);
      fflush(commandline.logfile);
      if (commandline.logfile != stderr)
        fprintf(stderr, "%s", errorStr);
        fflush(stderr);

      if (errno)
        perror(routineName);
      exit(level);
    }
  }
}
/* ------- end ---------------------------- Error.c ----------------- */
