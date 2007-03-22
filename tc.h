#ifndef _TC_H
#define _TC_H

extern char *prgname;
extern int opt_e;	/* specify SOR epsilon */
extern int opt_o;	/* specify SOR omega */
extern int opt_u;	/* use UMFPACK */
extern int opt_v;	/* verbose */

void bug(char *fmt, ...);
void warn(char *fmt, ...);
void warn_exit(char *fmt, ...);

#endif
