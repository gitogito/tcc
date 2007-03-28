#ifndef _TC_H_
#define _TC_H_

extern char *prgname;
extern int opt_e;	/* specify SOR epsilon */
extern int opt_o;	/* specify SOR omega */
extern int opt_u;	/* use UMFPACK */
extern int opt_v;	/* verbose */
extern int opt_y;	/* yydebug */

void bug(char *fmt, ...);
void warn(char *fmt, ...);
void warn_exit(char *fmt, ...);

#endif
