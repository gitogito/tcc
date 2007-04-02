#include <stdio.h>
#include <assert.h>
#include "sim.h"
#include "mem.h"
#include "tc.h"
#include "solvele.h"

Config *config_parser;

static int dir_array[NDIRS] = { DIR_LEFT, DIR_RIGHT, DIR_FRONT, DIR_BACK, DIR_BELOW, DIR_ABOVE };
static int dir_x[] = { DIR_LEFT, DIR_RIGHT };
static int dir_y[] = { DIR_FRONT, DIR_BACK };
static int dir_z[] = { DIR_BELOW, DIR_ABOVE };

static double *double_new(double v)
{
    double *p;

    p = EALLOC(double);
    *p = v;
    return p;
}

/* Coef */

static Coef get_coef(void)
{
    Coef coef;

    coef.index = -1;
    coef.value = 0.0;
    return coef;
}

/* Coefs */

static Coefs get_coefs(void)
{
    Coefs coefs;
    int idir, dir;

    for (idir = 0; idir < NELEMS(dir_array); ++idir) {
	dir = dir_array[idir];
	coefs.coef[dir] = get_coef();
    }
    coefs.coef0 = 0.0;
    coefs.cnst = 0.0;
    return coefs;
}

/* World */

static int world_to_index(World *self, iPoint ipoint)
{
    return ipoint.i + self->ni * (ipoint.j + self->nj * ipoint.k);
}

World *world_new(double x0, double y0, double z0,
	double xlen, double ylen, double zlen,
	double dx, double dy, double dz)
{
    World *self;

    if (xlen <= 0.0 || ylen <= 0.0 || zlen <= 0.0)
	warn_exit("length is negative for World");

    self = EALLOC(World);
    self->x0 = x0;
    self->y0 = y0;
    self->z0 = z0;
    self->xlen = xlen;
    self->ylen = ylen;
    self->zlen = zlen;
    self->dx = dx;
    self->dy = dy;
    self->dz = dz;
    self->ni = iround(xlen / dx) + 1;
    self->nj = iround(ylen / dy) + 1;
    self->nk = iround(zlen / dz) + 1;
    self->each = NULL;
    return self;
}

static iPoint *world_each(World *self)
{
    iPoint_ary *ipoint_ary;
    int i, j, k;

    if (self->each != NULL && self->each->index >= 0)
	return each_each(self->each);

    ipoint_ary = ipoint_ary_new();
    for (i = 0; i < self->ni; ++i) {
	for (j = 0; j < self->nj; ++j) {
	    for (k = 0; k < self->nk; ++k) {
		ipoint_ary_push(ipoint_ary, get_ipoint(i, j, k));
	    }
	}
    }
    self->each = each_new(ipoint_ary);
    return each_each(self->each);
}

static int world_inside_p(World *self, iPoint ipoint)
{
    if (ipoint.i < 0 || ipoint.i >= self->ni ||
	    ipoint.j < 0 || ipoint.j >= self->nj ||
	    ipoint.k < 0 || ipoint.k >= self->nk)
	return 0;
    else
	return 1;
}

/* Config */

static void config_parse(Config *self, char *fname)
{
    extern int yyparse();
    extern int yydebug;
    extern FILE *yyin;

    yydebug = opt_y;

    config_parser = self;
    if (fname == NULL) {
	yyin = stdin;
    } else {
	yyin = fopen(fname, "r");
	if (yyin == NULL)
	    warn_exit("can't open '%s'", fname);
    }
    yyparse();
    if (yyin != stdin)
	fclose(yyin);
}

static Config *config_new(char *fname)
{
    Config *self;

    self = EALLOC(Config);
    self->world = NULL;
    self->active_obj_ary = aryobj_new();
    self->fix_obj_ary = aryobj_new();
    self->heat_obj_ary = aryobj_new();
    self->lambda_obj_ary = aryobj_new();

    config_parse(self, fname);

    return self;
}

/* Sim */

int sim_active_p(Sim *self, iPoint ipoint)
{
    if (ipoint.i < 0 || ipoint.i >= self->world->ni ||
	    ipoint.j < 0 || ipoint.j >= self->world->nj ||
	    ipoint.k < 0 || ipoint.k >= self->world->nk)
	return 0;
    else
	return self->active_p_ary[ipoint.i][ipoint.j][ipoint.k];
}

static void sim_set_region_active(Sim *self)
{
    iPoint *p;
    AryObj *active_obj_ary;
    Obj *obj;
    int index;

    ALLOCATE_3D2(self->active_p_ary, int, self->world->ni, self->world->nj, self->world->nk, 0);
    active_obj_ary = self->config->active_obj_ary;
    for (index = 0; index < active_obj_ary->size; ++index) {
	obj = active_obj_ary->ptr[index];
	for (p = obj_each(obj); p != NULL; p = obj_each(obj)) {
	    if (!world_inside_p(self->world, *p))
		continue;
	    (self->active_p_ary)[p->i][p->j][p->k] = obj->uval.i;
	}
    }
}

static void sim_set_region_fix(Sim *self)
{
    iPoint *p;
    AryObj *obj_ary;
    Obj *obj;
    int i;

    ALLOCATE_3D2(self->fix_ary, double *, self->world->ni, self->world->nj, self->world->nk, NULL);
    obj_ary = self->config->fix_obj_ary;
    for (i = 0; i < obj_ary->size; ++i) {
	obj = obj_ary->ptr[i];
	for (p = obj_each(obj); p != NULL; p = obj_each(obj)) {
	    if (!world_inside_p(self->world, *p))
		continue;
	    self->fix_ary[p->i][p->j][p->k] = double_new(obj->uval.d);
	}
    }
}

static void sim_set_region_heat(Sim *self)
{
    iPoint *p;
    AryObj *obj_ary;
    int index;
    Obj *obj;
    int first;
    iPoint ipoint0;

    ALLOCATE_3D2(self->heat_ary, double *, self->world->ni, self->world->nj, self->world->nk, NULL);
    ALLOCATE_3D2(self->heat_ipoint_ary, iPoint *, self->world->ni, self->world->nj, self->world->nk, NULL);
    obj_ary = self->config->heat_obj_ary;
    for (index = 0; index < obj_ary->size; ++index) {
	obj = obj_ary->ptr[index];
	first = 1;
	for (p = obj_each(obj); p != NULL; p = obj_each(obj)) {
	    if (!world_inside_p(self->world, *p))
		continue;
	    if (first) {
		first = 0;
		ipoint0 = *p;
	    }
	    self->heat_ipoint_ary[p->i][p->j][p->k] = ipoint_new(ipoint0.i, ipoint0.j, ipoint0.k);
	    self->heat_ary[p->i][p->j][p->k] = double_new(obj->uval.d);
	}
    }
}

static void sim_set_region_lambda(Sim *self)
{
    int index;
    AryObj *obj_ary;
    Obj *obj;
    iPoint *p;

    ALLOCATE_3D2(self->lambda_ary, double, self->world->ni, self->world->nj, self->world->nk, 1.0);
    obj_ary = self->config->lambda_obj_ary;
    for (index = 0; index < obj_ary->size; ++index) {
	obj = obj_ary->ptr[index];
	for (p = obj_each(obj); p != NULL; p = obj_each(obj)) {
	    if (!world_inside_p(self->world, *p))
		continue;
	    self->lambda_ary[p->i][p->j][p->k] = obj->uval.d;
	}
    }
}

static void sim_set_region(Sim *self)
{
    self->world = self->config->world;
    sim_set_region_active(self);
    sim_set_region_fix(self);
    sim_set_region_heat(self);
    sim_set_region_lambda(self);
}

static void sim_set_matrix_const(Sim *self, int i, int j, int k)
{
    if (self->heat_ary[i][j][k] != NULL &&
	    ipoint_eq(*(self->heat_ipoint_ary[i][j][k]), get_ipoint(i, j, k)))
	self->coefs[i][j][k].cnst = *(self->heat_ary[i][j][k]);
}

static void sim_set_matrix_coef0(Sim *self, int i, int j, int k, double dx, double dy, double dz)
{
    int idir, dir;
    int ix, iy, iz;
    int dirx, diry, dirz;
    iPoint pp, ipoint_l;
    double l;

    for (idir = 0; idir < NELEMS(dir_array); ++idir) {
	dir = dir_array[idir];
	if (!sim_active_p(self, ipoint_add(get_ipoint(i, j, k), self->dir_to_ipoint[dir])))
	    continue;
	switch (dir) {
	case DIR_LEFT: case DIR_RIGHT:
	    for (iy = 0; iy < NELEMS(dir_y); ++iy) {
		diry = dir_y[iy];
		for (iz = 0; iz < NELEMS(dir_z); ++iz) {
		    dirz = dir_z[iz];
		    pp = ipoint_add(get_ipoint(i, j, k), self->dir_to_ipoint[diry]);
		    pp = ipoint_add(pp, self->dir_to_ipoint[dirz]);
		    if (sim_active_p(self, pp)) {
			ipoint_l = ipoint_offset(get_ipoint(i, j, k), dir, diry, dirz);
			l = self->lambda_ary[ipoint_l.i][ipoint_l.j][ipoint_l.k];
			self->coefs[i][j][k].coef0 += -l/dx*(dy/2)*(dz/2);
		    }
		}
	    }
	    break;
	case DIR_FRONT: case DIR_BACK:
	    for (iz = 0; iz < NELEMS(dir_z); ++iz) {
		dirz = dir_z[iz];
		for (ix = 0; ix < NELEMS(dir_x); ++ix) {
		    dirx = dir_x[ix];
		    pp = ipoint_add(get_ipoint(i, j, k), self->dir_to_ipoint[dirz]);
		    pp = ipoint_add(pp, self->dir_to_ipoint[dirx]);
		    if (sim_active_p(self, pp)) {
			ipoint_l = ipoint_offset(get_ipoint(i, j, k), dirx, dir, dirz);
			l = self->lambda_ary[ipoint_l.i][ipoint_l.j][ipoint_l.k];
			self->coefs[i][j][k].coef0 += -l/dy*(dz/2)*(dx/2);
		    }
		}
	    }
	    break;
	case DIR_BELOW: case DIR_ABOVE:
	    for (ix = 0; ix < NELEMS(dir_x); ++ix) {
		dirx = dir_x[ix];
		for (iy = 0; iy < NELEMS(dir_y); ++iy) {
		    diry = dir_y[iy];
		    pp = ipoint_add(get_ipoint(i, j, k), self->dir_to_ipoint[dirx]);
		    pp = ipoint_add(pp, self->dir_to_ipoint[diry]);
		    if (sim_active_p(self, pp)) {
			ipoint_l = ipoint_offset(get_ipoint(i, j, k), dirx, diry, dir);
			l = self->lambda_ary[ipoint_l.i][ipoint_l.j][ipoint_l.k];
			self->coefs[i][j][k].coef0 += -l/dz*(dx/2)*(dy/2);
		    }
		}
	    }
	    break;
	default:
	    bug("sim_set_matrix_coef0");
	}
    }
}

static void sim_set_matrix_coef(Sim *self, int i, int j, int k, double dx, double dy, double dz)
{
    int idir, dir;
    int dirx, diry, dirz;
    double value;
    int ix, iy, iz;
    iPoint pp, ipoint_l;
    double l;

    for (idir = 0; idir < NELEMS(dir_array); ++idir) {
	dir = dir_array[idir];
	if (!sim_active_p(self, ipoint_add(get_ipoint(i, j, k), self->dir_to_ipoint[dir]))) {
	    self->coefs[i][j][k].coef[dir].index = -1;
	    continue;
	}

	self->coefs[i][j][k].coef[dir].index =
	    world_to_index(self->world, ipoint_add(get_ipoint(i, j, k), self->dir_to_ipoint[dir]));
	value = 0.0;
	switch (dir) {
	case DIR_LEFT: case DIR_RIGHT:
	    for (iy = 0; iy < NELEMS(dir_y); ++iy) {
		diry = dir_y[iy];
		for (iz = 0; iz < NELEMS(dir_z); ++iz) {
		    dirz = dir_z[iz];
		    pp = ipoint_add(get_ipoint(i, j, k), self->dir_to_ipoint[diry]);
		    pp = ipoint_add(pp, self->dir_to_ipoint[dirz]);
		    if (sim_active_p(self, pp)) {
			ipoint_l = ipoint_offset(get_ipoint(i, j, k), dir, diry, dirz);
			l = self->lambda_ary[ipoint_l.i][ipoint_l.j][ipoint_l.k];
			value += l/dx*(dy/2)*(dz/2);
		    }
		}
	    }
	    break;
	case DIR_FRONT: case DIR_BACK:
	    for (iz = 0; iz < NELEMS(dir_z); ++iz) {
		dirz = dir_z[iz];
		for (ix = 0; ix < NELEMS(dir_x); ++ix) {
		    dirx = dir_x[ix];
		    pp = ipoint_add(get_ipoint(i, j, k), self->dir_to_ipoint[dirz]);
		    pp = ipoint_add(pp, self->dir_to_ipoint[dirx]);
		    if (sim_active_p(self, pp)) {
			ipoint_l = ipoint_offset(get_ipoint(i, j, k), dirx, dir, dirz);
			l = self->lambda_ary[ipoint_l.i][ipoint_l.j][ipoint_l.k];
			value += l/dy*(dz/2)*(dx/2);
		    }
		}
	    }
	    break;
	case DIR_BELOW: case DIR_ABOVE:
	    for (ix = 0; ix < NELEMS(dir_x); ++ix) {
		dirx = dir_x[ix];
		for (iy = 0; iy < NELEMS(dir_y); ++iy) {
		    diry = dir_y[iy];
		    pp = ipoint_add(get_ipoint(i, j, k), self->dir_to_ipoint[dirx]);
		    pp = ipoint_add(pp, self->dir_to_ipoint[diry]);
		    if (sim_active_p(self, pp)) {
			ipoint_l = ipoint_offset(get_ipoint(i, j, k), dirx, diry, dir);
			l = self->lambda_ary[ipoint_l.i][ipoint_l.j][ipoint_l.k];
			value += l/dz*(dx/2)*(dy/2);
		    }
		}
	    }
	    break;
	default:
	    bug("sim_set_matrix_coef");
	}
	self->coefs[i][j][k].coef[dir].value = value;
    }
}

static void sim_set_matrix(Sim *self)
{
    iPoint *p;
    double dx, dy, dz;

    ALLOCATE_3D(self->coefs, Coefs, self->world->ni, self->world->nj, self->world->nk);
    for (p = world_each(self->world); p != NULL; p = world_each(self->world)) {
	self->coefs[p->i][p->j][p->k] = get_coefs();
    }

    dx = self->world->dx;
    dy = self->world->dy;
    dz = self->world->dz;

    for (p = world_each(self->world); p != NULL; p = world_each(self->world)) {
	if (self->fix_ary[p->i][p->j][p->k] != NULL)
	    continue;

	sim_set_matrix_const(self, p->i, p->j, p->k);
	sim_set_matrix_coef0(self, p->i, p->j, p->k, dx, dy, dz);
	sim_set_matrix_coef(self, p->i, p->j, p->k, dx, dy, dz);
    }
}

Sim *sim_new(char *fname)
{
    Sim *self;

    self = EALLOC(Sim);
    if (opt_v)
	warn("configuring ...");
    self->config = config_new(fname);
    self->dir_to_ipoint[DIR_LEFT]  = get_ipoint(-1,  0,  0);
    self->dir_to_ipoint[DIR_RIGHT] = get_ipoint( 1,  0,  0);
    self->dir_to_ipoint[DIR_FRONT] = get_ipoint( 0, -1,  0);
    self->dir_to_ipoint[DIR_BACK]  = get_ipoint( 0,  1,  0);
    self->dir_to_ipoint[DIR_BELOW] = get_ipoint( 0,  0, -1);
    self->dir_to_ipoint[DIR_ABOVE] = get_ipoint( 0,  0,  1);

    if (opt_v)
	warn("setting region ...");
    sim_set_region(self);
    if (opt_v)
	warn("setting matrix ...");
    sim_set_matrix(self);

    return self;
}

Array3Dd sim_calc(Sim *self)
{
    int nindex;
    Solvele *solver;
    iPoint *p, *hfp;
    iPoint ipoint;
    int index, index2;
    int idir, dir;
    double c;
    double *sol;
    Array3Dd ary;

    nindex = self->world->ni * self->world->nj * self->world->nk;
    solver = solvele_new(nindex);
    for (p = world_each(self->world); p != NULL; p = world_each(self->world)) {
	if (self->fix_ary[p->i][p->j][p->k] != NULL) {
	    index = world_to_index(self->world, *p);
	    solvele_set_matrix(solver, index, index, 1.0);
	    solvele_set_vector(solver, index, *(self->fix_ary[p->i][p->j][p->k]));
	    continue;
	} else if (!self->active_p_ary[p->i][p->j][p->k]) {
	    index = world_to_index(self->world, *p);
	    solvele_set_matrix(solver, index, index, 1.0);
	    solvele_set_vector(solver, index, 0.0);
	    continue;
	}
	hfp = self->heat_ipoint_ary[p->i][p->j][p->k];
	if (hfp != NULL) {
	    if (!ipoint_eq(*p, *hfp)) {
		index = world_to_index(self->world, *p);
		solvele_set_matrix(solver, index, index, 1.0);
		solvele_set_matrix(solver, index, world_to_index(self->world, *hfp), -1.0);
		solvele_set_vector(solver, index, 0.0);
	    }
	    ipoint = *hfp;
	} else {
	    ipoint = *p;
	}
	index = world_to_index(self->world, ipoint);
	for (idir = 0; idir < NELEMS(dir_array); ++idir) {
	    dir = dir_array[idir];
	    index2 = self->coefs[p->i][p->j][p->k].coef[dir].index;
	    if (index2 < 0)
		continue;
	    c = self->coefs[p->i][p->j][p->k].coef[dir].value;
	    solvele_add_matrix(solver, index, index2, c);
	}
	solvele_add_matrix(solver, index, index, self->coefs[p->i][p->j][p->k].coef0);
	solvele_add_vector(solver, index, -self->coefs[p->i][p->j][p->k].cnst);
    }

    if (opt_v)
	warn("solving equations ...");
    sol = solvele_solve(solver, self->world->ni, self->world->nj, self->world->nk);

    ALLOCATE_3D(ary, double, self->world->ni, self->world->nj, self->world->nk);
    for (p = world_each(self->world); p != NULL; p = world_each(self->world)) {
	ary[p->i][p->j][p->k] = sol[world_to_index(self->world, *p)];
    }

    return ary;
}
