#include <stdio.h>
#include <assert.h>
#include "sim.h"
#include "mem.h"
#include "tc.h"
#include "solvele.h"

Config *config;
World *world;
Sim *sim;

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

static int world_to_index(iPoint *ipoint)
{
    return ipoint->i + world->ni * (ipoint->j + world->nj * ipoint->k);
}

World *world_new(double x0, double y0, double z0,
	double xlen, double ylen, double zlen,
	double dx, double dy, double dz)
{
    if (xlen <= 0.0 || ylen <= 0.0 || zlen <= 0.0)
	warn_exit("length is negative for World");

    world = EALLOC(World);
    world->x0 = x0;
    world->y0 = y0;
    world->z0 = z0;
    world->xlen = xlen;
    world->ylen = ylen;
    world->zlen = zlen;
    world->dx = dx;
    world->dy = dy;
    world->dz = dz;
    world->ni = iround(xlen / dx) + 1;
    world->nj = iround(ylen / dy) + 1;
    world->nk = iround(zlen / dz) + 1;
    world->each = NULL;
    return world;
}

static int world_each(iPoint **pp)
{
    iPoint_ary *ipoint_ary;
    IP_TYPE i, j, k;

    if (world->each != NULL && world->each->index >= 0)
	return each_each(world->each, pp);

    ipoint_ary = ipoint_ary_new();
    for (i = 0; i < world->ni; ++i) {
	for (j = 0; j < world->nj; ++j) {
	    for (k = 0; k < world->nk; ++k) {
		ipoint_ary_push(ipoint_ary, get_ipoint(i, j, k));
	    }
	}
    }
    world->each = each_new(ipoint_ary);
    return each_each(world->each, pp);
}

int world_inside_p(iPoint *ipoint)
{
    if (ipoint->i < 0 || ipoint->i >= world->ni ||
	    ipoint->j < 0 || ipoint->j >= world->nj ||
	    ipoint->k < 0 || ipoint->k >= world->nk)
	return 0;
    else
	return 1;
}

/* Config */

static void config_parse(char *fname)
{
    extern int yyparse();
    extern int yydebug;
    extern FILE *yyin;

    yydebug = opt_y;

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

static Config *config_new(void)
{
    Config *self;

    self = EALLOC(Config);
    self->active_obj_ary = aryobj_new();
    self->fix_obj_ary = aryobj_new();
    self->heat_obj_ary = aryobj_new();
    self->lambda_obj_ary = aryobj_new();
    return self;
}

/* Sim */

int sim_active_p(iPoint *ipoint)
{
    if (ipoint->i < 0 || ipoint->i >= world->ni ||
	    ipoint->j < 0 || ipoint->j >= world->nj ||
	    ipoint->k < 0 || ipoint->k >= world->nk)
	return 0;
    else
	return sim->active_p_ary[ipoint->i][ipoint->j][ipoint->k];
}

static void sim_set_region_active(void)
{
    iPoint *p;
    AryObj *active_obj_ary;
    Obj *obj;
    int index;

    ALLOCATE_3D2(sim->active_p_ary, int, world->ni, world->nj, world->nk, 0);
    active_obj_ary = config->active_obj_ary;
    for (index = 0; index < active_obj_ary->size; ++index) {
	obj = active_obj_ary->ptr[index];
	while (obj_each(obj, &p)) {
	    if (p == NULL)
		continue;
	    if (!world_inside_p(p))
		continue;
	    (sim->active_p_ary)[p->i][p->j][p->k] = obj->uval.i;
	}
    }
}

static void sim_set_region_fix(void)
{
    iPoint *p;
    AryObj *obj_ary;
    Obj *obj;
    int i;

    ALLOCATE_3D2(sim->fix_ary, double *, world->ni, world->nj, world->nk, NULL);
    obj_ary = config->fix_obj_ary;
    for (i = 0; i < obj_ary->size; ++i) {
	obj = obj_ary->ptr[i];
	while (obj_each(obj, &p)) {
	    if (p == NULL)
		continue;
	    if (!world_inside_p(p))
		continue;
	    sim->fix_ary[p->i][p->j][p->k] = double_new(obj->uval.d);
	}
    }
}

static void sim_set_region_heat(void)
{
    iPoint *p;
    AryObj *obj_ary;
    int index;
    Obj *obj;
    int first;
    iPoint ipoint0;

    ALLOCATE_3D2(sim->heat_ary, double *, world->ni, world->nj, world->nk, NULL);
    ALLOCATE_3D2(sim->heat_ipoint_ary, iPoint *, world->ni, world->nj, world->nk, NULL);
    obj_ary = config->heat_obj_ary;
    for (index = 0; index < obj_ary->size; ++index) {
	obj = obj_ary->ptr[index];
	first = 1;
	while (obj_each(obj, &p)) {
	    if (p == NULL)
		continue;
	    if (!sim_active_p(p))
		continue;
	    if (first) {
		first = 0;
		ipoint0 = *p;
	    }
	    sim->heat_ipoint_ary[p->i][p->j][p->k] = ipoint_new(ipoint0.i, ipoint0.j, ipoint0.k);
	    sim->heat_ary[p->i][p->j][p->k] = double_new(obj->uval.d);
	}
    }
}

static void sim_set_region_lambda(void)
{
    int index;
    AryObj *obj_ary;
    Obj *obj;
    iPoint *p;

    ALLOCATE_3D2(sim->lambda_ary, double, world->ni, world->nj, world->nk, 1.0);
    obj_ary = config->lambda_obj_ary;
    for (index = 0; index < obj_ary->size; ++index) {
	obj = obj_ary->ptr[index];
	while (obj_each(obj, &p)) {
	    if (p == NULL)
		continue;
	    if (!world_inside_p(p))
		continue;
	    sim->lambda_ary[p->i][p->j][p->k] = obj->uval.d;
	}
    }
}

static void sim_set_region(void)
{
    sim_set_region_active();
    sim_set_region_fix();
    sim_set_region_heat();
    sim_set_region_lambda();
}

static void sim_set_matrix_const(iPoint *p)
{
    if (sim->heat_ary[p->i][p->j][p->k] != NULL &&
	    ipoint_eq(sim->heat_ipoint_ary[p->i][p->j][p->k], p))
	sim->coefs[p->i][p->j][p->k].cnst = *(sim->heat_ary[p->i][p->j][p->k]);
}

static void sim_set_matrix_coef0(iPoint *p, double dx, double dy, double dz)
{
    int idir, dir;
    IP_TYPE ix, iy, iz;
    int dirx, diry, dirz;
    iPoint pp, ipoint_l;
    double l;

    for (idir = 0; idir < NELEMS(dir_array); ++idir) {
	dir = dir_array[idir];
	pp = ipoint_add(p, &sim->dir_to_ipoint[dir]);
	if (!sim_active_p(&pp))
	    continue;
	switch (dir) {
	case DIR_LEFT: case DIR_RIGHT:
	    for (iy = 0; iy < NELEMS(dir_y); ++iy) {
		diry = dir_y[iy];
		for (iz = 0; iz < NELEMS(dir_z); ++iz) {
		    dirz = dir_z[iz];
		    pp = ipoint_add(p, &sim->dir_to_ipoint[diry]);
		    pp = ipoint_add(&pp, &sim->dir_to_ipoint[dirz]);
		    if (sim_active_p(&pp)) {
			ipoint_l = ipoint_offset(p, dir, diry, dirz);
			l = sim->lambda_ary[ipoint_l.i][ipoint_l.j][ipoint_l.k];
			sim->coefs[p->i][p->j][p->k].coef0 += -l/dx*(dy/2)*(dz/2);
		    }
		}
	    }
	    break;
	case DIR_FRONT: case DIR_BACK:
	    for (iz = 0; iz < NELEMS(dir_z); ++iz) {
		dirz = dir_z[iz];
		for (ix = 0; ix < NELEMS(dir_x); ++ix) {
		    dirx = dir_x[ix];
		    pp = ipoint_add(p, &sim->dir_to_ipoint[dirz]);
		    pp = ipoint_add(&pp, &sim->dir_to_ipoint[dirx]);
		    if (sim_active_p(&pp)) {
			ipoint_l = ipoint_offset(p, dirx, dir, dirz);
			l = sim->lambda_ary[ipoint_l.i][ipoint_l.j][ipoint_l.k];
			sim->coefs[p->i][p->j][p->k].coef0 += -l/dy*(dz/2)*(dx/2);
		    }
		}
	    }
	    break;
	case DIR_BELOW: case DIR_ABOVE:
	    for (ix = 0; ix < NELEMS(dir_x); ++ix) {
		dirx = dir_x[ix];
		for (iy = 0; iy < NELEMS(dir_y); ++iy) {
		    diry = dir_y[iy];
		    pp = ipoint_add(p, &sim->dir_to_ipoint[dirx]);
		    pp = ipoint_add(&pp, &sim->dir_to_ipoint[diry]);
		    if (sim_active_p(&pp)) {
			ipoint_l = ipoint_offset(p, dirx, diry, dir);
			l = sim->lambda_ary[ipoint_l.i][ipoint_l.j][ipoint_l.k];
			sim->coefs[p->i][p->j][p->k].coef0 += -l/dz*(dx/2)*(dy/2);
		    }
		}
	    }
	    break;
	default:
	    bug("sim_set_matrix_coef0");
	}
    }
}

static void sim_set_matrix_coef(iPoint *p, double dx, double dy, double dz)
{
    int idir, dir;
    int dirx, diry, dirz;
    double value;
    IP_TYPE ix, iy, iz;
    iPoint pp, ipoint_l;
    double l;

    for (idir = 0; idir < NELEMS(dir_array); ++idir) {
	dir = dir_array[idir];
	pp = ipoint_add(p, &sim->dir_to_ipoint[dir]);
	if (!sim_active_p(&pp)) {
	    sim->coefs[p->i][p->j][p->k].coef[dir].index = -1;
	    continue;
	}

	pp = ipoint_add(p, &sim->dir_to_ipoint[dir]);
	sim->coefs[p->i][p->j][p->k].coef[dir].index = world_to_index(&pp);
	value = 0.0;
	switch (dir) {
	case DIR_LEFT: case DIR_RIGHT:
	    for (iy = 0; iy < NELEMS(dir_y); ++iy) {
		diry = dir_y[iy];
		for (iz = 0; iz < NELEMS(dir_z); ++iz) {
		    dirz = dir_z[iz];
		    pp = ipoint_add(p, &sim->dir_to_ipoint[diry]);
		    pp = ipoint_add(&pp, &sim->dir_to_ipoint[dirz]);
		    if (sim_active_p(&pp)) {
			ipoint_l = ipoint_offset(p, dir, diry, dirz);
			l = sim->lambda_ary[ipoint_l.i][ipoint_l.j][ipoint_l.k];
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
		    pp = ipoint_add(p, &sim->dir_to_ipoint[dirz]);
		    pp = ipoint_add(&pp, &sim->dir_to_ipoint[dirx]);
		    if (sim_active_p(&pp)) {
			ipoint_l = ipoint_offset(p, dirx, dir, dirz);
			l = sim->lambda_ary[ipoint_l.i][ipoint_l.j][ipoint_l.k];
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
		    pp = ipoint_add(p, &sim->dir_to_ipoint[dirx]);
		    pp = ipoint_add(&pp, &sim->dir_to_ipoint[diry]);
		    if (sim_active_p(&pp)) {
			ipoint_l = ipoint_offset(p, dirx, diry, dir);
			l = sim->lambda_ary[ipoint_l.i][ipoint_l.j][ipoint_l.k];
			value += l/dz*(dx/2)*(dy/2);
		    }
		}
	    }
	    break;
	default:
	    bug("sim_set_matrix_coef");
	}
	sim->coefs[p->i][p->j][p->k].coef[dir].value = value;
    }
}

static void sim_set_matrix(void)
{
    iPoint *p;
    double dx, dy, dz;

    ALLOCATE_3D(sim->coefs, Coefs, world->ni, world->nj, world->nk);
    while (world_each(&p)) {
	if (p == NULL)
	    continue;
	sim->coefs[p->i][p->j][p->k] = get_coefs();
    }

    dx = world->dx;
    dy = world->dy;
    dz = world->dz;

    while (world_each(&p)) {
	if (p == NULL)
	    continue;
	if (sim->fix_ary[p->i][p->j][p->k] != NULL)
	    continue;

	sim_set_matrix_const(p);
	sim_set_matrix_coef0(p, dx, dy, dz);
	sim_set_matrix_coef(p, dx, dy, dz);
    }
}

Sim *sim_new(char *fname)
{
    sim = EALLOC(Sim);
    world = NULL;

    sim->dir_to_ipoint[DIR_LEFT]  = get_ipoint(-1,  0,  0);
    sim->dir_to_ipoint[DIR_RIGHT] = get_ipoint( 1,  0,  0);
    sim->dir_to_ipoint[DIR_FRONT] = get_ipoint( 0, -1,  0);
    sim->dir_to_ipoint[DIR_BACK]  = get_ipoint( 0,  1,  0);
    sim->dir_to_ipoint[DIR_BELOW] = get_ipoint( 0,  0, -1);
    sim->dir_to_ipoint[DIR_ABOVE] = get_ipoint( 0,  0,  1);

    if (opt_v)
	warn("configuring ...");
    config = config_new();
    config_parse(fname);

    if (opt_v)
	warn("setting region ...");
    sim_set_region();

    if (opt_v)
	warn("setting matrix ...");
    sim_set_matrix();

    return sim;
}

Array3Dd sim_calc()
{
    int nindex;
    Solvele *solver;
    iPoint *p, *hfp;
    iPoint *ipoint;
    int index, index2;
    int idir, dir;
    double c;
    double *sol;
    Array3Dd ary;

    nindex = world->ni * world->nj * world->nk;
    solver = solvele_new(nindex);
    while (world_each(&p)) {
	if (p == NULL)
	    continue;
	if (sim->fix_ary[p->i][p->j][p->k] != NULL) {
	    index = world_to_index(p);
	    solvele_set_matrix(solver, index, index, 1.0);
	    solvele_set_vector(solver, index, *(sim->fix_ary[p->i][p->j][p->k]));
	    continue;
	} else if (!sim->active_p_ary[p->i][p->j][p->k]) {
	    index = world_to_index(p);
	    solvele_set_matrix(solver, index, index, 1.0);
	    solvele_set_vector(solver, index, 0.0);
	    continue;
	}
	hfp = sim->heat_ipoint_ary[p->i][p->j][p->k];
	if (hfp != NULL) {
	    if (!ipoint_eq(p, hfp)) {
		index = world_to_index(p);
		solvele_set_matrix(solver, index, index, 1.0);
		solvele_set_matrix(solver, index, world_to_index(hfp), -1.0);
		solvele_set_vector(solver, index, 0.0);
	    }
	    ipoint = hfp;
	} else {
	    ipoint = p;
	}
	index = world_to_index(ipoint);
	for (idir = 0; idir < NELEMS(dir_array); ++idir) {
	    dir = dir_array[idir];
	    index2 = sim->coefs[p->i][p->j][p->k].coef[dir].index;
	    if (index2 < 0)
		continue;
	    c = sim->coefs[p->i][p->j][p->k].coef[dir].value;
	    solvele_add_matrix(solver, index, index2, c);
	}
	solvele_add_matrix(solver, index, index, sim->coefs[p->i][p->j][p->k].coef0);
	solvele_add_vector(solver, index, -sim->coefs[p->i][p->j][p->k].cnst);
    }

    if (opt_v)
	warn("solving equations ...");
    sol = solvele_solve(solver, world->ni, world->nj, world->nk);

    ALLOCATE_3D(ary, double, world->ni, world->nj, world->nk);
    while (world_each(&p)) {
	if (p == NULL)
	    continue;
	ary[p->i][p->j][p->k] = sol[world_to_index(p)];
    }

    return ary;
}
