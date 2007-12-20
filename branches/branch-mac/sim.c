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

/* World */

int world_to_index(iPoint *ipoint)
{
    return sim->index_tbl[ipoint->i + world->ni * (ipoint->j + world->nj * ipoint->k)];
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
    self->fixheat_obj_ary = aryobj_new();
    self->heat_obj_ary = aryobj_new();
    self->lambda_obj_ary = aryobj_new();
    return self;
}

static void config_free(void)
{
    aryobj_free(config->lambda_obj_ary);
    aryobj_free(config->heat_obj_ary);
    aryobj_free(config->fixheat_obj_ary);
    aryobj_free(config->fix_obj_ary);
    aryobj_free(config->active_obj_ary);
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
    iPoint_ary *ipoint_ary;
    int idir, dir;
    iPoint pp;

    ALLOCATE_3D2(sim->active_p_ary, int, world->ni, world->nj, world->nk, 0);
    active_obj_ary = config->active_obj_ary;

    /* fill with 1 or 0 */
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

    /* check edge cells of noactive which contact an active cell */
    ipoint_ary = ipoint_ary_new();
    for (index = 0; index < active_obj_ary->size; ++index) {
	obj = active_obj_ary->ptr[index];
	while (obj_each(obj, &p)) {
	    if (p == NULL)
		continue;
	    if (!world_inside_p(p))
		continue;
	    if (!obj->uval.i) {
		for (idir = 0; idir < NELEMS(dir_array); ++idir) {
		    dir = dir_array[idir];
		    pp = ipoint_add(p, &(sim->dir_to_ipoint[dir]));
		    if (sim_active_p(&pp)) {
			ipoint_ary_push(ipoint_ary, *p);
			break;
		    }
		}
	    }
	}
    }

    /* set active on edge cells */
    for (index = 0; index < ipoint_ary->size; ++index) {
	p = &(ipoint_ary->ptr[index]);
	sim->active_p_ary[p->i][p->j][p->k] = 1;
    }

    ipoint_ary_free(ipoint_ary);
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
	    if (sim->fix_ary[p->i][p->j][p->k] == NULL)
		sim->fix_ary[p->i][p->j][p->k] = double_new(obj->uval.d);
	    else
		*(sim->fix_ary[p->i][p->j][p->k]) = obj->uval.d;
	}
    }
}

static void sim_free_region_fix(void)
{
    iPoint *p;
    AryObj *obj_ary;
    Obj *obj;
    int i;

    obj_ary = config->fix_obj_ary;
    for (i = 0; i < obj_ary->size; ++i) {
	obj = obj_ary->ptr[i];
	while (obj_each(obj, &p)) {
	    if (p == NULL)
		continue;
	    if (!world_inside_p(p))
		continue;
	    FREE(sim->fix_ary[p->i][p->j][p->k]);
	    /* set NULL for avoid to free twice */
	    sim->fix_ary[p->i][p->j][p->k] = NULL;
	}
    }
    FREE_3D(sim->fix_ary);
}

static void sim_set_region_fixheat(void)
{
    iPoint *p;
    AryObj *obj_ary;
    int index;
    Obj *obj;
    int first;
    iPoint ipoint0;

    ALLOCATE_3D2(sim->fixheat_ary, double *, world->ni, world->nj, world->nk, NULL);
    ALLOCATE_3D2(sim->fixheat_ipoint_ary, iPoint *, world->ni, world->nj, world->nk, NULL);
    obj_ary = config->fixheat_obj_ary;
    ipoint0 = get_ipoint(-1, -1, -1);	/* for shutting up compiler */
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
	    if (sim->fixheat_ipoint_ary[p->i][p->j][p->k] == NULL) {
		sim->fixheat_ipoint_ary[p->i][p->j][p->k] = ipoint_new(ipoint0.i, ipoint0.j, ipoint0.k);
		sim->fixheat_ary[p->i][p->j][p->k] = double_new(obj->uval.d);
	    } else {
		*(sim->fixheat_ipoint_ary[p->i][p->j][p->k]) = ipoint0;
		*(sim->fixheat_ary[p->i][p->j][p->k]) = obj->uval.d;
	    }
	}
    }
}

static void sim_free_region_fixheat(void)
{
    iPoint *p;
    AryObj *obj_ary;
    int index;
    Obj *obj;

    obj_ary = config->fixheat_obj_ary;
    for (index = 0; index < obj_ary->size; ++index) {
	obj = obj_ary->ptr[index];
	while (obj_each(obj, &p)) {
	    if (p == NULL)
		continue;
	    if (!sim_active_p(p))
		continue;
	    FREE(sim->fixheat_ary[p->i][p->j][p->k]);
	    /* set NULL for avoid to free twice */
	    sim->fixheat_ary[p->i][p->j][p->k] = NULL;
	    FREE(sim->fixheat_ipoint_ary[p->i][p->j][p->k]);
	    /* set NULL for avoid to free twice */
	    sim->fixheat_ipoint_ary[p->i][p->j][p->k] = NULL;
	}
    }
    FREE_3D(sim->fixheat_ipoint_ary);
    FREE_3D(sim->fixheat_ary);
}

static int exist_eighth_part_of_block_p(iPoint *p, int dir, int di, int dj, int dk)
{
    iPoint pp, pp2;
    int ci, cj, ck;

    pp = ipoint_add(p, &sim->dir_to_ipoint[dir]);
    if (!sim_active_p(&pp))
	return 0;

    pp.i = p->i + di;
    pp.j = p->j + dj;
    pp.k = p->k + dk;
    if (sim_active_p(&pp))
	return 1;

    /* reverse direction */
    pp2 = ipoint_add(p, &sim->dir_to_ipoint[dir]);
    pp = *p;

    switch (dir) {
    case DIR_LEFT:
    case DIR_RIGHT:
	ci = -1;
	cj = ck = 1;
	break;
    case DIR_FRONT:
    case DIR_BACK:
	cj = -1;
	ck = ci = 1;
	break;
    case DIR_BELOW:
    case DIR_ABOVE :
	ck = -1;
	ci = cj = 1;
	break;
    default:
	bug("exist_eighth_part_of_block_p");
	ci = cj = ck = -1;	/* for shutting up compiler */
    }

    pp.i = pp2.i + ci * di;
    pp.j = pp2.j + cj * dj;
    pp.k = pp2.k + ck * dk;
    if (sim_active_p(&pp))
	return 1;

    return 0;
}

static void sim_set_region_heat(void)
{
    Array3Di heat_p_ary;
    iPoint *p;
    AryObj *obj_ary;
    int index;
    Obj *obj;
    double heat_coef, sum_heat_coef;
    int idir, dir;
    int ix, iy, iz;
    IP_TYPE di, dj, dk;
    int dirx, diry, dirz;

    ALLOCATE_3D2(sim->heat_ary, double *, world->ni, world->nj, world->nk, NULL);
    ALLOCATE_3D(heat_p_ary, int, world->ni, world->nj, world->nk);
    obj_ary = config->heat_obj_ary;
    for (index = 0; index < obj_ary->size; ++index) {
	while (world_each(&p))
	    heat_p_ary[p->i][p->j][p->k] = 0;
	obj = obj_ary->ptr[index];
	while (obj_each(obj, &p)) {
	    if (p == NULL)
		continue;
	    if (!sim_active_p(p))
		continue;
	    heat_p_ary[p->i][p->j][p->k] = 1;
	}
	sum_heat_coef = 0;
	while (world_each(&p)) {
	    if (!heat_p_ary[p->i][p->j][p->k])
		continue;
	    heat_coef = 0.0;
	    for (idir = 0; idir < NELEMS(dir_array); ++idir) {
		dir = dir_array[idir];
		switch (dir) {
		case DIR_LEFT: case DIR_RIGHT:
		    for (iy = 0; iy < NELEMS(dir_y); ++iy) {
			diry = dir_y[iy];
			dj = sim->dir_to_ipoint[diry].j;
			for (iz = 0; iz < NELEMS(dir_z); ++iz) {
			    dirz = dir_z[iz];
			    dk = sim->dir_to_ipoint[dirz].k;
			    if (exist_eighth_part_of_block_p(p, dir, 0, dj, dk))
				heat_coef += 1.0 / (6 * 2 * 2);
			}
		    }
		    break;
		case DIR_FRONT: case DIR_BACK:
		    for (iz = 0; iz < NELEMS(dir_z); ++iz) {
			dirz = dir_z[iz];
			dk = sim->dir_to_ipoint[dirz].k;
			for (ix = 0; ix < NELEMS(dir_x); ++ix) {
			    dirx = dir_x[ix];
			    di = sim->dir_to_ipoint[dirx].i;
			    if (exist_eighth_part_of_block_p(p, idir, di, 0, dk))
				heat_coef += 1.0 / (6 * 2 * 2);
			}
		    }
		    break;
		case DIR_BELOW: case DIR_ABOVE:
		    for (ix = 0; ix < NELEMS(dir_x); ++ix) {
			dirx = dir_x[ix];
			di = sim->dir_to_ipoint[dirx].i;
			for (iy = 0; iy < NELEMS(dir_y); ++iy) {
			    diry = dir_y[iy];
			    dj = sim->dir_to_ipoint[diry].j;
			    if (exist_eighth_part_of_block_p(p, idir, di, dj, 0))
				heat_coef += 1.0 / (6 * 2 * 2);
			}
		    }
		    break;
		default:
		    bug("sim_set_region_heat");
		}
	    }
	    if (sim->heat_ary[p->i][p->j][p->k] == NULL)
		sim->heat_ary[p->i][p->j][p->k] = double_new(heat_coef * obj->uval.d);
	    else
		*(sim->heat_ary[p->i][p->j][p->k]) = heat_coef * obj->uval.d;
	    sum_heat_coef += heat_coef;
	}
	while (world_each(&p)) {
	    if (!heat_p_ary[p->i][p->j][p->k])
		continue;
	    *(sim->heat_ary[p->i][p->j][p->k]) /= sum_heat_coef;
	}
    }
    FREE_3D(heat_p_ary);
}

static void sim_free_region_heat(void)
{
    iPoint *p;
    AryObj *obj_ary;
    int index;
    Obj *obj;

    obj_ary = config->heat_obj_ary;
    for (index = 0; index < obj_ary->size; ++index) {
	obj = obj_ary->ptr[index];
	while (obj_each(obj, &p)) {
	    if (p == NULL)
		continue;
	    if (!sim_active_p(p))
		continue;
	    FREE(sim->heat_ary[p->i][p->j][p->k]);
	    /* set NULL for avoid to free twice */
	    sim->heat_ary[p->i][p->j][p->k] = NULL;
	}
    }
    FREE_3D(sim->heat_ary);
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

static void sim_free_region_lambda(void)
{
    FREE_3D(sim->lambda_ary);
}

static void sim_set_region(void)
{
    sim_set_region_active();
    sim_set_region_fix();
    sim_set_region_fixheat();
    sim_set_region_heat();
    sim_set_region_lambda();
}

static void sim_free_region(void)
{
    sim_free_region_lambda();
    sim_free_region_heat();
    sim_free_region_fixheat();
    sim_free_region_fix();
}

static void sim_add_matrix_coef0(iPoint *p0, iPoint *p, double dx, double dy, double dz)
{
    int idir, dir;
    int ix, iy, iz;
    int dirx, diry, dirz;
    IP_TYPE di, dj, dk;
    iPoint ipoint_l;
    double l;
    int index0, index;

    index0 = world_to_index(p0);
    index = world_to_index(p);
    for (idir = 0; idir < NELEMS(dir_array); ++idir) {
	dir = dir_array[idir];
	switch (dir) {
	case DIR_LEFT: case DIR_RIGHT:
	    for (iy = 0; iy < NELEMS(dir_y); ++iy) {
		diry = dir_y[iy];
		dj = sim->dir_to_ipoint[diry].j;
		for (iz = 0; iz < NELEMS(dir_z); ++iz) {
		    dirz = dir_z[iz];
		    dk = sim->dir_to_ipoint[dirz].k;
		    if (exist_eighth_part_of_block_p(p, dir, 0, dj, dk)) {
			ipoint_l = ipoint_offset(p, dir, diry, dirz);
			l = sim->lambda_ary[ipoint_l.i][ipoint_l.j][ipoint_l.k];
			solvele_add_matrix(sim->solver, index0, index0, -l/dx*(dy/2)*(dz/2));
		    }
		}
	    }
	    break;
	case DIR_FRONT: case DIR_BACK:
	    for (iz = 0; iz < NELEMS(dir_z); ++iz) {
		dirz = dir_z[iz];
		dk = sim->dir_to_ipoint[dirz].k;
		for (ix = 0; ix < NELEMS(dir_x); ++ix) {
		    dirx = dir_x[ix];
		    di = sim->dir_to_ipoint[dirx].i;
		    if (exist_eighth_part_of_block_p(p, idir, di, 0, dk)) {
			ipoint_l = ipoint_offset(p, dirx, dir, dirz);
			l = sim->lambda_ary[ipoint_l.i][ipoint_l.j][ipoint_l.k];
			solvele_add_matrix(sim->solver, index0, index0, -l/dy*(dz/2)*(dx/2));
		    }
		}
	    }
	    break;
	case DIR_BELOW: case DIR_ABOVE:
	    for (ix = 0; ix < NELEMS(dir_x); ++ix) {
		dirx = dir_x[ix];
		di = sim->dir_to_ipoint[dirx].i;
		for (iy = 0; iy < NELEMS(dir_y); ++iy) {
		    diry = dir_y[iy];
		    dj = sim->dir_to_ipoint[diry].j;
		    if (exist_eighth_part_of_block_p(p, idir, di, dj, 0)) {
			ipoint_l = ipoint_offset(p, dirx, diry, dir);
			l = sim->lambda_ary[ipoint_l.i][ipoint_l.j][ipoint_l.k];
			solvele_add_matrix(sim->solver, index0, index0, -l/dz*(dx/2)*(dy/2));
		    }
		}
	    }
	    break;
	default:
	    bug("sim_add_matrix_coef0");
	}
    }
}

static void sim_add_matrix_coef(iPoint *p0, iPoint *p, double dx, double dy, double dz)
{
    int idir, dir;
    int dirx, diry, dirz;
    double value;
    int ix, iy, iz;
    IP_TYPE di, dj, dk;
    iPoint pp, ipoint_l;
    double l;
    int index0, index, index2;

    index0 = world_to_index(p0);
    index = world_to_index(p);
    for (idir = 0; idir < NELEMS(dir_array); ++idir) {
	dir = dir_array[idir];
	pp = ipoint_add(p, &sim->dir_to_ipoint[dir]);
	if (!sim_active_p(&pp)) {
	    continue;
	}

	index2 = world_to_index(&pp);
	value = 0.0;
	switch (dir) {
	case DIR_LEFT: case DIR_RIGHT:
	    for (iy = 0; iy < NELEMS(dir_y); ++iy) {
		diry = dir_y[iy];
		dj = sim->dir_to_ipoint[diry].j;
		for (iz = 0; iz < NELEMS(dir_z); ++iz) {
		    dirz = dir_z[iz];
		    dk = sim->dir_to_ipoint[dirz].k;
		    if (exist_eighth_part_of_block_p(p, dir, 0, dj, dk)) {
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
		dk = sim->dir_to_ipoint[dirz].k;
		for (ix = 0; ix < NELEMS(dir_x); ++ix) {
		    dirx = dir_x[ix];
		    di = sim->dir_to_ipoint[dirx].i;
		    if (exist_eighth_part_of_block_p(p, idir, di, 0, dk)) {
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
		di = sim->dir_to_ipoint[dirx].i;
		for (iy = 0; iy < NELEMS(dir_y); ++iy) {
		    diry = dir_y[iy];
		    dj = sim->dir_to_ipoint[diry].j;
		    if (exist_eighth_part_of_block_p(p, idir, di, dj, 0)) {
			ipoint_l = ipoint_offset(p, dirx, diry, dir);
			l = sim->lambda_ary[ipoint_l.i][ipoint_l.j][ipoint_l.k];
			value += l/dz*(dx/2)*(dy/2);
		    }
		}
	    }
	    break;
	default:
	    bug("sim_add_matrix_coef");
	}
	solvele_add_matrix(sim->solver, index0, index2, value);
    }
}

static void sim_set_matrix(void)
{
    iPoint *p, *fh0p;
    double dx, dy, dz;
    int index, index_fh0p;
    double val;

    sim->solver = solvele_new(sim->index_size);

    dx = world->dx;
    dy = world->dy;
    dz = world->dz;

    while (world_each(&p)) {
	if (!sim_active_p(p))
	    continue;
	index = world_to_index(p);
	assert(index >= 0);
	if (!sim->active_p_ary[p->i][p->j][p->k]) {
	    solvele_set_matrix(sim->solver, index, index, 1.0);
	    solvele_set_vector(sim->solver, index, 0.0);
	} else if (sim->fix_ary[p->i][p->j][p->k] != NULL) {
	    solvele_set_matrix(sim->solver, index, index, 1.0);
	    val = *(sim->fix_ary[p->i][p->j][p->k]);
	    solvele_set_vector(sim->solver, index, val);
	} else if (sim->heat_ary[p->i][p->j][p->k] != NULL) {
	    sim_add_matrix_coef0(p, p, dx, dy, dz);
	    sim_add_matrix_coef(p, p, dx, dy, dz);
	    val = - *(sim->heat_ary[p->i][p->j][p->k]);
	    solvele_set_vector(sim->solver, index, val);
	} else if (sim->fixheat_ary[p->i][p->j][p->k] != NULL) {
	    fh0p = sim->fixheat_ipoint_ary[p->i][p->j][p->k];
	    if (ipoint_eq(p, fh0p)) {
		sim_add_matrix_coef0(p, p, dx, dy, dz);
		sim_add_matrix_coef(p, p, dx, dy, dz);
		val = - *(sim->fixheat_ary[p->i][p->j][p->k]);
		solvele_set_vector(sim->solver, index, val);
	    } else {
		index_fh0p = world_to_index(fh0p);
		solvele_set_matrix(sim->solver, index, index, -1.0);
		solvele_set_matrix(sim->solver, index, index_fh0p, 1.0);
		sim_add_matrix_coef0(fh0p, p, dx, dy, dz);
		sim_add_matrix_coef(fh0p, p, dx, dy, dz);
	    }
	} else {
	    sim_add_matrix_coef0(p, p, dx, dy, dz);
	    sim_add_matrix_coef(p, p, dx, dy, dz);
	}
    }
}

Sim *sim_new(char *fname, double eps_sor, double omega_sor)
{
    sim = EALLOC(Sim);
    world = NULL;
    sim->dir_to_ipoint[DIR_LEFT]  = get_ipoint(-1,  0,  0);
    sim->dir_to_ipoint[DIR_RIGHT] = get_ipoint( 1,  0,  0);
    sim->dir_to_ipoint[DIR_FRONT] = get_ipoint( 0, -1,  0);
    sim->dir_to_ipoint[DIR_BACK]  = get_ipoint( 0,  1,  0);
    sim->dir_to_ipoint[DIR_BELOW] = get_ipoint( 0,  0, -1);
    sim->dir_to_ipoint[DIR_ABOVE] = get_ipoint( 0,  0,  1);
    sim->fname = fname;
    sim->eps_sor = eps_sor;
    sim->omega_sor = omega_sor;

    return sim;
}

static double *sim_get_region(int region)
{
    int size;
    double *u;
    int i, j, k;
    iPoint ipoint;
    int index;
    double *p;

    size = sim->index_size;
    u = EALLOCN(double, size);
    for (k = 0; k < world->nk; ++k) {
	ipoint.k = k;
	for (j = 0; j < world->nj; ++j) {
	    ipoint.j = j;
	    for (i = 0; i < world->ni; ++i) {
		ipoint.i = i;
		if (!sim_active_p(&ipoint))
		    continue;
		index = world_to_index(&ipoint);
		switch (region) {
		case REGION_ACTIVE:
		    u[index] = 1.0;
		    break;
		case REGION_FIX:
		    p = sim->fix_ary[ipoint.i][ipoint.j][ipoint.k];
		    if (p != NULL)
			u[index] = *p;
		    else
			u[index] = -1.0;
		    break;
		case REGION_FIXHEAT:
		    p = sim->fixheat_ary[ipoint.i][ipoint.j][ipoint.k];
		    if (p != NULL)
			u[index] = *p;
		    else
			u[index] = 0.0;
		    break;
		case REGION_HEAT:
		    p = sim->heat_ary[ipoint.i][ipoint.j][ipoint.k];
		    if (p != NULL)
			u[index] = *p;
		    else
			u[index] = 0.0;
		    break;
		case REGION_LAMBDA:
		    u[index] = sim->lambda_ary[ipoint.i][ipoint.j][ipoint.k];
		    break;
		default:
		    bug("unknown region %d", region);
		    /* NOTREACHED */
		    u = NULL;	/* for shutting up compiler */
		}
	    }
	}
    }
    return u;
}

static void sim_set_index_tbl(void)
{
    int i, j, k;
    iPoint ipoint;
    int index0;

    sim->index_tbl = EALLOCN(int, world->ni * world->nj * world->nk);
    for (k = 0; k < world->nk; ++k) {
	for (j = 0; j < world->nj; ++j) {
	    for (i = 0; i < world->ni; ++i) {
		index0 = i + world->ni * (j + world->nj * k);
		sim->index_tbl[index0] = -1;
	    }
	}
    }
    sim->index_size = 0;
    for (k = 0; k < world->nk; ++k) {
	ipoint.k = k;
	for (j = 0; j < world->nj; ++j) {
	    ipoint.j = j;
	    for (i = 0; i < world->ni; ++i) {
		ipoint.i = i;
		if (sim_active_p(&ipoint)) {
		    index0 = i + world->ni * (j + world->nj * k);
		    sim->index_tbl[index0] = sim->index_size++;
		}
	    }
	}
    }
}

double *sim_calc()
{
    double *u;

    if (opt_v)
	warn("configuring ...");
    config = config_new();
    config_parse(sim->fname);
    warn("number of cells is %d x %d x %d", world->ni, world->nj, world->nk);

    if (opt_v)
	warn("setting region ...");
    sim_set_region();

    sim_set_index_tbl();

    if (opt_r == 0) {
	if (opt_v)
	    warn("setting matrix ...");
	sim_set_matrix();

	sim_free_region();
	config_free();

	if (opt_v)
	    warn("solving equations ...");

	u = solvele_solve(sim->solver, sim->eps_sor, sim->omega_sor,
		world->ni, world->nj, world->nk);
    } else {
	u = sim_get_region(opt_r);
	sim_free_region();
	config_free();
    }

    return u;
}
