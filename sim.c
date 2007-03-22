#include <stdio.h>
#include <assert.h>
#include "sim.h"
#include "mem.h"
#include "tc.h"
#include "solvele.h"

#define NELEMS(ary)	(sizeof(ary) / sizeof((ary)[0]))

Config *config_parser;

static int dir_array[NDIRS] = { DIR_LEFT, DIR_RIGHT, DIR_FRONT, DIR_BACK, DIR_BELOW, DIR_ABOVE };
static int dir_x[] = { DIR_LEFT, DIR_RIGHT };
static int dir_y[] = { DIR_FRONT, DIR_BACK };
static int dir_z[] = { DIR_BELOW, DIR_ABOVE };
static int axis_array[] = { AXIS_X, AXIS_Y, AXIS_Z };

static int iround(double x)
{
    return (int) (x + 0.5);
}

/* Point */

Point get_point(int i, int j, int k)
{
    Point point;

    point.i = i;
    point.j = j;
    point.k = k;
    return point;
}

static Point point_offset(Point point, int dirx, int diry, int dirz)
{
    int i, j, k;

    switch (dirx) {
    case DIR_LEFT:
	i = point.i - 1;
	break;
    case DIR_RIGHT:
	i = point.i;
	break;
    default:
	bug("unknown dir %d", dirx);
    }

    switch (diry) {
    case DIR_FRONT:
	j = point.j - 1;
	break;
    case DIR_BACK:
	j = point.j;
	break;
    default:
	bug("unknown dir %d", diry);
    }

    switch (dirz) {
    case DIR_BELOW:
	k = point.k - 1;
	break;
    case DIR_ABOVE:
	k = point.k;
	break;
    default:
	bug("unknown dir %d", dirz);
    }

    return get_point(i, j, k);
}

static int point_eq(Point point1, Point point2)
{
    if (point1.i == point2.i && point1.j == point2.j && point1.k == point2.k)
	return 1;
    else
	return 0;
}

static Point point_add(Point point1, Point point2)
{
    return get_point(
	    point1.i + point2.i,
	    point1.j + point2.j,
	    point1.k + point2.k
	    );
}

static Point *point_new(int i, int j, int k)
{
    Point *self;

    self = EALLOC(Point);
    self->i = i;
    self->j = j;
    self->k = k;
    return self;
}

/* Coef */

static Coef *coef_new(void)
{
    Coef *self;

    self = EALLOC(Coef);
    self->index = -1;
    self->value = 0.0;
    return self;
}


/* Coefs */

static Coefs *coefs_new(void)
{
    Coefs *self;
    int idir, dir;

    self = EALLOC(Coefs);
    for (idir = 0; idir < NELEMS(dir_array); ++idir) {
	dir = dir_array[idir];
	self->coef[dir] = coef_new();
    }
    self->coef0 = 0.0;
    self->cnst = 0.0;
    return self;
}

/* Each */

static Each *each_new(int size, Point *point_ary)
{
    Each *self;

    self = EALLOC(Each);
    self->size = size;
    self->point_ary = point_ary;
    self->index = 0;
    return self;
}

static Point *each_each(Each *self)
{
    if (self == NULL)
	bug("each_each");
    if (self->index < self->size) {
	return &self->point_ary[self->index++];
    } else {
	self->index = -1;
	return NULL;
    }
}

/* World */

static int world_to_index(World *self, Point point)
{
    return point.i + self->ni * (point.j + self->nj * point.k);
}

World *world_new(double x0, double y0, double z0,
	double xlen, double ylen, double zlen,
	double dx, double dy, double dz)
{
    World *self;

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

static Point *world_each_begin(World *self)
{
    int size;
    Point *point_ary;
    int n;
    int i, j, k;

    if (self->each != NULL && self->each->index >= 0)
	bug("world_each_begin");
    size = self->ni * self->nj * self->nk;
    point_ary = EALLOCN(Point, size);
    n = 0;
    for (i = 0; i < self->ni; ++i) {
	for (j = 0; j < self->nj; ++j) {
	    for (k = 0; k < self->nk; ++k) {
		point_ary[n++] = get_point(i, j, k);
	    }
	}
    }
    self->each = each_new(size, point_ary);
    return each_each(self->each);
}

static Point *world_each(World *self)
{
    return each_each(self->each);
}

/* Rect */

Rect *rect_new(World *world, double x, double y, double z, int axis, double len1, double len2)
{
    Rect *self;
    int ok, iaxis;

    ok = 0;
    for (iaxis = 0; iaxis < NELEMS(axis_array); ++iaxis) {
	if (axis == axis_array[iaxis]) {
	    ok = 1;
	    break;
	}
    }
    if (!ok)
	bug("unknown axis %d", axis);

    if (len1 <= 0.0 || len2 <= 0.0)
	warn_exit("length is negative for Rect");

    self = EALLOC(Rect);
    self->world = world;
    self->x = x;
    self->y = y;
    self->z = z;
    self->axis = axis;
    self->len1 = len1;
    self->len2 = len2;
    self->each = NULL;
    return self;
}

static Point *rect_each_begin(Rect *self)
{
    int xi, yi, zi;
    int size;
    int len1, len2;
    Point *point_ary;
    int i, j, k;

    if (self->each != NULL && self->each->index >= 0)
	bug("rect_each_begin");

    xi = iround((self->x - self->world->x0) / self->world->dx);
    yi = iround((self->y - self->world->y0) / self->world->dy);
    zi = iround((self->z - self->world->z0) / self->world->dz);
    size = 0;
    switch (self->axis) {
    case AXIS_X:
	len1 = iround(self->len1 / self->world->dy) + 1;
	len2 = iround(self->len2 / self->world->dz) + 1;
	point_ary = EALLOCN(Point, len1 * len2);
	for (j = yi; j < yi + len1; ++j) {
	    for (k = zi; k < zi + len2; ++k) {
		point_ary[size++] = get_point(xi, j, k);
	    }
	}
	break;
    case AXIS_Y:
	len1 = iround(self->len1 / self->world->dx) + 1;
	len2 = iround(self->len2 / self->world->dz) + 1;
	point_ary = EALLOCN(Point, len1 * len2);
	for (i = xi; i < xi + len1; ++i) {
	    for (k = zi; k < zi + len2; ++k) {
		point_ary[size++] = get_point(i, yi, k);
	    }
	}
	break;
    case AXIS_Z:
	len1 = iround(self->len1 / self->world->dx) + 1;
	len2 = iround(self->len2 / self->world->dy) + 1;
	point_ary = EALLOCN(Point, len1 * len2);
	for (i = xi; i < xi + len1; ++i) {
	    for (j = yi; j < yi + len2; ++j) {
		point_ary[size++] = get_point(i, j, zi);
	    }
	}
	break;
    default:
	bug("unknow axis %d for rect_each_begin", self->axis);
    }

    assert(size == len1 * len2);

    self->each = each_new(size, point_ary);
    return each_each(self->each);
}

static Point *rect_each(Rect *self)
{
    return each_each(self->each);
}

static void rect_offset(Rect *self)
{
    switch (self->axis) {
    case AXIS_X:
	self->len1 -= self->world->dy;
	self->len2 -= self->world->dz;
	break;
    case AXIS_Y:
	self->len1 -= self->world->dz;
	self->len2 -= self->world->dx;
	break;
    case AXIS_Z:
	self->len1 -= self->world->dx;
	self->len2 -= self->world->dy;
	break;
    default:
	bug("unknow axis %d", self->axis);
    }
}

/* Box */

Box *box_new(World *world, double x, double y, double z, double xlen, double ylen, double zlen)
{
    Box *self;

    if (xlen <= 0.0 || ylen <= 0.0 || zlen <= 0.0)
	warn_exit("length is negative for Rect");

    self = EALLOC(Box);
    self->world = world;
    self->rect = rect_new(world, x, y, z, AXIS_Z, xlen, ylen);
    self->sweeplen = zlen;
    self->each = NULL;
    return self;
}

static Point *box_each_begin(Box *self)
{
    Point *p;
    int leni;
    int size;
    Point *point_ary;
    int n;
    int k;

    if (self->each != NULL)
	bug("box_each_begin");
    p = rect_each_begin(self->rect);
    leni = iround(self->sweeplen / self->world->dz) + 1;
    size = self->rect->each->size * leni;
    point_ary = EALLOCN(Point, size);
    n = 0;
    do {
	for (k = p->k; k < p->k + leni; ++k) {
	    point_ary[n++] = get_point(p->i, p->j, k);
	}
	p = rect_each(self->rect);
    } while (p != NULL);
    assert(n == size);
    self->each = each_new(size, point_ary);
    return each_each(self->each);
}

static Point *box_each(Box *self)
{
    return each_each(self->each);
}

static void box_offset(Box *self)
{
    rect_offset(self->rect);
    self->sweeplen -= self->world->dz;
}

/* Obj */

Obj *obj_new(int objtype)
{
    Obj *self;

    self = EALLOC(Obj);
    self->objtype = objtype;
    return self;
}

static Point *obj_each_begin(Obj *self)
{
    Point *p;

    switch (self->objtype) {
    case OBJ_RECT:
	p = rect_each_begin(self->uobj.rect);
	break;
    case OBJ_BOX:
	p = box_each_begin(self->uobj.box);
	break;
    default:
	bug("unknown obj %d", self->objtype);
    }
    return p;
}

static Point *obj_each(Obj *self)
{
    Point *p;

    switch (self->objtype) {
    case OBJ_RECT:
	p = rect_each(self->uobj.rect);
	break;
    case OBJ_BOX:
	p = box_each(self->uobj.box);
	break;
    default:
	bug("unknown obj %d", self->objtype);
    }
    return p;
}

void obj_offset(Obj *self)
{
    switch (self->objtype) {
    case OBJ_RECT:
	rect_offset(self->uobj.rect);
	break;
    case OBJ_BOX:
	box_offset(self->uobj.box);
	break;
    default:
	bug("unknown obj %d", self->objtype);
    }
}

/* AryObj */

AryObj *aryobj_new(void)
{
    AryObj *self;

    self = EALLOC(AryObj);
    self->ptr = NULL;
    self->size = 0;
    return self;
}

void aryobj_push(AryObj *self, Obj *obj)
{
    ++self->size;
    self->ptr = erealloc(self->ptr, sizeof(Obj *) * self->size);
    self->ptr[self->size - 1] = obj;
}

/* Heatflow */

Heatflow *heatflow_new(int dir, double value)
{
    Heatflow *self;

    self = EALLOC(Heatflow);
    self->dir = dir;
    self->value = value;
    return self;
}

/* Config */

static void config_parse(Config *self)
{
    extern int yyparse();
    /* extern int yydebug; */

    /* yydebug = 1; */

    config_parser = self;
    yyparse();
}

static Config *config_new(void)
{
    Config *self;

    self = EALLOC(Config);
    self->world = NULL;
    self->active_obj_ary = aryobj_new();
    self->fix_obj_ary = aryobj_new();
    self->heatflow_obj_ary = aryobj_new();
    self->lambda_obj_ary = aryobj_new();

    config_parse(self);

    return self;
}

/* Sim */

int sim_active_p(Sim *self, Point point)
{
    if (point.i < 0 || point.i >= self->ni ||
	    point.j < 0 || point.j >= self->nj ||
	    point.k < 0 || point.k >= self->nk)
	return 0;
    else
	return self->active_p_ary[point.i][point.j][point.k];
}

static void sim_set_region_active(Sim *self)
{
    static int dir_xyz[3][2] = {
	{ DIR_LEFT, DIR_RIGHT },
	{ DIR_FRONT, DIR_BACK },
	{ DIR_BELOW, DIR_ABOVE }
    };
    Point *p;
    AryObj *active_obj_ary;
    Obj *obj;
    int index;
    int continue_p;
    int *dirs;

    ALLOCATE_3D2(self->active_p_ary, int, self->ni, self->nj, self->nk, 0);
    active_obj_ary = self->config->active_obj_ary;
    for (index = 0; index < active_obj_ary->size; ++index) {
	obj = active_obj_ary->ptr[index];
	for (p = obj_each_begin(obj); p != NULL; p = obj_each(obj)) {
	    (self->active_p_ary)[p->i][p->j][p->k] = obj->uval.i;
	}
    }

    /* remove active flag from the point having active alone plane in NDIR directions */
    do {
	continue_p = 0;
	for (p = world_each_begin(self->world); p != NULL; p = world_each(self->world)) {
	    if (!self->active_p_ary[p->i][p->j][p->k])
		continue;

	    for (index = 0; index < NELEMS(dir_xyz); ++index) {
		dirs = dir_xyz[index];
		if (!sim_active_p(self, point_add(*p, self->dir_to_point[dirs[0]])) &&
			!sim_active_p(self, point_add(*p, self->dir_to_point[dirs[1]]))) {
		    self->active_p_ary[p->i][p->j][p->k] = 0;
		    continue_p = 1;
		    warn("removed (%d, %d, %d)\n", p->i, p->j, p->k);
		}
	    }
	}
    } while (continue_p);
}

static void sim_set_region_fix(Sim *self)
{
    Point *p;
    AryObj *obj_ary;
    Obj *obj;
    int i;

    ALLOCATE_3D2(self->fix_ary, double, self->ni, self->nj, self->nk, -1.0);
    obj_ary = self->config->fix_obj_ary;
    for (i = 0; i < obj_ary->size; ++i) {
	obj = obj_ary->ptr[i];
	for (p = obj_each_begin(obj); p != NULL; p = obj_each(obj)) {
	    assert(obj->uval.d >= 0.0);
	    self->fix_ary[p->i][p->j][p->k] = obj->uval.d;
	}
    }
}

static void sim_set_region_heatflow(Sim *self)
{
    Point *p;
    int idir, dir;
    AryObj *obj_ary;
    int index;
    Obj *obj;
    int first;
    Point point0;

    ALLOCATE_3D(self->heatflow_ary, double *, self->ni, self->nj, self->nk);
    for (p = world_each_begin(self->world); p != NULL; p = world_each(self->world)) {
	self->heatflow_ary[p->i][p->j][p->k] = EALLOCN(double, NELEMS(dir_array));
	for (idir = 0; idir < NELEMS(dir_array); ++idir) {
	    dir = dir_array[idir];
	    self->heatflow_ary[p->i][p->j][p->k][dir] = -1.0;
	}
    }

    ALLOCATE_3D2(self->heatflow_point_ary, Point *, self->ni, self->nj, self->nk, NULL);
    obj_ary = self->config->heatflow_obj_ary;
    for (index = 0; index < obj_ary->size; ++index) {
	obj = obj_ary->ptr[index];
	first = 1;
	for (p = obj_each_begin(obj); p != NULL; p = obj_each(obj)) {
	    if (first) {
		first = 0;
		point0 = *p;
	    }
	    self->heatflow_point_ary[p->i][p->j][p->k] = point_new(point0.i, point0.j, point0.k);
	    self->heatflow_ary[p->i][p->j][p->k][obj->uval.h->dir] = obj->uval.h->value;
	}
    }
}

static void sim_set_region_lambda(Sim *self)
{
    int index;
    AryObj *obj_ary;
    Obj *obj;
    Point *p;

    ALLOCATE_3D2(self->lambda_ary, double, self->ni, self->nj, self->nk, 1.0);
    obj_ary = self->config->lambda_obj_ary;
    for (index = 0; index < obj_ary->size; ++index) {
	obj = obj_ary->ptr[index];
	for (p = obj_each_begin(obj); p != NULL; p = obj_each(obj)) {
	    self->lambda_ary[p->i][p->j][p->k] = obj->uval.d;
	}
    }
}

static void sim_set_region(Sim *self)
{
    self->world = self->config->world;
    self->ni = self->world->ni;
    self->nj = self->world->nj;
    self->nk = self->world->nk;

    sim_set_region_active(self);
    sim_set_region_fix(self);
    sim_set_region_heatflow(self);
    sim_set_region_lambda(self);
}

static void sim_set_matrix_const(Sim *self, int i, int j, int k, double dx, double dy, double dz)
{
    int idir, dir;

    for (idir = 0; idir < NELEMS(dir_array); ++idir) {
	dir = dir_array[idir];
	if (self->heatflow_ary[i][j][k][dir] >= 0.0 &&
		point_eq(*(self->heatflow_point_ary[i][j][k]), get_point(i, j, k)))
	    self->coefs[i][j][k]->cnst = self->heatflow_ary[i][j][k][dir];
    }
}

static void sim_set_matrix_coef0(Sim *self, int i, int j, int k, double dx, double dy, double dz)
{
    int idir, dir;
    int ix, iy, iz;
    int dirx, diry, dirz;
    Point pp, point_l;
    double l;

    for (idir = 0; idir < NELEMS(dir_array); ++idir) {
	dir = dir_array[idir];
	if (!sim_active_p(self, point_add(get_point(i, j, k), self->dir_to_point[dir])))
	    continue;
	switch (dir) {
	case DIR_LEFT: case DIR_RIGHT:
	    for (iy = 0; iy < NELEMS(dir_y); ++iy) {
		diry = dir_y[iy];
		for (iz = 0; iz < NELEMS(dir_z); ++iz) {
		    dirz = dir_z[iz];
		    pp = point_add(get_point(i, j, k), self->dir_to_point[diry]);
		    pp = point_add(pp, self->dir_to_point[dirz]);
		    if (sim_active_p(self, pp)) {
			point_l = point_offset(get_point(i, j, k), dir, diry, dirz);
			l = self->lambda_ary[point_l.i][point_l.j][point_l.k];
			self->coefs[i][j][k]->coef0 += -l/dx*(dy/2)*(dz/2);
		    }
		}
	    }
	    break;
	case DIR_FRONT: case DIR_BACK:
	    for (iz = 0; iz < NELEMS(dir_z); ++iz) {
		dirz = dir_z[iz];
		for (ix = 0; ix < NELEMS(dir_x); ++ix) {
		    dirx = dir_x[ix];
		    pp = point_add(get_point(i, j, k), self->dir_to_point[dirz]);
		    pp = point_add(pp, self->dir_to_point[dirx]);
		    if (sim_active_p(self, pp)) {
			point_l = point_offset(get_point(i, j, k), dirx, dir, dirz);
			l = self->lambda_ary[point_l.i][point_l.j][point_l.k];
			self->coefs[i][j][k]->coef0 += -l/dy*(dz/2)*(dx/2);
		    }
		}
	    }
	    break;
	case DIR_BELOW: case DIR_ABOVE:
	    for (ix = 0; ix < NELEMS(dir_x); ++ix) {
		dirx = dir_x[ix];
		for (iy = 0; iy < NELEMS(dir_y); ++iy) {
		    diry = dir_y[iy];
		    pp = point_add(get_point(i, j, k), self->dir_to_point[dirx]);
		    pp = point_add(pp, self->dir_to_point[diry]);
		    if (sim_active_p(self, pp)) {
			point_l = point_offset(get_point(i, j, k), dirx, diry, dir);
			l = self->lambda_ary[point_l.i][point_l.j][point_l.k];
			self->coefs[i][j][k]->coef0 += -l/dz*(dx/2)*(dy/2);
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
    Point pp, point_l;
    double l;

    for (idir = 0; idir < NELEMS(dir_array); ++idir) {
	dir = dir_array[idir];
	if (!sim_active_p(self, point_add(get_point(i, j, k), self->dir_to_point[dir]))) {
	    self->coefs[i][j][k]->coef[dir]->index = -1;
	    continue;
	}
	if (self->heatflow_ary[i][j][k][dir] >= 0.0)
	    warn_exit("heatflow(%g) comes from active cell at (%d, %d, %d), dir(%d)", self->heatflow_ary[i][j][k][dir], i, j, k, dir);

	self->coefs[i][j][k]->coef[dir]->index =
	    world_to_index(self->world, point_add(get_point(i, j, k), self->dir_to_point[dir]));
	value = 0.0;
	switch (dir) {
	case DIR_LEFT: case DIR_RIGHT:
	    for (iy = 0; iy < NELEMS(dir_y); ++iy) {
		diry = dir_y[iy];
		for (iz = 0; iz < NELEMS(dir_z); ++iz) {
		    dirz = dir_z[iz];
		    pp = point_add(get_point(i, j, k), self->dir_to_point[diry]);
		    pp = point_add(pp, self->dir_to_point[dirz]);
		    if (sim_active_p(self, pp)) {
			point_l = point_offset(get_point(i, j, k), dir, diry, dirz);
			l = self->lambda_ary[point_l.i][point_l.j][point_l.k];
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
		    pp = point_add(get_point(i, j, k), self->dir_to_point[dirz]);
		    pp = point_add(pp, self->dir_to_point[dirx]);
		    if (sim_active_p(self, pp)) {
			point_l = point_offset(get_point(i, j, k), dirx, dir, dirz);
			l = self->lambda_ary[point_l.i][point_l.j][point_l.k];
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
		    pp = point_add(get_point(i, j, k), self->dir_to_point[dirx]);
		    pp = point_add(pp, self->dir_to_point[diry]);
		    if (sim_active_p(self, pp)) {
			point_l = point_offset(get_point(i, j, k), dirx, diry, dir);
			l = self->lambda_ary[point_l.i][point_l.j][point_l.k];
			value += l/dz*(dx/2)*(dy/2);
		    }
		}
	    }
	    break;
	default:
	    bug("sim_set_matrix_coef");
	}
	self->coefs[i][j][k]->coef[dir]->value = value;
    }
}

static void sim_set_matrix(Sim *self)
{
    Point *p;
    double dx, dy, dz;

    ALLOCATE_3D2(self->u, double, self->ni, self->nj, self->nk, 0.0);
    for (p = world_each_begin(self->world); p != NULL; p = world_each(self->world)) {
	if (self->fix_ary[p->i][p->j][p->k] >= 0.0)
	    self->u[p->i][p->j][p->k] = self->fix_ary[p->i][p->j][p->k];
    }

    ALLOCATE_3D(self->coefs, Coefs *, self->ni, self->nj, self->nk);
    for (p = world_each_begin(self->world); p != NULL; p = world_each(self->world)) {
	self->coefs[p->i][p->j][p->k] = coefs_new();
    }

    dx = self->world->dx;
    dy = self->world->dy;
    dz = self->world->dz;

    for (p = world_each_begin(self->world); p != NULL; p = world_each(self->world)) {
	if (self->fix_ary[p->i][p->j][p->k] >= 0.0)
	    continue;

	sim_set_matrix_const(self, p->i, p->j, p->k, dx, dy, dz);
	sim_set_matrix_coef0(self, p->i, p->j, p->k, dx, dy, dz);
	sim_set_matrix_coef(self, p->i, p->j, p->k, dx, dy, dz);
    }
}

Sim *sim_new(void)
{
    Sim *self;

    self = EALLOC(Sim);
    if (opt_v)
	warn("configuring ...");
    self->config = config_new();
    self->dir_to_point[DIR_LEFT]  = get_point(-1,  0,  0);
    self->dir_to_point[DIR_RIGHT] = get_point( 1,  0,  0);
    self->dir_to_point[DIR_FRONT] = get_point( 0, -1,  0);
    self->dir_to_point[DIR_BACK]  = get_point( 0,  1,  0);
    self->dir_to_point[DIR_BELOW] = get_point( 0,  0, -1);
    self->dir_to_point[DIR_ABOVE] = get_point( 0,  0,  1);

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
    Point *p, *hfp;
    Point point;
    int index, index2;
    int idir, dir;
    double c;
    double *sol;
    Array3Dd ary;

    nindex = self->ni * self->nj * self->nk;
    solver = solvele_new(nindex);
    for (p = world_each_begin(self->world); p != NULL; p = world_each(self->world)) {
	hfp = self->heatflow_point_ary[p->i][p->j][p->k];
	if (self->fix_ary[p->i][p->j][p->k] >= 0.0 || !self->active_p_ary[p->i][p->j][p->k]) {
	    index = world_to_index(self->world, *p);
	    solvele_set_matrix(solver, index, index, 1.0);
	    solvele_set_vector(solver, index, self->u[p->i][p->j][p->k]);
	    continue;
	}
	if (hfp != NULL) {
	    if (!point_eq(*p, *hfp)) {
		index = world_to_index(self->world, *p);
		solvele_set_matrix(solver, index, index, 1.0);
		solvele_set_matrix(solver, index, world_to_index(self->world, *hfp), -1.0);
		solvele_set_vector(solver, index, 0.0);
	    }
	    point = *hfp;
	} else {
	    point = *p;
	}
	for (idir = 0; idir < NELEMS(dir_array); ++idir) {
	    dir = dir_array[idir];
	    index = world_to_index(self->world, point);
	    index2 = self->coefs[p->i][p->j][p->k]->coef[dir]->index;
	    if (index2 < 0)
		continue;
	    c = self->coefs[p->i][p->j][p->k]->coef[dir]->value;
	    solvele_add_matrix(solver, index, index2, c);
	}
	solvele_add_matrix(solver, index, index, self->coefs[p->i][p->j][p->k]->coef0);
	solvele_add_vector(solver, index, -self->coefs[p->i][p->j][p->k]->cnst);
    }

    if (opt_v)
	warn("solving equations ...");
    sol = solvele_solve(solver, self->ni, self->nj, self->nk);

    ALLOCATE_3D(ary, double, self->ni, self->nj, self->nk);
    for (p = world_each_begin(self->world); p != NULL; p = world_each(self->world)) {
	ary[p->i][p->j][p->k] = sol[world_to_index(self->world, *p)];
    }

    return ary;
}
