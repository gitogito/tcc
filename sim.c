#include <stdio.h>
#include <stdlib.h>
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
    if (x < 0.0)
	return (int) (x - 0.5);
    else
	return (int) (x + 0.5);
}

static double min_d(double a, double b)
{
    if (a < b)
        return a;
    else
        return b;
}

static double max_d(double a, double b)
{
    if (a > b)
        return a;
    else
        return b;
}

/* Vector2d */

static Vector2d vector2d_sub(Vector2d va, Vector2d vb)
{
    Vector2d v;

    v.x = va.x - vb.x;
    v.y = va.y - vb.y;
    return v;
}

static double vector2d_outer_prod(Vector2d va, Vector2d vb)
{
    return va.x * vb.y - va.y * vb.x;
}

static int vector2d_counter_clock_p(Vector2d va, Vector2d vb, Vector2d vc)
{
    return vector2d_outer_prod(vector2d_sub(va, vb), vector2d_sub(vc, vb)) <= 0.0;
}

/* Vector2d_ary */

Vector2d_ary *vector2d_ary_new(void)
{
    Vector2d_ary *self;

    self = EALLOC(Vector2d_ary);
    self->size = 0;
    self->ptr = NULL;
    return self;
}

void vector2d_ary_push(Vector2d_ary *self, Vector2d vector2d)
{
    ++(self->size);
    self->ptr = erealloc(self->ptr, sizeof(Vector2d) * self->size);
    self->ptr[self->size - 1] = vector2d;
}

void vector2d_ary_delete_at(Vector2d_ary *self, int index)
{
    int i;

    for (i = index + 1; i < self->size; ++i)
	self->ptr[i - 1] = self->ptr[i];
    --(self->size);
}

void vector2d_ary_rotate(Vector2d_ary *self, int index)
{
    Vector2d *ptr;
    int i;

    ptr = EALLOCN(Vector2d, self->size);
    for (i = 0; i < self->size; ++i)
	ptr[i] = self->ptr[(i + index) % self->size];
    self->ptr = ptr;
}

/* iPoint */

iPoint get_ipoint(int i, int j, int k)
{
    iPoint ipoint;

    ipoint.i = i;
    ipoint.j = j;
    ipoint.k = k;
    return ipoint;
}

static iPoint ipoint_offset(iPoint ipoint, int dirx, int diry, int dirz)
{
    int i, j, k;

    switch (dirx) {
    case DIR_LEFT:
	i = ipoint.i - 1;
	break;
    case DIR_RIGHT:
	i = ipoint.i;
	break;
    default:
	bug("unknown dir %d", dirx);
    }

    switch (diry) {
    case DIR_FRONT:
	j = ipoint.j - 1;
	break;
    case DIR_BACK:
	j = ipoint.j;
	break;
    default:
	bug("unknown dir %d", diry);
    }

    switch (dirz) {
    case DIR_BELOW:
	k = ipoint.k - 1;
	break;
    case DIR_ABOVE:
	k = ipoint.k;
	break;
    default:
	bug("unknown dir %d", dirz);
    }

    return get_ipoint(i, j, k);
}

static int ipoint_eq(iPoint ipoint1, iPoint ipoint2)
{
    if (ipoint1.i == ipoint2.i && ipoint1.j == ipoint2.j && ipoint1.k == ipoint2.k)
	return 1;
    else
	return 0;
}

static iPoint ipoint_add(iPoint ipoint1, iPoint ipoint2)
{
    return get_ipoint(
	    ipoint1.i + ipoint2.i,
	    ipoint1.j + ipoint2.j,
	    ipoint1.k + ipoint2.k
	    );
}

static iPoint *ipoint_new(int i, int j, int k)
{
    iPoint *self;

    self = EALLOC(iPoint);
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

static Each *each_new(int size, iPoint *ipoint_ary)
{
    Each *self;

    self = EALLOC(Each);
    self->size = size;
    self->ipoint_ary = ipoint_ary;
    self->index = 0;
    return self;
}

static iPoint *each_each(Each *self)
{
    if (self == NULL)
	bug("each_each");
    if (self->index < self->size) {
	return &self->ipoint_ary[self->index++];
    } else {
	self->index = -1;
	return NULL;
    }
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

static iPoint *world_each_begin(World *self)
{
    int size;
    iPoint *ipoint_ary;
    int n;
    int i, j, k;

    if (self->each != NULL && self->each->index >= 0)
	bug("world_each_begin");
    size = self->ni * self->nj * self->nk;
    ipoint_ary = EALLOCN(iPoint, size);
    n = 0;
    for (i = 0; i < self->ni; ++i) {
	for (j = 0; j < self->nj; ++j) {
	    for (k = 0; k < self->nk; ++k) {
		ipoint_ary[n++] = get_ipoint(i, j, k);
	    }
	}
    }
    self->each = each_new(size, ipoint_ary);
    return each_each(self->each);
}

static iPoint *world_each(World *self)
{
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

/* Sweep */

Sweep *sweep_new(World *world, int axis, double len, Obj *obj)
{
    Sweep *self;

    self = EALLOC(Sweep);
    self->world = world;
    self->axis = axis;
    self->len = len;
    self->obj = obj;
    self->each = NULL;
    return self;
}

static iPoint *obj_each_begin(Obj *self);
static iPoint *obj_each(Obj *self);
static int obj_each_size(Obj *obj);

static iPoint *sweep_each_begin(Sweep *self)
{
    iPoint *p;
    int leni;
    int size;
    iPoint *ipoint_ary;
    int n;
    int i, j, k;

    if (self->each != NULL && self->each->index >= 0)
	bug("sweep_each_begin");
    p = obj_each_begin(self->obj);
    switch (self->axis) {
    case AXIS_X:
	leni = iround(self->len / self->world->dx) + 1;
	size = obj_each_size(self->obj) * leni;
	ipoint_ary = EALLOCN(iPoint, size);
	n = 0;
	do {
	    for (i = p->i; i < p->i + leni; ++i) {
		ipoint_ary[n++] = get_ipoint(i, p->j, p->k);
	    }
	    p = obj_each(self->obj);
	} while (p != NULL);
	break;
    case AXIS_Y:
	leni = iround(self->len / self->world->dy) + 1;
	size = obj_each_size(self->obj) * leni;
	ipoint_ary = EALLOCN(iPoint, size);
	n = 0;
	do {
	    for (j = p->j; j < p->j + leni; ++j) {
		ipoint_ary[n++] = get_ipoint(p->i, j, p->k);
	    }
	    p = obj_each(self->obj);
	} while (p != NULL);
	break;
    case AXIS_Z:
	leni = iround(self->len / self->world->dz) + 1;
	size = obj_each_size(self->obj) * leni;
	ipoint_ary = EALLOCN(iPoint, size);
	n = 0;
	do {
	    for (k = p->k; k < p->k + leni; ++k) {
		ipoint_ary[n++] = get_ipoint(p->i, p->j, k);
	    }
	    p = obj_each(self->obj);
	} while (p != NULL);
	break;
    default:
	bug("unknow axis %d", self->axis);
    }
    assert(n == size);
    self->each = each_new(size, ipoint_ary);
    return each_each(self->each);
}

static iPoint *sweep_each(Sweep *self)
{
    return each_each(self->each);
}

static void sweep_offset(Sweep *self)
{
    switch (self->axis) {
    case AXIS_X:
	self->len -= self->world->dx;
	break;
    case AXIS_Y:
	self->len -= self->world->dy;
	break;
    case AXIS_Z:
	self->len -= self->world->dz;
	break;
    default:
	bug("unknow axis %d", self->axis);
    }
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

static iPoint *rect_each_begin(Rect *self)
{
    int xi, yi, zi;
    int size;
    int len1, len2;
    iPoint *ipoint_ary;
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
	ipoint_ary = EALLOCN(iPoint, len1 * len2);
	for (j = yi; j < yi + len1; ++j) {
	    for (k = zi; k < zi + len2; ++k) {
		ipoint_ary[size++] = get_ipoint(xi, j, k);
	    }
	}
	break;
    case AXIS_Y:
	len1 = iround(self->len1 / self->world->dx) + 1;
	len2 = iround(self->len2 / self->world->dz) + 1;
	ipoint_ary = EALLOCN(iPoint, len1 * len2);
	for (i = xi; i < xi + len1; ++i) {
	    for (k = zi; k < zi + len2; ++k) {
		ipoint_ary[size++] = get_ipoint(i, yi, k);
	    }
	}
	break;
    case AXIS_Z:
	len1 = iround(self->len1 / self->world->dx) + 1;
	len2 = iround(self->len2 / self->world->dy) + 1;
	ipoint_ary = EALLOCN(iPoint, len1 * len2);
	for (i = xi; i < xi + len1; ++i) {
	    for (j = yi; j < yi + len2; ++j) {
		ipoint_ary[size++] = get_ipoint(i, j, zi);
	    }
	}
	break;
    default:
	bug("unknow axis %d for rect_each_begin", self->axis);
    }

    assert(size == len1 * len2);

    self->each = each_new(size, ipoint_ary);
    return each_each(self->each);
}

static iPoint *rect_each(Rect *self)
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

/* Triangle_z */

/*
 *       (x1 + dx2, y1 + dy2)
 *         **
 *       **  ***
 *     **       ***
 *   *****************
 *   <---------------> dx
 * (x1, y1)
 */
static Triangle_z *triangle_z_new(World *world, double x1, double y1, double dx, double dx2, double dy2)
{
    Triangle_z *self;

    self = EALLOC(Triangle_z);
    self->world = world;
    self->x1 = x1;
    self->y1 = y1;
    self->dx = dx;
    self->dx2 = dx2;
    self->dy2 = dy2;
    self->each = NULL;
    return self;
}

static iPoint *triangle_z_each_begin(Triangle_z *self)
{
    int i1, j1, i2, j2, ix, jx;
    double a12, b12, ax2, bx2;
    int size;
    int i, j, i12, ix2, n;
    int istart, iend, jstart, jend;
    iPoint *ipoint_ary;

    if (self->each != NULL && self->each->index >= 0)
	bug("triangle_z_each_begin");

    i1 = iround((self->x1 - self->world->x0) / self->world->dx);
    j1 = iround((self->y1 - self->world->y0) / self->world->dy);
    i2 = iround((self->x1 + self->dx2 - self->world->x0) / self->world->dx);
    j2 = iround((self->y1 + self->dy2 - self->world->y0) / self->world->dy);
    ix = (((self->x1+self->dx) - self->world->x0) / self->world->dx);
    jx = j1;

    a12 = (double) (i1 - i2) / (j1 - j2);
    b12 = (double) (i2*j1 - i1*j2) / (j1 - j2);
    ax2 = (double) (ix - i2) / (jx - j2);
    bx2 = (double) (i2*jx - ix*j2) / (jx - j2);

    if (j1 < j2) {
        jstart = j1;
        jend = j2;
    } else {
        jstart = j2;
        jend = j1;
    }

    size = 0;
    for (j = jstart; j <= jend; ++j) {
        i12 = iround(a12 * j + b12);
        ix2 = iround(ax2 * j + bx2);
        if (i12 < ix2) {
            istart = i12;
            iend = ix2;
        } else {
            istart = ix2;
            iend = i12;
        }
        for (i = istart; i <= iend; ++i)
            ++size;
    }

    ipoint_ary = EALLOCN(iPoint, size);
    n = 0;
    for (j = jstart; j <= jend; ++j) {
        i12 = iround(a12 * j + b12);
        ix2 = iround(ax2 * j + bx2);
        if (i12 < ix2) {
            istart = i12;
            iend = ix2;
        } else {
            istart = ix2;
            iend = i12;
        }
        for (i = istart; i <= iend; ++i)
            ipoint_ary[n++] = get_ipoint(i, j, 0);
    }
    assert(size == n);
    self->each = each_new(size, ipoint_ary);

    return each_each(self->each);
}

static iPoint *triangle_z_each(Triangle_z *self)
{
    return each_each(self->each);
}

/* Triangle */

Triangle *triangle_new(World *world, double x1, double y1, double z1,
	int axis, double du2, double dv2, double du3, double dv3)
{
    Triangle *self;

    self = EALLOC(Triangle);
    self->world = world;
    self->axis = axis;
    switch (self->axis) {
    case AXIS_X:
        self->u1 = y1;
        self->v1 = z1;
        self->wi = iround((x1 - self->world->x0) / self->world->dx);
        break;
    case AXIS_Y:
        self->u1 = x1;
        self->v1 = z1;
        self->wi = iround((y1 - self->world->y0) / self->world->dy);
        break;
    case AXIS_Z:
        self->u1 = x1;
        self->v1 = y1;
        self->wi = iround((z1 - self->world->z0) / self->world->dz);
        break;
    default:
        bug("unknow axis %d", self->axis);
    }
    self->du2 = du2;
    self->dv2 = dv2;
    self->du3 = du3;
    self->dv3 = dv3;
    self->tr1 = NULL;
    self->tr2 = NULL;
    self->each = NULL;
    return self;
}

static iPoint *triangle_each2(Triangle *self, iPoint *p)
{
    int i2, j2, k2;

    assert(p->k == 0);
    switch (self->axis) {
    case AXIS_X:
        i2 = self->wi;
        j2 = p->i;
        k2 = p->j;
        break;
    case AXIS_Y:
        i2 = p->j;
        j2 = self->wi;
        k2 = p->i;
        break;
    case AXIS_Z:
        i2 = p->i;
        j2 = p->j;
        k2 = self->wi;
        break;
    default:
        bug("unknow axis %d", self->axis);
    }
    return ipoint_new(i2, j2, k2);
}

static iPoint *triangle_each_begin(Triangle *self)
{
    double u1, v1, u2, v2, u3, v3;
    double ua, va;
    double ux, vx;
    double ub, vb;
    double a, b, dx;
    iPoint *p1, *p2;
    iPoint *ipoint_ary;
    int size, n;

    if (self->each != NULL && self->each->index >= 0)
	bug("triangle_each_begin");

    u1 = self->u1;
    v1 = self->v1;
    u2 = self->u1 + self->du2;
    v2 = self->v1 + self->dv2;
    u3 = self->u1 + self->du3;
    v3 = self->v1 + self->dv3;

    if (v1 == v2) {
        self->tr1 = triangle_z_new(self->world, u1, v1, u2 - u1, u3 - u1, v3 - v1);
    } else if (v2 == v3) {
        self->tr1 = triangle_z_new(self->world, u2, v2, u3 - u2, u1 - u2, v1 - v2);
    } else if (v3 == v1) {
        self->tr1 = triangle_z_new(self->world, u3, v3, u1 - u3, u2 - u3, v2 - v3);
    } else {
        /*
         *           * pa
         *          * *
         *         *   *
         *        *     *
         *       *  dx   *
         *      *<------->* px
         *     *      ****
         *    *   ****
         *   *****
         *  *
         *  pb
         */
        if (v1 < max_d(v2, v3) && v1 > min_d(v2, v3)) {
            ua = u2;
            va = v2;
            ux = u1;
            vx = v1;
            ub = u3;
            vb = v3;
        } else if (v2 < max_d(v3, v1) && v2 > min_d(v3, v1)) {
            ua = u3;
            va = v3;
            ux = u2;
            vx = v2;
            ub = u1;
            vb = v1;
        } else if (v3 < max_d(v1, v2) && v3 > min_d(v1, v2)) {
            ua = u1;
            va = v1;
            ux = u3;
            vx = v3;
            ub = u2;
            vb = v2;
        } else {
            bug("not reached");
        }
        a = (ua - ub) / (va - vb);
        b = (ub*va - ua*vb) / (va - vb);
        dx = a * vx + b - ux;
        self->tr1 = triangle_z_new(self->world, ux, vx, dx, ua - ux, va - vx);
        self->tr2 = triangle_z_new(self->world, ux, vx, dx, ub - ux, vb - vx);
    }
    assert(self->tr1 != NULL);

    p1 = triangle_z_each_begin(self->tr1);
    if (self->tr2 == NULL) {
        size = self->tr1->each->size;
    } else {
        p2 = triangle_z_each_begin(self->tr2);
        size = self->tr1->each->size + self->tr2->each->size;
    }
    ipoint_ary = EALLOCN(iPoint, size);
    n = 0;
    for (; p1 != NULL; p1 = triangle_z_each(self->tr1))
        ipoint_ary[n++] = *triangle_each2(self, p1);
    if (self->tr2 != NULL) {
        for (; p2 != NULL; p2 = triangle_z_each(self->tr2))
            ipoint_ary[n++] = *triangle_each2(self, p2);
    }
    assert(size == n);

    self->each = each_new(size, ipoint_ary);

    return each_each(self->each);
}

static iPoint *triangle_each(Triangle *self)
{
    return each_each(self->each);
}

/* Polygon */

Polygon *polygon_new(World *world, double x1, double y1, double z1,
	int axis, Vector2d_ary *dudv_ary)
{
    Polygon *self;
    Vector2d v0, v;
    int index;

    self = EALLOC(Polygon);
    self->world = world;
    self->axis = axis;
    self->vector2d_ary = vector2d_ary_new();
    switch (self->axis) {
    case AXIS_X:
	v0.x = y1;
	v0.y = z1;
	self->w = x1;
	break;
    case AXIS_Y:
	v0.x = x1;
	v0.y = z1;
	self->w = y1;
	break;
    case AXIS_Z:
	v0.x = x1;
	v0.y = y1;
	self->w = z1;
	break;
    default:
	bug("unknown axis %d", axis);
    }
    vector2d_ary_push(self->vector2d_ary, v0);
    for (index = 0; index < dudv_ary->size; ++index) {
	v.x = dudv_ary->ptr[index].x + v0.x;
	v.y = dudv_ary->ptr[index].y + v0.y;
	vector2d_ary_push(self->vector2d_ary, v);
    }
    self->each = NULL;
    return self;
}

static iPoint *polygon_each_begin(Polygon *self)
{
    Vector2d_ary *ary;
    int index;
    int size;
    iPoint *ipoint_ary;
    Vector2d va, vb, vc;
    double x1, y1, z1;
    Triangle *tr;
    iPoint *p;

    ary = vector2d_ary_new();
    for (index = 0; index < self->vector2d_ary->size; ++index) {
	vector2d_ary_push(ary, self->vector2d_ary->ptr[index]);
    }

    size = 0;
    ipoint_ary = NULL;
    while (ary->size >= 3) {
	va = ary->ptr[0];
	vb = ary->ptr[1];
	vc = ary->ptr[2];
	if (!vector2d_counter_clock_p(va, vb, vc))
	    warn_exit("points of polygon must be placed counterclockwise");
	switch (self->axis) {
	case AXIS_X:
	    x1 = va.x;
	    y1 = va.y;
	    z1 = self->w;
	    break;
	case AXIS_Y:
	    x1 = va.x;
	    y1 = self->w;
	    z1 = va.y;
	    break;
	case AXIS_Z:
	    x1 = va.x;
	    y1 = va.y;
	    z1 = self->w;
	    break;
	default:
	    bug("unknow axis %d", self->axis);
	}
	tr = triangle_new(self->world, x1, y1, z1, self->axis,
		vb.x - x1, vb.y - y1, vc.x - x1, vc.y - y1);
	for (p = triangle_each_begin(tr); p != NULL; p = triangle_each(tr)) {
	    ++size;
	    ipoint_ary = erealloc(ipoint_ary, sizeof(iPoint) * size);
	    ipoint_ary[size - 1] = *p;
	}
	vector2d_ary_delete_at(ary, 1);
	vector2d_ary_rotate(ary, 1);
    }
    self->each = each_new(size, ipoint_ary);
    return each_each(self->each);
}

static iPoint *polygon_each(Polygon *self)
{
    return each_each(self->each);
}

/* Box */

Box *box_new(World *world, double x, double y, double z, double xlen, double ylen, double zlen)
{
    Box *self;
    Obj *obj;

    if (xlen <= 0.0 || ylen <= 0.0 || zlen <= 0.0)
	warn_exit("length is negative for Rect");

    self = EALLOC(Box);
    self->world = world;
    obj = obj_new(OBJ_RECT);
    obj->uobj.rect = rect_new(world, x, y, z, AXIS_Z, xlen, ylen);
    self->sweep = sweep_new(self->world, AXIS_Z, zlen, obj);
    return self;
}

static iPoint *box_each_begin(Box *self)
{
    return sweep_each_begin(self->sweep);
}

static iPoint *box_each(Box *self)
{
    return sweep_each(self->sweep);
}

static void box_offset(Box *self)
{
    sweep_offset(self->sweep);
}

/* Obj */

Obj *obj_new(int objtype)
{
    Obj *self;

    self = EALLOC(Obj);
    self->objtype = objtype;
    return self;
}

static iPoint *obj_each_begin(Obj *self)
{
    iPoint *p;

    switch (self->objtype) {
    case OBJ_RECT:
	p = rect_each_begin(self->uobj.rect);
	break;
    case OBJ_TRIANGLE:
	p = triangle_each_begin(self->uobj.triangle);
	break;
    case OBJ_POLYGON:
	p = polygon_each_begin(self->uobj.polygon);
	break;
    case OBJ_BOX:
	p = box_each_begin(self->uobj.box);
	break;
    case OBJ_SWEEP:
	p = sweep_each_begin(self->uobj.sweep);
	break;
    default:
	bug("unknown obj %d", self->objtype);
    }
    return p;
}

static iPoint *obj_each(Obj *self)
{
    iPoint *p;

    switch (self->objtype) {
    case OBJ_RECT:
	p = rect_each(self->uobj.rect);
	break;
    case OBJ_TRIANGLE:
	p = triangle_each(self->uobj.triangle);
	break;
    case OBJ_POLYGON:
	p = polygon_each(self->uobj.polygon);
	break;
    case OBJ_BOX:
	p = box_each(self->uobj.box);
	break;
    case OBJ_SWEEP:
	p = sweep_each(self->uobj.sweep);
	break;
    default:
	bug("unknown obj %d", self->objtype);
    }
    return p;
}

static int obj_each_size(Obj *self)
{
    int size;

    switch (self->objtype) {
    case OBJ_RECT:
	size = self->uobj.rect->each->size;
	break;
    case OBJ_TRIANGLE:
	size = self->uobj.triangle->each->size;
	break;
    case OBJ_POLYGON:
	size = self->uobj.polygon->each->size;
	break;
    case OBJ_BOX:
	size = self->uobj.box->sweep->each->size;
	break;
    case OBJ_SWEEP:
	size = self->uobj.sweep->each->size;
	break;
    default:
	bug("unknown obj %d", self->objtype);
    }
    return size;
}

void obj_offset(Obj *self)
{
    switch (self->objtype) {
    case OBJ_RECT:
	rect_offset(self->uobj.rect);
	break;
    case OBJ_TRIANGLE:
        warn_exit("not yet obj_offset for triangle");
	/* triangle_offset(self->uobj.triangle); */
	break;
    case OBJ_POLYGON:
        warn_exit("not yet obj_offset for polygon");
	/* polygon_offset(self->uobj.polygon); */
	break;
    case OBJ_BOX:
	box_offset(self->uobj.box);
	break;
    case OBJ_SWEEP:
	sweep_offset(self->uobj.sweep);
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
    int ok, idir;

    ok = 0;
    for (idir = 0; idir < NELEMS(dir_array); ++idir) {
	if (dir == dir_array[idir]) {
	    ok = 1;
	    break;
	}
    }
    if (!ok)
	bug("unknown dir %d", dir);

    if (value <= 0.0)
	warn_exit("value is negative for Heatflow");

    self = EALLOC(Heatflow);
    self->dir = dir;
    self->value = value;
    return self;
}

/* Config */

static void config_parse(Config *self)
{
    extern int yyparse();
    extern int yydebug;

    yydebug = opt_y;

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

int sim_active_p(Sim *self, iPoint ipoint)
{
    if (ipoint.i < 0 || ipoint.i >= self->ni ||
	    ipoint.j < 0 || ipoint.j >= self->nj ||
	    ipoint.k < 0 || ipoint.k >= self->nk)
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

    ALLOCATE_3D2(self->active_p_ary, int, self->ni, self->nj, self->nk, 0);
    active_obj_ary = self->config->active_obj_ary;
    for (index = 0; index < active_obj_ary->size; ++index) {
	obj = active_obj_ary->ptr[index];
	for (p = obj_each_begin(obj); p != NULL; p = obj_each(obj)) {
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

    ALLOCATE_3D2(self->fix_ary, double, self->ni, self->nj, self->nk, -1.0);
    obj_ary = self->config->fix_obj_ary;
    for (i = 0; i < obj_ary->size; ++i) {
	obj = obj_ary->ptr[i];
	for (p = obj_each_begin(obj); p != NULL; p = obj_each(obj)) {
	    assert(obj->uval.d >= 0.0);
	    if (!world_inside_p(self->world, *p))
		continue;
	    self->fix_ary[p->i][p->j][p->k] = obj->uval.d;
	}
    }
}

static void sim_set_region_heatflow(Sim *self)
{
    iPoint *p;
    int idir, dir;
    AryObj *obj_ary;
    int index;
    Obj *obj;
    int first;
    iPoint ipoint0;

    ALLOCATE_3D(self->heatflow_ary, double *, self->ni, self->nj, self->nk);
    for (p = world_each_begin(self->world); p != NULL; p = world_each(self->world)) {
	self->heatflow_ary[p->i][p->j][p->k] = EALLOCN(double, NELEMS(dir_array));
	for (idir = 0; idir < NELEMS(dir_array); ++idir) {
	    dir = dir_array[idir];
	    self->heatflow_ary[p->i][p->j][p->k][dir] = -1.0;
	}
    }

    ALLOCATE_3D2(self->heatflow_ipoint_ary, iPoint *, self->ni, self->nj, self->nk, NULL);
    obj_ary = self->config->heatflow_obj_ary;
    for (index = 0; index < obj_ary->size; ++index) {
	obj = obj_ary->ptr[index];
	first = 1;
	for (p = obj_each_begin(obj); p != NULL; p = obj_each(obj)) {
	    if (!world_inside_p(self->world, *p))
		continue;
	    if (first) {
		first = 0;
		ipoint0 = *p;
	    }
	    self->heatflow_ipoint_ary[p->i][p->j][p->k] = ipoint_new(ipoint0.i, ipoint0.j, ipoint0.k);
	    self->heatflow_ary[p->i][p->j][p->k][obj->uval.h->dir] = obj->uval.h->value;
	}
    }
}

static void sim_set_region_lambda(Sim *self)
{
    int index;
    AryObj *obj_ary;
    Obj *obj;
    iPoint *p;

    ALLOCATE_3D2(self->lambda_ary, double, self->ni, self->nj, self->nk, 1.0);
    obj_ary = self->config->lambda_obj_ary;
    for (index = 0; index < obj_ary->size; ++index) {
	obj = obj_ary->ptr[index];
	for (p = obj_each_begin(obj); p != NULL; p = obj_each(obj)) {
	    if (!world_inside_p(self->world, *p))
		continue;
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

static void sim_set_matrix_const(Sim *self, int i, int j, int k)
{
    int idir, dir;

    for (idir = 0; idir < NELEMS(dir_array); ++idir) {
	dir = dir_array[idir];
	if (self->heatflow_ary[i][j][k][dir] >= 0.0 &&
		ipoint_eq(*(self->heatflow_ipoint_ary[i][j][k]), get_ipoint(i, j, k)))
	    self->coefs[i][j][k]->cnst = self->heatflow_ary[i][j][k][dir];
    }
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
		    pp = ipoint_add(get_ipoint(i, j, k), self->dir_to_ipoint[dirz]);
		    pp = ipoint_add(pp, self->dir_to_ipoint[dirx]);
		    if (sim_active_p(self, pp)) {
			ipoint_l = ipoint_offset(get_ipoint(i, j, k), dirx, dir, dirz);
			l = self->lambda_ary[ipoint_l.i][ipoint_l.j][ipoint_l.k];
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
		    pp = ipoint_add(get_ipoint(i, j, k), self->dir_to_ipoint[dirx]);
		    pp = ipoint_add(pp, self->dir_to_ipoint[diry]);
		    if (sim_active_p(self, pp)) {
			ipoint_l = ipoint_offset(get_ipoint(i, j, k), dirx, diry, dir);
			l = self->lambda_ary[ipoint_l.i][ipoint_l.j][ipoint_l.k];
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
    iPoint pp, ipoint_l;
    double l;

    for (idir = 0; idir < NELEMS(dir_array); ++idir) {
	dir = dir_array[idir];
	if (!sim_active_p(self, ipoint_add(get_ipoint(i, j, k), self->dir_to_ipoint[dir]))) {
	    self->coefs[i][j][k]->coef[dir]->index = -1;
	    continue;
	}
	if (self->heatflow_ary[i][j][k][dir] >= 0.0)
	    warn_exit("heatflow(%g) comes from active cell at (%d, %d, %d), dir(%d)", self->heatflow_ary[i][j][k][dir], i, j, k, dir);

	self->coefs[i][j][k]->coef[dir]->index =
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
	self->coefs[i][j][k]->coef[dir]->value = value;
    }
}

static void sim_set_matrix(Sim *self)
{
    iPoint *p;
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

	sim_set_matrix_const(self, p->i, p->j, p->k);
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

    nindex = self->ni * self->nj * self->nk;
    solver = solvele_new(nindex);
    for (p = world_each_begin(self->world); p != NULL; p = world_each(self->world)) {
	hfp = self->heatflow_ipoint_ary[p->i][p->j][p->k];
	if (self->fix_ary[p->i][p->j][p->k] >= 0.0 || !self->active_p_ary[p->i][p->j][p->k]) {
	    index = world_to_index(self->world, *p);
	    solvele_set_matrix(solver, index, index, 1.0);
	    solvele_set_vector(solver, index, self->u[p->i][p->j][p->k]);
	    continue;
	}
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
	for (idir = 0; idir < NELEMS(dir_array); ++idir) {
	    dir = dir_array[idir];
	    index = world_to_index(self->world, ipoint);
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
