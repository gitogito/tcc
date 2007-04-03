#include <assert.h>
#include <math.h>
#include "sim.h"
#include "mem.h"
#include "tc.h"

int iround(double x)
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

static int alloc_size(int size)
{
    if (size < 10)
	return 10;
    else
	return (int) (size * 1.2);
}

/* iPoint_ary */

iPoint_ary * ipoint_ary_new(void)
{
    iPoint_ary *self;

    self = EALLOC(iPoint_ary);
    self->size = 0;
    self->alloc_size = 0;
    self->ptr = NULL;
    return self;
}

void ipoint_ary_push(iPoint_ary *self, iPoint ipoint)
{
    ++(self->size);
    if (self->size > self->alloc_size) {
	self->alloc_size = alloc_size(self->size);
	self->ptr = (iPoint *) erealloc(self->ptr, sizeof(iPoint) * self->alloc_size);
    }
    self->ptr[self->size - 1] = ipoint;
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

static int vector2d_inner_triangle_p(Vector2d vp, Vector2d va, Vector2d vb, Vector2d vc)
{
    double val_ab, val_bc, val_ca;

    val_ab = vector2d_outer_prod(vector2d_sub(vb, va), vector2d_sub(vp, va));
    val_bc = vector2d_outer_prod(vector2d_sub(vc, vb), vector2d_sub(vp, vb));
    val_ca = vector2d_outer_prod(vector2d_sub(va, vc), vector2d_sub(vp, vc));
    if ((val_ab >= 0.0 && val_bc >= 0.0 && val_ca >= 0.0) ||
	    (val_ab <= 0.0 && val_bc <= 0.0 && val_ca <= 0.0))
	return 1;
    else
	return 0;
}

/* Vector2d_ary */

Vector2d_ary *vector2d_ary_new(void)
{
    Vector2d_ary *self;

    self = EALLOC(Vector2d_ary);
    self->size = 0;
    self->alloc_size = 0;
    self->ptr = NULL;
    return self;
}

void vector2d_ary_push(Vector2d_ary *self, Vector2d vector2d)
{
    ++(self->size);
    if (self->size > self->alloc_size) {
	self->alloc_size = alloc_size(self->size);
	self->ptr = (Vector2d *) erealloc(self->ptr, sizeof(Vector2d) * self->alloc_size);
    }
    self->ptr[self->size - 1] = vector2d;
}

static void vector2d_ary_delete_at(Vector2d_ary *self, int index)
{
    int i;

    for (i = index + 1; i < self->size; ++i)
	self->ptr[i - 1] = self->ptr[i];
    --(self->size);
}

static void vector2d_ary_rotate(Vector2d_ary *self, int index)
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

iPoint ipoint_offset(iPoint ipoint, int dirx, int diry, int dirz)
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

int ipoint_eq(iPoint ipoint1, iPoint ipoint2)
{
    if (ipoint1.i == ipoint2.i && ipoint1.j == ipoint2.j && ipoint1.k == ipoint2.k)
	return 1;
    else
	return 0;
}

iPoint ipoint_add(iPoint ipoint1, iPoint ipoint2)
{
    return get_ipoint(
	    ipoint1.i + ipoint2.i,
	    ipoint1.j + ipoint2.j,
	    ipoint1.k + ipoint2.k
	    );
}

iPoint *ipoint_new(int i, int j, int k)
{
    iPoint *self;

    self = EALLOC(iPoint);
    self->i = i;
    self->j = j;
    self->k = k;
    return self;
}

/* Each */

Each *each_new(iPoint_ary *ipoint_ary)
{
    Each *self;

    self = EALLOC(Each);
    self->ipoint_ary = ipoint_ary;
    self->index = 0;
    return self;
}

int each_each(Each *self, iPoint **pp)
{
    if (self == NULL)
	bug("each_each");
    if (self->index < self->ipoint_ary->size) {
	*pp = &(self->ipoint_ary->ptr[self->index++]);
	return 1;
    } else {
	self->index = -1;
	return 0;
    }
}

/* Sweep */

static int obj_dim(Obj *obj);

Sweep *sweep_new(int axis, double len, Obj *obj)
{
    Sweep *self;

    if (len <= 0.0)
	warn_exit("length is negative for Sweep");
    if (obj_dim(obj) != 2)
	warn_exit("sweep can't take a %dD object", obj_dim(obj));

    self = EALLOC(Sweep);
    self->axis = axis;
    self->len = len;
    self->obj = obj;
    self->each = NULL;
    return self;
}

static int sweep_each(Sweep *self, iPoint **pp)
{
    iPoint *p;
    int leni;
    iPoint_ary *ipoint_ary;
    int i, j, k;
    int ret;

    if (self->each != NULL && self->each->index >= 0)
	return each_each(self->each, pp);

    ret = obj_each(self->obj, &p);
    if (!ret)
	return ret;
    ipoint_ary = ipoint_ary_new();
    switch (self->axis) {
    case AXIS_X:
	leni = iround(self->len / sim->world->dx) + 1;
	do {
	    if (p != NULL) {
		for (i = p->i; i < p->i + leni; ++i) {
		    ipoint_ary_push(ipoint_ary, get_ipoint(i, p->j, p->k));
		}
	    }
	    ret = obj_each(self->obj, &p);
	} while (ret);
	break;
    case AXIS_Y:
	leni = iround(self->len / sim->world->dy) + 1;
	do {
	    if (p != NULL) {
		for (j = p->j; j < p->j + leni; ++j) {
		    ipoint_ary_push(ipoint_ary, get_ipoint(p->i, j, p->k));
		}
	    }
	    ret = obj_each(self->obj, &p);
	} while (ret);
	break;
    case AXIS_Z:
	leni = iround(self->len / sim->world->dz) + 1;
	do {
	    if (p != NULL) {
		for (k = p->k; k < p->k + leni; ++k) {
		    ipoint_ary_push(ipoint_ary, get_ipoint(p->i, p->j, k));
		}
	    }
	    ret = obj_each(self->obj, &p);
	} while (ret);
	break;
    default:
	bug("unknow axis %d", self->axis);
    }
    self->each = each_new(ipoint_ary);
    return each_each(self->each, pp);
}

static void sweep_offset(Sweep *self)
{
    switch (self->axis) {
    case AXIS_X:
	self->len -= sim->world->dx;
	obj_offset(self->obj);
	break;
    case AXIS_Y:
	self->len -= sim->world->dy;
	obj_offset(self->obj);
	break;
    case AXIS_Z:
	self->len -= sim->world->dz;
	obj_offset(self->obj);
	break;
    default:
	bug("unknow axis %d", self->axis);
    }
    if (self->len <= 0.0)
	warn_exit("length of sweep becomes negative");
}

/* Edge */

static void obj_edge(Obj *self);

Edge *edge_new(Obj *obj)
{
    Edge *self;

    self = EALLOC(Edge);
    self->obj = obj;
    obj_edge(self->obj);
    return self;
}

int edge_each(Edge *self, iPoint **pp)
{
    return obj_each(self->obj, pp);
}

static int edge_dim(Edge *self)
{
    return obj_dim(self->obj);
}

/* Rect */

Rect *rect_new(double x, double y, double z, int axis, double len1, double len2)
{
    static int axis_array[] = { AXIS_X, AXIS_Y, AXIS_Z };

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
    self->x = x;
    self->y = y;
    self->z = z;
    self->axis = axis;
    self->len1 = len1;
    self->len2 = len2;
    self->each = NULL;
    self->edge = 0;
    return self;
}

static int rect_each(Rect *self, iPoint **pp)
{
    int xi, yi, zi;
    int len1, len2;
    iPoint_ary *ipoint_ary;
    int i, j, k;

    if (self->each != NULL && self->each->index >= 0)
	return each_each(self->each, pp);

    xi = iround((self->x - sim->world->x0) / sim->world->dx);
    yi = iround((self->y - sim->world->y0) / sim->world->dy);
    zi = iround((self->z - sim->world->z0) / sim->world->dz);
    ipoint_ary = ipoint_ary_new();
    switch (self->axis) {
    case AXIS_X:
	len1 = iround(self->len1 / sim->world->dy) + 1;
	len2 = iround(self->len2 / sim->world->dz) + 1;
	if (self->edge) {
	    for (j = yi; j < yi + len1; ++j) {
		ipoint_ary_push(ipoint_ary, get_ipoint(xi, j, zi           ));
		ipoint_ary_push(ipoint_ary, get_ipoint(xi, j, zi + len2 - 1));
	    }
	    for (k = zi + 1; k < zi + len2 - 1; ++k) {
		ipoint_ary_push(ipoint_ary, get_ipoint(xi, yi           , k));
		ipoint_ary_push(ipoint_ary, get_ipoint(xi, yi + len1 - 1, k));
	    }
	} else {
	    for (j = yi; j < yi + len1; ++j) {
		for (k = zi; k < zi + len2; ++k) {
		    ipoint_ary_push(ipoint_ary, get_ipoint(xi, j, k));
		}
	    }
	}
	break;
    case AXIS_Y:
	len1 = iround(self->len1 / sim->world->dx) + 1;
	len2 = iround(self->len2 / sim->world->dz) + 1;
	if (self->edge) {
	    for (i = xi; i < xi + len1; ++i) {
		ipoint_ary_push(ipoint_ary, get_ipoint(i, yi, zi           ));
		ipoint_ary_push(ipoint_ary, get_ipoint(i, yi, zi + len2 - 1));
	    }
	    for (k = zi + 1; k < zi + len2 - 1; ++k) {
		ipoint_ary_push(ipoint_ary, get_ipoint(xi           , yi, k));
		ipoint_ary_push(ipoint_ary, get_ipoint(xi + len1 - 1, yi, k));
	    }
	} else {
	    for (i = xi; i < xi + len1; ++i) {
		for (k = zi; k < zi + len2; ++k) {
		    ipoint_ary_push(ipoint_ary, get_ipoint(i, yi, k));
		}
	    }
	}
	break;
    case AXIS_Z:
	len1 = iround(self->len1 / sim->world->dx) + 1;
	len2 = iround(self->len2 / sim->world->dy) + 1;
	if (self->edge) {
	    for (i = xi; i < xi + len1; ++i) {
		ipoint_ary_push(ipoint_ary, get_ipoint(i, yi           , zi));
		ipoint_ary_push(ipoint_ary, get_ipoint(i, yi + len2 - 1, zi));
	    }
	    for (j = yi + 1; j < yi + len2 - 1; ++j) {
		ipoint_ary_push(ipoint_ary, get_ipoint(xi           , j, zi));
		ipoint_ary_push(ipoint_ary, get_ipoint(xi + len1 - 1, j, zi));
	    }
	} else {
	    for (i = xi; i < xi + len1; ++i) {
		for (j = yi; j < yi + len2; ++j) {
		    ipoint_ary_push(ipoint_ary, get_ipoint(i, j, zi));
		}
	    }
	}
	break;
    default:
	bug("unknow axis %d for rect_each", self->axis);
    }

    self->each = each_new(ipoint_ary);
    return each_each(self->each, pp);
}

static void rect_offset(Rect *self)
{
    switch (self->axis) {
    case AXIS_X:
	self->len1 -= sim->world->dy;
	self->len2 -= sim->world->dz;
	break;
    case AXIS_Y:
	self->len1 -= sim->world->dz;
	self->len2 -= sim->world->dx;
	break;
    case AXIS_Z:
	self->len1 -= sim->world->dx;
	self->len2 -= sim->world->dy;
	break;
    default:
	bug("unknow axis %d", self->axis);
    }
    if (self->len1 <= 0.0 || self->len2 <= 0.0)
	warn_exit("length of rect becomes negative");
}

static void rect_edge(Rect *self)
{
    self->edge = 1;
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
static Triangle_z *triangle_z_new(double x1, double y1, double dx, double dx2, double dy2)
{
    Triangle_z *self;

    self = EALLOC(Triangle_z);
    self->x1 = x1;
    self->y1 = y1;
    self->dx = dx;
    self->dx2 = dx2;
    self->dy2 = dy2;
    self->each = NULL;
    return self;
}

static int triangle_z_each(Triangle_z *self, iPoint **pp)
{
    int i1, j1, i2, j2, ix, jx;
    double a12, b12, ax2, bx2;
    int i, j, i12, ix2;
    int istart, iend, jstart, jend;
    iPoint_ary *ipoint_ary;

    if (self->each != NULL && self->each->index >= 0)
	return each_each(self->each, pp);

    i1 = iround((self->x1 - sim->world->x0) / sim->world->dx);
    j1 = iround((self->y1 - sim->world->y0) / sim->world->dy);
    i2 = iround((self->x1 + self->dx2 - sim->world->x0) / sim->world->dx);
    j2 = iround((self->y1 + self->dy2 - sim->world->y0) / sim->world->dy);
    ix = iround(((self->x1+self->dx) - sim->world->x0) / sim->world->dx);
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

    ipoint_ary = ipoint_ary_new();
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
	    ipoint_ary_push(ipoint_ary, get_ipoint(i, j, 0));
    }
    self->each = each_new(ipoint_ary);

    return each_each(self->each, pp);
}

/* Triangle */

Triangle *triangle_new(double x1, double y1, double z1,
	int axis, double du2, double dv2, double du3, double dv3)
{
    Triangle *self;

    self = EALLOC(Triangle);
    self->axis = axis;
    switch (self->axis) {
    case AXIS_X:
	self->u1 = y1;
	self->v1 = z1;
	self->wi = iround((x1 - sim->world->x0) / sim->world->dx);
	break;
    case AXIS_Y:
	self->u1 = x1;
	self->v1 = z1;
	self->wi = iround((y1 - sim->world->y0) / sim->world->dy);
	break;
    case AXIS_Z:
	self->u1 = x1;
	self->v1 = y1;
	self->wi = iround((z1 - sim->world->z0) / sim->world->dz);
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

static int triangle_each(Triangle *self, iPoint **pp)
{
    double u1, v1, u2, v2, u3, v3;
    double ua, va;
    double ux, vx;
    double ub, vb;
    double a, b, dx;
    iPoint *p1, *p2;
    iPoint_ary *ipoint_ary;
    int ret1, ret2;

    if (self->each != NULL && self->each->index >= 0)
	return each_each(self->each, pp);

    u1 = self->u1;
    v1 = self->v1;
    u2 = self->u1 + self->du2;
    v2 = self->v1 + self->dv2;
    u3 = self->u1 + self->du3;
    v3 = self->v1 + self->dv3;

    if (v1 == v2) {
	self->tr1 = triangle_z_new(u1, v1, u2 - u1, u3 - u1, v3 - v1);
    } else if (v2 == v3) {
	self->tr1 = triangle_z_new(u2, v2, u3 - u2, u1 - u2, v1 - v2);
    } else if (v3 == v1) {
	self->tr1 = triangle_z_new(u3, v3, u1 - u3, u2 - u3, v2 - v3);
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
	self->tr1 = triangle_z_new(ux, vx, dx, ua - ux, va - vx);
	self->tr2 = triangle_z_new(ux, vx, dx, ub - ux, vb - vx);
    }
    assert(self->tr1 != NULL);

    ret1 = triangle_z_each(self->tr1, &p1);
    if (self->tr2 != NULL) {
	ret2 = triangle_z_each(self->tr2, &p2);
    }
    ipoint_ary = ipoint_ary_new();
    for (; ret1; ret1 = triangle_z_each(self->tr1, &p1)) {
	if (p1 == NULL)
	    continue;
	ipoint_ary_push(ipoint_ary, *triangle_each2(self, p1));
    }
    if (self->tr2 != NULL) {
	for (; ret2; ret2 = triangle_z_each(self->tr2, &p2)) {
	    if (p1 == NULL)
		continue;
	    ipoint_ary_push(ipoint_ary, *triangle_each2(self, p2));
	}
    }
    self->each = each_new(ipoint_ary);

    return each_each(self->each, pp);
}

static void triangle_offset(Triangle *self)
{
    switch (self->axis) {
    case AXIS_X:
	self->du2 -= sim->world->dy;
	self->du2 -= sim->world->dz;
	self->du3 -= sim->world->dy;
	self->du3 -= sim->world->dz;
	break;
    case AXIS_Y:
	self->du2 -= sim->world->dz;
	self->dv2 -= sim->world->dx;
	self->du3 -= sim->world->dz;
	self->dv3 -= sim->world->dx;
	break;
    case AXIS_Z:
	self->du2 -= sim->world->dx;
	self->dv2 -= sim->world->dy;
	self->du3 -= sim->world->dx;
	self->dv3 -= sim->world->dy;
	break;
    default:
	bug("unknow axis %d", self->axis);
    }
}

/* Ellipse */

Ellipse *ellipse_new(double x, double y, double z, int axis, double ru, double rv)
{
    Ellipse *self;

    if (ru <= 0.0 || rv <= 0.0)
	warn_exit("length is negative for Ellipse");

    self = EALLOC(Ellipse);
    self->x = x;
    self->y = y;
    self->z = z;
    self->axis = axis;
    self->ru = ru;
    self->rv = rv;
    self->each = NULL;
    self->edge = 0;
    return self;
}

static void ellipse_ipoint_ary_add(Ellipse *self, iPoint_ary *ipoint_ary, int axis,
	int wc, int u1, int u2, int v)
{
    int u;
    iPoint ipoint;

    if (self->edge) {
	switch (axis) {
	case AXIS_X:
	    ipoint = get_ipoint(wc, u1, v);
	    break;
	case AXIS_Y:
	    ipoint = get_ipoint(u1, wc, v);
	    break;
	case AXIS_Z:
	    ipoint = get_ipoint(u1, v, wc);
	    break;
	default:
	    bug("unknown axis %d", axis);
	}
	ipoint_ary_push(ipoint_ary, ipoint);

	switch (axis) {
	case AXIS_X:
	    ipoint = get_ipoint(wc, u2, v);
	    break;
	case AXIS_Y:
	    ipoint = get_ipoint(u2, wc, v);
	    break;
	case AXIS_Z:
	    ipoint = get_ipoint(u2, v, wc);
	    break;
	default:
	    bug("unknown axis %d", axis);
	}
	ipoint_ary_push(ipoint_ary, ipoint);
    } else {
	for (u = u1; u <= u2; ++u) {
	    switch (axis) {
	    case AXIS_X:
		ipoint = get_ipoint(wc, u, v);
		break;
	    case AXIS_Y:
		ipoint = get_ipoint(u, wc, v);
		break;
	    case AXIS_Z:
		ipoint = get_ipoint(u, v, wc);
		break;
	    default:
		bug("unknown axis %d", axis);
	    }
	    ipoint_ary_push(ipoint_ary, ipoint);
	}
    }
}

static int ellipse_each(Ellipse *self, iPoint **pp)
{
    int uc, vc, wc, ru, rv;
    int ui, vi, ri, u1, v1;
    iPoint_ary *ipoint_ary;

    if (self->each != NULL && self->each->index >= 0)
	return each_each(self->each, pp);

    switch (self->axis) {
    case AXIS_X:
	uc = iround((self->y - sim->world->y0) / sim->world->dy);
	vc = iround((self->z - sim->world->z0) / sim->world->dz);
	wc = iround((self->x - sim->world->x0) / sim->world->dx);
	ru = iround(self->ru / sim->world->dy);
	rv = iround(self->rv / sim->world->dz);
	break;
    case AXIS_Y:
	uc = iround((self->x - sim->world->x0) / sim->world->dx);
	vc = iround((self->z - sim->world->z0) / sim->world->dz);
	wc = iround((self->y - sim->world->y0) / sim->world->dy);
	ru = iround(self->ru / sim->world->dx);
	rv = iround(self->rv / sim->world->dz);
	break;
    case AXIS_Z:
	uc = iround((self->x - sim->world->x0) / sim->world->dx);
	vc = iround((self->y - sim->world->y0) / sim->world->dy);
	wc = iround((self->z - sim->world->z0) / sim->world->dz);
	ru = iround(self->ru / sim->world->dx);
	rv = iround(self->rv / sim->world->dy);
	break;
    default:
	bug("unknown axis %d", self->axis);
    }

    if (ru == 0 && rv == 0)
	return 0;

    ipoint_ary = ipoint_ary_new();
    if (ru > rv) {
	ui = ri = ru;  vi = 0;
	while (ui >= vi) {
	    u1 = (int)((long)ui * rv / ru);
	    v1 = (int)((long)vi * rv / ru);
	    ellipse_ipoint_ary_add(self, ipoint_ary, self->axis, wc, uc - ui, uc + ui, vc - v1);
	    ellipse_ipoint_ary_add(self, ipoint_ary, self->axis, wc, uc - ui, uc + ui, vc + v1);
	    ellipse_ipoint_ary_add(self, ipoint_ary, self->axis, wc, uc - vi, uc + vi, vc - u1);
	    ellipse_ipoint_ary_add(self, ipoint_ary, self->axis, wc, uc - vi, uc + vi, vc + u1);
	    if ((ri -= (vi++ << 1) + 1) <= 0)
		ri += (ui-- - 1) << 1;
	}
    } else {
	ui = ri = rv;  vi = 0;
	while (ui >= vi) {
	    u1 = (int)((long)ui * ru / rv);
	    v1 = (int)((long)vi * ru / rv);
	    ellipse_ipoint_ary_add(self, ipoint_ary, self->axis, wc, uc - u1, uc + u1, vc - vi);
	    ellipse_ipoint_ary_add(self, ipoint_ary, self->axis, wc, uc - u1, uc + u1, vc + vi);
	    ellipse_ipoint_ary_add(self, ipoint_ary, self->axis, wc, uc - v1, uc + v1, vc - ui);
	    ellipse_ipoint_ary_add(self, ipoint_ary, self->axis, wc, uc - v1, uc + v1, vc + ui);
	    if ((ri -= (vi++ << 1) - 1) < 0)
		ri += (ui-- - 1) << 1;
	}
    }
    self->each = each_new(ipoint_ary);
    return each_each(self->each, pp);
}

static void ellipse_offset(Ellipse *self)
{
    switch (self->axis) {
    case AXIS_X:
	self->ru -= sim->world->dy;
	self->rv -= sim->world->dz;
	break;
    case AXIS_Y:
	self->ru -= sim->world->dz;
	self->rv -= sim->world->dx;
	break;
    case AXIS_Z:
	self->ru -= sim->world->dx;
	self->rv -= sim->world->dy;
	break;
    default:
	bug("unknow axis %d", self->axis);
    }
    if (self->ru <= 0.0 || self->rv <= 0.0)
	warn_exit("length of ellipse becomes negative");
}

static void ellipse_edge(Ellipse *self)
{
    self->edge = 1;
}

/* Circle */

Circle *circle_new(double x, double y, double z, int axis, double r)
{
    Circle *self;

    if (r <= 0.0)
	warn_exit("length is negative for Circle");

    self = EALLOC(Circle);
    self->ellipse = ellipse_new(x, y, z, axis, r, r);
    return self;
}

static int circle_each(Circle *self, iPoint **pp)
{
    return ellipse_each(self->ellipse, pp);
}

static void circle_offset(Circle *self)
{
    ellipse_offset(self->ellipse);
}

static void circle_edge(Circle *self)
{
    ellipse_edge(self->ellipse);
}

/* Polygon */

Polygon *polygon_new(double x1, double y1, double z1,
	int axis, Vector2d_ary *dudv_ary)
{
    Polygon *self;
    Vector2d v0, v;
    int index;

    self = EALLOC(Polygon);
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

static int polygon_each(Polygon *self, iPoint **pp)
{
    Vector2d_ary *ary;
    int index;
    iPoint_ary *ipoint_ary;
    Vector2d va, vb, vc;
    int ok;
    double x1, y1, z1;
    Triangle *tr;
    iPoint *p;

    if (self->each != NULL && self->each->index >= 0)
	return each_each(self->each, pp);

    ary = vector2d_ary_new();
    for (index = 0; index < self->vector2d_ary->size; ++index) {
	vector2d_ary_push(ary, self->vector2d_ary->ptr[index]);
    }

    ipoint_ary = ipoint_ary_new();
    while (ary->size >= 3) {
	for (index = 0; index < ary->size; ++index) {
	}
	va = ary->ptr[0];
	vb = ary->ptr[1];
	vc = ary->ptr[2];
	if (vector2d_counter_clock_p(va, vb, vc)) {
	    ok = 1;
	    for (index = 3; index < ary->size; ++index) {
		if (vector2d_inner_triangle_p(ary->ptr[index], va, vb, vc)) {
		    ok = 0;
		    break;
		}
	    }
	    if (ok) {
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
		tr = triangle_new(x1, y1, z1, self->axis,
			vb.x - x1, vb.y - y1, vc.x - x1, vc.y - y1);
		while (triangle_each(tr, &p)) {
		    if (p == NULL)
			continue;
		    ipoint_ary_push(ipoint_ary, *p);
		}
		vector2d_ary_delete_at(ary, 1);
	    }
	}
	vector2d_ary_rotate(ary, 1);
    }
    self->each = each_new(ipoint_ary);
    return each_each(self->each, pp);
}

static void polygon_offset(Polygon *self)
{
    int index;

    switch (self->axis) {
    case AXIS_X:
	for (index = 0; index < self->vector2d_ary->size; ++index) {
	    self->vector2d_ary->ptr[index].x -= sim->world->dy;
	    self->vector2d_ary->ptr[index].y -= sim->world->dz;
	}
	break;
    case AXIS_Y:
	for (index = 0; index < self->vector2d_ary->size; ++index) {
	    self->vector2d_ary->ptr[index].x -= sim->world->dx;
	    self->vector2d_ary->ptr[index].y -= sim->world->dz;
	}
	break;
    case AXIS_Z:
	for (index = 0; index < self->vector2d_ary->size; ++index) {
	    self->vector2d_ary->ptr[index].x -= sim->world->dx;
	    self->vector2d_ary->ptr[index].y -= sim->world->dy;
	}
	break;
    default:
	bug("unknow axis %d", self->axis);
    }
}

/* Box */

Box *box_new(double x, double y, double z, double xlen, double ylen, double zlen)
{
    Box *self;
    Obj *obj;

    if (xlen <= 0.0 || ylen <= 0.0 || zlen <= 0.0)
	warn_exit("length is negative for Box");

    self = EALLOC(Box);
    obj = obj_new(OBJ_RECT);
    obj->uobj.rect = rect_new(x, y, z, AXIS_Z, xlen, ylen);
    self->sweep = sweep_new(AXIS_Z, zlen, obj);
    return self;
}

static int box_each(Box *self, iPoint **pp)
{
    return sweep_each(self->sweep, pp);
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

int obj_each(Obj *self, iPoint **pp)
{
    switch (self->objtype) {
    case OBJ_RECT:
	return rect_each(self->uobj.rect, pp);
	break;
    case OBJ_TRIANGLE:
	return triangle_each(self->uobj.triangle, pp);
	break;
    case OBJ_ELLIPSE:
	return ellipse_each(self->uobj.ellipse, pp);
	break;
    case OBJ_CIRCLE:
	return circle_each(self->uobj.circle, pp);
	break;
    case OBJ_POLYGON:
	return polygon_each(self->uobj.polygon, pp);
	break;
    case OBJ_BOX:
	return box_each(self->uobj.box, pp);
	break;
    case OBJ_SWEEP:
	return sweep_each(self->uobj.sweep, pp);
	break;
    case OBJ_EDGE:
	return edge_each(self->uobj.edge, pp);
	break;
    default:
	bug("unknown obj %d", self->objtype);
    }
    return 0;	/* NOTREACHED */
}

void obj_offset(Obj *self)
{
    switch (self->objtype) {
    case OBJ_RECT:
	rect_offset(self->uobj.rect);
	break;
    case OBJ_TRIANGLE:
	triangle_offset(self->uobj.triangle);
	break;
    case OBJ_ELLIPSE:
	ellipse_offset(self->uobj.ellipse);
	break;
    case OBJ_CIRCLE:
	circle_offset(self->uobj.circle);
	break;
    case OBJ_POLYGON:
	polygon_offset(self->uobj.polygon);
	break;
    case OBJ_BOX:
	box_offset(self->uobj.box);
	break;
    case OBJ_SWEEP:
	sweep_offset(self->uobj.sweep);
	break;
    case OBJ_EDGE:
	warn_exit("edge_offset is not implemented");
	break;
    default:
	bug("unknown obj %d", self->objtype);
    }
}

static void obj_edge(Obj *self)
{
    switch (self->objtype) {
    case OBJ_RECT:
	rect_edge(self->uobj.rect);
	break;
    case OBJ_TRIANGLE:
	warn_exit("triangle_edge is not implemented");
	break;
    case OBJ_ELLIPSE:
	ellipse_edge(self->uobj.ellipse);
	break;
    case OBJ_CIRCLE:
	circle_edge(self->uobj.circle);
	break;
    case OBJ_POLYGON:
	warn_exit("polygon_edge is not implemented");
	break;
    case OBJ_BOX:
	warn_exit("box_edge is not implemented");
	break;
    case OBJ_SWEEP:
	warn_exit("sweep_edge is not implemented");
	break;
    case OBJ_EDGE:
	warn_exit("edge_edge is not implemented");
	break;
    default:
	bug("unknown obj %d", self->objtype);
    }
}

static int obj_dim(Obj *self)
{
    switch (self->objtype) {
    case OBJ_RECT:
	return 2;
	break;
    case OBJ_TRIANGLE:
	return 2;
	break;
    case OBJ_ELLIPSE:
	return 2;
	break;
    case OBJ_CIRCLE:
	return 2;
	break;
    case OBJ_POLYGON:
	return 2;
	break;
    case OBJ_BOX:
	return 3;
	break;
    case OBJ_SWEEP:
	return 3;
	break;
    case OBJ_EDGE:
	return edge_dim(self->uobj.edge);
	break;
    default:
	bug("unknown obj %d", self->objtype);
    }
    /* NOTREACHED */
    return -1;
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
    self->ptr = (Obj **) erealloc(self->ptr, sizeof(Obj *) * self->size);
    self->ptr[self->size - 1] = obj;
}
