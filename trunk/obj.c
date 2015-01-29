#include <assert.h>
#include <math.h>
#include "sim.h"
#include "mem.h"
#include "tcc.h"

IP_TYPE iround(double x)
{
    if (x < 0.0)
	return (IP_TYPE) (x - 0.5);
    else
	return (IP_TYPE) (x + 0.5);
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

iPoint_ary *ipoint_ary_new(void)
{
    iPoint_ary *self;

    self = EALLOC(iPoint_ary);
    self->size = 0;
    self->alloc_size = 0;
    self->ptr = NULL;
    return self;
}

void ipoint_ary_free(iPoint_ary *self)
{
    if (self == NULL)
	return;
    FREE(self->ptr);
    FREE(self);
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

static int vector2d_clock_p(Vector2d va, Vector2d vb, Vector2d vc)
{
    return !vector2d_counter_clock_p(va, vb, vc);
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

static void vector2d_ary_free(Vector2d_ary *self)
{
    if (self == NULL)
	return;
    FREE(self->ptr);
    FREE(self);
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

iPoint get_ipoint(IP_TYPE i, IP_TYPE j, IP_TYPE k)
{
    iPoint ipoint;

    ipoint.i = i;
    ipoint.j = j;
    ipoint.k = k;
    return ipoint;
}

iPoint ipoint_offset(iPoint *ipoint, int dirx, int diry, int dirz)
{
    IP_TYPE i, j, k;

    i = j = k = -1;

    switch (dirx) {
    case DIR_LEFT:
	i = ipoint->i - 1;
	break;
    case DIR_RIGHT:
	i = ipoint->i;
	break;
    default:
	bug("unknown dir %d", dirx);
    }

    switch (diry) {
    case DIR_FRONT:
	j = ipoint->j - 1;
	break;
    case DIR_BACK:
	j = ipoint->j;
	break;
    default:
	bug("unknown dir %d", diry);
    }

    switch (dirz) {
    case DIR_BELOW:
	k = ipoint->k - 1;
	break;
    case DIR_ABOVE:
	k = ipoint->k;
	break;
    default:
	bug("unknown dir %d", dirz);
    }

    assert(i >= 0 && j >= 0 && k >= 0);

    return get_ipoint(i, j, k);
}

int ipoint_eq(iPoint *ipoint1, iPoint *ipoint2)
{
    if (ipoint1->i == ipoint2->i && ipoint1->j == ipoint2->j && ipoint1->k == ipoint2->k)
	return 1;
    else
	return 0;
}

iPoint ipoint_add(iPoint *ipoint1, iPoint *ipoint2)
{
    return get_ipoint(
	    ipoint1->i + ipoint2->i,
	    ipoint1->j + ipoint2->j,
	    ipoint1->k + ipoint2->k
	    );
}

iPoint *ipoint_new(IP_TYPE i, IP_TYPE j, IP_TYPE k)
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

static void each_free(Each *self)
{
    if (self == NULL)
	return;
    ipoint_ary_free(self->ipoint_ary);
    FREE(self);
}

int each_each(Each *self, iPoint **pp)
{
    if (self == NULL)
	bug("each_each");
    if (self->index < self->ipoint_ary->size) {
	*pp = &(self->ipoint_ary->ptr[self->index++]);
	return 1;
    } else {
	self->index = 0;
	return 0;
    }
}

/* Sweep */

static int obj_dim(Obj *obj);

Sweep *sweep_new(int axis, double len, Obj *obj)
{
    Sweep *self;

    if (len < 0.0)
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

void obj_free(Obj *self);

static void sweep_free(Sweep *self)
{
    if (self == NULL)
	return;
    obj_free(self->obj);
    each_free(self->each);
    FREE(self);
}

static int sweep_each(Sweep *self, iPoint **pp)
{
    iPoint *p;
    IP_TYPE leni;
    iPoint_ary *ipoint_ary;
    IP_TYPE i, j, k;
    int ret;

    if (self->each != NULL && self->each->index >= 0)
	return each_each(self->each, pp);

    ret = obj_each(self->obj, &p);
    if (!ret)
	return ret;
    ipoint_ary = ipoint_ary_new();
    switch (self->axis) {
    case AXIS_X:
	leni = iround(self->len / world->dx) + 1;
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
	leni = iround(self->len / world->dy) + 1;
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
	leni = iround(self->len / world->dz) + 1;
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
	self->len -= world->dx;
	obj_offset(self->obj);
	break;
    case AXIS_Y:
	self->len -= world->dy;
	obj_offset(self->obj);
	break;
    case AXIS_Z:
	self->len -= world->dz;
	obj_offset(self->obj);
	break;
    default:
	bug("unknow axis %d", self->axis);
    }
    if (self->len < 0.0)
	warn_exit("length of sweep becomes negative");
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

    if (len1 < 0.0 || len2 < 0.0)
	warn_exit("length is negative for Rect");

    self = EALLOC(Rect);
    self->x = x;
    self->y = y;
    self->z = z;
    self->axis = axis;
    self->len1 = len1;
    self->len2 = len2;
    self->each = NULL;
    return self;
}

static void rect_free(Rect *self)
{
    if (self == NULL)
	return;
    each_free(self->each);
    FREE(self);
}

static int rect_each(Rect *self, iPoint **pp)
{
    IP_TYPE xi, yi, zi;
    IP_TYPE len1, len2;
    iPoint_ary *ipoint_ary;
    IP_TYPE i, j, k;

    if (self->each != NULL && self->each->index >= 0)
	return each_each(self->each, pp);

    xi = iround((self->x - world->x0) / world->dx);
    yi = iround((self->y - world->y0) / world->dy);
    zi = iround((self->z - world->z0) / world->dz);
    ipoint_ary = ipoint_ary_new();
    switch (self->axis) {
    case AXIS_X:
	len1 = iround(self->len1 / world->dy) + 1;
	len2 = iround(self->len2 / world->dz) + 1;
	for (j = yi; j < yi + len1; ++j) {
	    for (k = zi; k < zi + len2; ++k) {
		ipoint_ary_push(ipoint_ary, get_ipoint(xi, j, k));
	    }
	}
	break;
    case AXIS_Y:
	len1 = iround(self->len1 / world->dx) + 1;
	len2 = iround(self->len2 / world->dz) + 1;
	for (i = xi; i < xi + len1; ++i) {
	    for (k = zi; k < zi + len2; ++k) {
		ipoint_ary_push(ipoint_ary, get_ipoint(i, yi, k));
	    }
	}
	break;
    case AXIS_Z:
	len1 = iround(self->len1 / world->dx) + 1;
	len2 = iround(self->len2 / world->dy) + 1;
	for (i = xi; i < xi + len1; ++i) {
	    for (j = yi; j < yi + len2; ++j) {
		ipoint_ary_push(ipoint_ary, get_ipoint(i, j, zi));
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
	self->len1 -= world->dy;
	self->len2 -= world->dz;
	break;
    case AXIS_Y:
	self->len1 -= world->dz;
	self->len2 -= world->dx;
	break;
    case AXIS_Z:
	self->len1 -= world->dx;
	self->len2 -= world->dy;
	break;
    default:
	bug("unknow axis %d", self->axis);
    }
    if (self->len1 < 0.0 || self->len2 < 0.0)
	warn_exit("length of rect becomes negative");
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

static void triangle_z_free(Triangle_z *self)
{
    if (self == NULL)
	return;
    each_free(self->each);
    FREE(self);
}

static void line_set_ipoint_ary(iPoint_ary *ipoint_ary,
	int x1, int y1, int x2, int y2);

static int triangle_z_each(Triangle_z *self, iPoint **pp)
{
    IP_TYPE i1, j1, i2, j2, ix, jx;
    double a12, b12, ax2, bx2;
    IP_TYPE i, j, i12, ix2;
    IP_TYPE istart, iend, jstart, jend;
    iPoint_ary *ipoint_ary;

    if (self->each != NULL && self->each->index >= 0)
	return each_each(self->each, pp);

    i1 = iround((self->x1 - world->x0) / world->dx);
    j1 = iround((self->y1 - world->y0) / world->dy);
    i2 = iround((self->x1 + self->dx2 - world->x0) / world->dx);
    j2 = iround((self->y1 + self->dy2 - world->y0) / world->dy);
    ix = iround(((self->x1+self->dx) - world->x0) / world->dx);
    jx = j1;

    ipoint_ary = ipoint_ary_new();

    if (j1 == j2) {
	/* do nothing */
    } else {
	if (j1 < j2) {
	    jstart = j1;
	    jend = j2;
	} else {	/* j1 > j2 */
	    jstart = j2;
	    jend = j1;
	}

	a12 = (double) (i1 - i2) / (j1 - j2);
	b12 = (double) (i2*j1 - i1*j2) / (j1 - j2);
	ax2 = (double) (ix - i2) / (jx - j2);
	bx2 = (double) (i2*jx - ix*j2) / (jx - j2);

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
    }

    /* three lines of edges */
    line_set_ipoint_ary(ipoint_ary, i1, j1, i2, j2);
    line_set_ipoint_ary(ipoint_ary, i2, j2, ix, jx);
    line_set_ipoint_ary(ipoint_ary, ix, jx, i1, j1);

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
	self->wi = iround((x1 - world->x0) / world->dx);
	break;
    case AXIS_Y:
	self->u1 = x1;
	self->v1 = z1;
	self->wi = iround((y1 - world->y0) / world->dy);
	break;
    case AXIS_Z:
	self->u1 = x1;
	self->v1 = y1;
	self->wi = iround((z1 - world->z0) / world->dz);
	break;
    default:
	bug("unknow axis %d", self->axis);
    }
    self->du2 = du2;
    self->dv2 = dv2;
    self->du3 = du3;
    self->dv3 = dv3;
    self->each = NULL;
    return self;
}

static void triangle_free(Triangle *self)
{
    if (self == NULL)
	return;
    each_free(self->each);
    FREE(self);
}

static iPoint *triangle_each2(Triangle *self, iPoint *p)
{
    IP_TYPE i2, j2, k2;

    i2 = j2 = k2 = -1;
    assert(p->k == 0);
    switch (self->axis) {
    case AXIS_X:
	i2 = self->wi;
	j2 = p->i;
	k2 = p->j;
	break;
    case AXIS_Y:
	i2 = p->i;
	j2 = self->wi;
	k2 = p->j;
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
    Triangle_z *tr1, *tr2;

    if (self->each != NULL && self->each->index >= 0)
	return each_each(self->each, pp);

    u1 = self->u1;
    v1 = self->v1;
    u2 = self->u1 + self->du2;
    v2 = self->v1 + self->dv2;
    u3 = self->u1 + self->du3;
    v3 = self->v1 + self->dv3;

    tr2 = NULL;	/* for shutting up compiler */
    if (v1 == v2) {
	tr1 = triangle_z_new(u1, v1, u2 - u1, u3 - u1, v3 - v1);
    } else if (v2 == v3) {
	tr1 = triangle_z_new(u2, v2, u3 - u2, u1 - u2, v1 - v2);
    } else if (v3 == v1) {
	tr1 = triangle_z_new(u3, v3, u1 - u3, u2 - u3, v2 - v3);
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
	    ua = va = ux = vx = ub = vb = -1.0;	/* for shutting up compiler */
	}
	a = (ua - ub) / (va - vb);
	b = (ub*va - ua*vb) / (va - vb);
	dx = a * vx + b - ux;
	tr1 = triangle_z_new(ux, vx, dx, ua - ux, va - vx);
	tr2 = triangle_z_new(ux, vx, dx, ub - ux, vb - vx);
    }
    assert(tr1 != NULL);

    ret1 = triangle_z_each(tr1, &p1);
    if (tr2 != NULL) {
	ret2 = triangle_z_each(tr2, &p2);
    } else {
	ret2 = 0;
    }
    ipoint_ary = ipoint_ary_new();
    for (; ret1; ret1 = triangle_z_each(tr1, &p1)) {
	if (p1 == NULL)
	    continue;
	ipoint_ary_push(ipoint_ary, *triangle_each2(self, p1));
    }
    if (tr2 != NULL) {
	for (; ret2; ret2 = triangle_z_each(tr2, &p2)) {
	    if (p1 == NULL)
		continue;
	    ipoint_ary_push(ipoint_ary, *triangle_each2(self, p2));
	}
    }
    triangle_z_free(tr1);
    triangle_z_free(tr2);

    self->each = each_new(ipoint_ary);

    return each_each(self->each, pp);
}

static void triangle_offset(Triangle *self)
{
    switch (self->axis) {
    case AXIS_X:
	self->du2 -= world->dy;
	self->du2 -= world->dz;
	self->du3 -= world->dy;
	self->du3 -= world->dz;
	break;
    case AXIS_Y:
	self->du2 -= world->dx;
	self->dv2 -= world->dz;
	self->du3 -= world->dx;
	self->dv3 -= world->dz;
	break;
    case AXIS_Z:
	self->du2 -= world->dx;
	self->dv2 -= world->dy;
	self->du3 -= world->dx;
	self->dv3 -= world->dy;
	break;
    default:
	bug("unknow axis %d", self->axis);
    }
}

/* Ellipse */

Ellipse *ellipse_new(double x, double y, double z, int axis, double ru, double rv)
{
    Ellipse *self;

    if (ru < 0.0 || rv < 0.0)
	warn_exit("length is negative for Ellipse");

    self = EALLOC(Ellipse);
    self->x = x;
    self->y = y;
    self->z = z;
    self->axis = axis;
    self->ru = ru;
    self->rv = rv;
    self->each = NULL;
    return self;
}

static void ellipse_free(Ellipse *self)
{
    if (self == NULL)
	return;
    each_free(self->each);
    FREE(self);
}

static void ellipse_ipoint_ary_add(Ellipse *self, iPoint_ary *ipoint_ary, int axis,
	IP_TYPE wc, IP_TYPE u1, IP_TYPE u2, IP_TYPE v)
{
    IP_TYPE u;
    iPoint ipoint;

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

static int ellipse_each(Ellipse *self, iPoint **pp)
{
    IP_TYPE uc, vc, wc, ru, rv;
    IP_TYPE ui, vi, ri, u1, v1;
    iPoint_ary *ipoint_ary;

    if (self->each != NULL && self->each->index >= 0)
	return each_each(self->each, pp);

    switch (self->axis) {
    case AXIS_X:
	uc = iround((self->y - world->y0) / world->dy);
	vc = iround((self->z - world->z0) / world->dz);
	wc = iround((self->x - world->x0) / world->dx);
	ru = iround(self->ru / world->dy);
	rv = iround(self->rv / world->dz);
	break;
    case AXIS_Y:
	uc = iround((self->x - world->x0) / world->dx);
	vc = iround((self->z - world->z0) / world->dz);
	wc = iround((self->y - world->y0) / world->dy);
	ru = iround(self->ru / world->dx);
	rv = iround(self->rv / world->dz);
	break;
    case AXIS_Z:
	uc = iround((self->x - world->x0) / world->dx);
	vc = iround((self->y - world->y0) / world->dy);
	wc = iround((self->z - world->z0) / world->dz);
	ru = iround(self->ru / world->dx);
	rv = iround(self->rv / world->dy);
	break;
    default:
	bug("unknown axis %d", self->axis);
	uc = vc = wc = ru = rv = -1;	/* for shutting up compiler */
    }

    if (ru == 0 && rv == 0)
	return 0;

    ipoint_ary = ipoint_ary_new();
    if (ru > rv) {
	ui = ri = ru;  vi = 0;
	while (ui >= vi) {
	    u1 = (IP_TYPE)((long)ui * rv / ru);
	    v1 = (IP_TYPE)((long)vi * rv / ru);
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
	    u1 = (IP_TYPE)((long)ui * ru / rv);
	    v1 = (IP_TYPE)((long)vi * ru / rv);
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
	self->ru -= world->dy;
	self->rv -= world->dz;
	break;
    case AXIS_Y:
	self->ru -= world->dz;
	self->rv -= world->dx;
	break;
    case AXIS_Z:
	self->ru -= world->dx;
	self->rv -= world->dy;
	break;
    default:
	bug("unknow axis %d", self->axis);
    }
    if (self->ru < 0.0 || self->rv < 0.0)
	warn_exit("length of ellipse becomes negative");
}

/* Ellipseperi */

Ellipseperi *ellipseperi_new(double x, double y, double z, int axis, double ru, double rv, double angle_st, double angle_en)
{
    Ellipseperi *self;

    if (ru < 0.0 || rv < 0.0)
	warn_exit("length is negative for Ellipseperi");

    self = EALLOC(Ellipseperi);
    self->x = x;
    self->y = y;
    self->z = z;
    self->axis = axis;
    self->ru = ru;
    self->rv = rv;
    self->angle_st = angle_st;
    self->angle_en = angle_en;
    self->each = NULL;
    return self;
}

static void ellipseperi_free(Ellipseperi *self)
{
    if (self == NULL)
	return;
    each_free(self->each);
    FREE(self);
}

static void ellipseperi_ipoint_ary_add(Ellipseperi *self, iPoint_ary *ipoint_ary, int axis,
	IP_TYPE wc, IP_TYPE u, IP_TYPE v)
{
    iPoint ipoint;

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

static double renormalize_angle(double angle, double angle0)
{
    if (angle < angle0)
        return renormalize_angle(angle + 360.0, angle0);
    else if (angle < angle0 + 360.0)
        return angle;
    else
        return renormalize_angle(angle - 360.0, angle0);
}

static void ellipseperi_ipoint_ary_add_sub(Ellipseperi *self, iPoint_ary *ipoint_ary, int axis,
	IP_TYPE wc, IP_TYPE uc, IP_TYPE vc, IP_TYPE du, IP_TYPE dv)
{
    double angle;

    angle = atan2(dv, du) / M_PI * 180.0;
    angle = renormalize_angle(angle, self->angle_st);
    if (angle >= self->angle_st && angle <= self->angle_en)
        ellipseperi_ipoint_ary_add(self, ipoint_ary, self->axis, wc, uc + du, vc + dv);
}

static int ellipseperi_each(Ellipseperi *self, iPoint **pp)
{
    IP_TYPE uc, vc, wc, ru, rv;
    IP_TYPE ui, vi, ri, u1, v1;
    iPoint_ary *ipoint_ary;

    if (self->each != NULL && self->each->index >= 0)
	return each_each(self->each, pp);

    switch (self->axis) {
    case AXIS_X:
	uc = iround((self->y - world->y0) / world->dy);
	vc = iround((self->z - world->z0) / world->dz);
	wc = iround((self->x - world->x0) / world->dx);
	ru = iround(self->ru / world->dy);
	rv = iround(self->rv / world->dz);
	break;
    case AXIS_Y:
	uc = iround((self->x - world->x0) / world->dx);
	vc = iround((self->z - world->z0) / world->dz);
	wc = iround((self->y - world->y0) / world->dy);
	ru = iround(self->ru / world->dx);
	rv = iround(self->rv / world->dz);
	break;
    case AXIS_Z:
	uc = iround((self->x - world->x0) / world->dx);
	vc = iround((self->y - world->y0) / world->dy);
	wc = iround((self->z - world->z0) / world->dz);
	ru = iround(self->ru / world->dx);
	rv = iround(self->rv / world->dy);
	break;
    default:
	bug("unknown axis %d", self->axis);
	uc = vc = wc = ru = rv = -1;	/* for shutting up compiler */
    }

    if (ru == 0 && rv == 0)
	return 0;

    ipoint_ary = ipoint_ary_new();
    if (ru > rv) {
        ui = ri = ru;  vi = 0;
        while (ui >= vi) {
            u1 = (IP_TYPE)((long)ui * rv / ru);
            v1 = (IP_TYPE)((long)vi * rv / ru);
            ellipseperi_ipoint_ary_add_sub(self, ipoint_ary, self->axis, wc, uc, vc, -ui, -v1);
            ellipseperi_ipoint_ary_add_sub(self, ipoint_ary, self->axis, wc, uc, vc,  ui, -v1);
            ellipseperi_ipoint_ary_add_sub(self, ipoint_ary, self->axis, wc, uc, vc, -ui,  v1);
            ellipseperi_ipoint_ary_add_sub(self, ipoint_ary, self->axis, wc, uc, vc,  ui,  v1);
            ellipseperi_ipoint_ary_add_sub(self, ipoint_ary, self->axis, wc, uc, vc, -vi, -u1);
            ellipseperi_ipoint_ary_add_sub(self, ipoint_ary, self->axis, wc, uc, vc,  vi, -u1);
            ellipseperi_ipoint_ary_add_sub(self, ipoint_ary, self->axis, wc, uc, vc, -vi,  u1);
            ellipseperi_ipoint_ary_add_sub(self, ipoint_ary, self->axis, wc, uc, vc,  vi,  u1);
            if ((ri -= (vi++ << 1) + 1) <= 0)
                ri += (ui-- - 1) << 1;
        }
    } else {
        ui = ri = rv;  vi = 0;
        while (ui >= vi) {
            u1 = (IP_TYPE)((long)ui * ru / rv);
            v1 = (IP_TYPE)((long)vi * ru / rv);
            ellipseperi_ipoint_ary_add_sub(self, ipoint_ary, self->axis, wc, uc, vc, -u1, -vi);
            ellipseperi_ipoint_ary_add_sub(self, ipoint_ary, self->axis, wc, uc, vc,  u1, -vi);
            ellipseperi_ipoint_ary_add_sub(self, ipoint_ary, self->axis, wc, uc, vc, -u1,  vi);
            ellipseperi_ipoint_ary_add_sub(self, ipoint_ary, self->axis, wc, uc, vc,  u1,  vi);
            ellipseperi_ipoint_ary_add_sub(self, ipoint_ary, self->axis, wc, uc, vc, -v1, -ui);
            ellipseperi_ipoint_ary_add_sub(self, ipoint_ary, self->axis, wc, uc, vc,  v1, -ui);
            ellipseperi_ipoint_ary_add_sub(self, ipoint_ary, self->axis, wc, uc, vc, -v1,  ui);
            ellipseperi_ipoint_ary_add_sub(self, ipoint_ary, self->axis, wc, uc, vc,  v1,  ui);
            if ((ri -= (vi++ << 1) - 1) < 0)
                ri += (ui-- - 1) << 1;
        }
    }
    self->each = each_new(ipoint_ary);
    return each_each(self->each, pp);
}

static void ellipseperi_offset(Ellipseperi *self)
{
    switch (self->axis) {
    case AXIS_X:
	self->ru -= world->dy;
	self->rv -= world->dz;
	break;
    case AXIS_Y:
	self->ru -= world->dz;
	self->rv -= world->dx;
	break;
    case AXIS_Z:
	self->ru -= world->dx;
	self->rv -= world->dy;
	break;
    default:
	bug("unknow axis %d", self->axis);
    }
    if (self->ru < 0.0 || self->rv < 0.0)
	warn_exit("length of ellipseperi becomes negative");
}

/* Circle */

Circle *circle_new(double x, double y, double z, int axis, double r)
{
    Circle *self;

    if (r < 0.0)
	warn_exit("length is negative for Circle");

    self = EALLOC(Circle);
    self->ellipse = ellipse_new(x, y, z, axis, r, r);
    return self;
}

static void circle_free(Circle *self)
{
    if (self == NULL)
	return;
    ellipse_free(self->ellipse);
    FREE(self);
}

static int circle_each(Circle *self, iPoint **pp)
{
    return ellipse_each(self->ellipse, pp);
}

static void circle_offset(Circle *self)
{
    ellipse_offset(self->ellipse);
}

/* Circleperi */

Circleperi *circleperi_new(double x, double y, double z, int axis, double r, double angle_st, double angle_en)
{
    Circleperi *self;

    if (r < 0.0)
	warn_exit("length is negative for Circleperi");

    self = EALLOC(Circleperi);
    self->ellipseperi = ellipseperi_new(x, y, z, axis, r, r, angle_st, angle_en);
    return self;
}

static void circleperi_free(Circleperi *self)
{
    if (self == NULL)
	return;
    ellipseperi_free(self->ellipseperi);
    FREE(self);
}

static int circleperi_each(Circleperi *self, iPoint **pp)
{
    return ellipseperi_each(self->ellipseperi, pp);
}

static void circleperi_offset(Circleperi *self)
{
    ellipseperi_offset(self->ellipseperi);
}

/* Polygon */

Polygon *polygon_new(double x1, double y1, double z1,
	int axis, Vector2d_ary *dudv_ary)
{
    Polygon *self;
    Vector2d v0, v;
    int index;

    v0.x = v0.y = -1.0;	/* for shutting up compiler */

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
    self->rotate_p_func = vector2d_counter_clock_p;
    self->each = NULL;
    return self;
}

Polygon *polygon_new2(double x1, double y1, double z1,
	int axis, Vector2d_ary *uv_ary)
{
    Vector2d_ary *dudv_ary;
    double u, v;
    int i;
    Vector2d vec2d;

    dudv_ary = vector2d_ary_new();
    switch (axis) {
    case AXIS_X:
	u = y1;
	v = z1;
	break;
    case AXIS_Y:
	u = x1;
	v = z1;
	break;
    case AXIS_Z:
	u = x1;
	v = y1;
	break;
    default:
	bug("unknow axis %d", axis);
	u = -1.0;	/* for shutting up compiler */
	v = -1.0;	/* for shutting up compiler */
    }
    for (i = 0; i < uv_ary->size; ++i) {
	vec2d.x = uv_ary->ptr[i].x - u;
	vec2d.y = uv_ary->ptr[i].y - v;
	vector2d_ary_push(dudv_ary, vec2d);
    }
    return polygon_new(x1, y1, z1, axis, dudv_ary);
}

static void polygon_free(Polygon *self)
{
    if (self == NULL)
	return;
    vector2d_ary_free(self->vector2d_ary);
    each_free(self->each);
    FREE(self);
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
    int n_fail;

    if (self->each != NULL && self->each->index >= 0)
	return each_each(self->each, pp);

    ary = vector2d_ary_new();
    for (index = 0; index < self->vector2d_ary->size; ++index) {
	vector2d_ary_push(ary, self->vector2d_ary->ptr[index]);
    }

    ipoint_ary = ipoint_ary_new();
    n_fail = 0;
    while (ary->size >= 3) {
	va = ary->ptr[0];
	vb = ary->ptr[1];
	vc = ary->ptr[2];
	ok = 1;
	if (!self->rotate_p_func(va, vb, vc)) {
	    ok = 0;
	} else {
	    for (index = 3; index < ary->size; ++index) {
		if (vector2d_inner_triangle_p(ary->ptr[index], va, vb, vc)) {
		    ok = 0;
		    break;
		}
	    }
	    if (ok) {
		switch (self->axis) {
		case AXIS_X:
		    x1 = self->w;
		    y1 = va.x;
		    z1 = va.y;
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
		    x1 = y1 = z1 = -1.0;	/* for shutting up compiler */
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
	if (ok) {
	    n_fail = 0;
	} else {
	    ++n_fail;
	    if (n_fail >= ary->size) {
		if (self->rotate_p_func == vector2d_counter_clock_p) {
		    self->rotate_p_func = vector2d_clock_p;
		    return polygon_each(self, pp);
		} else {
		    warn_exit("invalid polygon");
		}
	    }
	    vector2d_ary_rotate(ary, 1);
	}
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
	    self->vector2d_ary->ptr[index].x -= world->dy;
	    self->vector2d_ary->ptr[index].y -= world->dz;
	}
	break;
    case AXIS_Y:
	for (index = 0; index < self->vector2d_ary->size; ++index) {
	    self->vector2d_ary->ptr[index].x -= world->dx;
	    self->vector2d_ary->ptr[index].y -= world->dz;
	}
	break;
    case AXIS_Z:
	for (index = 0; index < self->vector2d_ary->size; ++index) {
	    self->vector2d_ary->ptr[index].x -= world->dx;
	    self->vector2d_ary->ptr[index].y -= world->dy;
	}
	break;
    default:
	bug("unknow axis %d", self->axis);
    }
}

/* Line */

Line *line_new(double x1, double y1, double z1,
	int axis, Vector2d_ary *dudv_ary)
{
    Line *self;
    Vector2d v0, v;
    int index;

    v0.x = v0.y = -1.0;	/* for shutting up compiler */

    self = EALLOC(Line);
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

Line *line_new2(double x1, double y1, double z1,
	int axis, Vector2d_ary *uv_ary)
{
    Vector2d_ary *dudv_ary;
    double u, v;
    int i;
    Vector2d vec2d;

    dudv_ary = vector2d_ary_new();
    switch (axis) {
    case AXIS_X:
	u = y1;
	v = z1;
	break;
    case AXIS_Y:
	u = x1;
	v = z1;
	break;
    case AXIS_Z:
	u = x1;
	v = y1;
	break;
    default:
	bug("unknow axis %d", axis);
	u = -1.0;	/* for shutting up compiler */
	v = -1.0;	/* for shutting up compiler */
    }
    for (i = 0; i < uv_ary->size; ++i) {
	vec2d.x = uv_ary->ptr[i].x - u;
	vec2d.y = uv_ary->ptr[i].y - v;
	vector2d_ary_push(dudv_ary, vec2d);
    }
    return line_new(x1, y1, z1, axis, dudv_ary);
}

static void line_free(Line *self)
{
    if (self == NULL)
	return;
    vector2d_ary_free(self->vector2d_ary);
    each_free(self->each);
    FREE(self);
}

static void line_set_ipoint_ary(iPoint_ary *ipoint_ary,
	int x1, int y1, int x2, int y2)
{
    int dx, dy, s, step;

    dx = abs(x2 - x1);  dy = abs(y2 - y1);
    if (dx > dy) {
        if (x1 > x2) {
            step = (y1 > y2) ? 1 : -1;
            s = x1;  x1 = x2;  x2 = s;  y1 = y2;
        } else step = (y1 < y2) ? 1: -1;
        ipoint_ary_push(ipoint_ary, get_ipoint(x1, y1, 0));
        s = dx >> 1;
        while (++x1 <= x2) {
            if ((s -= dy) < 0) {
                s += dx;  y1 += step;
            };
	    ipoint_ary_push(ipoint_ary, get_ipoint(x1, y1, 0));
        }
    } else {
        if (y1 > y2) {
            step = (x1 > x2) ? 1 : -1;
            s = y1;  y1 = y2;  y2 = s;  x1 = x2;
        } else step = (x1 < x2) ? 1 : -1;
	ipoint_ary_push(ipoint_ary, get_ipoint(x1, y1, 0));
        s = dy >> 1;
        while (++y1 <= y2) {
            if ((s -= dx) < 0) {
                s += dy;  x1 += step;
            }
	    ipoint_ary_push(ipoint_ary, get_ipoint(x1, y1, 0));
        }
    }
}

static int line_each(Line *self, iPoint **pp)
{
    iPoint_ary *ipoint_ary;
    double x1, y1, x2, y2;
    int ix1, iy1, ix2, iy2;
    int index;

    if (self->each != NULL && self->each->index >= 0)
	return each_each(self->each, pp);

    ipoint_ary = ipoint_ary_new();
    for (index = 0; index < self->vector2d_ary->size - 1; ++index) {
	x1 = self->vector2d_ary->ptr[index].x;
	y1 = self->vector2d_ary->ptr[index].y;
	x2 = self->vector2d_ary->ptr[index + 1].x;
	y2 = self->vector2d_ary->ptr[index + 1].y;
	switch (self->axis) {
	case AXIS_X:
	    ix1 = iround((x1 - world->y0) / world->dy);
	    iy1 = iround((y1 - world->z0) / world->dz);
	    ix2 = iround((x2 - world->y0) / world->dy);
	    iy2 = iround((y2 - world->z0) / world->dz);
	    break;
	case AXIS_Y:
	    ix1 = iround((x1 - world->x0) / world->dx);
	    iy1 = iround((y1 - world->z0) / world->dz);
	    ix2 = iround((x2 - world->x0) / world->dx);
	    iy2 = iround((y2 - world->z0) / world->dz);
	    break;
	case AXIS_Z:
	    ix1 = iround((x1 - world->x0) / world->dx);
	    iy1 = iround((y1 - world->y0) / world->dy);
	    ix2 = iround((x2 - world->x0) / world->dx);
	    iy2 = iround((y2 - world->y0) / world->dy);
	    break;
	default:
	    bug("unknown axis %d", self->axis);
	    ix1 = iy1 = ix2 = iy2 = -1;	/* for shutting compiler */
	}
	line_set_ipoint_ary(ipoint_ary, ix1, iy1, ix2, iy2);
    }

    for (index = 0; index < ipoint_ary->size; ++index) {
	switch (self->axis) {
	case AXIS_X:
	    ipoint_ary->ptr[index].k = ipoint_ary->ptr[index].j;
	    ipoint_ary->ptr[index].j = ipoint_ary->ptr[index].i;
	    ipoint_ary->ptr[index].i = iround((self->w - world->x0) / world->dx);
	    break;
	case AXIS_Y:
	    ipoint_ary->ptr[index].k = ipoint_ary->ptr[index].j;
	    ipoint_ary->ptr[index].j = iround((self->w - world->y0) / world->dy);
	    break;
	case AXIS_Z:
	    ipoint_ary->ptr[index].k = iround((self->w - world->z0) / world->dz);
	    break;
	default:
	    bug("unknown axis %d", self->axis);
	}
    }

    self->each = each_new(ipoint_ary);
    return each_each(self->each, pp);
}

/* Box */

Box *box_new(double x, double y, double z, double xlen, double ylen, double zlen)
{
    Box *self;
    Obj *obj;

    if (xlen < 0.0 || ylen < 0.0 || zlen < 0.0)
	warn_exit("length is negative for Box");

    self = EALLOC(Box);
    obj = obj_new(OBJ_RECT);
    obj->uobj.rect = rect_new(x, y, z, AXIS_Z, xlen, ylen);
    self->sweep = sweep_new(AXIS_Z, zlen, obj);
    return self;
}

static void box_free(Box *self)
{
    if (self == NULL)
	return;
    sweep_free(self->sweep);
    FREE(self);
}

static int box_each(Box *self, iPoint **pp)
{
    return sweep_each(self->sweep, pp);
}

static void box_offset(Box *self)
{
    sweep_offset(self->sweep);
}

/* ObjAry */

ObjAry *objary_new(AryObj *aryobj)
{
    ObjAry *self;

    self = EALLOC(ObjAry);
    self->aryobj = aryobj;
    self->each_obj_index = 0;
    return self;
}

static void objary_free(ObjAry *self)
{
    if (self == NULL)
	return;
    aryobj_free(self->aryobj);
    FREE(self);
}

static int objary_each(ObjAry *self, iPoint **pp)
{
    int each_ret;

    while (self->each_obj_index < self->aryobj->size) {
	each_ret = obj_each(self->aryobj->ptr[self->each_obj_index], pp);
	if (each_ret)
	    return each_ret;
	++(self->each_obj_index);
    }
    self->each_obj_index = 0;
    return 0;
}

static void objary_offset(ObjAry *self)
{
    int index;

    for (index = 0; index < self->aryobj->size; ++index)
	obj_offset(self->aryobj->ptr[index]);
}

/* Obj */

Obj *obj_new(int objtype)
{
    Obj *self;

    self = EALLOC(Obj);
    self->objtype = objtype;
    return self;
}

void obj_free(Obj *self)
{
    if (self == NULL)
	return;
    switch (self->objtype) {
    case OBJ_RECT:
	rect_free(self->uobj.rect);
	break;
    case OBJ_TRIANGLE:
	triangle_free(self->uobj.triangle);
	break;
    case OBJ_ELLIPSE:
	ellipse_free(self->uobj.ellipse);
	break;
    case OBJ_ELLIPSEPERI:
	ellipseperi_free(self->uobj.ellipseperi);
	break;
    case OBJ_CIRCLE:
	circle_free(self->uobj.circle);
	break;
    case OBJ_CIRCLEPERI:
	circleperi_free(self->uobj.circleperi);
	break;
    case OBJ_POLYGON:
	polygon_free(self->uobj.polygon);
	break;
    case OBJ_LINE:
	line_free(self->uobj.line);
	break;
    case OBJ_BOX:
	box_free(self->uobj.box);
	break;
    case OBJ_SWEEP:
	sweep_free(self->uobj.sweep);
	break;
    case OBJ_OBJARY:
	objary_free(self->uobj.objary);
	break;
    default:
	bug("unknown obj %d", self->objtype);
    }
    return;
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
    case OBJ_ELLIPSEPERI:
	return ellipseperi_each(self->uobj.ellipseperi, pp);
	break;
    case OBJ_CIRCLE:
	return circle_each(self->uobj.circle, pp);
	break;
    case OBJ_CIRCLEPERI:
	return circleperi_each(self->uobj.circleperi, pp);
	break;
    case OBJ_POLYGON:
	return polygon_each(self->uobj.polygon, pp);
	break;
    case OBJ_LINE:
	return line_each(self->uobj.line, pp);
	break;
    case OBJ_BOX:
	return box_each(self->uobj.box, pp);
	break;
    case OBJ_SWEEP:
	return sweep_each(self->uobj.sweep, pp);
	break;
    case OBJ_OBJARY:
	return objary_each(self->uobj.objary, pp);
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
    case OBJ_ELLIPSEPERI:
	ellipseperi_offset(self->uobj.ellipseperi);
	break;
    case OBJ_CIRCLE:
	circle_offset(self->uobj.circle);
	break;
    case OBJ_CIRCLEPERI:
	circleperi_offset(self->uobj.circleperi);
	break;
    case OBJ_POLYGON:
	polygon_offset(self->uobj.polygon);
	break;
    case OBJ_LINE:
	warn_exit("line_offset is not implemented");
	break;
    case OBJ_BOX:
	box_offset(self->uobj.box);
	break;
    case OBJ_SWEEP:
	sweep_offset(self->uobj.sweep);
	break;
    case OBJ_OBJARY:
	objary_offset(self->uobj.objary);
	break;
    default:
	bug("unknown obj %d", self->objtype);
    }
}

static int obj_dim(Obj *self)
{
    int dim, index;
    AryObj *aryobj;

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
    case OBJ_ELLIPSEPERI:
	return 2;
	break;
    case OBJ_CIRCLE:
	return 2;
	break;
    case OBJ_CIRCLEPERI:
	return 2;
	break;
    case OBJ_POLYGON:
	return 2;
	break;
    case OBJ_LINE:
	return 2;
	break;
    case OBJ_BOX:
	return 3;
	break;
    case OBJ_SWEEP:
	return 3;
	break;
    case OBJ_OBJARY:
	dim = -1;
	aryobj = self->uobj.objary->aryobj;
	for (index = 0; index < aryobj->size; ++index) {
	    if (dim < 0)
		dim = obj_dim(aryobj->ptr[index]);
	    else if (dim != obj_dim(aryobj->ptr[index]))
		return -1;
	}
	return dim;
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
    self->alloc_size = 0;
    self->size = 0;
    return self;
}

void aryobj_free(AryObj *self)
{
    int index;

    if (self == NULL)
	return;
    for (index = 0; index < self->size; ++index)
	obj_free(self->ptr[index]);
    FREE(self->ptr);
    FREE(self);
}

void aryobj_push(AryObj *self, Obj *obj)
{
    ++(self->size);
    if (self->size > self->alloc_size) {
	self->alloc_size = alloc_size(self->size);
	self->ptr = (Obj **) erealloc(self->ptr, sizeof(Obj *) * self->alloc_size);
    }
    self->ptr[self->size - 1] = obj;
}
