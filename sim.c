#include <stdio.h>
#include <assert.h>
#include "sim.h"
#include "mem.h"
#include "tc.h"

static int iround(double x)
{
    return (int) (x + 0.5);
}

/* Point */

static Point get_point(int i, int j, int k)
{
    Point point;

    point.i = i;
    point.j = j;
    point.k = k;
    return point;
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
    if (self->index < self->size)
	return &self->point_ary[self->index++];
    else
	return NULL;
}

/* World */

static World *world_new(double x, double y, double z,
	double xlen, double ylen, double zlen,
	double dx, double dy, double dz)
{
    World *self;

    self = EALLOC(World);
    self->x = x;
    self->y = y;
    self->z = z;
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

    if (self->each != NULL)
	bug("world_each_begin");
    size = self->ni * self->nj * self->nk;
    point_ary = EALLOCN(Point, size);
    n = 0;
    for (k = 0; k < self->nk; ++k) {
	for (j = 0; j < self->nj; ++j) {
	    for (i = 0; i < self->ni; ++i) {
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

static Rect *rect_new(World *world, double x, double y, double z, int axis, double len1, double len2)
{
    Rect *self;

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

    if (self->each != NULL)
	bug("rect_each_begin");

    xi = iround((self->x - self->world->x) / self->world->dx);
    yi = iround((self->y - self->world->y) / self->world->dy);
    zi = iround((self->z - self->world->z) / self->world->dz);
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

/* Box */

static Box *box_new(World *world, double x, double y, double z, double xlen, double ylen, double zlen)
{
    Box *self;

    self = EALLOC(Box);
    self->world = world;
    self->rect = rect_new(world, x, y, z, AXIS_Z, xlen, ylen);
    self->sweeplen = zlen;
    self->each = NULL;
    return self;
}

static Point *box_each_begin(Box *self)
{
    Point *pnt;
    int leni;
    int size;
    Point *point_ary;
    int n;
    int k;

    if (self->each != NULL)
	bug("box_each_begin");
    pnt = rect_each_begin(self->rect);
    leni = iround(self->sweeplen / self->world->dz) + 1;
    size = self->rect->each->size * leni;
    point_ary = EALLOCN(Point, size);
    n = 0;
    do {
	for (k = pnt->k; k < pnt->k + leni; ++k) {
	    point_ary[n++] = *pnt;
	}
	pnt = rect_each(self->rect);
    } while (pnt != NULL);
    assert(n == size);
    self->each = each_new(size, point_ary);
    return each_each(self->each);
}

static Point *box_each(Box *self)
{
    return each_each(self->each);
}

/* Obj */

static Obj *obj_new(int objtype, int valtype)
{
    Obj *self;

    self = EALLOC(Obj);
    self->objtype = objtype;
    self->valtype = valtype;
    return self;
}

static Point *obj_each_begin(Obj *self)
{
    Point *pnt;

    switch (self->objtype) {
    case OBJ_RECT:
	pnt = rect_each_begin(self->uobj.rect);
	break;
    case OBJ_BOX:
	pnt = box_each_begin(self->uobj.box);
	break;
    default:
	bug("unknown obj %d", self->objtype);
    }
    return pnt;
}

static Point *obj_each(Obj *self)
{
    Point *pnt;

    switch (self->objtype) {
    case OBJ_RECT:
	pnt = rect_each(self->uobj.rect);
	break;
    case OBJ_BOX:
	pnt = box_each(self->uobj.box);
	break;
    default:
	bug("unknown obj %d", self->objtype);
    }
    return pnt;
}

/* AryObj */

static AryObj *aryobj_new(void)
{
    AryObj *self;

    self = EALLOC(AryObj);
    self->ptr = NULL;
    self->size = 0;
    return self;
}

static void aryobj_push(AryObj *self, Obj *obj)
{
    ++self->size;
    self->ptr = erealloc(self->ptr, sizeof(Obj) * self->size);
    self->ptr[self->size - 1] = obj;
}

/* Config */

static Config *config_new(void)
{
    Config *self;
    Obj *obj;

    self = EALLOC(Config);
    self->world = world_new(0.0, 0.0, 0.0,
	    1.0, 1.0, 1.0,
	    0.5, 0.5, 0.5);
    self->active_obj_ary = aryobj_new();
    obj = obj_new(OBJ_BOX, OBJVAL_I);
    obj->uobj.box = box_new(self->world, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    obj->uval.i = 1;
    aryobj_push(self->active_obj_ary, obj);

    self->fix_obj_ary = aryobj_new();
    obj = obj_new(OBJ_RECT, OBJVAL_D);
    obj->uobj.rect = rect_new(self->world, 0.0, 0.0, 0.0, AXIS_X, 1.0, 1.0);
    obj->uval.d = 0.0;
    aryobj_push(self->fix_obj_ary, obj);
    obj = obj_new(OBJ_RECT, OBJVAL_D);
    obj->uobj.rect = rect_new(self->world, 1.0, 0.0, 0.0, AXIS_X, 1.0, 1.0);
    obj->uval.d = 1.0;
    aryobj_push(self->fix_obj_ary, obj);

    self->heatflow_obj_ary = aryobj_new();

    self->lambda_obj_ary = aryobj_new();

    return self;
}

/* Sim */

static void sim_set_region_active(Sim *self)
{
    Point *pnt;
    AryObj *active_obj_ary;
    Obj *obj;
    int i;

    ALLOCATE_3D(self->active_p_ary, int, self->ni, self->nj, self->nk);
    for (pnt = world_each_begin(self->world); pnt != NULL; pnt = world_each(self->world)) {
	(self->active_p_ary)[pnt->i][pnt->j][pnt->k] = -1;
	printf("%d, %d, %d\n", pnt->i, pnt->j, pnt->k);
    }
    active_obj_ary = self->config->active_obj_ary;
    for (i = 0; i < active_obj_ary->size; ++i) {
	obj = active_obj_ary->ptr[i];
	for (pnt = obj_each_begin(obj); pnt != NULL; pnt = obj_each(obj)) {
	    (self->active_p_ary)[pnt->i][pnt->j][pnt->k] = obj->uval.i;
	}
    }
    /*
    # remove active flag from the point having active alone plane in NDIR directions
    begin
      continue_p = false
      @world.each do |i, j, k|
        next unless @active_p_ary[i, j, k]
        pnt = Point.new(i, j, k)
        [
          [:LEFT, :RIGHT],
          [:FRONT, :BACK],
          [:BELOW, :ABOVE]
        ].each do |dirs|
          if not active_p(pnt + @dir_to_point[dirs[0]]) and
            not active_p(pnt + @dir_to_point[dirs[1]])
            @active_p_ary[i, j, k] = false
            continue_p = true
            STDERR.puts "removed (#{i}, #{j}, #{k})"
          end
        end
      end
    end while continue_p
    */
    fprintf(stderr, "XXX sim_set_region_active\n");
}

static void sim_set_region_fix(Sim *self)
{
    Point *pnt;
    AryObj *obj_ary;
    Obj *obj;
    int i;

    ALLOCATE_3D(self->fix_ary, double, self->ni, self->nj, self->nk);
    for (pnt = world_each_begin(self->world); pnt != NULL; pnt = world_each(self->world)) {
	self->fix_ary[pnt->i][pnt->j][pnt->k] = -1.0;
    }
    obj_ary = self->config->fix_obj_ary;
    for (i = 0; i < obj_ary->size; ++i) {
	obj = obj_ary->ptr[i];
	for (pnt = obj_each_begin(obj); pnt != NULL; pnt = obj_each(obj)) {
	    self->fix_ary[pnt->i][pnt->j][pnt->k] = obj->uval.d;
	}
    }
}

static void sim_set_region_heatflow(Sim *self)
{
    fprintf(stderr, "XXX sim_set_region_heatflow\n");
}

static void sim_set_region_lambda(Sim *self)
{
    fprintf(stderr, "XXX sim_set_region_lambda\n");
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

Sim *sim_new(void)
{
    Sim *self;

    self = EALLOC(Sim);
    self->config = config_new();
    sim_set_region(self);

    return self;
}

Array3Dd sim_calc(Sim *sim)
{
    Array3Dd ary;

    ALLOCATE_3D(ary, double, 1, 1, 1);
    return ary;
}
