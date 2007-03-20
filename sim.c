#include <stdio.h>
#include <assert.h>
#include "sim.h"
#include "mem.h"
#include "tc.h"
#include "solvele.h"

#define NELEMS(ary)	(sizeof(ary) / sizeof((ary)[0]))

enum {
    DIR_LEFT = 0,
    DIR_RIGHT,
    DIR_FRONT,
    DIR_BACK,
    DIR_BELOW,
    DIR_ABOVE
};

static int dir_array[NDIRS] = { DIR_LEFT, DIR_RIGHT, DIR_FRONT, DIR_BACK, DIR_BELOW, DIR_ABOVE };
static int dir_x[] = { DIR_LEFT, DIR_RIGHT };
static int dir_y[] = { DIR_FRONT, DIR_BACK };
static int dir_z[] = { DIR_BELOW, DIR_ABOVE };

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

static Point point_add(Point ponint1, Point ponint2)
{
    return get_point(
	    ponint1.i + ponint2.i,
	    ponint1.j + ponint2.j,
	    ponint1.k + ponint2.k
	    );
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
    int dir;

    self = EALLOC(Coefs);
    for (dir = 0; dir < NELEMS(dir_array); ++dir) {
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

    if (self->each != NULL && self->each->index >= 0)
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

    if (self->each != NULL && self->each->index >= 0)
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
	    point_ary[n++] = *p;
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
	    0.1, 0.1, 0.1);
    self->active_obj_ary = aryobj_new();
    {
	obj = obj_new(OBJ_BOX, OBJVAL_I);
	obj->uobj.box = box_new(self->world, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
	obj->uval.i = 1;
	aryobj_push(self->active_obj_ary, obj);
    }

    self->fix_obj_ary = aryobj_new();
    {
	obj = obj_new(OBJ_RECT, OBJVAL_D);
	obj->uobj.rect = rect_new(self->world, 0.0, 0.0, 0.0, AXIS_X, 0.5, 0.5);
	obj->uval.d = 0.0;
	aryobj_push(self->fix_obj_ary, obj);
    }
    {
	obj = obj_new(OBJ_RECT, OBJVAL_D);
	obj->uobj.rect = rect_new(self->world, 1.0, 0.5, 0.5, AXIS_X, 0.5, 0.5);
	obj->uval.d = 1.0;
	aryobj_push(self->fix_obj_ary, obj);
    }

    self->heatflow_obj_ary = aryobj_new();

    self->lambda_obj_ary = aryobj_new();

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
    Point *p;
    AryObj *active_obj_ary;
    Obj *obj;
    int i;

    ALLOCATE_3D2(self->active_p_ary, int, self->ni, self->nj, self->nk, -1);
    active_obj_ary = self->config->active_obj_ary;
    for (i = 0; i < active_obj_ary->size; ++i) {
	obj = active_obj_ary->ptr[i];
	for (p = obj_each_begin(obj); p != NULL; p = obj_each(obj)) {
	    (self->active_p_ary)[p->i][p->j][p->k] = obj->uval.i;
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

static void sim_set_matrix_const(Sim *self, int i, int j, int k, Point point)
{
    int dir;

    static int XXX = 1;
    if (XXX) {
	fprintf(stderr, "XXX sim_set_matrix_const\n");
	XXX = 0;
    }

    for (dir = 0; dir < NELEMS(dir_array); ++dir) {
	/*
	 * next if @heatflow_ary[i, j, k][dir].nil?
	 * @coefs[*point.to_a].const += @heatflow_ary[i, j, k][dir]
	 */
    }
}

static void sim_set_matrix_coef0(Sim *self,
	int i, int j, int k, Point point, double dx, double dy, double dz)
{
    int dir;
    int ix, iy, iz;
    Point pp;
    double l;

    static int XXX = 1;
    if (XXX) {
	fprintf(stderr, "XXX sim_set_matrix_coef0\n");
	XXX = 0;
    }

    for (dir = 0; dir < NELEMS(dir_array); ++dir) {
	/* next unless @heatflow_ary[i, j, k][dir].nil? */
	if (!sim_active_p(self, point_add(point, self->dir_to_point[dir])))
	    continue;
	switch (dir) {
	case DIR_LEFT: case DIR_RIGHT:
	    for (iy = 0; iy < NELEMS(dir_y); ++iy) {
		for (iz = 0; iz < NELEMS(dir_z); ++iz) {
		    pp = point_add(point, self->dir_to_point[dir_y[iy]]);
		    pp = point_add(pp, self->dir_to_point[dir_z[iz]]);
		    if (sim_active_p(self, pp)) {
			/*
			   pnt_l = point.offset(dir, diry, dirz);
			   l = @lambda_ary[*pnt_l.to_a];
			   */
			l = 1.0;
			self->coefs[point.i][point.j][point.k]->coef0 += -l/dx*(dy/2)*(dz/2);
		    }
		}
	    }
	    break;
	case DIR_FRONT: case DIR_BACK:
	    for (iz = 0; iz < NELEMS(dir_z); ++iz) {
		for (ix = 0; ix < NELEMS(dir_x); ++ix) {
		    pp = point_add(point, self->dir_to_point[dir_z[iz]]);
		    pp = point_add(pp, self->dir_to_point[dir_x[ix]]);
		    if (sim_active_p(self, pp)) {
			/*
			   pnt_l = point.offset(dirx, dir, dirz);
			   l = @lambda_ary[*pnt_l.to_a];
			   */
			l = 1.0;
			self->coefs[point.i][point.j][point.k]->coef0 += -l/dy*(dz/2)*(dx/2);
		    }
		}
	    }
	    break;
	case DIR_BELOW: case DIR_ABOVE:
	    for (ix = 0; ix < NELEMS(dir_x); ++ix) {
		for (iy = 0; iy < NELEMS(dir_y); ++iy) {
		    pp = point_add(point, self->dir_to_point[dir_x[ix]]);
		    pp = point_add(pp, self->dir_to_point[dir_y[iy]]);
		    if (sim_active_p(self, pp)) {
			/*
			   pnt_l = point.offset(dirx, diry, dir);
			   l = @lambda_ary[*pnt_l.to_a];
			   */
			l = 1.0;
			self->coefs[point.i][point.j][point.k]->coef0 += -l/dz*(dx/2)*(dy/2);
		    }
		}
	    }
	    break;
	default:
	    warn_exit("bug in sim_set_matrix_coef0");
	}
    }
}

static void sim_set_matrix_coef(Sim *self,
	int i, int j, int k, Point point, double dx, double dy, double dz)
{
    int dir;
    double value;
    int ix, iy, iz;
    Point pp;
    double l;

    static int XXX = 1;
    if (XXX) {
	fprintf(stderr, "XXX sim_set_matrix_coef\n");
	XXX = 0;
    }

    for (dir = 0; dir < NELEMS(dir_array); ++dir) {
	if (!sim_active_p(self, point_add(point, self->dir_to_point[dir]))) {
	    self->coefs[point.i][point.j][point.k]->coef[dir]->index = -1;
	    continue;
	}
	/*
	unless @heatflow_ary[i, j, k][dir].nil?
	  raise "heatflow comes from active cell at (#{i}, #{j}, #{k})"
	end
	*/
	self->coefs[point.i][point.j][point.k]->coef[dir]->index =
	    world_to_index(self->world, point_add(point, self->dir_to_point[dir]));
	value = 0.0;
	switch (dir) {
	case DIR_LEFT: case DIR_RIGHT:
	    for (iy = 0; iy < NELEMS(dir_y); ++iy) {
		for (iz = 0; iz < NELEMS(dir_z); ++iz) {
		    pp = point_add(point, self->dir_to_point[dir_y[iy]]);
		    pp = point_add(pp, self->dir_to_point[dir_z[iz]]);
		    if (sim_active_p(self, pp)) {
			/*
			   pnt_l = point.offset(dir, diry, dirz);
			   l = @lambda_ary[*pnt_l.to_a];
			   */
			l = 1.0;
			value += l/dx*(dy/2)*(dz/2);
		    }
		}
	    }
	    break;
	case DIR_FRONT: case DIR_BACK:
	    for (iz = 0; iz < NELEMS(dir_z); ++iz) {
		for (ix = 0; ix < NELEMS(dir_x); ++ix) {
		    pp = point_add(point, self->dir_to_point[dir_z[iz]]);
		    pp = point_add(pp, self->dir_to_point[dir_x[ix]]);
		    if (sim_active_p(self, pp)) {
			/*
			   pnt_l = point.offset(dirx, dir, dirz);
			   l = @lambda_ary[*pnt_l.to_a];
			   */
			l = 1.0;
			value += l/dy*(dz/2)*(dx/2);
		    }
		}
	    }
	    break;
	case DIR_BELOW: case DIR_ABOVE:
	    for (ix = 0; ix < NELEMS(dir_x); ++ix) {
		for (iy = 0; iy < NELEMS(dir_y); ++iy) {
		    pp = point_add(point, self->dir_to_point[dir_x[ix]]);
		    pp = point_add(pp, self->dir_to_point[dir_y[iy]]);
		    if (sim_active_p(self, pp)) {
			/*
			   pnt_l = point.offset(dirx, diry, dir);
			   l = @lambda_ary[*pnt_l.to_a];
			   */
			l = 1.0;
			value += l/dz*(dx/2)*(dy/2);
		    }
		}
	    }
	    break;
	default:
	    warn_exit("bug in sim_set_matrix_coef");
	}
	self->coefs[point.i][point.j][point.k]->coef[dir]->value = value;
    }
}

static void sim_set_matrix(Sim *self)
{
    Point *p;
    double dx, dy, dz;
    Point point;

    static int XXX = 1;
    if (XXX) {
	fprintf(stderr, "XXX sim_set_matrix\n");
	XXX = 0;
    }

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

	/*
	  point = @heatflow_point_ary[i, j, k]
	  if point.nil?
	    point = Point.new(i, j, k)
	  elsif point != Point.new(i, j, k)
	    next
	  end
	  */

	point = get_point(p->i, p->j, p->k);

	sim_set_matrix_const(self, p->i, p->j, p->k, point);
	sim_set_matrix_coef0(self, p->i, p->j, p->k, point, dx, dy, dz);
	sim_set_matrix_coef(self, p->i, p->j, p->k, point, dx, dy, dz);
    }
}

Sim *sim_new(void)
{
    Sim *self;

    self = EALLOC(Sim);
    self->config = config_new();
    self->dir_to_point[DIR_LEFT]  = get_point(-1,  0,  0);
    self->dir_to_point[DIR_RIGHT] = get_point( 1,  0,  0);
    self->dir_to_point[DIR_FRONT] = get_point( 0, -1,  0);
    self->dir_to_point[DIR_BACK]  = get_point( 0,  1,  0);
    self->dir_to_point[DIR_BELOW] = get_point( 0,  0, -1);
    self->dir_to_point[DIR_ABOVE] = get_point( 0,  0,  1);

    sim_set_region(self);
    sim_set_matrix(self);

    return self;
}

Array3Dd sim_calc(Sim *self)
{
    int nindex;
    Solvele *solver;
    Point *p;
    int index, index2;
    int dir;
    double c;
    double *sol;
    Array3Dd ary;

    nindex = self->ni * self->nj * self->nk;
    solver = solvele_new(nindex);
    for (p = world_each_begin(self->world); p != NULL; p = world_each(self->world)) {
	index = world_to_index(self->world, *p);
	if (self->fix_ary[p->i][p->j][p->k] >= 0.0 || !self->active_p_ary[p->i][p->j][p->k]) {
	    solvele_set_matrix(solver, index, index, 1.0);
	    solvele_set_vector(solver, index, self->u[p->i][p->j][p->k]);
	} else if (0) {
	    /*
	      elsif not @heatflow_point_ary[i, j, k].nil? and (point != @heatflow_point_ary[i, j, k])
		solver.set_matrix(index, index, 1.0);
		solver.set_matrix(index, @world.to_index(@heatflow_point_ary[i, j, k]), -1.0);
		solver.set_vector(index, 0.0);
		*/
	} else {
	    for (dir = 0; dir < NELEMS(dir_array); ++dir) {
		index2 = self->coefs[p->i][p->j][p->k]->coef[dir]->index;
		if (index2 < 0)
		    continue;
		c = self->coefs[p->i][p->j][p->k]->coef[dir]->value;
		solvele_set_matrix(solver, index, index2, c);
	    }
	    solvele_set_matrix(solver, index, index, self->coefs[p->i][p->j][p->k]->coef0);
	    solvele_set_vector(solver, index, -self->coefs[p->i][p->j][p->k]->cnst);
	}
    }

    sol = solvele_solve(solver, 0);

    ALLOCATE_3D(ary, double, self->ni, self->nj, self->nk);
    index = 0;
    for (p = world_each_begin(self->world); p != NULL; p = world_each(self->world)) {
	ary[p->i][p->j][p->k] = sol[index++];
    }

    return ary;
}
