#ifndef _SIM_H_
#define _SIM_H_

#define NDIRS   (2 * 3) /* 2 * (number of dimensions) */

#define NELEMS(ary)	(sizeof(ary) / sizeof((ary)[0]))

int iround(double x);

enum {
    DIR_LEFT = 0,
    DIR_RIGHT,
    DIR_FRONT,
    DIR_BACK,
    DIR_BELOW,
    DIR_ABOVE
};

typedef double *** Array3Dd;
typedef int *** Array3Di;

typedef struct Point {
    double x;
    double y;
    double z;
} Point;

typedef struct Vector2d {
    double x;
    double y;
} Vector2d;

typedef struct Vector2d_ary {
    int size;
    Vector2d *ptr;
} Vector2d_ary;

Vector2d_ary *vector2d_ary_new(void);
void vector2d_ary_push(Vector2d_ary *self, Vector2d vector2d);

typedef struct iPoint {
    int i;
    int j;
    int k;
} iPoint;

iPoint get_ipoint(int i, int j, int k);
iPoint ipoint_offset(iPoint ipoint, int dirx, int diry, int dirz);
int ipoint_eq(iPoint ipoint1, iPoint ipoint2);
iPoint ipoint_add(iPoint ipoint1, iPoint ipoint2);
iPoint *ipoint_new(int i, int j, int k);

typedef iPoint *** Array3Dp;

typedef struct Coef {
    int index;
    double value;
} Coef;

typedef struct Coefs {
    Coef *coef[NDIRS];
    double coef0;
    double cnst;
} Coefs;

typedef Coefs *** Array3Dc;

typedef struct Obj Obj;

typedef struct Each {
    int size;
    iPoint *ipoint_ary;
    int index;
    Obj *obj;
} Each;

Each *each_new(int size, iPoint *ipoint_ary);
iPoint *each_each(Each *self);

typedef struct World {
    double x0;
    double y0;
    double z0;
    double xlen;
    double ylen;
    double zlen;
    double dx;
    double dy;
    double dz;
    int ni;
    int nj;
    int nk;
    Each *each;
} World;

World *world_new(double x, double y, double z,
	double xlen, double ylen, double zlen,
	double dx, double dy, double dz);

enum {
    AXIS_X,
    AXIS_Y,
    AXIS_Z
};

typedef struct Sweep {
    World *world;
    int axis;
    double len;
    Obj *obj;
    Each *each;
} Sweep;

Sweep *sweep_new(World *world, int axis, double len, Obj *obj);

typedef struct Edge {
    World *world;
    Obj *obj;
} Edge;

Edge *edge_new(World *world, Obj *obj);

typedef struct Rect {
    World *world;
    double x;
    double y;
    double z;
    int axis;
    double len1;
    double len2;
    Each *each;
} Rect;

Rect *rect_new(World *world, double x, double y, double z, int axis, double len1, double len2);

typedef struct Triangle_z {
    World *world;
    double x1;
    double y1;
    double dx;
    double dx2;
    double dy2;
    Each *each;
} Triangle_z;

typedef struct Triangle {
    World *world;
    int axis;
    double u1;
    double v1;
    int wi;
    double du2;
    double dv2;
    double du3;
    double dv3;
    Triangle_z *tr1;
    Triangle_z *tr2;
    Each *each;
} Triangle;

Triangle *triangle_new(World *world, double x1, double y1, double z1,
	int axis, double du2, double dv2, double du3, double dv3);

typedef struct Ellipse {
    World *world;
    double x;
    double y;
    double z;
    int axis;
    double ru;
    double rv;
    Each *each;
    int edge;
} Ellipse;

Ellipse *ellipse_new(World *world, double x, double y, double z, int axis, double ru, double rv);

typedef struct Circle {
    Ellipse *ellipse;
} Circle;

Circle *circle_new(World *world, double x, double y, double z, int axis, double r);

typedef struct Polygon {
    World *world;
    int axis;
    double w;
    Vector2d_ary *vector2d_ary;
    Each *each;
} Polygon;

Polygon *polygon_new(World *world, double x1, double y1, double z1,
	int axis, Vector2d_ary *dudv_ary);

typedef struct Box {
    World *world;
    Sweep *sweep;
} Box;

Box *box_new(World *world, double x, double y, double z, double xlen, double ylen, double zlen);

enum {
    OBJ_RECT,
    OBJ_TRIANGLE,
    OBJ_ELLIPSE,
    OBJ_CIRCLE,
    OBJ_POLYGON,
    OBJ_BOX,
    OBJ_SWEEP,
    OBJ_EDGE,
};

union uobj {
    Rect *rect;
    Triangle *triangle;
    Ellipse *ellipse;
    Circle *circle;
    Polygon *polygon;
    Box *box;
    Sweep *sweep;
    Edge *edge;
};

union uval {
    int i;
    double d;
};

struct Obj {
    int objtype;
    union uobj uobj;
    union uval uval;
};

Obj *obj_new(int objtype);
iPoint *obj_each_begin(Obj *self);
iPoint *obj_each(Obj *self);
void obj_offset(Obj *self);

typedef struct AryObj {
    Obj **ptr;
    int size;
} AryObj;

AryObj *aryobj_new(void);
void aryobj_push(AryObj *self, Obj *obj);

typedef struct Config {
    World *world;
    AryObj *active_obj_ary;
    AryObj *fix_obj_ary;
    AryObj *heat_obj_ary;
    AryObj *lambda_obj_ary;
} Config;

extern Config *config_parser;

typedef struct Sim {
    Config *config;
    iPoint dir_to_ipoint[NDIRS];
    World *world;
    int ni;
    int nj;
    int nk;
    Array3Di active_p_ary;
    Array3Dd *fix_ary;
    Array3Dd *heat_ary;
    Array3Dp *heat_ipoint_ary;
    Array3Dd lambda_ary;
    Array3Dd u;
    Array3Dc *coefs;
} Sim;

int sim_active_p(Sim *self, iPoint ipoint);
Sim *sim_new(void);
Array3Dd sim_calc(Sim *sim);

#endif
