#ifndef _SIM_H_
#define _SIM_H_

#include <stdio.h>
#include "solvele.h"

#define NDIRS   	(2 * 3) /* 2 * (number of dimensions) */

#define NELEMS(ary)	(sizeof(ary) / sizeof((ary)[0]))

typedef short IP_TYPE;

IP_TYPE iround(double x);

enum {
    DIR_LEFT = 0,
    DIR_RIGHT,
    DIR_FRONT,
    DIR_BACK,
    DIR_BELOW,
    DIR_ABOVE
};

enum {
    REGION_ACTIVE = 1,	/* opt_r = 0 for calculation */
    REGION_FIX,
    REGION_FIXHEAT,
    REGION_HEAT,
    REGION_LAMBDA,
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
    int alloc_size;
    Vector2d *ptr;
} Vector2d_ary;

Vector2d_ary *vector2d_ary_new(void);
void vector2d_ary_push(Vector2d_ary *self, Vector2d vector2d);

typedef struct iPoint {
    IP_TYPE i;
    IP_TYPE j;
    IP_TYPE k;
} iPoint;

iPoint get_ipoint(IP_TYPE i, IP_TYPE j, IP_TYPE k);
iPoint ipoint_offset(iPoint *ipoint, int dirx, int diry, int dirz);
int ipoint_eq(iPoint *ipoint1, iPoint *ipoint2);
iPoint ipoint_add(iPoint *ipoint1, iPoint *ipoint2);
iPoint *ipoint_new(IP_TYPE i, IP_TYPE j, IP_TYPE k);

typedef iPoint *** Array3Dp;

typedef struct Obj Obj;

typedef struct iPoint_ary {
    int size;
    int alloc_size;
    iPoint *ptr;
} iPoint_ary;

iPoint_ary *ipoint_ary_new(void);
void ipoint_ary_push(iPoint_ary *self, iPoint ipoint);

typedef struct Each {
    iPoint_ary *ipoint_ary;
    int index;
    Obj *obj;
} Each;

Each *each_new(iPoint_ary *ipoint_ary);
int each_each(Each *self, iPoint **pp);

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
    IP_TYPE ni;
    IP_TYPE nj;
    IP_TYPE nk;
    Each *each;
} World;

extern World *world;

int world_to_index(iPoint *ipoint);
World *world_new(double x, double y, double z,
	double xlen, double ylen, double zlen,
	double dx, double dy, double dz);
int world_inside_p(iPoint *ipoint);

enum {
    AXIS_X,
    AXIS_Y,
    AXIS_Z
};

typedef struct Sweep {
    int axis;
    double len;
    Obj *obj;
    Each *each;
} Sweep;

Sweep *sweep_new(int axis, double len, Obj *obj);

typedef struct Edge {
    Obj *obj;
} Edge;

Edge *edge_new(Obj *obj);

typedef struct Rect {
    double x;
    double y;
    double z;
    int axis;
    double len1;
    double len2;
    Each *each;
    int edge;
} Rect;

Rect *rect_new(double x, double y, double z, int axis, double len1, double len2);

typedef struct Triangle_z {
    double x1;
    double y1;
    double dx;
    double dx2;
    double dy2;
    Each *each;
} Triangle_z;

typedef struct Triangle {
    int axis;
    double u1;
    double v1;
    IP_TYPE wi;
    double du2;
    double dv2;
    double du3;
    double dv3;
    Triangle_z *tr1;
    Triangle_z *tr2;
    Each *each;
} Triangle;

Triangle *triangle_new(double x1, double y1, double z1,
	int axis, double du2, double dv2, double du3, double dv3);

typedef struct Ellipse {
    double x;
    double y;
    double z;
    int axis;
    double ru;
    double rv;
    Each *each;
    int edge;
} Ellipse;

Ellipse *ellipse_new(double x, double y, double z, int axis, double ru, double rv);

typedef struct Circle {
    Ellipse *ellipse;
} Circle;

Circle *circle_new(double x, double y, double z, int axis, double r);

typedef struct Polygon {
    int axis;
    double w;
    Vector2d_ary *vector2d_ary;
    int (*rotate_p_func)(Vector2d va, Vector2d vb, Vector2d vc);
    Each *each;
} Polygon;

Polygon *polygon_new(double x1, double y1, double z1,
	int axis, Vector2d_ary *dudv_ary);

typedef struct Box {
    Sweep *sweep;
} Box;

Box *box_new(double x, double y, double z, double xlen, double ylen, double zlen);

typedef struct AryObj {
    Obj **ptr;
    int alloc_size;
    int size;
} AryObj;

AryObj *aryobj_new(void);
void aryobj_free(AryObj *self);
void aryobj_push(AryObj *self, Obj *obj);

typedef struct ObjAry {
    AryObj *aryobj;
    int each_obj_index;
} ObjAry;

ObjAry *objary_new(AryObj *aryobj);

enum {
    OBJ_RECT,
    OBJ_TRIANGLE,
    OBJ_ELLIPSE,
    OBJ_CIRCLE,
    OBJ_POLYGON,
    OBJ_BOX,
    OBJ_SWEEP,
    OBJ_EDGE,
    OBJ_OBJARY,
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
    ObjAry *objary;
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
int obj_each(Obj *self, iPoint **pp);
void obj_offset(Obj *self);

typedef struct Config {
    AryObj *active_obj_ary;
    AryObj *fix_obj_ary;
    AryObj *fixheat_obj_ary;
    AryObj *heat_obj_ary;
    AryObj *lambda_obj_ary;
} Config;

extern Config *config;

typedef struct Sim {
    iPoint dir_to_ipoint[NDIRS];
    Array3Di active_p_ary;
    Array3Dd *fix_ary;
    Array3Dd *fixheat_ary;
    Array3Dp *fixheat_ipoint_ary;
    Array3Dd *heat_ary;
    Array3Dd lambda_ary;
    Solvele *solver;
    char *fname;
} Sim;

extern Sim *sim;

int sim_active_p(iPoint *ipoint);
Sim *sim_new(char *fname);
double *sim_calc(void);

#endif
