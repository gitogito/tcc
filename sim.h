#ifndef _SIM_H_
#define _SIM_H_

#define NDIRS   (2 * 3) /* 2 * (number of dimensions) */

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
    int i;
    int j;
    int k;
} Point;

typedef Point *** Array3Dp;

Point get_point(int i, int j, int k);

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
    Point *point_ary;
    int index;
    Obj *obj;
} Each;

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

typedef struct Heatflow {
    int dir;
    double value;
} Heatflow;

Heatflow *heatflow_new(int dir, double value);

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
    double x2;
    double y2;
    Each *each;
} Triangle_z;

typedef struct Triangle {
    World *world;
    int axis;
    double u1;
    double v1;
    int wi;
    double u2;
    double v2;
    double u3;
    double v3;
    Triangle_z *tr1;
    Triangle_z *tr2;
    Each *each;
} Triangle;

Triangle *triangle_new(World *world, double x1, double y1, double z1,
	int axis, double u2, double v2, double u3, double v3);

typedef struct Box {
    World *world;
    Sweep *sweep;
} Box;

Box *box_new(World *world, double x, double y, double z, double xlen, double ylen, double zlen);

enum {
    OBJ_RECT,
    OBJ_TRIANGLE,
    OBJ_BOX,
    OBJ_SWEEP
};

union uobj {
    Rect *rect;
    Triangle *triangle;
    Box *box;
    Sweep *sweep;
};

union uval {
    int i;
    double d;
    Heatflow *h;
};

struct Obj {
    int objtype;
    union uobj uobj;
    union uval uval;
};

Obj *obj_new(int objtype);
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
    AryObj *heatflow_obj_ary;
    AryObj *lambda_obj_ary;
} Config;

extern Config *config_parser;

typedef struct Sim {
    Config *config;
    Point dir_to_point[NDIRS];
    World *world;
    int ni;
    int nj;
    int nk;
    Array3Di active_p_ary;
    Array3Dd fix_ary;
    Array3Dd *heatflow_ary;
    Array3Dp *heatflow_point_ary;
    Array3Dd lambda_ary;
    Array3Dd u;
    Array3Dc *coefs;
} Sim;

int sim_active_p(Sim *self, Point point);
Sim *sim_new(void);
Array3Dd sim_calc(Sim *sim);

#endif
