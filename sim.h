#ifndef _SIM_H_
#define _SIM_H_

typedef double *** Array3Dd;
typedef int *** Array3Di;

typedef struct Point {
    int i;
    int j;
    int k;
} Point;

typedef struct Each {
    int size;
    Point *point_ary;
    int index;
} Each;

typedef struct World {
    double x;
    double y;
    double z;
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

enum {
    AXIS_X,
    AXIS_Y,
    AXIS_Z
};

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

typedef struct Box {
    World *world;
    Rect *rect;
    double sweeplen;
    Each * each;
} Box;

enum {
    OBJ_RECT,
    OBJ_BOX
};

union uobj {
    Rect *rect;
    Box *box;
};

enum {
    OBJVAL_I,
    OBJVAL_D,
};

union uval {
    int i;
    double d;
} uval;

typedef struct Obj {
    int objtype;
    union uobj uobj;
    int valtype;
    union uval uval;
} Obj;

typedef struct AryObj {
    Obj **ptr;
    int size;
} AryObj;

typedef struct Config {
    World *world;
    AryObj *active_obj_ary;
    AryObj *fix_obj_ary;
    AryObj *heatflow_obj_ary;
    AryObj *lambda_obj_ary;
} Config;

typedef struct Sim {
    Config *config;
    World *world;
    int ni;
    int nj;
    int nk;
    Array3Di active_p_ary;
    Array3Dd fix_ary;
} Sim;

Sim *sim_new(void);
Array3Dd sim_calc(Sim *sim);

#endif
