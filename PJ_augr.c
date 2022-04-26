/******************************************************************************
 * Copyright (c) 2016, Marcin Ciura
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *****************************************************************************/

#define PJ_LIB__

#include <assert.h>
#include <stdarg.h>
#include "projects.h"

#define EPSILON 1e-7
#define ANGLE_DELTA 1e-7

PROJ_HEAD(augr, "pseudo-Authagraph")
    "\n\tMisc\n\tlat_1= lon_1= alpha= xshift= yscale= bb= cc= roundness=";

typedef struct {
    double x, y, z;
} V;

static double dot(V* u, V* v) {
    return u->x * v->x + u->y * v->y + u->z * v->z;
}

static double length(V* v) {
    return sqrt(v->x * v->x + v->y * v->y + v->z * v->z);
}

static V vadd(V* u, V* v) {
    V w = {u->x + v->x, u->y + v->y, u->z + v->z};
    return w;
}

static V scale(double a, V* v) {
    V w = {a * v->x, a * v->y, a * v->z};
    return w;
}

static V vnormalize(V* v) {
    V w = scale(1. / length(v), v);
    return w;
}

static V cross(V* u, V* v) {
    V w = {
        u->y * v->z - u->z * v->y,
        u->z * v->x - u->x * v->z,
        u->x * v->y - u->y * v->x
    };
    return w;
}

static double great_circle_distance(V* u, V* v) {
    V w = cross(u, v);
    return atan2(length(&w), dot(u, v));
}

static double edge_point_great_circle_distance(V* e, V* v) {
    return fabs(M_PI / 2 - great_circle_distance(e, v));
}

static double great_circle_distance_ratio(V* u, V* v, V* w, double duw) {
    double duv = great_circle_distance(u, v);
    double dvw = great_circle_distance(v, w);
    double numerator = duv;
    double smallest_error = fabs(duw - (duv + dvw));
    double error = fabs(duw - (duv + M_PI - dvw));

    if (error < smallest_error) {
        smallest_error = error;
    }
    duv = M_PI - duv;
    error = fabs(duw - (duv + dvw));
    if (error < smallest_error) {
        smallest_error = error;
        numerator = duv;
    }
    error = fabs(duw - (duv + M_PI - dvw));
    if (error < smallest_error) {
        smallest_error = error;
        numerator = duv;
    }
    return numerator / duw;
}

static V to_geocentric(double lat, double lon) {
    V v;
    v.x = cos(lat) * cos(lon);
    v.y = cos(lat) * sin(lon);
    v.z = sin(lat);
    return v;
}

static LP to_geodetic(V* v) {
    LP lp;
    lp.lam = atan2(v->y, v->x);
    lp.phi = atan2(v->z, hypot(v->x, v->y));
    return lp;
}

typedef struct {
    double r;
    V v;
} Q;

static double norm(Q* q) {
    return sqrt(
        q->r * q->r + q->v.x * q->v.x + q->v.y * q->v.y + q->v.z * q->v.z);
}

static Q conjugate(Q* q) {
    Q r = {q->r, scale(-1, &q->v)};
    return r;
}

static Q normalize(Q* q) {
    double n = norm(q);
    Q r = {q->r / n, scale(1. / n, &q->v)};
    return r;
}

static Q mul(Q* p, Q* q) {
    V v1 = cross(&p->v, &q->v);
    V v2 = scale(p->r, &q->v);
    V v3 = scale(q->r, &p->v);
    V v4 = vadd(&v2, &v3);
    Q r = {p->r * q->r - dot(&p->v, &q->v), vadd(&v1, &v4)};
    return r;
}

static Q rotation_from_vectors(V* u, V* v) {
    Q r = {dot(u, v) + 1., cross(u, v)};
    return normalize(&r);
}

static V rotate(Q* rot, Q* v, Q* rot_1) {
    Q tmp = mul(rot, v);
    return mul(&tmp, rot_1).v;
}

struct pj_opaque {
    double phi1;
    double lam1;
    double alpha;
    double xshift;
    double yscale;
    double bb;
    double cc;
    double roundness;
    Q vertices[4];
    V faces[4][3];
    XY triangles[4][3];
    V edge_01;
    V edge_02;
    V edge_03;
    V edge_12;
    V edge_13;
    V edge_23;
    double edge_01_dot_2;
    double edge_01_dot_3;
    double edge_02_dot_1;
    double edge_02_dot_3;
    double edge_03_dot_1;
    double edge_03_dot_2;
    double edge_12_dot_0;
    double edge_12_dot_3;
    double edge_13_dot_0;
    double edge_13_dot_2;
    double edge_23_dot_0;
    double edge_23_dot_1;
    double distance_01;
    double distance_02;
    double distance_03;
    double distance_12;
    double distance_13;
    double distance_23;
};

struct xy_face {
    XY xy;
    int face;
};

static int find_face(V* v, struct pj_opaque* PO) {
    if (dot(&PO->edge_01, v) * PO->edge_01_dot_3 <= EPSILON &&
        dot(&PO->edge_02, v) * PO->edge_02_dot_3 <= EPSILON &&
        dot(&PO->edge_12, v) * PO->edge_12_dot_3 <= EPSILON) {
      return 3;
    }
    if (dot(&PO->edge_01, v) * PO->edge_01_dot_2 <= EPSILON &&
        dot(&PO->edge_03, v) * PO->edge_03_dot_2 <= EPSILON &&
        dot(&PO->edge_13, v) * PO->edge_13_dot_2 <= EPSILON) {
      return 2;
    }
    if (dot(&PO->edge_02, v) * PO->edge_02_dot_1 <= EPSILON &&
        dot(&PO->edge_03, v) * PO->edge_03_dot_1 <= EPSILON &&
        dot(&PO->edge_23, v) * PO->edge_23_dot_1 <= EPSILON) {
      return 1;
    }
    if (dot(&PO->edge_12, v) * PO->edge_12_dot_0 <= EPSILON &&
        dot(&PO->edge_13, v) * PO->edge_13_dot_0 <= EPSILON &&
        dot(&PO->edge_23, v) * PO->edge_23_dot_0 <= EPSILON) {
      return 0;
    }
    assert(0);
}

static XY intersection(XY* p1, XY* p2, XY* p3, XY* p4) {
    double denominator, d12, d34;
    XY p;

    if (p2->x != p2->x) {
        return *p1;
    }
    if (p4->x != p4->x) {
        return *p3;
    }
    denominator =
        (p1->x - p2->x) * (p3->y - p4->y) -
        (p1->y - p2->y) * (p3->x - p4->x);
    if (denominator == 0) {
        return *p1;
    }
    d12 = p1->x * p2->y - p1->y * p2->x;
    d34 = p3->x * p4->y - p3->y * p4->x;
    p.x = (d12 * (p3->x - p4->x) - (p1->x - p2->x) * d34) / denominator;
    p.y = (d12 * (p3->y - p4->y) - (p1->y - p2->y) * d34) / denominator;
    return p;
}

static XY incenter(XY* p1, XY* p2, XY* p3) {
    double d12 = hypot(p2->x - p1->x, p2->y - p1->y);
    double d13 = hypot(p3->x - p1->x, p3->y - p1->y);
    double d23 = hypot(p3->x - p2->x, p3->y - p2->y);
    double perimeter = d12 + d13 + d23;
    XY p;
    if (perimeter == 0) {
        return *p1;
    }
    p.x = (d12 * p3->x + d13 * p2->x + d23 * p1->x) / perimeter;
    p.y = (d12 * p3->y + d13 * p2->y + d23 * p1->y) / perimeter;
    return p;
}

static struct xy_face project_internal(V* v, struct pj_opaque* PO) {
    struct xy_face result;
    V *v0, *v1, *v2;
    V *edge01, *edge02, *edge12;
    V tmp, tmp1;
    double d01, d02, d12;
    double w01, w12, w20;
    XY *triangle;
    XY p01, p12, p20;
    XY p0112, p0120, p1220;

    switch (result.face = find_face(v, PO)) {
    case 0:
        v0 = &PO->vertices[1].v;
        v1 = &PO->vertices[2].v;
        v2 = &PO->vertices[3].v;
        d01 = PO->distance_12;
        d02 = PO->distance_13;
        d12 = PO->distance_23;
        edge01 = &PO->edge_12;
        edge02 = &PO->edge_13;
        edge12 = &PO->edge_23;
        break;
    case 1:
        v0 = &PO->vertices[0].v;
        v1 = &PO->vertices[2].v;
        v2 = &PO->vertices[3].v;
        d01 = PO->distance_02;
        d02 = PO->distance_03;
        d12 = PO->distance_23;
        edge01 = &PO->edge_02;
        edge02 = &PO->edge_03;
        edge12 = &PO->edge_23;
        break;
    case 2:
        v0 = &PO->vertices[0].v;
        v1 = &PO->vertices[1].v;
        v2 = &PO->vertices[3].v;
        d01 = PO->distance_01;
        d02 = PO->distance_03;
        d12 = PO->distance_13;
        edge01 = &PO->edge_01;
        edge02 = &PO->edge_03;
        edge12 = &PO->edge_13;
        break;
    case 3:
        v0 = &PO->vertices[0].v;
        v1 = &PO->vertices[1].v;
        v2 = &PO->vertices[2].v;
        d01 = PO->distance_01;
        d02 = PO->distance_02;
        d12 = PO->distance_12;
        edge01 = &PO->edge_01;
        edge02 = &PO->edge_02;
        edge12 = &PO->edge_12;
        break;
    default:
        assert(0);
    }

    tmp = cross(v2, v);
    tmp1 = cross(edge01, &tmp);
    tmp = vnormalize(&tmp1);
    w01 = great_circle_distance_ratio(v0, &tmp, v1, d01);

    tmp = cross(v0, v);
    tmp1 = cross(edge12, &tmp);
    tmp = vnormalize(&tmp1);
    w12 = great_circle_distance_ratio(v1, &tmp, v2, d12);

    tmp = cross(v1, v);
    tmp1 = cross(edge02, &tmp);
    tmp = vnormalize(&tmp1);
    w20 = great_circle_distance_ratio(v2, &tmp, v0, d02);

    triangle = PO->triangles[result.face];
    p01.x = triangle[0].x + w01 * (triangle[1].x - triangle[0].x);
    p01.y = triangle[0].y + w01 * (triangle[1].y - triangle[0].y);
    p12.x = triangle[1].x + w12 * (triangle[2].x - triangle[1].x);
    p12.y = triangle[1].y + w12 * (triangle[2].y - triangle[1].y);
    p20.x = triangle[2].x + w20 * (triangle[0].x - triangle[2].x);
    p20.y = triangle[2].y + w20 * (triangle[0].y - triangle[2].y);

    p0112 = intersection(&triangle[2], &p01, &triangle[0], &p12);
    p0120 = intersection(&triangle[2], &p01, &triangle[1], &p20);
    p1220 = intersection(&triangle[0], &p12, &triangle[1], &p20);

    result.xy = incenter(&p0112, &p0120, &p1220);
    return result;
}

static XY fix_xy(XY* xy, struct pj_opaque* PO) {
    xy->x -= PO->xshift;
    if (xy->x < 0.) {
        xy->x += 4.;
    } else if (xy->x > 4.) {
        xy->x -= 4.;
    }
    xy->y *= PO->yscale;
    return *xy;
}

static void fix_new_xy_face(
        struct xy_face* xy_face, struct xy_face* new_xy_face,
        struct pj_opaque* PO) {
    if (abs(xy_face->face - new_xy_face->face) == 2) {
        if (xy_face->face + new_xy_face->face == 4) {
            new_xy_face->xy.x = 2. * PO->triangles[1][0].x - new_xy_face->xy.x;
            new_xy_face->xy.y = 2. * PO->triangles[1][0].y - new_xy_face->xy.y;
        } else {
            new_xy_face->xy.x = 4. - new_xy_face->xy.x;
            new_xy_face->xy.y = -new_xy_face->xy.y;
        }
    }
    if (fabs(xy_face->xy.x - new_xy_face->xy.x) > 2.) {
        if (fabs(xy_face->xy.x + 4. - new_xy_face->xy.x) <
            fabs(xy_face->xy.x - 4. - new_xy_face->xy.x)) {
            new_xy_face->xy.x -= 4.;
        } else {
            new_xy_face->xy.x += 4.;
        }
    }
}

static struct xy_face project_near_point(
        V* v, struct xy_face* xy_face, struct pj_opaque* PO) {
    struct xy_face new_xy_face = project_internal(v, PO);

    fix_new_xy_face(xy_face, &new_xy_face, PO);
    return new_xy_face;
}

static double roundness_border(
        V* v, V* va, V* vb, V* edge, struct pj_opaque* PO, Q* vc) {
    double edge_length = great_circle_distance(va, vb);
    double cast_length;
    V tmp = cross(v, edge);
    V tmp1 = cross(edge, &tmp);

    vc->r = 0.;
    vc->v = vnormalize(&tmp1);
    cast_length = great_circle_distance(va, &vc->v);
    return PO->roundness * cast_length * (edge_length - cast_length) / edge_length;
}

static Q perpendicular_rotation(V* v, V* edge, double angle) {
    V tmp = cross(v, edge);
    Q result;

    tmp = vnormalize(&tmp);
    result.r = cos(angle / 2);
    result.v = scale(sin(angle / 2), &tmp);
    return result;
}

static void interpolate(V* v, struct xy_face* xy_face, struct pj_opaque* PO) {
    V *v0, *v1, *v2, *va, *vb;
    V *edge_01, *edge_02, *edge_12, *edge;
    Q vc;
    double h01, h02, h12, minh, border;

    switch (xy_face->face) {
    case 0:
        v0 = &PO->vertices[1].v;
        v1 = &PO->vertices[2].v;
        v2 = &PO->vertices[3].v;
        edge_01 = &PO->edge_12;
        edge_02 = &PO->edge_13;
        edge_12 = &PO->edge_23;
        break;
    case 1:
        v0 = &PO->vertices[0].v;
        v1 = &PO->vertices[2].v;
        v2 = &PO->vertices[3].v;
        edge_01 = &PO->edge_02;
        edge_02 = &PO->edge_03;
        edge_12 = &PO->edge_23;
        break;
    case 2:
        v0 = &PO->vertices[0].v;
        v1 = &PO->vertices[1].v;
        v2 = &PO->vertices[3].v;
        edge_01 = &PO->edge_01;
        edge_02 = &PO->edge_03;
        edge_12 = &PO->edge_13;
        break;
    case 3:
        v0 = &PO->vertices[0].v;
        v1 = &PO->vertices[1].v;
        v2 = &PO->vertices[2].v;
        edge_01 = &PO->edge_01;
        edge_02 = &PO->edge_02;
        edge_12 = &PO->edge_12;
        break;
    default:
        assert(0);
    }
    h01 = edge_point_great_circle_distance(edge_01, v);
    h02 = edge_point_great_circle_distance(edge_02, v);
    h12 = edge_point_great_circle_distance(edge_12, v);

    minh = h01;
    edge = edge_01;
    va = v0;
    vb = v1;
    if (h02 < minh) {
        minh = h02;
        edge = edge_02;
        vb = v2;
    }
    if (h12 < minh) {
        minh = h12;
        edge = edge_12;
        va = v1;
        vb = v2;
    }
    border = roundness_border(v, va, vb, edge, PO, &vc);
    if (minh < border) {
        Q rot, rot_1;
        V tmpv;
        XY xys[4], xya, xyb;
        double s = minh / border;
        double t = (1. + sin(M_PI / 2 * s)) / 2.;

        rot  = perpendicular_rotation(&vc.v, edge, border);
        rot_1 = conjugate(&rot);
        tmpv = rotate(&rot, &vc, &rot_1);
        if (great_circle_distance(&tmpv, v) < border) {
            minh = -minh;
            t = 1. - t;
        }
        xys[0] = project_near_point(&tmpv, xy_face, PO).xy;
        tmpv = rotate(&rot_1, &vc, &rot);
        xys[3] = project_near_point(&tmpv, xy_face, PO).xy;

        rot  = perpendicular_rotation(&vc.v, edge, border - ANGLE_DELTA);
        rot_1 = conjugate(&rot);
        tmpv = rotate(&rot, &vc, &rot_1);
        xys[1] = project_near_point(&tmpv, xy_face, PO).xy;
        tmpv = rotate(&rot_1, &vc, &rot);
        xys[2] = project_near_point(&tmpv, xy_face, PO).xy;

        xya.x = xys[0].x +
                (xys[1].x - xys[0].x) * (border + minh) / ANGLE_DELTA / 2.;
        xya.y = xys[0].y +
                (xys[1].y - xys[0].y) * (border + minh) / ANGLE_DELTA / 2.;
        xyb.x = xys[3].x +
                (xys[2].x - xys[3].x) * (border - minh) / ANGLE_DELTA / 2.;
        xyb.y = xys[3].y +
                (xys[2].y - xys[3].y) * (border - minh) / ANGLE_DELTA / 2.;

        xy_face->xy.x = (1. - t) * xya.x + t * xyb.x;
        xy_face->xy.y = (1. - t) * xya.y + t * xyb.y;
    }
}

static XY s_forward (LP lp, PJ *P) {
    struct pj_opaque* PO = P->opaque;
    V v = to_geocentric(lp.phi, lp.lam);
    struct xy_face xy_face = project_internal(&v, PO);

    if (PO->roundness != 0.) {
        interpolate(&v, &xy_face, PO);
    }
    return fix_xy(&xy_face.xy, PO);
}

static int on_the_left(XY* xy1, XY* xy2, XY* xya) {
    double dx = xy2->x - xy1->x;
    double dy = xy2->y - xy1->y;

    return (dy * (xya->x - xy1->x) - dx * (xya->y - xy1->y) < 0.);
}

static int find_face_on_plane(XY* xy, struct pj_opaque* PO) {
    if (on_the_left(&PO->triangles[0][2], &PO->triangles[0][1], xy)) {
        xy->x += 4.;
        return 1;
    }
    if (on_the_left(&PO->triangles[0][0], &PO->triangles[0][1], xy)) {
        return 0;
    }
    if (on_the_left(&PO->triangles[2][1], &PO->triangles[2][0], xy)) {
        return 3;
    }
    if (on_the_left(&PO->triangles[1][2], &PO->triangles[1][0], xy)) {
        return 2;
    }
    if (on_the_left(&PO->triangles[1][2], &PO->triangles[1][1], xy)) {
        return 1;
    }
    assert(0);
}

static XY diff2(XY* u, XY* v) {
    XY w = {u->x - v->x, u->y - v->y};
    return w;
}

static double cross2(XY* u, XY* v) {
    return u->x * v->y - u->y * v->x;
}

static LP gnomonic(struct xy_face* xy_face, struct pj_opaque* PO) {
    XY u, v, w;
    double n;
    double b[3];
    XY *planar_vertices;
    V *vertices;
    V image, tmp0, tmp1, tmp2;

    xy_face->xy.x += PO->xshift;
    if (xy_face->xy.x > 4.) {
        xy_face->xy.x -= 4.;
    }
    xy_face->xy.y /= PO->yscale;
    xy_face->face = find_face_on_plane(&xy_face->xy, PO);

    /* https://www.cs.ubc.ca/~heidrich/Papers/JGT.05.pdf */
    planar_vertices = PO->triangles[xy_face->face];
    u = diff2(&planar_vertices[1], &planar_vertices[0]);
    v = diff2(&planar_vertices[2], &planar_vertices[0]);
    w = diff2(&xy_face->xy, &planar_vertices[0]);
    n = cross2(&u, &v);
    b[2] = cross2(&u, &w) / n;
    b[1] = cross2(&w, &v) / n;
    b[0] = 1. - b[1] - b[2];
    vertices = PO->faces[xy_face->face];
    tmp0 = scale(b[0], &vertices[0]);
    tmp1 = scale(b[1], &vertices[1]);
    tmp2 = scale(b[2], &vertices[2]);
    image = vadd(&tmp0, &tmp1);
    image = vadd(&image, &tmp2);
    return to_geodetic(&image);
}

static XY round_project(
        LP* lp, struct xy_face* xy_face, struct pj_opaque* PO) {
    LP lp0 = *lp;
    V v;
    struct xy_face new_xy_face;

    assert(lp0.phi == lp0.phi && lp0.lam == lp0.lam);
    v = to_geocentric(lp0.phi, lp0.lam);
    if (PO->roundness == 0.) {
        new_xy_face = project_near_point(&v, xy_face, PO);
    } else {
        new_xy_face = project_internal(&v, PO);
        interpolate(&v, &new_xy_face, PO);
        fix_new_xy_face(xy_face, &new_xy_face, PO);
    }
    assert(new_xy_face.xy.x == new_xy_face.xy.x &&
           new_xy_face.xy.y == new_xy_face.xy.y);
    new_xy_face.xy.x -= xy_face->xy.x;
    new_xy_face.xy.y -= xy_face->xy.y;
    return new_xy_face.xy;
}

static void solve2(double A[2][2], XY* xy, LP* s) {
    double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
    s->lam = (xy->x * A[1][1] - A[0][1] * xy->y) / det;
    s->phi = (A[0][0] * xy->y - xy->x * A[1][0]) / det;
}

static LP sub2(LP* lp0, LP* lp1) {
    LP lp;
    lp.lam = lp0->lam - lp1->lam;
    lp.phi = lp0->phi - lp1->phi;
    if (lp.phi > M_PI / 2) {
        lp.phi = M_PI - lp.phi;
        lp.lam += M_PI;
    } else if (lp.phi < -M_PI / 2) {
        lp.phi = -M_PI - lp.phi;
        lp.lam += M_PI;
    }
    if (lp.lam > M_PI) {
        lp.lam -= 2 * M_PI;
    } else if (lp.lam < -M_PI) {
        lp.lam += 2 * M_PI;
    }
    return lp;
}

static double dot_lp2(LP* lp0, LP* lp1) {
    return lp0->lam * lp1->lam + lp0->phi * lp1->phi;
}

static double norm_xy2(XY* xy) {
    return sqrt(xy->x * xy->x + xy->y * xy->y);
}

static void approximate_jacobian(
        LP* lp0, XY* xy0, double A[2][2],
        struct xy_face* xy_face, struct pj_opaque* PO) {
    LP lp;
    XY xy;
    double delta = ANGLE_DELTA * cos(lp0->phi);

    lp.lam = lp0->lam + ANGLE_DELTA;
    lp.phi = lp0->phi;
    xy = round_project(&lp, xy_face, PO);
    A[0][0] = (xy.x - xy0->x) / ANGLE_DELTA;
    A[1][0] = (xy.y - xy0->y) / ANGLE_DELTA;

    lp.lam = lp0->lam;
    lp.phi = lp0->phi + delta;
    xy = round_project(&lp, xy_face, PO);
    A[0][1] = (xy.x - xy0->x) / delta;
    A[1][1] = (xy.y - xy0->y) / delta;
}

static void update_jacobian(double df, LP* s, double d, double V[2]) {
    double tmp = df - V[0] * s->lam - V[1] * s->phi;
    V[0] -= tmp * s->lam / d;
    V[1] -= tmp * s->phi / d;
}

static LP broyden(XY* xy, struct pj_opaque* PO) {
    struct xy_face xy_face;
    LP lp0;
    XY xy0;
    double g0;
    double A[2][2];
    int first;
    int iter;

    xy_face.xy = *xy;
    lp0 = gnomonic(&xy_face, PO);
    xy0 = round_project(&lp0, &xy_face, PO);
    g0 = norm_xy2(&xy0);
    approximate_jacobian(&lp0, &xy0, A, &xy_face, PO);
    first = 1;
    for (iter = 0; iter < 100; iter++) {
        LP lp1, s;
        XY xy1;
        double g1, d;

        solve2(A, &xy0, &s);
        if (s.lam < -1.) {
            s.lam = -1.;
        } else if (s.lam > +1.) {
            s.lam = +1.;
        }
        if (s.phi < -1.) {
            s.phi = -1.;
        } else if (s.phi > +1.) {
            s.phi = +1.;
        }
        lp1 = sub2(&lp0, &s);
        xy1 = round_project(&lp1, &xy_face, PO);
#if 0
        printf(
            "%d %lf %lf -> %lf %lf %lf %lf [%lf %lf %lf %lf] [%lf %lf]\n",
            iter, lp0.phi, lp0.lam, xy0.x, xy0.y, xy1.x, xy1.y,
            A[0][0], A[0][1], A[1][0], A[1][1], s.lam, s.phi);
#endif
        g1 = norm_xy2(&xy1);
        if (g1 < EPSILON) {
            break;
        }
        if (g1 > g0 && !first) {
            approximate_jacobian(&lp0, &xy0, A, &xy_face, PO);
            first = 1;
            continue;
        }
        d = dot_lp2(&s, &s);
        update_jacobian(xy1.x + xy0.x, &s, d, A[0]);
        update_jacobian(xy1.y + xy0.y, &s, d, A[1]);

        lp0 = lp1;
        xy0 = xy1;
        g0 = g1;
        first = 0;
    }
#if 0
    printf("iteration %d\n", iter);
#endif
    return lp0;
}

static LP s_inverse (XY xy, PJ *P) {
    return broyden(&xy, P->opaque);
}

static void *freeup_new (PJ *P) {
    if (0==P)
        return 0;
    if (0==P->opaque)
        return pj_dealloc (P);
    pj_dealloc (P->opaque);
    return pj_dealloc(P);
}

static void freeup (PJ *P) {
    freeup_new (P);
}

static void setup(PJ *P) {
    static V north_pole = {0, 0, 1};
    static V first_meridian = {1, 0, 0};

    struct pj_opaque *PO = P->opaque;
    Q rot, rot_1, tmp, tmp1;
    V v0, vtmp;
    double side_a, side_b, side_c, semiperimeter, altitude, foot_b;
    XY projected[6];
    int i, j;

    tmp.r = 0.;
    tmp.v.x = +1.;
    tmp.v.y = +PO->bb;
    tmp.v.z = +PO->cc;
    PO->vertices[0] = normalize(&tmp);
    tmp.v.x = -1.;
    tmp.v.z = -PO->cc;
    PO->vertices[1] = normalize(&tmp);
    tmp.v.x = +1.;
    tmp.v.y = -PO->bb;
    PO->vertices[2] = normalize(&tmp);
    tmp.v.x = -1.;
    tmp.v.z = +PO->cc;
    PO->vertices[3] = normalize(&tmp);

    rot = rotation_from_vectors(&PO->vertices[0].v, &north_pole);
    rot_1 = conjugate(&rot);
    for (i = 0; i < 4; i++) {
        PO->vertices[i].v = rotate(&rot, &PO->vertices[i], &rot_1);
    }

    v0 = to_geocentric(PO->phi1, PO->lam1);
    tmp = rotation_from_vectors(&PO->vertices[0].v, &v0);
    vtmp.x = cos(PO->alpha);
    vtmp.y = sin(PO->alpha);
    vtmp.z = 0.;
    tmp1 = rotation_from_vectors(&first_meridian, &vtmp);
    rot = mul(&tmp, &tmp1);
    rot_1 = conjugate(&rot);
    for (i = 0; i < 4; i++) {
        PO->vertices[i].v = rotate(&rot, &PO->vertices[i], &rot_1);
    }

    PO->edge_01 = cross(&PO->vertices[0].v, &PO->vertices[1].v);
    PO->edge_02 = cross(&PO->vertices[0].v, &PO->vertices[2].v);
    PO->edge_03 = cross(&PO->vertices[0].v, &PO->vertices[3].v);
    PO->edge_12 = cross(&PO->vertices[1].v, &PO->vertices[2].v);
    PO->edge_13 = cross(&PO->vertices[1].v, &PO->vertices[3].v);
    PO->edge_23 = cross(&PO->vertices[2].v, &PO->vertices[3].v);

    PO->edge_01_dot_2 = dot(&PO->edge_01, &PO->vertices[2].v);
    PO->edge_01_dot_3 = dot(&PO->edge_01, &PO->vertices[3].v);
    PO->edge_02_dot_1 = dot(&PO->edge_02, &PO->vertices[1].v);
    PO->edge_02_dot_3 = dot(&PO->edge_02, &PO->vertices[3].v);
    PO->edge_03_dot_1 = dot(&PO->edge_03, &PO->vertices[1].v);
    PO->edge_03_dot_2 = dot(&PO->edge_03, &PO->vertices[2].v);
    PO->edge_12_dot_0 = dot(&PO->edge_12, &PO->vertices[0].v);
    PO->edge_12_dot_3 = dot(&PO->edge_12, &PO->vertices[3].v);
    PO->edge_13_dot_0 = dot(&PO->edge_13, &PO->vertices[0].v);
    PO->edge_13_dot_2 = dot(&PO->edge_13, &PO->vertices[2].v);
    PO->edge_23_dot_0 = dot(&PO->edge_23, &PO->vertices[0].v);
    PO->edge_23_dot_1 = dot(&PO->edge_23, &PO->vertices[1].v);

    PO->distance_01 = great_circle_distance(
        &PO->vertices[0].v, &PO->vertices[1].v);
    PO->distance_02 = great_circle_distance(
        &PO->vertices[0].v, &PO->vertices[2].v);
    PO->distance_03 = great_circle_distance(
        &PO->vertices[0].v, &PO->vertices[3].v);
    PO->distance_12 = great_circle_distance(
        &PO->vertices[1].v, &PO->vertices[2].v);
    PO->distance_13 = great_circle_distance(
        &PO->vertices[1].v, &PO->vertices[3].v);
    PO->distance_23 = great_circle_distance(
        &PO->vertices[2].v, &PO->vertices[3].v);

    side_a = sqrt(PO->bb * PO->bb + PO->cc * PO->cc);
    side_b = sqrt(1. + PO->cc * PO->cc);
    side_c = sqrt(1. + PO->bb * PO->bb);
    semiperimeter = (side_a + side_b + side_c) / 2.;
    altitude = 4. / (side_a * side_a) * sqrt(
        semiperimeter *
        (semiperimeter - side_a) *
        (semiperimeter - side_b) *
        (semiperimeter - side_c));
    if (side_b >= side_c) {
        foot_b = sqrt(
            4. * (side_b / side_a) * (side_b / side_a) - altitude * altitude);
    } else {
        foot_b = 2. - sqrt(
            4. * (side_c / side_a) * (side_c / side_a) - altitude * altitude);
    }
    projected[0].x = 0.;
    projected[0].y = 0.;
    projected[1].x = foot_b;
    projected[1].y = altitude;
    projected[2].x = 2.;
    projected[2].y = 0.;
    projected[3].x = 2. + foot_b;
    projected[3].y = altitude;
    projected[4].x = 4.;
    projected[4].y = 0.;
    projected[5].x = 4. + foot_b;
    projected[5].y = altitude;
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 3; j++) {
            static int vertex_indices[6] = {3, 2, 1, 0, 3, 2};
            static int indices[4][3] =
                {{2, 1, 0}, {3, 5, 4}, {3, 2, 4}, {3, 2, 1}};

            PO->faces[i][j] = PO->vertices[vertex_indices[indices[i][j]]].v;
            PO->triangles[i][j] = projected[indices[i][j]];
        }
    }

    P->fwd = s_forward;
    P->inv = s_inverse;
}

PJ *PROJECTION(augr) {
    struct pj_opaque *PO = pj_calloc(1, sizeof (struct pj_opaque));

    if (0 == PO)
        return freeup_new(P);
    P->opaque = PO;
    PO->phi1 = pj_param(P->ctx, P->params, "tlat_1").i ?
        pj_param(P->ctx, P->params, "rlat_1").f : +73.10012877 * DEG_TO_RAD;
    PO->lam1 = pj_param(P->ctx, P->params, "tlon_1").i ?
        pj_param(P->ctx, P->params, "rlon_1").f : +120.95842104 * DEG_TO_RAD;
    PO->alpha = pj_param(P->ctx, P->params, "talpha").i ?
        pj_param(P->ctx, P->params, "ralpha").f : -23.55970851 * DEG_TO_RAD;
    PO->xshift = pj_param(P->ctx, P->params, "txshift").i ?
        pj_param(P->ctx, P->params, "dxshift").f : 1.6;
    PO->yscale = pj_param(P->ctx, P->params, "tyscale").i ?
        pj_param(P->ctx, P->params, "dyscale").f : 1.03107608;
    PO->bb = pj_param(P->ctx, P->params, "tbb").i ?
        pj_param(P->ctx, P->params, "dbb").f : 0.76372191;
    PO->cc = pj_param(P->ctx, P->params, "tcc").i ?
        pj_param(P->ctx, P->params, "dcc").f : 0.96292104;
    PO->roundness = tan((pj_param(P->ctx, P->params, "troundness").i ?
        pj_param(P->ctx, P->params, "rroundness").f : 20.0) * DEG_TO_RAD);
    setup(P);
    return P;
}

int pj_augr_selftest (void) {
    return 0;
}
