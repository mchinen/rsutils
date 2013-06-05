/*
 *  rs_geometry.h
 *
 *  Created by Administrator on 12/21/10.
 *  Copyright 2010 Roughsoft. All rights reserved.
 *
 */

#ifndef __FILTER_GEOMETRY_H__
#define __FILTER_GEOMETRY_H__

//get M_PI
#define _USE_MATH_DEFINES
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

enum RSCubeOverlapResult {
    eRSCubeOverlapResultNone       = 0x0,
    //hit the top side of a
    eRSCubeOverlapResultTop        = 0x0001,
    //hit the bottom side of A
    eRSCubeOverlapResultBottom     = 0x0002,
    //hit to left side of A
    eRSCubeOverlapResultLeft       = 0x0004,
    //hit the right side of A
    eRSCubeOverlapResultRight      = 0x0008,
    eRSCubeOverlapResultEverywhere = 0x000F
};

// this is useful when a result from filter_cube_overlaps is needed to be used from
// the perspective of the b cube, as the returned value is from the a cube.
void rs_overlap_result_flip(int *inoutResult);

typedef struct RSVec3 {
    float x;
    float y;
    float z;
} RSVec3;

typedef struct RSCube {
    RSVec3 pos;
    RSVec3 size;
} RSCube;

// gets a radians based angle for the result based on its value
// the radian result is arbitrary with top being 0, left is PI/2. 
float rs_overlap_result_radians(int result);

//non-zero if intersects, zero otherwise.
int rs_cube_overlaps(RSCube *a, RSCube *b, int *outResult);

// returns a normalized vector for an cube overlap result code
void rs_overlap_result_normal(int result, RSVec3* outNormal);

void rs_vec_zero(RSVec3 *a);
void rs_vec_add(RSVec3 *a, RSVec3 *b, RSVec3 *out);
void rs_vec_scale(RSVec3 *a, float scalar, RSVec3 *out);
//subtraction of a and b.
//@param out can point to a or b (or a different RSVec3) and the value will be overwritten
void rs_vec_subtract(RSVec3 *a, RSVec3 *b, RSVec3 *out);

//euclidean distance.
float rs_vec_length(RSVec3 *a);
float rs_vec_dist(RSVec3 *a, RSVec3 *b);

//standard dot product
float rs_vec_dot(RSVec3 *a, RSVec3 *b);

void rs_vec_cross(RSVec3 *a, RSVec3 *b, RSVec3 *out);

//normalizes a vector to a distance of 1.0
//returns the original euclidean distance of in.
//in case it was zero then out is not written to.
float rs_vec_norm(RSVec3 *in, RSVec3 *out);

void rs_vec_abs(RSVec3 *in, RSVec3 *out);

void rs_vec_neg(RSVec3 *in, RSVec3 *out);

void rs_vec_rotatexy(RSVec3 *in, float radAngle, RSVec3 *out);

float rs_vec_anglexy(RSVec3 *a, RSVec3 *b, RSVec3 *ref);

void rs_cube_rotatexy_bounds(RSCube* in, float radangle, RSCube* out);

/// does collision detection for rotated cubes
/// outResult contains the side of the cube a that was hit
int rs_cube_overlaps_rotated(RSCube *a, float aAngleRad,
                                 RSCube *b, float bAngleRad,
                                 int *outResult);
/// projects an input vector in onto an axis vector
/// axis does not need to be normalized, but can be.
/// result is in normal coordinate distance scale.
float rs_vec_project_axis_vec(RSVec3 *in, RSVec3* axis);

/// takes an axis aligned cube and rotates the corners by radangle
void rs_cube_rotated_points(RSCube* in, float radangle, RSVec3 points[4]);

#ifdef __cplusplus
};
#endif

#endif
